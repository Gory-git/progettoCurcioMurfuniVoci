section .note.GNU-stack noalloc noexec nowrite progbits

section .data
align 16
one_float:    dd 1.0, 1.0, 1.0, 1.0    ; Vector of four 1.0 values for SSE
zero_float:   dd 0.0, 0.0, 0.0, 0.0    ; Vector of four 0.0 values for SSE

section .text
[BITS 32]
global computeScores
extern copy_vector
extern alloc_vector
extern memset

; SSE-optimized implementation of computeScores function
; VECTOR computeScores(MATRIX tranMat, float alfaB, int maxBias, VECTOR d, int numPages)
%define TYPE_SIZE 4

computeScores:
    ; Standard prologue
    push ebp
    mov ebp, esp
    sub esp, 48                    ; Space for local variables
    push ebx                       ; Save callee-saved registers
    push esi
    push edi

    ; Setup local variables
    ; [ebp-4]  = b (loop counter)
    ; [ebp-8]  = i (inner loop counter)
    ; [ebp-12] = j (innermost loop counter)
    ; [ebp-16] = alfaB (scalar)
    ; [ebp-20] = unoAlfaB (1-alfaB) (scalar)
    ; [ebp-24] = aligned_count (for vectorized operations)
    ; [ebp-28] = base_row_addr (address of current row)
    ; [ebp-32] = tmp
    ; [ebp-36] = tmp
    ; [ebp-40] = tmp
    ; [ebp-44] = tmp
    
    ; Get parameters
    ; [ebp+8]  = tranMat
    ; [ebp+12] = alfaB (float)
    ; [ebp+16] = maxBias
    ; [ebp+20] = d
    ; [ebp+24] = numPages

    ; Create a scalar copy of alfaB for safer access
    movss xmm0, [ebp+12]
    movss [ebp-16], xmm0

    ; Calculate unoAlfaB = 1.0 - alfaB 
    movss xmm0, [one_float]        ; Load 1.0
    movss xmm1, [ebp-16]           ; Load alfaB
    subss xmm0, xmm1               ; 1.0 - alfaB
    movss [ebp-20], xmm0           ; Store to unoAlfaB local
    
    ; Call copy_vector to initialize ret
    push dword [ebp+24]            ; numPages
    push dword [ebp+20]            ; d
    call copy_vector
    add esp, 8
    mov edi, eax                   ; EDI = ret

    ; Allocate somma
    push dword [ebp+24]            ; numPages
    call alloc_vector
    add esp, 4
    mov esi, eax                   ; ESI = somma

    ; Clear somma with memset
    mov ecx, [ebp+24]              ; numPages
    imul ecx, TYPE_SIZE            ; numPages * TYPE_SIZE
    push ecx                       ; size
    push dword 0                   ; value
    push esi                       ; somma
    call memset
    add esp, 12

    ; Calculate aligned count for vectorized operations
    mov eax, [ebp+24]              ; numPages
    mov ecx, eax                   ; Copy numPages
    shr ecx, 2                     ; Divide by 4 (number of complete SSE blocks)
    shl ecx, 2                     ; Multiply by 4 again (aligned count)
    mov [ebp-24], ecx              ; Store aligned_count

    ; Outer loop - iterate maxBias times
    mov dword [ebp-4], 0           ; b = 0
outer_loop:
    mov eax, [ebp-4]               ; Load b
    cmp eax, [ebp+16]              ; Compare with maxBias
    jge end_outer_loop             ; Exit if b >= maxBias

    ; Inner loop - for each page
    mov dword [ebp-8], 0           ; i = 0
inner_loop:
    mov eax, [ebp-8]               ; Load i
    cmp eax, [ebp+24]              ; Compare with numPages
    jge end_inner_loop             ; Exit if i >= numPages

    ; Calculate somma[i] = unoAlfaB * d[i] using SSE
    mov eax, [ebp-8]               ; i
    mov ebx, [ebp+20]              ; d
    
    ; Load values and multiply
    movss xmm0, [ebp-20]           ; unoAlfaB
    movss xmm1, [ebx + eax*4]      ; d[i]
    mulss xmm0, xmm1               ; unoAlfaB * d[i]
    
    ; Store to somma[i]
    movss [esi + eax*4], xmm0      ; somma[i] = result

    ; Calculate row address for tranMat[i]
    mov eax, [ebp-8]               ; i
    mov ebx, [ebp+24]              ; numPages
    imul eax, ebx                  ; i * numPages
    shl eax, 2                     ; i * numPages * 4 (sizeof float)
    add eax, [ebp+8]               ; tranMat + (i * numPages * 4)
    mov [ebp-28], eax              ; Store base_row_addr

    ; Zero accumulator for dot product
    xorps xmm7, xmm7               ; Clear accumulator register

    ; Vectorized loop - process 4 elements at a time
    mov dword [ebp-12], 0          ; j = 0
vectorized_loop:
    mov ecx, [ebp-12]              ; Load j
    cmp ecx, [ebp-24]              ; Compare with aligned_count
    jge remainder_loop             ; Process remaining elements

    ; Process 4 elements at once using SSE
    mov eax, [ebp-28]              ; Load row address
    movups xmm0, [eax + ecx*4]     ; Load 4 elements from tranMat[i,j...j+3]
    movups xmm1, [edi + ecx*4]     ; Load 4 elements from ret[j...j+3]
    
    ; Perform vectorized multiplication
    mulps xmm0, xmm1               ; tranMat[i,j...j+3] * ret[j...j+3]
    
    ; Manual horizontal add using SSE/SSE2 instructions
    ; Copy xmm0 to xmm2
    movaps xmm2, xmm0              ; xmm2 = [a, b, c, d]
    
    ; Shuffle xmm2 to get [b, a, d, c]
    shufps xmm2, xmm0, 0xB1        ; 0xB1 = 10110001 binary (swap within pairs)
    
    ; Add to get [a+b, b+a, c+d, d+c] which is effectively [a+b, a+b, c+d, c+d]
    addps xmm0, xmm2               ; xmm0 = [a+b, a+b, c+d, c+d]
    
    ; Copy xmm0 to xmm2 again
    movaps xmm2, xmm0              ; xmm2 = [a+b, a+b, c+d, c+d]
    
    ; Shuffle to get [c+d, c+d, a+b, a+b]
    shufps xmm2, xmm0, 0x4E        ; 0x4E = 01001110 binary (swap upper/lower halves)
    
    ; Add to get [a+b+c+d, a+b+c+d, c+d+a+b, c+d+a+b]
    addps xmm0, xmm2               ; xmm0 now has the sum in all elements
    
    ; Add lowest element to accumulator
    addss xmm7, xmm0               ; Add the sum to the accumulator
    
    ; Advance to next 4 elements
    add dword [ebp-12], 4
    jmp vectorized_loop

remainder_loop:
    ; Process remaining elements (not aligned to 4)
    mov ecx, [ebp-12]              ; Load j
    cmp ecx, [ebp+24]              ; Compare with numPages
    jge end_innermost_loop         ; Exit if j >= numPages

    ; Calculate tranMat[i,j] * ret[j]
    mov eax, [ebp-28]              ; Load row address
    
    ; Perform scalar multiplication and accumulation
    movss xmm0, [eax + ecx*4]      ; tranMat[i,j]
    movss xmm1, [edi + ecx*4]      ; ret[j]
    mulss xmm0, xmm1               ; tranMat[i,j] * ret[j]
    addss xmm7, xmm0               ; Add to accumulated result

    ; Next j
    inc dword [ebp-12]
    jmp remainder_loop

end_innermost_loop:
    ; Multiply accumulated value by alfaB
    movss xmm0, [ebp-16]           ; alfaB
    mulss xmm7, xmm0               ; alfaB * accumulated
    
    ; Add somma[i]
    mov eax, [ebp-8]               ; i
    movss xmm0, [esi + eax*4]      ; somma[i]
    addss xmm7, xmm0               ; Add to result
    
    ; Store to ret[i]
    movss [edi + eax*4], xmm7      ; ret[i] = result

    ; Next i
    inc dword [ebp-8]
    jmp inner_loop

end_inner_loop:
    ; Next b
    inc dword [ebp-4]
    jmp outer_loop

end_outer_loop:
    ; Return ret
    mov eax, edi

    ; Epilogue
    pop edi
    pop esi
    pop ebx
    mov esp, ebp
    pop ebp
    ret