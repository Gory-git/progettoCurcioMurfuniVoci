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
    ; [ebp-28] = row_offset (row offset calculation)
    ; [ebp-32] = tmp index
    ; [ebp-36] = tranMat base
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
    
    ; Save tranMat base
    mov eax, [ebp+8]
    mov [ebp-36], eax
    
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

    ; Calculate somma[i] = unoAlfaB * d[i]
    mov eax, [ebp-8]               ; i
    mov ebx, [ebp+20]              ; d
    
    ; Load values and multiply
    movss xmm0, [ebp-20]           ; unoAlfaB
    movss xmm1, [ebx + eax*4]      ; d[i]
    mulss xmm0, xmm1               ; unoAlfaB * d[i]
    
    ; Store to somma[i]
    movss [esi + eax*4], xmm0      ; somma[i] = result

    ; Zero accumulator for dot product
    xorps xmm7, xmm7               ; Clear accumulator register for row calculation

    ; Innermost loop - calculate row*vector product
    mov dword [ebp-12], 0          ; j = 0
innermost_loop:
    mov ecx, [ebp-12]              ; Load j
    cmp ecx, [ebp+24]              ; Compare with numPages
    jge end_innermost_loop         ; Exit if j >= numPages

    ; Calculate the index for tranMat[i,j]
    mov eax, [ebp-8]               ; i
    mov ebx, [ebp+24]              ; numPages
    imul eax, ebx                  ; i * numPages
    add eax, ecx                   ; i * numPages + j
    mov [ebp-32], eax              ; Store index

    ; Calculate tranMat[i,j] * ret[j] * alfaB
    mov edx, [ebp-36]              ; tranMat base
    mov eax, [ebp-32]              ; Load index
    movss xmm0, [edx + eax*4]      ; tranMat[i*numPages + j]
    movss xmm1, [edi + ecx*4]      ; ret[j]
    mulss xmm0, xmm1               ; tranMat[i,j] * ret[j]
    
    ; Multiply by alfaB
    movss xmm1, [ebp-16]           ; alfaB
    mulss xmm0, xmm1               ; alfaB * tranMat[i,j] * ret[j]
    
    ; Add to accumulated result
    addss xmm7, xmm0               ; accumulate

    ; Next j
    inc dword [ebp-12]
    jmp innermost_loop

end_innermost_loop:
    ; Add somma[i] to accumulated value and store to ret[i]
    mov eax, [ebp-8]               ; i
    movss xmm0, [esi + eax*4]      ; somma[i]
    addss xmm0, xmm7               ; somma[i] + row calculation
    movss [edi + eax*4], xmm0      ; ret[i] = result

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