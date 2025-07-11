section .note.GNU-stack noalloc noexec nowrite progbits

section .text
[BITS 32]
global computeScores
extern copy_vector
extern alloc_vector
extern memset

; Simple implementation of computeScores function
; VECTOR computeScores(MATRIX tranMat, float alfaB, int maxBias, VECTOR d, int numPages)
%define TYPE_SIZE 4

computeScores:
    ; Prologue
    push ebp
    mov ebp, esp
    sub esp, 24                    ; Allocate space for local variables
    push ebx                       ; Save callee-saved registers
    push esi
    push edi

    ; Setup local variables
    ; [ebp-4]  = b (loop counter)
    ; [ebp-8]  = i
    ; [ebp-12] = j
    ; [ebp-16] = accumulated value (float)
    ; [ebp-20] = unoAlfaB (1-alfaB)
    
    ; Get parameters
    ; [ebp+8]  = tranMat
    ; [ebp+12] = alfaB (float)
    ; [ebp+16] = maxBias
    ; [ebp+20] = d
    ; [ebp+24] = numPages

    ; Call copy_vector to initialize ret
    push dword [ebp+24]        ; numPages
    push dword [ebp+20]        ; d
    call copy_vector
    add esp, 8
    mov edi, eax               ; EDI = ret

    ; Allocate somma
    push dword [ebp+24]        ; numPages
    call alloc_vector
    add esp, 4
    mov esi, eax               ; ESI = somma

    ; Clear somma with memset
    mov ecx, [ebp+24]          ; numPages
    imul ecx, TYPE_SIZE        ; numPages * TYPE_SIZE
    push ecx                   ; size
    push dword 0               ; value
    push esi                   ; somma
    call memset
    add esp, 12

    ; Calculate unoAlfaB = 1.0 - alfaB
    fld1                       ; Load 1.0 onto FPU stack
    fld dword [ebp+12]         ; Load alfaB
    fsubp st1, st0             ; 1.0 - alfaB
    fstp dword [ebp-20]        ; Store to unoAlfaB local

    ; Outer loop - iterate maxBias times
    mov dword [ebp-4], 0       ; b = 0
outer_loop:
    mov eax, [ebp-4]           ; Load b
    cmp eax, [ebp+16]          ; Compare with maxBias
    jge end_outer_loop         ; Exit if b >= maxBias

    ; Inner loop - for each page
    mov dword [ebp-8], 0       ; i = 0
inner_loop:
    mov eax, [ebp-8]           ; Load i
    cmp eax, [ebp+24]          ; Compare with numPages
    jge end_inner_loop         ; Exit if i >= numPages

    ; Calculate somma[i] = unoAlfaB * d[i]
    mov eax, [ebp-8]           ; i
    mov ebx, [ebp+20]          ; d
    fld dword [ebp-20]         ; unoAlfaB
    fld dword [ebx + eax*4]    ; d[i]
    fmulp st1, st0             ; unoAlfaB * d[i]
    fstp dword [esi + eax*4]   ; somma[i] = result

    ; Initialize accumulated value for innermost loop
    fldz                       ; Load 0.0 for accumulation
    fstp dword [ebp-16]        ; Store to local

    ; Innermost loop - for each page (j)
    mov dword [ebp-12], 0      ; j = 0
innermost_loop:
    mov ecx, [ebp-12]          ; Load j
    cmp ecx, [ebp+24]          ; Compare with numPages
    jge end_innermost_loop     ; Exit if j >= numPages

    ; Calculate tranMat[i,j] * ret[j]
    mov eax, [ebp-8]           ; i
    mov ebx, [ebp+24]          ; numPages
    imul eax, ebx              ; i * numPages
    add eax, ecx               ; i * numPages + j
    mov ebx, [ebp+8]           ; tranMat
    fld dword [ebx + eax*4]    ; tranMat[i,j]
    
    mov ebx, edi               ; ret
    fld dword [ebx + ecx*4]    ; ret[j]
    
    fmulp st1, st0             ; tranMat[i,j] * ret[j]
    fld dword [ebp-16]         ; Load accumulated value
    faddp st1, st0             ; Add to accumulated
    fstp dword [ebp-16]        ; Store back to local

    ; Next j
    inc dword [ebp-12]
    jmp innermost_loop

end_innermost_loop:
    ; Multiply accumulated value by alfaB
    fld dword [ebp-16]         ; Load accumulated
    fld dword [ebp+12]         ; alfaB
    fmulp st1, st0             ; alfaB * accumulated
    
    ; Add somma[i]
    mov eax, [ebp-8]           ; i
    fld dword [esi + eax*4]    ; somma[i]
    faddp st1, st0             ; Add to result
    
    ; Store to ret[i]
    fstp dword [edi + eax*4]   ; ret[i] = result

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