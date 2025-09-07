; ====================================================================
; VECTOR computeScores(MATRIX tranMat, float alfaB, int maxBias, VECTOR d, int numPages)
; Versione x86-32 con SSE:
; - SIMD completa (mulps/addps)
; - Unrolling ×2 su blocchi da 8
; - Riduzione orizzontale
; - Variabili stabili in .bss (no offset numerici)
; ====================================================================

section .note.GNU-stack noalloc noexec nowrite progbits

section .data
align 16
one_float:      dd 1.0, 1.0, 1.0, 1.0
zero_float:     dd 0.0, 0.0, 0.0, 0.0

section .bss
align 16
alpha:          resd 1
oneMinusAlpha:  resd 1
tranBase:       resd 1
dPtr:           resd 1
retPtr:         resd 1
sommaPtr:       resd 1
numPagesVar:    resd 1
maxBiasVar:     resd 1
biasLoopVar:    resd 1
pageLoopVar:    resd 1
innerLoopVar:   resd 1
vecBound4:      resd 1
vecBound8:      resd 1

section .text
[BITS 32]
global computeScores
extern copy_vector
extern alloc_vector
extern dealloc_vector

computeScores:
    push    ebp
    mov     ebp, esp
    push    ebx
    push    esi
    push    edi

    ; --- Validazione parametri ---
    mov     eax, [ebp+8]      ; tranMat
    test    eax, eax
    jz      .error_exit
    mov     eax, [ebp+20]     ; d
    test    eax, eax
    jz      .error_exit
    mov     eax, [ebp+24]     ; numPages
    test    eax, eax
    jle     .error_exit

    ; --- Salvataggio parametri in BSS ---
    mov     eax, [ebp+12]     ; alfaB
    mov     [alpha], eax
    movss   xmm0, [one_float]
    subss   xmm0, [alpha]
    movss   [oneMinusAlpha], xmm0

    mov     eax, [ebp+8]
    mov     [tranBase], eax
    mov     eax, [ebp+20]
    mov     [dPtr], eax
    mov     eax, [ebp+24]
    mov     [numPagesVar], eax
    mov     eax, [ebp+16]
    mov     [maxBiasVar], eax

    ; --- Pre-calcolo bound ---
    mov     eax, [numPagesVar]
    mov     ebx, eax
    and     ebx, 0FFFFFFFCh
    mov     [vecBound4], ebx
    and     eax, 0FFFFFFF8h
    mov     [vecBound8], eax

    ; --- ret = copy_vector(d, numPages) ---
    push    dword [numPagesVar]
    push    dword [dPtr]
    call    copy_vector
    add     esp, 8
    test    eax, eax
    jz      .error_exit
    mov     [retPtr], eax

    ; --- somma = alloc_vector(numPages) ---
    push    dword [numPagesVar]
    call    alloc_vector
    add     esp, 4
    test    eax, eax
    jz      .error_cleanup_ret
    mov     [sommaPtr], eax

    ; --- Bias loop ---
    mov     dword [biasLoopVar], 0

.bias_loop:
    mov     eax, [biasLoopVar]
    cmp     eax, [maxBiasVar]
    jge     .done

    ; --- Page loop ---
    mov     dword [pageLoopVar], 0

.page_loop:
    mov     eax, [pageLoopVar]
    cmp     eax, [numPagesVar]
    jge     .end_page_loop

    ; somma[i] = (1-alfa) * d[i]
    mov     esi, [sommaPtr]
    mov     edi, [dPtr]
    mov     ecx, [pageLoopVar]
    movss   xmm0, [oneMinusAlpha]
    movss   xmm1, [edi + ecx*4]
    mulss   xmm0, xmm1
    movss   [esi + ecx*4], xmm0

    ; --- Dot product riga·ret ---
    xorps   xmm7, xmm7
    xorps   xmm6, xmm6
    mov     edi, [retPtr]
    mov     dword [innerLoopVar], 0

.dot_vec8:
    mov     ecx, [innerLoopVar]
    cmp     ecx, [vecBound8]
    jge     .dot_vec4

    ; base = tranMat + (i*numPages + j)
    mov     eax, [pageLoopVar]
    mov     ebx, [numPagesVar]
    imul    eax, ebx
    add     eax, ecx
    mov     edx, [tranBase]
    lea     edx, [edx + eax*4]

    movups  xmm0, [edx]
    movups  xmm1, [edi + ecx*4]
    mulps   xmm0, xmm1
    addps   xmm7, xmm0

    movups  xmm2, [edx+16]
    movups  xmm3, [edi + ecx*4 + 16]
    mulps   xmm2, xmm3
    addps   xmm6, xmm2

    add     ecx, 8
    mov     [innerLoopVar], ecx
    jmp     .dot_vec8

.dot_vec4:
    mov     ecx, [innerLoopVar]
    cmp     ecx, [vecBound4]
    jge     .dot_reduce

    mov     eax, [pageLoopVar]
    mov     ebx, [numPagesVar]
    imul    eax, ebx
    add     eax, ecx
    mov     edx, [tranBase]
    lea     edx, [edx + eax*4]

    movups  xmm0, [edx]
    movups  xmm1, [edi + ecx*4]
    mulps   xmm0, xmm1
    addps   xmm7, xmm0

    add     ecx, 4
    mov     [innerLoopVar], ecx

.dot_reduce:
    addps   xmm7, xmm6
    movaps  xmm0, xmm7
    shufps  xmm1, xmm7, 0b11101110
    addps   xmm0, xmm1
    movhlps xmm1, xmm0
    addss   xmm0, xmm1

; --- Tail scalare ---
.dot_tail:
    mov     ecx, [innerLoopVar]
    cmp     ecx, [numPagesVar]
    jge     .dot_done

    mov     eax, [pageLoopVar]
    mov     ebx, [numPagesVar]
    imul    eax, ebx
    add     eax, ecx
    mov     edx, [tranBase]

    movss   xmm4, [edx + eax*4]
    movss   xmm5, [edi + ecx*4]
    movaps  xmm6, xmm4
    mulss   xmm6, xmm5
    addss   xmm0, xmm6

    inc     ecx
    mov     [innerLoopVar], ecx
    jmp     .dot_tail

.dot_done:
    ; alfa * dot + somma[i]
    movss   xmm1, [alpha]
    mulss   xmm0, xmm1
    mov     esi, [sommaPtr]
    mov     edi, [retPtr]
    mov     eax, [pageLoopVar]
    movss   xmm6, [esi + eax*4]
    addss   xmm0, xmm6
    movss   [edi + eax*4], xmm0

    inc     dword [pageLoopVar]
    jmp     .page_loop

.end_page_loop:
    inc     dword [biasLoopVar]
    jmp     .bias_loop

.done:
    ; free somma
    push    dword [sommaPtr]
    call    dealloc_vector
    add     esp, 4

    mov     eax, [retPtr]
    jmp     .exit

.error_cleanup_ret:
    push    dword [retPtr]
    call    dealloc_vector
    add     esp, 4

.error_exit:
    xor     eax, eax

.exit:
    pop     edi
    pop     esi
    pop     ebx
    mov     esp, ebp
    pop     ebp
    ret
