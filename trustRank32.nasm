section .note.GNU-stack noalloc noexec nowrite progbits

[BITS 32]

%define TYPE_SIZE 4

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

section .text
global computeScores
extern copy_vector
extern alloc_vector
extern dealloc_vector

; --------------------------------------------------------------------
; intf:
;   [ebp+8]   = tranMat (float*)
;   [ebp+12]  = alfaB   (float)
;   [ebp+16]  = maxBias (int)
;   [ebp+20]  = d       (float*)
;   [ebp+24]  = numPages(int)
; return eax = ret vector (float*)
; --------------------------------------------------------------------

computeScores:
    push    ebp
    mov     ebp, esp
    push    ebx
    push    esi
    push    edi

    ; --- Validazione parametri ---
    mov     eax, [ebp+8]       ; tranMat
    test    eax, eax
    jz      .error_exit
    mov     eax, [ebp+20]      ; d
    test    eax, eax
    jz      .error_exit
    mov     eax, [ebp+24]      ; numPages
    test    eax, eax
    jle     .error_exit

    ; --- Salvataggi base ---
    mov     eax, [ebp+12]      ; alpha
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

    ; --- ret = copy_vector(d, n) ---
    push    dword [numPagesVar]
    push    dword [dPtr]
    call    copy_vector
    add     esp, 8
    test    eax, eax
    jz      .error_exit
    mov     [retPtr], eax

    ; --- somma = alloc_vector(n) ---
    push    dword [numPagesVar]
    call    alloc_vector
    add     esp, 4
    test    eax, eax
    jz      .error_cleanup_ret
    mov     [sommaPtr], eax

    ; bias loop
    mov     dword [biasLoopVar], 0

.bias_loop:
    mov     eax, [biasLoopVar]
    cmp     eax, [maxBiasVar]
    jge     .done

    ; ------------------------------------------------------------
    ; somma = (1 - alpha) * d      (calcolata una volta per bias)
    ; ------------------------------------------------------------
    mov     esi, [dPtr]          ; src d
    mov     edi, [sommaPtr]      ; dst somma
    mov     ecx, [numPagesVar]   ; count
    movss   xmm5, [oneMinusAlpha]
    shufps  xmm5, xmm5, 0x00     ; broadcast

    ; vettoriale 4 a 4
.somma_vec4:
    cmp     ecx, 4
    jb      .somma_tail

    movups  xmm0, [esi]          ; d (unaligned: d può non essere allineato)
    mulps   xmm0, xmm5
    ; Se somma è sicuramente 16B aligned, usare movaps:
    movaps  [edi], xmm0
    add     esi, 16
    add     edi, 16
    sub     ecx, 4
    jmp     .somma_vec4

.somma_tail:
    test    ecx, ecx
    jz      .page_loop_start
.somma_tail_loop:
    movss   xmm0, [esi]
    mulss   xmm0, xmm5
    movss   [edi], xmm0
    add     esi, 4
    add     edi, 4
    dec     ecx
    jnz     .somma_tail_loop

    ; ------------------------------------------------------------
    ; Loop sulle righe (pagine)
    ; ret[i] = somma[i] + alpha * dot(tranMat[i,*], ret)
    ; ------------------------------------------------------------
.page_loop_start:
    mov     dword [pageLoopVar], 0

.page_loop:
    mov     eax, [pageLoopVar]
    cmp     eax, [numPagesVar]
    jge     .end_page_loop

    ; rowPtr = tranBase + i*n*4
    mov     ebx, [numPagesVar]
    mov     edx, [tranBase]
    mov     ecx, eax             ; ecx = i
    imul    ecx, ebx             ; ecx = i*n
    shl     ecx, 2               ; *4 bytes
    lea     edx, [edx + ecx]     ; edx = rowPtr

    ; ptr cursori
    mov     esi, [retPtr]        ; ret base (allineato)
    mov     edi, edx             ; row ptr

    ; accumulatori
    xorps   xmm7, xmm7
    xorps   xmm6, xmm6

    ; prefetch riga
    prefetcht0 [edi]

    ; --- vettoriale 8 a 8 (unrolling ×2) ---
    mov     ecx, [numPagesVar]
    mov     ebx, ecx
    and     ebx, 0xFFFFFFF8      ; blocchi da 8
    test    ebx, ebx
    jz      .vec4_check
    xor     eax, eax             ; processed = 0

.vec8_loop:
    cmp     eax, ebx
    jge     .vec4_check

    movups  xmm0, [edi]          ; 1° blocco 4
    movaps  xmm1, [esi]
    mulps   xmm0, xmm1
    addps   xmm7, xmm0

    movups  xmm2, [edi+16]       ; 2° blocco 4
    movaps  xmm3, [esi+16]
    mulps   xmm2, xmm3
    addps   xmm6, xmm2

    add     edi, 32
    add     esi, 32
    add     eax, 8
    jmp     .vec8_loop

.vec4_check:
    mov     ecx, [numPagesVar]
    sub     ecx, eax             ; restanti dopo blocchi da 8
    cmp     ecx, 4
    jb      .tail_scalar

    ; un blocco da 4
    movups  xmm0, [edi]
    movaps  xmm1, [esi]
    mulps   xmm0, xmm1
    addps   xmm7, xmm0
    add     edi, 16
    add     esi, 16
    sub     ecx, 4
    add     eax, 4

.tail_scalar:
    ; riduzione parziale vettoriale
    addps   xmm7, xmm6
    movaps  xmm0, xmm7
    movhlps xmm1, xmm0           ; (a2,a3,*,*)
    addps   xmm0, xmm1           ; (a0+a2, a1+a3, ...)
    movaps  xmm1, xmm0
    shufps  xmm1, xmm1, 0x55     ; dup lane1
    addss   xmm0, xmm1           ; xmm0.low = somma 4-lanesc

    ; tail scalare rimanente
    mov     ecx, [numPagesVar]
    sub     ecx, eax             ; quanti ancora
    jz      .dot_ready

.tail_loop:
    movss   xmm2, [edi]
    movss   xmm3, [esi]
    mulss   xmm2, xmm3
    addss   xmm0, xmm2
    add     edi, 4
    add     esi, 4
    dec     ecx
    jnz     .tail_loop

.dot_ready:
    ; ret[i] = somma[i] + alpha * dot
    movss   xmm1, [alpha]
    mulss   xmm0, xmm1
    mov     edx, [sommaPtr]
    mov     eax, [pageLoopVar]
    shl     eax, 2
    addss   xmm0, [edx + eax]
    mov     edx, [retPtr]
    movss   [edx + eax], xmm0

    inc     dword [pageLoopVar]
    jmp     .page_loop

.end_page_loop:
    inc     dword [biasLoopVar]
    jmp     .bias_loop

.done:
    ; free somma e ritorna ret
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
