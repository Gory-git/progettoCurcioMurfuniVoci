section .note.GNU-stack noalloc noexec nowrite progbits

BITS 64
default rel

extern alloc_vector
extern dealloc_vector
extern copy_vector

section .data
align 32
const_one:      dq 1.0

section .bss
align 32
alpha:          resq 1
oneMinusAlpha:  resq 1
tranMatPtr:     resq 1
dPtr:           resq 1
maxBiasVar:     resd 1

section .text
global computeScores

; SysV AMD64:
;   rdi = tranMat (double*)
;   xmm0= alpha   (double)
;   rsi = maxBias (int)
;   rdx = d       (double*)
;   ecx = numPages(int)
; return rax = ret vector (double*)

%define REG_numPages    r12
%define REG_retA        r13
%define REG_retB        r14
%define REG_somma       r15
%define REG_biasLoop    r8d
%define REG_pageLoop    r9
%define REG_innerLoop   r10
%define REG_rowBase     rbx

%define YMM_alpha           ymm15
%define YMM_oneMinusAlpha   ymm14
%define YMM_acc0            ymm13
%define YMM_acc1            ymm12
%define YMM_tmp1            ymm11
%define YMM_tmp2            ymm10

computeScores:
    push    rbp
    mov     rbp, rsp
    push    rbx
    push    r12
    push    r13
    push    r14
    push    r15

    ; Validazione parametri
    test    rdi, rdi
    jz      .error_exit
    test    rdx, rdx
    jz      .error_exit
    test    ecx, ecx
    jle     .error_exit

    mov     r12d, ecx
    mov     [tranMatPtr], rdi
    mov     [dPtr], rdx
    mov     [maxBiasVar], esi

    ; alpha e (1-alpha)
    movq    [alpha], xmm0
    movsd   xmm1, [const_one]
    subsd   xmm1, [alpha]
    movq    [oneMinusAlpha], xmm1

    vbroadcastsd YMM_alpha,         [alpha]
    vbroadcastsd YMM_oneMinusAlpha, [oneMinusAlpha]

    ; retA = copy_vector(d, n)
    vzeroupper
    mov     rdi, [dPtr]
    mov     esi, r12d
    call    copy_vector
    mov     REG_retA, rax

    ; retB = alloc_vector(n)
    mov     edi, r12d
    call    alloc_vector
    mov     REG_retB, rax

    ; somma = alloc_vector(n)
    mov     edi, r12d
    call    alloc_vector
    mov     REG_somma, rax

    ; bias loop
    xor     REG_biasLoop, REG_biasLoop

.bias_loop:
    cmp     REG_biasLoop, dword [maxBiasVar]
    jge     .done

    ; ------------------------------------------------------------
    ; somma = (1 - alpha) * d        (una volta per bias)
    ; ------------------------------------------------------------
    mov     rsi, [dPtr]
    xor     REG_pageLoop, REG_pageLoop
    mov     rax, REG_numPages
    and     rax, -4

.somma_vec:
    cmp     REG_pageLoop, rax
    jge     .somma_tail

    vmovupd YMM_tmp1, [rsi + REG_pageLoop*8]
    vmulpd  YMM_tmp1, YMM_tmp1, YMM_oneMinusAlpha
    ; somma/ret allocati devono essere allineati a 32B:
    vmovapd [REG_somma + REG_pageLoop*8], YMM_tmp1
    add     REG_pageLoop, 4
    jmp     .somma_vec

.somma_tail:
    cmp     REG_pageLoop, REG_numPages
    jge     .page_loop_start
    mov     rdx, REG_pageLoop
    shl     rdx, 3
    movsd   xmm0, [rsi + rdx]
    mulsd   xmm0, [oneMinusAlpha]
    movsd   [REG_somma + rdx], xmm0
    inc     REG_pageLoop
    jmp     .somma_tail

; ------------------------------------------------------------
; Loop sulle righe
; ------------------------------------------------------------
.page_loop_start:
    mov     rdx, [tranMatPtr]
    xor     REG_pageLoop, REG_pageLoop

.page_loop:
    cmp     REG_pageLoop, REG_numPages
    jge     .end_iter

    ; rowBase = tranMat + i*n*8
    mov     rax, REG_pageLoop
    imul    rax, REG_numPages
    shl     rax, 3
    lea     REG_rowBase, [rdx + rax]

    prefetcht0 [REG_rowBase]

    vxorpd  YMM_acc0, YMM_acc0, YMM_acc0
    vxorpd  YMM_acc1, YMM_acc1, YMM_acc1
    xorpd   xmm7, xmm7

    xor     REG_innerLoop, REG_innerLoop
    mov     rcx, REG_numPages
    and     rcx, -8              ; blocchi da 8 (2×4)

.dot_vec_unroll:
    cmp     REG_innerLoop, rcx
    jge     .dot_tail

    ; primo blocco (4 doppi)
    vmovupd YMM_tmp1, [REG_rowBase + REG_innerLoop*8]  ; tranMat: unaligned
    vmovapd YMM_tmp2, [REG_retA    + REG_innerLoop*8]  ; retA:    aligned
    vfmadd231pd YMM_acc0, YMM_tmp1, YMM_tmp2

    add     REG_innerLoop, 4

    ; secondo blocco (4 doppi)
    vmovupd YMM_tmp1, [REG_rowBase + REG_innerLoop*8]
    vmovapd YMM_tmp2, [REG_retA    + REG_innerLoop*8]
    vfmadd231pd YMM_acc1, YMM_tmp1, YMM_tmp2

    add     REG_innerLoop, 4
    jmp     .dot_vec_unroll

.dot_tail:
    cmp     REG_innerLoop, REG_numPages
    jge     .dot_reduce
    mov     rax, REG_innerLoop
    shl     rax, 3
    movsd   xmm2, [REG_rowBase + rax]
    movsd   xmm3, [REG_retA    + rax]
    mulsd   xmm2, xmm3
    addsd   xmm7, xmm2
    inc     REG_innerLoop
    jmp     .dot_tail

.dot_reduce:
    vaddpd  YMM_acc0, YMM_acc0, YMM_acc1
    vextractf128 xmm0, YMM_acc0, 0
    vextractf128 xmm1, YMM_acc0, 1
    vaddpd  xmm0, xmm0, xmm1
    vhaddpd xmm0, xmm0, xmm0
    addsd   xmm0, xmm7

    ; retB[i] = somma[i] + alpha * dot
    mulsd   xmm0, [alpha]
    mov     rax, REG_pageLoop
    shl     rax, 3
    addsd   xmm0, [REG_somma + rax]
    movsd   [REG_retB + rax], xmm0

    inc     REG_pageLoop
    jmp     .page_loop

.end_iter:
    ; swap retA <-> retB
    mov     rax, REG_retA
    mov     REG_retA, REG_retB
    mov     REG_retB, rax

    inc     REG_biasLoop
    jmp     .bias_loop

.done:
    vzeroupper
    mov     rdi, REG_retB
    call    dealloc_vector
    mov     rdi, REG_somma
    call    dealloc_vector
    mov     rax, REG_retA
    jmp     .exit

.error_exit:
    xor     rax, rax

.exit:
    pop     r15
    pop     r14
    pop     r13
    pop     r12
    pop     rbx
    mov     rsp, rbp
    pop     rbp
    ret
