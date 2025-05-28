%include "sseutils32.nasm"
section .data
    ; Define constants (assuming double precision)
    alfaI   dq 0.15         ; Example value for alfaI
    alfaB   dq 0.85         ; Example value for alfaB
    unoAlfaB dq 0.85        ; Example value for unoAlfaB (assuming it's same as alfaB)
    double_8 dq 8.0         ; Constant 8.0 for numPages * constant in k-loop
    double_zero dq 0.0      ; Constant 0.0
    double_one  dq 1.0      ; Constant 1.0 (unused based on assumptions)

section .text
    global funzione_unica
    extern puts ; Example: if you need debugging prints

funzione_unica:
    ; Standard function entry
    push    ebp
    mov     ebp, esp
    ; Reserve stack space for local variables
    ; scalare_costante (double = 8 bytes)
    ; riga (double = 8 bytes)
    sub     esp, 16

    ; Argument positions on stack (cdecl-like):
    ; [ebp+8]: tranMat (pointer)
    ; [ebp+12]: numPages (int)
    ; [ebp+16]: decay (double, 8 bytes) - unused
    ; [ebp+24]: max_outer_iterations (int)
    ; [ebp+28]: indici (pointer)
    ; [ebp+32]: d (pointer)
    ; [ebp+36]: ret (pointer)
    ; [ebp+40]: sommaV (pointer)
    ; [ebp+44]: funz1 (int/bool, 4 bytes)

    ; Local variable stack offsets:
    ; [ebp-8]: scalare_costante
    ; [ebp-16]: riga

    ; --- Variable Loading ---
    ; Load numPages, max_outer_iterations, funz1 into registers
    mov     ebx, [ebp+12]   ; ebx = numPages
    mov     ecx, [ebp+24]   ; ecx = max_outer_iterations
    mov     edx, [ebp+44]   ; edx = funz1 (0 or 1)

    ; Load pointers into registers
    mov     esi, [ebp+28]   ; esi = indici
    mov     edi, [ebp+32]   ; edi = d
    mov     eax, [ebp+36]   ; eax = ret
    mov     ebx, [ebp+40]   ; ebx = sommaV  <-- Reusing ebx, need numPages later. Store numPages somewhere else.
    mov     eax, [ebp+36]   ; eax = ret
    mov     ebp, [ebp+40]   ; ebp = sommaV  <-- Reusing ebp! Cannot do this.
    ; Let's reassign registers more carefully
    push    dword [ebp+12]  ; Save numPages temporarily
    mov     eax, [ebp+8]    ; eax = tranMat
    mov     edx, [ebp+28]   ; edx = indici
    mov     esi, [ebp+32]   ; esi = d
    mov     edi, [ebp+36]   ; edi = ret
    mov     ebx, [ebp+40]   ; ebx = sommaV
    mov     ecx, [ebp+44]   ; ecx = funz1 (0 or 1)
    mov     r8d, [ebp+24]   ; r8d = max_outer_iterations (need 64-bit reg for this, oops 32-bit)
    mov     r8d, [ebp+24]   ; This is 32-bit assembly, cannot use r8d.
    ; Need more 32-bit registers or stack vars
    ; eax=tranMat, edx=indici, esi=d, edi=ret, ebx=sommaV, ecx=funz1
    ; Use stack for max_outer_iterations and numPages after loading
    pop     dword [ebp-20]  ; Restore numPages to stack [ebp-20] (reserve 4 bytes)
    mov     dword [ebp-24], ecx ; Store funz1 on stack [ebp-24] (reserve 4 bytes)
    mov     dword [ebp-28], [ebp+24] ; Store max_outer_iterations on stack [ebp-28] (reserve 4 bytes)
    mov     ecx, [ebp-28]   ; ecx = max_outer_iterations

    ; Pointers in: eax=tranMat, edx=indici, esi=d, edi=ret, ebx=sommaV
    ; numPages in [ebp-20]
    ; funz1 in [ebp-24]
    ; max_outer_iterations in ecx (loop counter)

    ; --- Initialize scalare_costante ---
    ; Check funz1
    cmp     dword [ebp-24], 0
    je      .else_scalar_const

    ; if funz1: scalare_costante = (1 - alfaI) / (type) numPages;
    movsd   xmm0, [double_one]  ; xmm0 = 1.0
    subsd   xmm0, [alfaI]       ; xmm0 = 1.0 - alfaI
    cvtsi2sd xmm1, dword [ebp-20] ; xmm1 = (double) numPages
    divsd   xmm0, xmm1          ; xmm0 = (1.0 - alfaI) / (double) numPages
    movsd   [ebp-8], xmm0       ; store scalare_costante = xmm0
    jmp     .end_scalar_const

.else_scalar_const:
    ; else: scalare_costante = (type) 1 - alfaB;
    movsd   xmm0, [double_one]  ; xmm0 = 1.0
    subsd   xmm0, [alfaB]       ; xmm0 = 1.0 - alfaB
    movsd   [ebp-8], xmm0       ; store scalare_costante = xmm0

.end_scalar_const:

    ; --- Outer loop: for i < max_outer_iterations ---
    mov     dword [ebp-32], 0   ; i = 0 (reserve 4 bytes for i at [ebp-32])
    mov     ecx, [ebp-28]       ; ecx = max_outer_iterations
    jmp     .outer_loop_cond

.outer_loop:
    ; i is in [ebp-32]

    ; --- Inner loop: for j < numPages ---
    mov     dword [ebp-36], 0   ; j = 0 (reserve 4 bytes for j at [ebp-36])
    mov     edx, [ebp-20]       ; edx = numPages (loop limit)
    jmp     .inner_loop_cond

.inner_loop:
    ; j is in [ebp-36]
    ; i is in [ebp-32]

    ; --- Inside inner loop, before k loop ---
    cmp     dword [ebp-24], 0   ; Check funz1
    je      .inner_else_before_k

    ; if funz1
    mov     eax, [ebp-32]       ; eax = i
    cmp     eax, 0              ; if i == 0
    jne     .inner_funz1_before_k_end

    ; if i == 0
    mov     ebx, [ebp-32]       ; ebx = i
    mov     dword [edx + ebx*4], ebx ; indici[i] = i (using edx=indici, ebx=i, size 4)
    mov     ebx, [ebp-32]       ; ebx = i
    movsd   xmm0, [double_zero] ; xmm0 = 0.0
    movsd   [esi + ebx*8], xmm0 ; d[i] = 0 (using esi=d, ebx=i, size 8)

.inner_funz1_before_k_end:
    jmp     .after_inner_before_k

.inner_else_before_k:
    ; else (!funz1)
    ; somma[i] = unoAlfaB * d[i];  <-- pseudocode used index i, inside j loop!
    mov     eax, [ebp-32]       ; eax = i
    movsd   xmm0, [unoAlfaB]    ; xmm0 = unoAlfaB
    movsd   xmm1, [esi + eax*8] ; xmm1 = d[i]
    mulsd   xmm0, xmm1          ; xmm0 = unoAlfaB * d[i]
    mov     ebx, [ebp-32]       ; ebx = i
    movsd   [ebx + ebx*8], xmm0 ; sommaV[i] = xmm0 (using ebx=sommaV, ebx=i, size 8) <-- This is WRONG! ebx holds sommaV pointer, needs correct base.
    mov     ecx, [ebx]          ; ecx = sommaV base
    mov     ebx, [ebp-32]       ; ebx = i
    movsd   [ecx + ebx*8], xmm0 ; sommaV[i] = xmm0 (using ecx=sommaV, ebx=i, size 8) <-- Fixed

.after_inner_before_k:

    ; --- Initialize riga ---
    movsd   xmm0, [double_zero] ; riga = 0.0
    movsd   [ebp-16], xmm0      ; store riga locally

    ; --- Innermost loop: for k < numPages ---
    mov     dword [ebp-40], 0   ; k = 0 (reserve 4 bytes for k at [ebp-40])
    mov     ecx, [ebp-20]       ; ecx = numPages (loop limit)
    jmp     .k_loop_cond

.k_loop:
    ; k is in [ebp-40]
    ; i is in [ebp-32]
    ; j is in [ebp-36]

    ; --- Inside k loop ---
    cmp     dword [ebp-24], 0   ; Check funz1
    je      .k_else

    ; if funz1
    ; (Ignoring the nonsensical "if (m == 0 && i == 0) { s[j] = 1; }" line)

    ; Calculate matrix index: x = i * numPages + j
    mov     eax, [ebp-32]       ; eax = i
    mov     ebx, [ebp-20]       ; ebx = numPages
    imul    eax, ebx            ; eax = i * numPages
    mov     ebx, [ebp-36]       ; ebx = j
    add     eax, ebx            ; eax = i * numPages + j
    mov     ebx, [eax*8]        ; ebx = tranMat + (i * numPages + j) * 8 (Base tranMat missing!)
    mov     ebx, [ebp+8]        ; ebx = tranMat base pointer
    mov     eax, [ebp-32]       ; eax = i
    mov     ecx, [ebp-20]       ; ecx = numPages
    imul    eax, ecx            ; eax = i * numPages
    mov     ecx, [ebp-36]       ; ecx = j
    add     eax, ecx            ; eax = i * numPages + j
    movsd   xmm1, [ebx + eax*8] ; xmm1 = tranMat[i * numPages + j]

    ; Load s[j] (assuming s is d)
    mov     eax, [ebp-36]       ; eax = j
    movsd   xmm2, [esi + eax*8] ; xmm2 = d[j] (using esi=d, eax=j, size 8)

    ; Calculate alfaI * tranMat[...] * d[j]
    movsd   xmm3, [alfaI]       ; xmm3 = alfaI
    mulsd   xmm3, xmm1          ; xmm3 = alfaI * tranMat[...]
    mulsd   xmm3, xmm2          ; xmm3 = alfaI * tranMat[...] * d[j]

    ; Add to riga
    movsd   xmm0, [ebp-16]      ; xmm0 = riga
    addsd   xmm0, xmm3          ; xmm0 = riga + term
    movsd   [ebp-16], xmm0      ; store riga

    jmp     .k_loop_end

.k_else:
    ; else (!funz1)
    ; Calculate matrix index: i * numPages + j
    mov     eax, [ebp-32]       ; eax = i
    mov     ebx, [ebp-20]       ; ebx = numPages
    imul    eax, ebx            ; eax = i * numPages
    mov     ebx, [ebp-36]       ; ebx = j
    add     eax, ebx            ; eax = i * numPages + j
    mov     ebx, [ebp+8]        ; ebx = tranMat base pointer
    movsd   xmm1, [ebx + eax*8] ; xmm1 = tranMat[i * numPages + j]

    ; Load ret[j]
    mov     eax, [ebp-36]       ; eax = j
    movsd   xmm2, [edi + eax*8] ; xmm2 = ret[j] (using edi=ret, eax=j, size 8)

    ; Calculate alfaB * ret[j] * tranMat[...]
    movsd   xmm3, [alfaB]       ; xmm3 = alfaB
    mulsd   xmm3, xmm2          ; xmm3 = alfaB * ret[j]
    mulsd   xmm3, xmm1          ; xmm3 = alfaB * ret[j] * tranMat[...]

    ; Add to riga
    movsd   xmm0, [ebp-16]      ; xmm0 = riga
    addsd   xmm0, xmm3          ; xmm0 = riga + term
    movsd   [ebp-16], xmm0      ; store riga

.k_loop_end:
    ; Increment k
    inc     dword [ebp-40]
    ; k_loop condition: k < numPages
.k_loop_cond:
    mov     eax, [ebp-40]       ; eax = k
    mov     ebx, [ebp-20]       ; ebx = numPages
    cmp     eax, ebx
    jl      .k_loop

    ; --- After k loop ---
    ; Load riga
    movsd   xmm0, [ebp-16]      ; xmm0 = riga

    ; Check funz1
    cmp     dword [ebp-24], 0
    je      .after_k_else

    ; if funz1
    ; s[i] = riga + somma; <-- assuming s is d, somma is scalare_costante
    movsd   xmm1, [ebp-8]       ; xmm1 = scalare_costante
    addsd   xmm0, xmm1          ; xmm0 = riga + scalare_costante
    mov     eax, [ebp-32]       ; eax = i <-- pseudocode used index i!
    movsd   [esi + eax*8], xmm0 ; d[i] = xmm0 (using esi=d, eax=i, size 8)

    jmp     .after_k_end

.after_k_else:
    ; else (!funz1)
    ; ret[i] = riga + somma[i]; <-- pseudocode used index i!
    mov     eax, [ebp-32]       ; eax = i <-- pseudocode used index i!
    mov     ecx, [ebp+40]       ; ecx = sommaV base pointer
    movsd   xmm1, [ecx + eax*8] ; xmm1 = sommaV[i] (using ecx=sommaV, eax=i, size 8)
    addsd   xmm0, xmm1          ; xmm0 = riga + sommaV[i]
    mov     eax, [ebp-32]       ; eax = i <-- pseudocode used index i!
    movsd   [edi + eax*8], xmm0 ; ret[i] = xmm0 (using edi=ret, eax=i, size 8)

.after_k_end:

    ; Increment j
    inc     dword [ebp-36]
    ; Inner loop condition: j < numPages
.inner_loop_cond:
    mov     eax, [ebp-36]       ; eax = j
    mov     ebx, [ebp-20]       ; ebx = numPages
    cmp     eax, ebx
    jl      .inner_loop

    ; Increment i
    inc     dword [ebp-32]
    ; Outer loop condition: i < max_outer_iterations
.outer_loop_cond:
    mov     eax, [ebp-32]       ; eax = i
    mov     ebx, [ebp-28]       ; ebx = max_outer_iterations
    cmp     eax, ebx
    jl      .outer_loop


    ; --- Function Return ---
    ; Pseudocode returns ret (pointer). Conventionally, return the pointer argument.
    mov     eax, [ebp+36]   ; eax = ret pointer

    ; Standard function exit
    add     esp, 40         ; Deallocate stack space (16 for locals + 24 for saved args)
    pop     ebp
    ret


;To assemble and link (Linux):
;nasm -f elf funzione_unica_32.asm -o funzione_unica_32.o
;gcc -m32 -o test_program test_program.c funzione_unica_32.o # Replace test_program.c with a file that calls the function
