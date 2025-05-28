section .data
align 8
alfaI:          dq 0.0     ; Definisci queste costanti (valori di esempio)
alfaB:          dq 0.0
unoAlfaB:       dq 0.0
one_double:     dq 1.0
zero_double:    dq 0.0

section .text
global funzione_unica

; VECTOR funzione_unica(MATRIX tranMatInv, int numPages, type decay, int max_outer_iterations, int* indici, VECTOR d, VECTOR ret, VECTOR somma, bool funz1, MATRIX tranMatParam)
; Parametri nei registri (System V ABI):
; rdi: tranMatInv (ptr)
; rsi: numPages (int)
; rdx: decay (type - ignorato)
; rcx: max_outer_iterations (int)
; r8:  indici (ptr)
; r9:  d (ptr)
; Parametri sullo stack (dal 7° argomento in poi):
; [rbp+16]: ret (ptr)
; [rbp+24]: somma (ptr)
; [rbp+32]: funz1 (bool, probabilmente 4 byte sullo stack, letta come int)
; [rbp+40]: tranMatParam (ptr)

funzione_unica:
    push    rbp            ; Salva base pointer
    mov     rbp, rsp       ; Imposta stack frame

    ; Salva registri callee-saved
    push    rbx
    push    r12
    push    r13
    push    r14
    push    r15
    ; rsp è ora rbp_entry - 8 (per rbp) - 5*8 (per rbx, r12-r15) = rbp_entry - 48
    ; rsp è (qualcosa) + 8 (da call) - 48 = (qualcosa) - 40.
    ; Per allineare rsp a 16 byte prima di allocare spazio locale (se necessario per chiamate future):
    ; Se (qualcosa) era 16N, rsp è 16N-40. Non è allineato.
    ; Se questa funzione chiamasse altre funzioni passando argomenti SSE sullo stack, l'allineamento di rsp sarebbe cruciale.
    ; Non ci sono chiamate in questo codice, quindi l'allineamento stretto di rsp *prima* di sub rsp
    ; non è critico per il funzionamento *interno*, ma è buona norma per la robustezza.
    ; Lo spazio locale di 80 byte renderà rsp = (qualcosa) - 40 - 80 = (qualcosa) - 120.
    ; Se (qualcosa) è 16N, allora 16N - 120 è 16M - 8. Quindi (qualcosa) - 120 non è allineato.
    ; La ABI richiede che RSP sia 16-byte allineato *prima* di una CALL.
    ; All'ingresso, RSP % 16 == 8. Dopo push RBP, RSP % 16 == 0.
    ; Dopo 5 push di registri callee-saved, RSP % 16 == 0 - 5*8 = 0 - 40 = -40.
    ; -40 % 16 = -40 + 3*16 = -40 + 48 = 8.
    ; Quindi, dopo i push dei callee-saved, RSP % 16 == 8.
    ; Per mantenere l'allineamento per future chiamate, lo spazio allocato (80 byte in questo caso)
    ; dovrebbe essere tale che 80 % 16 == 8, oppure allocare un po' di più.
    ; 80 % 16 == 0. Quindi (RSP % 16 == 8) - (80 % 16 == 0) = 8. RSP rimane non allineato.
    ; Per correggere, allocare 80 + 8 = 88 bytes, o 80 - 8 = 72 (se basta).
    ; Dato che i commenti originali allocano 80, manteniamo 80 per ora.
    ; Se ci fossero chiamate `call` qui, avremmo bisogno di allineare rsp.

    sub     rsp, 80       ; Alloca 80 byte per le variabili locali

    ; Mappa variabili locali (offset da rbp):
    ; [rbp-8]:  scalare_costante (double)
    ; [rbp-16]: funz1_val (int)
    ; [rbp-24]: numPages_val (int)
    ; [rbp-32]: max_outer_iterations_val (int)
    ; [rbp-40]: ret_ptr (ptr)
    ; [rbp-48]: indici_ptr (ptr)
    ; [rbp-56]: d_ptr (ptr)
    ; [rbp-64]: somma_ptr (ptr)
    ; [rbp-72]: tranMatInv_ptr (ptr)
    ; [rbp-80]: tranMatParam_ptr (ptr)

    ; Copia parametri nelle variabili locali sullo stack
    mov     [rbp - 72], rdi ; tranMatInv_ptr
    mov     dword [rbp - 24], esi ; numPages_val (esi è 32-bit)
                                ; rdx (decay) è ignorato
    mov     dword [rbp - 32], ecx ; max_outer_iterations_val (ecx è 32-bit)
    mov     [rbp - 48], r8  ; indici_ptr
    mov     [rbp - 56], r9  ; d_ptr

    mov     rax, [rbp + 16] ; ret (parametro dallo stack)
    mov     [rbp - 40], rax ; ret_ptr (locale)
    mov     rax, [rbp + 24] ; somma (parametro dallo stack)
    mov     [rbp - 64], rax ; somma_ptr (locale)
    mov     eax, [rbp + 32] ; funz1 (parametro dallo stack, letto come int 32-bit)
    mov     [rbp - 16], eax ; funz1_val (locale)
    mov     rax, [rbp + 40] ; tranMatParam (parametro dallo stack)
    mov     [rbp - 80], rax ; tranMatParam_ptr (locale)

    ; r12 usato per contatore loop i
    ; r13 usato per contatore loop j
    ; r14 usato per contatore loop k
    ; xmm15 (o un xmm libero) per riga (double)

    ; Inizializza scalare_costante se funz1
    mov     eax, [rbp - 16] ; Carica funz1_val
    test    eax, eax        ; Controlla se funz1 è true (non-zero)
    jz      .L_calc_scalare_end_64

    ; scalare_costante = (1 - alfaI) / (type) numPages;
    movsd   xmm0, [rel one_double]   ; Carica 1.0
    movsd   xmm1, [rel alfaI]        ; Carica alfaI
    subsd   xmm0, xmm1               ; 1.0 - alfaI

    mov     eax, [rbp - 24]          ; Carica numPages_val (int)
    cvtsi2sd xmm1, eax               ; Converte numPages (da eax) in double
    divsd   xmm0, xmm1               ; (1 - alfaI) / (double)numPages
    movsd   [rbp - 8], xmm0          ; Salva scalare_costante

.L_calc_scalare_end_64:

    ; Loop esterno: for i < max_outer_iterations
    mov     r12, 0        ; i = 0 (r12 come i)
.L_outer_loop_i_64:
    mov     eax, [rbp - 32] ; Carica max_outer_iterations_val (int)
    cdqe                    ; Estende eax a rax (per confronto con r12q)
    cmp     r12, rax        ; Confronta i (r12) con max_outer_iterations_val (rax)
    jge     .L_outer_loop_end_i_64 ; Se i >= max_outer_iterations, fine loop

    ; Loop medio: for j < numPages
    mov     r13, 0        ; j = 0 (r13 come j)
.L_middle_loop_j_64:
    mov     eax, [rbp - 24] ; Carica numPages_val (int)
    cdqe                    ; Estende eax a rax (per confronto con r13q)
    cmp     r13, rax        ; Confronta j (r13) con numPages_val (rax)
    jge     .L_middle_loop_end_j_64 ; Se j >= numPages, fine loop

    ; Azioni dentro il loop j (dipende da funz1)
    mov     eax, [rbp - 16] ; Carica funz1_val
    test    eax, eax        ; Controlla se funz1 è true
    jnz     .L_funz1_branch_64 ; Se funz1 è true

    ; Ramo Else (!funz1)
    ; somma[i] = unoAlfaB * d[i]; (Condizionato da somma e d non NULL)
    mov     rbx, [rbp - 64] ; Carica somma_ptr
    test    rbx, rbx        ; Controlla se somma è NULL
    jz      .L_skip_somma_d_access_64 ; Salta se somma è NULL

    mov     rcx, [rbp - 56] ; Carica d_ptr
    test    rcx, rcx        ; Controlla se d è NULL
    jz      .L_skip_somma_d_access_64 ; Salta se d è NULL

    ; Sia somma che d non sono NULL
    ; r12 contiene i
    mov     rax, r12        ; Copia i in rax per indirizzamento
    lea     rdx, [rcx + rax*8] ; Calcola indirizzo di d[i] (d_ptr + i * 8)
    movsd   xmm1, [rdx]     ; Carica d[i]

    movsd   xmm0, [rel unoAlfaB] ; Carica unoAlfaB
    mulsd   xmm0, xmm1      ; unoAlfaB * d[i]

    lea     rdx, [rbx + rax*8] ; Calcola indirizzo di somma[i] (somma_ptr + i * 8)
    movsd   [rdx], xmm0     ; Salva risultato in somma[i]

.L_skip_somma_d_access_64:
    jmp     .L_after_funz1_branch_64 ; Vai al calcolo di riga

.L_funz1_branch_64:
    ; Se funz1 è true
    ; if i == 0
    cmp     r12, 0          ; Confronta i (r12) con 0
    jnz     .L_skip_i_zero_init_64 ; Se i != 0, salta init

    ; Blocco i == 0
    ; indici[i] = i; (Condizionato da indici non NULL)
    mov     rbx, [rbp - 48] ; Carica indici_ptr
    test    rbx, rbx        ; Controlla se indici è NULL
    jz      .L_skip_indici_access_64

    mov     rax, r12        ; Copia i (che è 0) in rax
    lea     rcx, [rbx + rax*4] ; Calcola indirizzo di indici[i] (indici_ptr + i * 4)
    mov     [rcx], eax      ; Salva i (eax contiene la parte bassa di rax, che è 0) in indici[i] (mov dword)

.L_skip_indici_access_64:
    ; d[i] = 0; (Condizionato da d non NULL)
    mov     rbx, [rbp - 56] ; Carica d_ptr
    test    rbx, rbx        ; Controlla se d è NULL
    jz      .L_skip_d_zero_access_64

    mov     rax, r12        ; Copia i (che è 0) in rax
    lea     rcx, [rbx + rax*8] ; Calcola indirizzo di d[i] (d_ptr + i * 8)
    movsd   xmm0, [rel zero_double] ; Carica 0.0
    movsd   [rcx], xmm0     ; Salva 0.0 in d[i]

.L_skip_d_zero_access_64:
.L_skip_i_zero_init_64:

.L_after_funz1_branch_64:
    ; type riga = 0
    movsd   xmm15, [rel zero_double] ; Carica 0.0 in xmm15 (usato per riga)

    ; Loop più interno: for k < numPages
    mov     r14, 0        ; k = 0 (r14 come k)
.L_innermost_loop_k_64:
    mov     eax, [rbp - 24] ; Carica numPages_val (int)
    cdqe                    ; Estende eax a rax
    cmp     r14, rax        ; Confronta k (r14) con numPages_val (rax)
    jge     .L_innermost_loop_end_k_64 ; Se k >= numPages, fine loop

    ; Calcolo dentro il loop k (dipende da funz1)
    ; xmm15 contiene riga corrente
    mov     eax, [rbp - 16] ; Carica funz1_val
    test    eax, eax        ; Controlla se funz1 è true
    jnz     .L_funz1_calc_k_loop_64

    ; Ramo Else (!funz1 calcolo)
    ; riga = riga + alfaB * ret[j] * tranMatParam[i * numPages + j];
    movsd   xmm1, [rel alfaB]      ; Carica alfaB in xmm1

    mov     rbx, [rbp - 40]        ; Carica ret_ptr
    mov     rax, r13               ; Copia j (r13) in rax
    lea     rdx, [rbx + rax*8]     ; Indirizzo di ret[j]
    movsd   xmm2, [rdx]            ; Carica ret[j] in xmm2

    mov     rbx, [rbp - 80]        ; Carica tranMatParam_ptr
    movsxd  rax, dword [rbp - 24]  ; Carica numPages (int) e estende a 64-bit in rax
    mov     rdx, r12               ; Copia i (r12) in rdx
    imul    rdx, rax               ; rdx = i * numPages
    add     rdx, r13               ; rdx = (i * numPages) + j
    lea     rcx, [rbx + rdx*8]     ; Indirizzo di tranMatParam[offset]
    movsd   xmm3, [rcx]            ; Carica tranMatParam[offset] in xmm3

    mulsd   xmm1, xmm2             ; xmm1 = alfaB * ret[j]
    mulsd   xmm1, xmm3             ; xmm1 = (alfaB * ret[j]) * tranMatParam[offset]
    addsd   xmm15, xmm1            ; riga (xmm15) = riga + xmm1
    jmp     .L_end_k_calc_64

.L_funz1_calc_k_loop_64:
    ; Calcolo funz1
    ; riga = riga + alfaI * ret[j] * tranMatInv[i * numPages + j];
    movsd   xmm1, [rel alfaI]      ; Carica alfaI in xmm1

    mov     rbx, [rbp - 40]        ; Carica ret_ptr (s è ret)
    mov     rax, r13               ; Copia j (r13) in rax
    lea     rdx, [rbx + rax*8]     ; Indirizzo di ret[j] (s[j])
    movsd   xmm2, [rdx]            ; Carica ret[j] (s[j]) in xmm2

    mov     rbx, [rbp - 72]        ; Carica tranMatInv_ptr
    movsxd  rax, dword [rbp - 24]  ; Carica numPages (int) e estende a 64-bit in rax
    mov     rdx, r12               ; Copia i (r12) in rdx
    imul    rdx, rax               ; rdx = i * numPages
    add     rdx, r13               ; rdx = (i * numPages) + j
    lea     rcx, [rbx + rdx*8]     ; Indirizzo di tranMatInv[offset]
    movsd   xmm3, [rcx]            ; Carica tranMatInv[offset] in xmm3

    mulsd   xmm1, xmm2             ; xmm1 = alfaI * ret[j]
    mulsd   xmm1, xmm3             ; xmm1 = (alfaI * ret[j]) * tranMatInv[offset]
    addsd   xmm15, xmm1            ; riga (xmm15) = riga + xmm1

.L_end_k_calc_64:
    inc     r14                    ; Incrementa k (r14)
    jmp     .L_innermost_loop_k_64

.L_innermost_loop_end_k_64:
    ; xmm15 contiene il valore finale di 'riga' per questa combinazione (i,j)
    mov     eax, [rbp - 16] ; Carica funz1_val
    test    eax, eax        ; Controlla se funz1 è true
    jnz     .L_funz1_assign_ret_64

    ; Ramo Else (!funz1 assegnazione)
    ; ret[i] = riga + somma[i]; (Condizionato da somma non NULL)
    mov     rbx, [rbp - 40]   ; Carica ret_ptr
    mov     rcx, [rbp - 64]   ; Carica somma_ptr
    test    rcx, rcx          ; Controlla se somma è NULL
    jz      .L_final_skip_ret_assign_64 ; Salta se somma è NULL

    ; xmm15 contiene riga
    mov     rax, r12          ; Copia i (r12) in rax
    lea     rdx, [rcx + rax*8] ; Indirizzo di somma[i]
    movsd   xmm1, [rdx]       ; Carica somma[i] in xmm1
    addsd   xmm15, xmm1       ; riga (xmm15) = riga + somma[i]

    lea     rdx, [rbx + rax*8] ; Indirizzo di ret[i]
    movsd   [rdx], xmm15      ; Salva in ret[i]
    jmp     .L_final_skip_ret_assign_64

.L_funz1_assign_ret_64:
    ; Assegnazione funz1
    ; ret[i] = riga + scalare_costante;
    mov     rbx, [rbp - 40]   ; Carica ret_ptr
    ; xmm15 contiene riga
    movsd   xmm1, [rbp - 8]   ; Carica scalare_costante in xmm1
    addsd   xmm15, xmm1       ; riga (xmm15) = riga + scalare_costante

    mov     rax, r12          ; Copia i (r12) in rax
    lea     rdx, [rbx + rax*8] ; Indirizzo di ret[i]
    movsd   [rdx], xmm15      ; Salva in ret[i]

.L_final_skip_ret_assign_64:
    inc     r13               ; Incrementa j (r13)
    jmp     .L_middle_loop_j_64

.L_middle_loop_end_j_64:
    inc     r12               ; Incrementa i (r12)
    jmp     .L_outer_loop_i_64

.L_outer_loop_end_i_64:
    ; Ritorna ret (il puntatore)
    mov     rax, [rbp - 40] ; Carica ret_ptr in rax (valore di ritorno)

    ; Epilogo
    add     rsp, 80         ; Dealloca spazio variabili locali
    pop     r15
    pop     r14
    pop     r13
    pop     r12
    pop     rbx
    pop     rbp
    ret