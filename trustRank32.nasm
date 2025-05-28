section .data
align 8
alfaI:          dq 0.0     ; Define these constants (example values)
alfaB:          dq 0.0
unoAlfaB:       dq 0.0
one_double:     dq 1.0
zero_double:    dq 0.0

section .text
global funzione_unica

; VECTOR funzione_unica(MATRIX tranMatInv, int numPages, type decay, int max_outer_iterations, int* indici, VECTOR d, VECTOR ret, VECTOR somma, bool funz1, MATRIX tranMatParam)

; Parametri sullo stack (da ebp+8 in su, ordine cdecl):
; [ebp+8]:  tranMatInv (ptr)
; [ebp+12]: numPages (int)
; [ebp+16]: decay (type - non usato)
; [ebp+20]: max_outer_iterations (int)
; [ebp+24]: indici (ptr)
; [ebp+28]: d (ptr)
; [ebp+32]: ret (ptr)
; [ebp+36]: somma (ptr)
; [ebp+40]: funz1 (int, 0 o 1)
; [ebp+44]: tranMatParam (ptr)

funzione_unica:
    push    ebp
    mov     ebp, esp

    ; Salvataggio registri callee-saved e altri usati localmente
    push    ebx
    push    esi
    push    edi
    push    eax
    push    ecx
    push    edx
    push    eax            ; Dummy push come nell'originale (7 push totali = 28 bytes)
                           ; esp ora è ebp - 28

    ; Allocazione spazio per variabili locali
    ; Necessari 60 bytes per le variabili.
    ; Per allineamento a 16 byte di ESP (ebp - 28 - X deve essere multiplo di 16):
    ; Se ebp è ...0, ebp-28 è ...4. Serve X = ...4.
    ; 60 byte necessari. Il successivo valore che termina con 4 è 68.
    sub     esp, 68        ; Alloca 68 bytes (60 usati + 8 padding)
                           ; esp ora è ebp - 28 - 68 = ebp - 96 (allineato a 16-byte)

    ; Mappa delle variabili locali (offset da ebp):
    ; [ebp-36]: scalare_costante (double, 8 bytes) ; ebp - 28 - 8
    ; [ebp-40]: i (int, 4 bytes)                   ; ebp - 36 - 4
    ; [ebp-44]: j (int, 4 bytes)
    ; [ebp-48]: k (int, 4 bytes)
    ; [ebp-56]: riga (double, 8 bytes)
    ; [ebp-60]: funz1_val (int, 4 bytes)
    ; [ebp-64]: numPages_val (int, 4 bytes)
    ; [ebp-68]: max_outer_iterations_val (int, 4 bytes)
    ; [ebp-72]: ret_ptr (ptr, 4 bytes)
    ; [ebp-76]: indici_ptr (ptr, 4 bytes)
    ; [ebp-80]: d_ptr (ptr, 4 bytes)
    ; [ebp-84]: somma_ptr (ptr, 4 bytes)
    ; [ebp-88]: tranMatInv_ptr (ptr, 4 bytes)
    ; [ebp-92]: tranMatParam_ptr (ptr, 4 bytes)

    ; Copia parametri nelle variabili locali sullo stack
    mov     eax, [ebp + 24]  ; indici (parametro)
    mov     [ebp - 76], eax  ; indici_ptr (locale)
    mov     eax, [ebp + 28]  ; d
    mov     [ebp - 80], eax  ; d_ptr
    mov     eax, [ebp + 32]  ; ret
    mov     [ebp - 72], eax  ; ret_ptr
    mov     eax, [ebp + 36]  ; somma
    mov     [ebp - 84], eax  ; somma_ptr
    mov     eax, [ebp + 8]   ; tranMatInv
    mov     [ebp - 88], eax  ; tranMatInv_ptr
    mov     eax, [ebp + 44]  ; tranMatParam
    mov     [ebp - 92], eax  ; tranMatParam_ptr

    mov     eax, [ebp + 12]  ; numPages
    mov     [ebp - 64], eax  ; numPages_val
    mov     eax, [ebp + 20]  ; max_outer_iterations
    mov     [ebp - 68], eax  ; max_outer_iterations_val
    mov     eax, [ebp + 40]  ; funz1
    mov     [ebp - 60], eax  ; funz1_val

    ; Inizializza scalare_costante se funz1
    mov     eax, [ebp - 60] ; Carica funz1_val
    test    eax, eax        ; Controlla se funz1 è true (non-zero)
    jz      L_calc_scalare_end

    ; scalare_costante = (1 - alfaI) / (type) numPages;
    movsd   xmm0, [one_double]   ; Carica 1.0
    movsd   xmm1, [alfaI]        ; Carica alfaI
    subsd   xmm0, xmm1           ; 1.0 - alfaI

    mov     eax, [ebp - 64]      ; Carica numPages_val
    cvtsi2sd xmm1, eax           ; Converte numPages in double
    divsd   xmm0, xmm1           ; (1 - alfaI) / (double)numPages
    movsd   [ebp - 36], xmm0     ; Salva scalare_costante

L_calc_scalare_end:

    ; Loop esterno: for i < max_outer_iterations
    mov     dword [ebp - 40], 0   ; i = 0
L_outer_loop_i:
    mov     eax, [ebp - 40]     ; Carica i
    cmp     eax, [ebp - 68]     ; Confronta i con max_outer_iterations_val
    jge     L_outer_loop_end_i  ; Se i >= max_outer_iterations, fine loop

    ; Loop medio: for j < numPages
    mov     dword [ebp - 44], 0   ; j = 0
L_middle_loop_j:
    mov     ebx, [ebp - 44]     ; Carica j
    cmp     ebx, [ebp - 64]     ; Confronta j con numPages_val
    jge     L_middle_loop_end_j ; Se j >= numPages, fine loop

    ; Azioni dentro il loop j (dipende da funz1)
    mov     eax, [ebp - 60] ; Carica funz1_val
    test    eax, eax        ; Controlla se funz1 è true
    jnz     L_funz1_branch  ; Se funz1 è true

    ; Ramo Else (!funz1)
    ; somma[i] = unoAlfaB * d[i]; (Condizionato da somma e d non NULL)
    mov     ecx, [ebp - 84] ; Carica somma_ptr
    test    ecx, ecx        ; Controlla se somma è NULL
    jz      L_skip_somma_d_access ; Salta se somma è NULL

    mov     edx, [ebp - 80] ; Carica d_ptr
    test    edx, edx        ; Controlla se d è NULL
    jz      L_skip_somma_d_access ; Salta se d è NULL

    ; Sia somma che d non sono NULL
    mov     eax, [ebp - 40]   ; Carica i (indice loop esterno)
    mov     esi, eax          ; Usa esi per i (indice)
    imul    esi, 8            ; Calcola offset in byte per double (i * 8)

    movsd   xmm0, [unoAlfaB]  ; Carica unoAlfaB
                              ; edx contiene già d_ptr
    movsd   xmm1, [edx + esi] ; Carica d[i]
    mulsd   xmm0, xmm1        ; unoAlfaB * d[i]

                              ; ecx contiene già somma_ptr
    movsd   [ecx + esi], xmm0 ; Salva risultato in somma[i]

L_skip_somma_d_access:
    jmp     L_after_funz1_branch ; Vai al calcolo di riga

L_funz1_branch:
    ; Se funz1 è true
    ; if i == 0
    mov     eax, [ebp - 40] ; Carica i
    test    eax, eax        ; Controlla se i == 0
    jnz     L_skip_i_zero_init ; Se i != 0, salta init

    ; Blocco i == 0
    ; indici[i] = i; (Condizionato da indici non NULL)
    mov     ecx, [ebp - 76] ; Carica indici_ptr
    test    ecx, ecx        ; Controlla se indici è NULL
    jz      L_skip_indici_access_f1

    mov     eax, [ebp - 40]   ; Carica i (che è 0)
    mov     esi, eax          ; Usa esi per i (indice)
    imul    esi, 4            ; Calcola offset in byte per int (i * 4)

    mov     [ecx + esi], eax  ; Salva i in indici[i]

L_skip_indici_access_f1:

    ; d[i] = 0; (Condizionato da d non NULL)
    mov     ecx, [ebp - 80] ; Carica d_ptr
    test    ecx, ecx        ; Controlla se d è NULL
    jz      L_skip_d_zero_access_f1

    mov     eax, [ebp - 40]   ; Carica i (che è 0)
    mov     esi, eax          ; Usa esi per i (indice)
    imul    esi, 8            ; Calcola offset in byte per double (i * 8)

    movsd   xmm0, [zero_double] ; Carica 0.0
    movsd   [ecx + esi], xmm0   ; Salva 0.0 in d[i]

L_skip_d_zero_access_f1:
L_skip_i_zero_init:

L_after_funz1_branch:
    ; type riga = 0
    movsd   xmm0, [zero_double] ; Carica 0.0
    movsd   [ebp - 56], xmm0    ; Salva in riga

    ; Loop più interno: for k < numPages
    mov     dword [ebp - 48], 0   ; k = 0
L_innermost_loop_k:
    mov     edx, [ebp - 48]     ; Carica k
    cmp     edx, [ebp - 64]     ; Confronta k con numPages_val
    jge     L_innermost_loop_end_k ; Se k >= numPages, fine loop

    ; Calcolo dentro il loop k (dipende da funz1)
    ; *** ATTENZIONE: La logica originale AT&T qui usa gli indici i e j, NON k.
    ; *** Questo significa che lo stesso termine viene sommato a 'riga' numPages volte.
    mov     eax, [ebp - 60] ; Carica funz1_val
    test    eax, eax        ; Controlla se funz1 è true
    jnz     L_funz1_calc_k_loop

    ; Ramo Else (!funz1 calcolo)
    ; riga = riga + alfaB * ret[j] * tranMatParam[i * numPages + j];
    movsd   xmm0, [ebp - 56]  ; Carica riga corrente
    movsd   xmm1, [alfaB]     ; Carica alfaB

    mov     esi, [ebp - 72]   ; Carica ret_ptr
    mov     eax, [ebp - 44]   ; Carica j (indice loop medio)
    imul    eax, 8            ; Offset per ret[j] (j * 8)
    movsd   xmm2, [esi + eax] ; Carica ret[j]

    mov     esi, [ebp - 92]   ; Carica tranMatParam_ptr
    mov     edi, [ebp - 40]   ; Carica i (indice loop esterno)
    mov     ecx, [ebp - 64]   ; Carica numPages_val
    imul    edi, ecx          ; edi = i * numPages
    add     edi, [ebp - 44]   ; edi = (i * numPages) + j
    imul    edi, 8            ; Offset per tranMatParam[i*numPages+j]
    movsd   xmm3, [esi + edi] ; Carica tranMatParam[i*numPages+j]

    mulsd   xmm1, xmm2        ; xmm1 = alfaB * ret[j]
    mulsd   xmm1, xmm3        ; xmm1 = (alfaB * ret[j]) * tranMatParam[i*numPages+j]
    addsd   xmm0, xmm1        ; riga = riga + xmm1
    movsd   [ebp - 56], xmm0  ; Salva riga aggiornata
    jmp     L_end_k_calc

L_funz1_calc_k_loop:
    ; Calcolo funz1
    ; riga = riga + alfaI * ret[j] * tranMatInv[i * numPages + j]; (s è ret)
    movsd   xmm0, [ebp - 56]  ; Carica riga corrente
    movsd   xmm1, [alfaI]     ; Carica alfaI

    mov     esi, [ebp - 72]   ; Carica ret_ptr (s è ret)
    mov     eax, [ebp - 44]   ; Carica j (indice loop medio)
    imul    eax, 8            ; Offset per ret[j] (j * 8)
    movsd   xmm2, [esi + eax] ; Carica ret[j] (o s[j])

    mov     esi, [ebp - 88]   ; Carica tranMatInv_ptr
    mov     edi, [ebp - 40]   ; Carica i (indice loop esterno)
    mov     ecx, [ebp - 64]   ; Carica numPages_val
    imul    edi, ecx          ; edi = i * numPages
    add     edi, [ebp - 44]   ; edi = (i * numPages) + j
    imul    edi, 8            ; Offset per tranMatInv[i*numPages+j]
    movsd   xmm3, [esi + edi] ; Carica tranMatInv[i*numPages+j]

    mulsd   xmm1, xmm2        ; xmm1 = alfaI * ret[j]
    mulsd   xmm1, xmm3        ; xmm1 = (alfaI * ret[j]) * tranMatInv[i*numPages+j]
    addsd   xmm0, xmm1        ; riga = riga + xmm1
    movsd   [ebp - 56], xmm0  ; Salva riga aggiornata

L_end_k_calc:
    inc     dword [ebp - 48]    ; Incrementa k
    jmp     L_innermost_loop_k  ; Torna all'inizio del loop k

L_innermost_loop_end_k:
    ; Dopo il loop k (ancora dentro il loop j)
    ; *** ATTENZIONE: Il loop k sopra ha aggiunto lo stesso termine (calcolato con i e j) numPages volte.
    ; *** ATTENZIONE: La seguente assegnazione a ret[i] avviene dentro il loop j.
    ; *** Questo significa che ret[i] sarà sovrascritto numPages volte, e il suo valore finale
    ; *** sarà basato sul calcolo quando j = numPages - 1.
    ; *** Questo comportamento è altamente sospetto e probabilmente un bug nella logica originale.

    mov     eax, [ebp - 60] ; Carica funz1_val
    test    eax, eax        ; Controlla se funz1 è true
    jnz     L_funz1_assign_ret

    ; Ramo Else (!funz1 assegnazione)
    ; ret[i] = riga + somma[i]; (Condizionato da somma non NULL)
    mov     ecx, [ebp - 72]   ; Carica ret_ptr
    mov     edx, [ebp - 84]   ; Carica somma_ptr
    test    edx, edx          ; Controlla se somma è NULL
    jz      L_skip_ret_assign_val ; Salta se somma è NULL

    movsd   xmm0, [ebp - 56]  ; Carica riga (che è numPages * termine_ij)

    mov     eax, [ebp - 40]   ; Carica i (indice loop esterno)
    mov     esi, eax          ; Usa esi per i
    imul    esi, 8            ; Offset per ret[i] e somma[i] (i * 8)

    movsd   xmm1, [edx + esi] ; Carica somma[i]
    addsd   xmm0, xmm1        ; riga + somma[i]

                              ; ecx contiene ancora ret_ptr
    movsd   [ecx + esi], xmm0 ; Salva in ret[i]
    jmp     L_final_skip_ret_assign

L_funz1_assign_ret:
    ; Assegnazione funz1
    ; ret[i] = riga + scalare_costante;
    mov     ecx, [ebp - 72]   ; Carica ret_ptr
    movsd   xmm0, [ebp - 56]  ; Carica riga (che è numPages * termine_ij)
    movsd   xmm1, [ebp - 36]  ; Carica scalare_costante
    addsd   xmm0, xmm1        ; riga + scalare_costante

    mov     eax, [ebp - 40]   ; Carica i (indice loop esterno)
    mov     esi, eax          ; Usa esi per i
    imul    esi, 8            ; Offset per ret[i] (i * 8)

                              ; ecx contiene ancora ret_ptr
    movsd   [ecx + esi], xmm0 ; Salva in ret[i]

L_skip_ret_assign_val:      ; Saltato qui se somma è NULL nel caso !funz1
L_final_skip_ret_assign:    ; Punto di skip comune

    inc     dword [ebp - 44]    ; Incrementa j
    jmp     L_middle_loop_j     ; Torna all'inizio del loop j

L_middle_loop_end_j:
    inc     dword [ebp - 40]    ; Incrementa i
    jmp     L_outer_loop_i      ; Torna all'inizio del loop i

L_outer_loop_end_i:

    ; Ritorna ret (il puntatore)
    mov     eax, [ebp - 72] ; Carica ret_ptr in eax (come da originale, ret viene ritornato)

    ; Epilogo
    add     esp, 68         ; Dealloca spazio variabili locali
    pop     eax             ; Ripristina dummy eax
    pop     edx
    pop     ecx
    pop     eax
    pop     edi
    pop     esi
    pop     ebx
    pop     ebp
    ret