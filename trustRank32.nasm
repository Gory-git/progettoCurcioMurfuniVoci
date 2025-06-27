; Sposta le dichiarazioni extern e global all'inizio
extern alfaI        ; float (dichiarato in C)
extern alfaB        ; float (dichiarato in C)
extern unoAlfaB     ; float (dichiarato in C)
global funzione_unica
global one_float    ; costante float 1.0f

section .note.GNU-stack noalloc noexec nowrite progbits

section .data
    ; Correzione: rimuovi 'f' dal valore float
    one_float dd 1.0

section .text
    ; global funzione_unica ; Già dichiarato all'inizio

; VECTOR funzione_unica(MATRIX tranMat, int numPages, type decay, int max_outer_iterations, int* indici, VECTOR d, VECTOR ret, void* somma_param, bool funz1, MATRIX tranMatParam);
; Argomenti sullo stack (da [ebp+8] in su - convenzione cdecl):
; [ebp+8]  = Matrix (tranMat o tranMatInv - la matrice effettiva da usare)
; [ebp+12] = numPages (int) -> useremo ESI
; [ebp+16] = decay (float) -> useremo XMM3 temporaneamente nel loop k, or loaded for scalar multiplication
; [ebp+20] = max_outer_iterations (int) -> useremo EDI
; [ebp+24] = indici (int*) - base address of the indici vector
; [ebp+28] = d (VECTOR) - base address of the d vector
; [ebp+32] = ret (VECTOR) - base address of the ret vector (vector being updated)
; [ebp+36] = somma_param (VECTOR if funz1=false, float* if funz1=true) - base address of somma vector or pointer to scalar somma
; [ebp+40] = funz1 (bool) - boolean flag
; [ebp+44] = tranMatParam (MATRIX) - Not used in this assembly implementation, [ebp+8] is the correct matrix

funzione_unica:
    push ebp             ; Salva il vecchio EBP
    mov ebp, esp         ; Imposta il nuovo frame pointer
    sub esp, 32          ; Alloca spazio sul stack (ad esempio per allineamento o variabili locali non usate esplicitamente qui)

    ; Carica i valori costanti nei registri appropriati
    mov esi, [ebp+12]    ; ESI = numPages
    mov edi, [ebp+20]    ; EDI = max_outer_iterations

    ; Loop esterno per il numero massimo di iterazioni (outer_iter)
    xor ecx, ecx         ; ECX = outer_iter = 0

.outer_iterations_loop:
    cmp ecx, edi         ; Confronta outer_iter con max_outer_iterations
    jge .done            ; Se outer_iter >= max_outer_iterations, esci dalla funzione

    ; Loop interno sui singoli elementi del vettore (i)
    ; Calcoliamo ret[i] per ogni i da 0 a numPages-1
    xor eax, eax         ; EAX = i = 0

.i_loop:
    cmp eax, esi         ; Confronta i con numPages
    jge .next_outer_iteration ; Se i >= numPages, una iterazione esterna è completa

    ; --- Inizia il calcolo di ret[i] per l'elemento corrente i ---

    ; Controlla il flag funz1 ([ebp+40]) per scegliere il percorso (selectSeed o computeScores)
    mov edx, [ebp+40]    ; EDX = funz1 flag (usa EDX come registro temporaneo)
    cmp edx, 0
    je .funz0_prep       ; Se funz1 è falso, salta al setup specifico per funz0

    ; --- Percorso funz1 = true (corrisponde a selectSeed) ---
    ; Gestisce le inizializzazioni che avvengono solo nella prima iterazione esterna (outer_iter == 0)
    cmp ecx, 0           ; Confronta outer_iter (ECX) con 0
    jne .skip_init_funz1 ; Se outer_iter > 0, salta le inizializzazioni

    ; If outer_iter == 0: indici[i] = i; d[i] = 0; ret[i] = 1.0;
    ; NOTE: The indici list should be sorted after selectSeed,
    ; but this function only sets the initial values of the passed indici vector.
    push edx             ; Save EDX temporarily
    mov edx, [ebp+24]    ; EDX = base address of the indici vector
    mov [edx + eax*4], eax ; indici[i] = i  (using i = EAX)  ; Assuming int is 4 bytes
    pop edx              ; Restore EDX

    push edx             ; Save EDX temporarily
    mov edx, [ebp+28]    ; EDX = base address of the d vector
    ; Assuming type (float) is 4 bytes and int is 4 bytes, zero out the 4 bytes
    mov dword [edx + eax*4], 0 ; d[i] = 0 (using i = EAX)
    pop edx              ; Restore EDX

    push edx             ; Save EDX temporarily
    mov edx, [ebp+32]    ; EDX = base address of the ret vector (s in selectSeed)
    movss xmm0, [one_float] ; Load the float value 1.0f in XMM0
    movss [edx + eax*4], xmm0 ; ret[i] = 1.0f (using i = EAX)
    pop edx              ; Restore EDX


.skip_init_funz1:
    ; Initialize the accumulator for the matrix-vector multiplication part
    xorps xmm0, xmm0     ; XMM0 = riga = 0.0f (float)

    ; Innermost loop for matrix-vector multiplication (k index)
    ; Calculate sum_k (Matrix[i][k] * Vector[k])
    xor ebx, ebx         ; EBX = k = 0

.k_loop_funz1:
    cmp ebx, esi         ; Compare k with numPages
    jge .riga_done_funz1 ; If k >= numPages, the sum for ret[i] is complete

    ; Calculate the linear index for matrix access: i * numPages + k
    ; i is in EAX, numPages is in ESI, k is in EBX
    push eax             ; Save i (EAX) because EAX is needed elsewhere
    mov edx, eax         ; EDX = i (save i temporarily in EDX)
    imul edx, esi        ; EDX = i * numPages
    add edx, ebx         ; EDX = i * numPages + k  <-- Using EBX (k) for the column/vector index

    ; Load the matrix element: tranMatInv[i * numPages + k]
    ; Matrix base address is [ebp+8]
    mov eax, [ebp+8]     ; EAX = base address of tranMatInv
    movss xmm1, [eax + edx*4] ; XMM1 = tranMatInv[i * numPages + k] (float matrix element)
    pop eax              ; Restore i into EAX

    ; Load the vector element: ret[k]
    ; ret vector base address is [ebp+32]
    push eax             ; Save i (EAX) again
    mov eax, [ebp+32]    ; EAX = base address of the ret vector
    movss xmm2, [eax + ebx*4] ; XMM2 = ret[k] (float vector element) <-- Using EBX (k)
    pop eax              ; Restore i into EAX

    ; Calculate the product term: Matrix[...] * ret[...]
    mulss xmm1, xmm2     ; XMM1 = tranMatInv[...] * ret[...]

    ; Multiply by the decay factor (alfaI)
    ; Decay factor is [ebp+16]
    movss xmm3, [ebp+16] ; XMM3 = decay (alfaI) (load the float from the stack)
    mulss xmm1, xmm3     ; XMM1 = decay * product_term

    ; Accumulate the term into XMM0 (riga)
    addss xmm0, xmm1     ; XMM0 = XMM0 + decay * term

    inc ebx              ; Increment k
    jmp .k_loop_funz1    ; Continue k loop

.riga_done_funz1:
    ; XMM0 now contains the sum of products (matrix-vector part)

    ; Add the scalar somma term (from the float pointer at [ebp+36])
    ; somma_param at [ebp+36] is float* in funz1 mode
    mov edx, [ebp+36]    ; EDX = pointer to the scalar somma (&somma)
    movss xmm4, [edx]    ; XMM4 = value of the scalar somma (*somma_param)
    addss xmm0, xmm4     ; XMM0 = riga + scalar_somma

    ; Store the final result into ret[i]
    mov edx, [ebp+32]    ; EDX = base address of the ret vector
    movss [edx + eax*4], xmm0 ; ret[i] = final result (using i = EAX)

    ; Jump to the end of ret[i] calculation to proceed to the next i element
    jmp .next_i          ; Jump to i increment (join point)

    ; --- Path funz0 = false (corresponds to computeScores) ---
.funz0_prep:
    ; Calculate the additive part: additive_term[i] = (1-alfaB) * d[i]
    ; somma_param at [ebp+36] is VECTOR additive_term in funz0 mode (used to store the *intermediate* results of d * (1-alfaB))
    ; the d vector is at [ebp+28]
    ; unoAlfaB is a global variable (declared extern at start and defined in C)

    push edx             ; Save EDX - This push/pop pair is needed here
    mov edx, [ebp+28]    ; EDX = base address of the d vector
    movss xmm1, [edx + eax*4] ; XMM1 = d[i] (using i = EAX)
    pop edx              ; Restore EDX

    movss xmm2, [unoAlfaB] ; Load the global variable unoAlfaB (1 - alfaB)
    mulss xmm1, xmm2     ; XMM1 = d[i] * unoAlfaB  (This is the value of additive_term[i])

    ; Store the calculated additive_term[i] into the vector pointed by somma_param
    push edx             ; Save EDX - This push/pop pair is needed here
    mov edx, [ebp+36]    ; EDX = base address of the somma/additive_term vector (somma_param)
    movss [edx + eax*4], xmm1 ; somma[i] = d[i] * unoAlfaB (using i = EAX)
    pop edx              ; Restore EDX

    ; Initialize the accumulator for the matrix-vector multiplication part
    xorps xmm0, xmm0     ; XMM0 = riga = 0.0f (float)

    ; Initialize the innermost loop for matrix-vector multiplication (k index)
    xor ebx, ebx         ; EBX = k = 0

.k_loop_funz0:
    cmp ebx, esi         ; Compare k with numPages
    jge .riga_done_funz0 ; If k >= numPages, the sum for ret[i] is complete

    ; Calculate the linear index for matrix access: i * numPages + k
    ; i is in EAX, numPages is in ESI, k is in EBX
    push eax             ; Save i (EAX)
    mov edx, eax         ; EDX = i (save i temporarily in EDX)
    imul edx, esi        ; EDX = i * numPages
    add edx, ebx         ; EDX = i * numPages + k  <-- Using EBX (k)

    ; Load the matrix element: tranMat[i * numPages + k]
    ; Matrix base address is [ebp+8]
    mov eax, [ebp+8]     ; EAX = base address of tranMat
    movss xmm1, [eax + edx*4] ; XMM1 = tranMat[i * numPages + k] (float matrix element)
    pop eax              ; Restore i

    ; Load the vector element: ret[k]
    ; ret vector base address is [ebp+32]. This is the vector from the *previous* iteration.
    push eax             ; Save i (EAX)
    mov eax, [ebp+32]     ; EAX = base address of the ret vector (old values)
    movss xmm2, [eax + ebx*4] ; XMM2 = ret[k] (float vector element) <-- Using EBX (k)
    pop eax              ; Restore i

    ; Calculate the product term: Matrix[...] * ret[...]
    mulss xmm1, xmm2     ; XMM1 = tranMat[...] * ret[...]

    ; Multiply by the decay factor (alfaB)
    ; Decay factor is [ebp+16]
    movss xmm3, [ebp+16] ; XMM3 = decay (alfaB) (load the float from the stack)
    mulss xmm1, xmm3     ; XMM1 = decay * product_term

    ; Accumulate the term into XMM0 (riga)
    addss xmm0, xmm1     ; XMM0 = XMM0 + decay * term

    inc ebx              ; Increment k
    jmp .k_loop_funz0    ; Continue k loop

.riga_done_funz0:
    ; XMM0 has sum_k(alfaB * tranMat[i][k] * ret_old[k])

    ; Add the additive_term[i] vector element
    ; additive_term vector is pointed by somma_param [ebp+36]
    ; i index is EAX
    mov edx, [ebp+36]    ; EDX = base address of the additive_term vector (somma_param)
    movss xmm4, [edx + eax*4] ; XMM4 = somma[i] (using i = EAX)
    addss xmm0, xmm4     ; XMM0 = riga + additive_term[i]

    ; Store the final result for ret[i] into the *ret* vector (in-place update as per C)
    mov edx, [ebp+32]    ; EDX = base address of ret vector
    movss [edx + eax*4], xmm0 ; ret[i] = final result (using i = EAX)

    jmp .next_i ; Continue to next i

    ; --- End calculation of ret[i] for the current i element ---

.next_i: ; Join point for both paths (funz1 and funz0)
    inc eax              ; Increment i
    jmp .i_loop          ; Continue i loop (i=0..numPages-1)

.next_outer_iteration:
    ; After iterating over all i (0..numPages-1), one outer iteration is complete.
    ; The ret vector has been updated in-place.
    ; The additive_term vector (somma) is already calculated for all i in the first i-loop iteration of this outer iteration.
    inc ecx              ; Increment outer_iter
    jmp .outer_iterations_loop ; Continue outer loop

.done:
    ; The function is finished. The return value is the pointer to the updated ret vector.
    mov eax, [ebp+32]    ; Load the pointer to ret into EAX for the return value

    leave                ; Restore EBP and SP to their values before function entry
    ret                  ; Return from subroutine