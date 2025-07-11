; ----------------------------------------------------------
; VECTOR computeScores(MATRIX tranMat, type alfaB, int maxBias, VECTOR d, int numPages)
; ----------------------------------------------------------
; Translated from C++ to x86-64 assembly with AVX instructions
; Assumes 'type' is double (double precision floating point)
;
; Parameters in registers according to x64 calling convention:
; tranMat (MATRIX) -> rdi (1st parameter)
; alfaB (type/double) -> xmm0 (1st float parameter)
; maxBias (int) -> esi (2nd integer parameter)
; d (VECTOR) -> rdx (3rd parameter)
; numPages (int) -> ecx (4th parameter)
;
; Local registers usage:
; rax - general purpose, also used for return value
; rbx - used for indexing
; r8 - loop counter
; r9 - temporary storage
; r10 - temporary storage
; r11 - temporary storage
; r12 - storing numPages
; r13 - storing ret vector pointer
; r14 - storing somma vector pointer
; ymm0 - alfaB
; ymm1 - unoAlfaB (1-alfaB)
; ymm2-ymm7 - temporary calculations
; ----------------------------------------------------------

%include "sseutils64.nasm"

; External function declarations
extern copy_vector
extern alloc_vector
extern dealloc_vector

section .text
global computeScores

computeScores:
    start                       ; Macro from sseutils64.nasm - pushes rbp and sets up stack frame
    
    ; Save non-volatile registers we'll use
    push rbx
    push r12
    push r13
    push r14
    push r15
    
    ; Store parameters
    mov r12d, ecx              ; r12d = numPages
    
    ; Allocate ret vector: ret = copy_vector(d, numPages)
    mov rdi, rdx               ; First parameter (src = d)
    mov rsi, rcx               ; Second parameter (size = numPages)
    call copy_vector           ; Call copy_vector(d, numPages)
    mov r13, rax               ; r13 = pointer to ret vector
    
    ; Calculate unoAlfaB = 1 - alfaB
    vmovsd xmm1, [rel one]     ; xmm1 = 1.0
    vsubsd xmm1, xmm1, xmm0    ; xmm1 = 1.0 - alfaB
    
    ; Allocate somma vector
    mov rdi, r12               ; First parameter (n = numPages)
    call alloc_vector          ; Call alloc_vector(numPages)
    mov r14, rax               ; r14 = pointer to somma vector
    
    ; Initialize somma to zeros
    mov r8, 0                  ; Initialize loop counter
    vxorpd ymm10, ymm10, ymm10 ; Zero out ymm10 register
    
.zero_loop:
    cmp r8, r12                ; Compare counter with numPages
    jge .zero_end              ; If counter >= numPages, end loop
    
    mov rbx, r8                ; rbx = counter
    imul rbx, 8                ; rbx = counter * sizeof(double)
    vmovsd qword [r14+rbx], xmm10  ; somma[counter] = 0
    
    inc r8                     ; Increment counter
    jmp .zero_loop             ; Repeat loop
    
.zero_end:
    ; Begin bias loop (for (int b = 0; b < maxBias; b++))
    mov r15d, 0                ; r15d = bias counter (b)
    
.bias_loop:
    cmp r15d, esi              ; Compare b with maxBias
    jge .bias_end              ; If b >= maxBias, end loop
    
    ; Begin pages loop (for (int i = 0; i < numPages; i++))
    mov r8d, 0                 ; r8d = page counter (i)
    
.pages_loop:
    cmp r8d, r12d              ; Compare i with numPages
    jge .pages_end             ; If i >= numPages, end loop
    
    ; Calculate somma[i] = unoAlfaB * d[i]
    mov rbx, r8                ; rbx = i
    imul rbx, 8                ; rbx = i * sizeof(double)
    vmovsd xmm2, [rdx+rbx]     ; xmm2 = d[i]
    vmulsd xmm2, xmm2, xmm1    ; xmm2 = unoAlfaB * d[i]
    vmovsd [r14+rbx], xmm2     ; somma[i] = unoAlfaB * d[i]
    
    ; Initialize riga = 0
    vxorpd xmm3, xmm3, xmm3    ; xmm3 = riga = 0
    
    ; Begin inner loop (for (int j = 0; j < numPages; j++))
    mov r9d, 0                 ; r9d = inner loop counter (j)
    
.inner_loop:
    cmp r9d, r12d              ; Compare j with numPages
    jge .inner_end             ; If j >= numPages, end loop
    
    ; Calculate part of riga: riga += alfaB * ret[j] * tranMat[i*numPages + j]
    mov rbx, r9                ; rbx = j
    imul rbx, 8                ; rbx = j * sizeof(double)
    vmovsd xmm4, [r13+rbx]     ; xmm4 = ret[j]
    vmulsd xmm4, xmm4, xmm0    ; xmm4 = alfaB * ret[j]
    
    ; Calculate index for tranMat[i*numPages + j]
    mov r10, r8                ; r10 = i
    imul r10, r12              ; r10 = i * numPages
    add r10, r9                ; r10 = i * numPages + j
    imul r10, 8                ; r10 = (i * numPages + j) * sizeof(double)
    
    vmovsd xmm5, [rdi+r10]     ; xmm5 = tranMat[i*numPages + j]
    vmulsd xmm4, xmm4, xmm5    ; xmm4 = alfaB * ret[j] * tranMat[i*numPages + j]
    vaddsd xmm3, xmm3, xmm4    ; riga += alfaB * ret[j] * tranMat[i*numPages + j]
    
    inc r9d                    ; Increment j
    jmp .inner_loop            ; Repeat inner loop
    
.inner_end:
    ; Calculate ret[i] = riga + somma[i]
    mov rbx, r8                ; rbx = i
    imul rbx, 8                ; rbx = i * sizeof(double)
    vaddsd xmm3, xmm3, [r14+rbx] ; xmm3 = riga + somma[i]
    vmovsd [r13+rbx], xmm3     ; ret[i] = riga + somma[i]
    
    inc r8d                    ; Increment i
    jmp .pages_loop            ; Repeat pages loop
    
.pages_end:
    inc r15d                   ; Increment b
    jmp .bias_loop             ; Repeat bias loop
    
.bias_end:
    ; Free somma vector
    mov rdi, r14               ; Parameter for dealloc_vector
    call dealloc_vector        ; Call dealloc_vector(somma)
    
    ; Return ret vector
    mov rax, r13               ; Return value = ret vector pointer
    
    ; Restore non-volatile registers
    pop r15
    pop r14
    pop r13
    pop r12
    pop rbx
    
    stop                       ; Macro from sseutils64.nasm - restores stack frame and returns

section .data
    align 8
    one:    dq 1.0             ; Constant 1.0 for calculations
    zero:   dq 0.0             ; Constant 0.0 for initialization