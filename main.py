from functions import (mobius, euler, system_solution, legendre_symbol, jacobi_symbol,
                       pollards_rho, logarithm, cipolla, solovay_strassen,
                       n_generation, encryption, decryption, key_generation,
                       elgamal_encryption, elgamal_decryption)


print(f"Mobius result: {mobius(100)}")
print(f"Euler result: {euler(100)}")

system_of_linear_comp_eq = [(2, 3), (3, 5), (7, 11)]
print(f"System solution: {system_solution(system_of_linear_comp_eq)}")


print(f"Legendre symbol: {legendre_symbol(219, 383)}")
print(f"Jacobi symbol: {jacobi_symbol(219, 383)}")

print(f"Pollards_rho factorisation: {pollards_rho(8051)}")

print(f"Pollards_rho algorithm for logarithms: {logarithm(2, 5, 3)}")

print(f"Cipolla's algorithm for sqrt: {cipolla(10, 13)}")

print(f"Solovay–Strassen primality test: {solovay_strassen(11)}")

print("RSA")
message = 10235
n, e, d = n_generation()
enc = encryption(message, e, n)
print(f"Encrypted value: {enc}")
dec = decryption(enc, d, n)
print(f"Decrypted value: {dec}")

print("ElGamal encryption")
message_elgamal = 217
p = 4397
g = 2
x, y = key_generation(p, g)
a, b = elgamal_encryption(message_elgamal, p, g, y)
print(f"Encrypted value: {a, b}")
elg_dec = elgamal_decryption(a, b, x, p)
print(f"Decrypted value: {elg_dec}")




