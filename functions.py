import math
import random
import numpy as np
from decimal import *
from fractions import Fraction


# utils
def is_prime(n):
    for i in range(2, int(n/2) + 1):
        if np.mod(n, i) == 0:
            return False
    return True


def gcd(p, q):
    while q != 0:
        p, q = q, np.mod(p, q)
    return p


def find_N(M, m):
    N = 1
    while True:
        if Decimal(M*N % m) == 1:
            break
        N += 1
    return N


def find_X(a, p):
    for i in range(1, p):
        if Decimal(Decimal(math.pow(i, 2)) % p) == a:
            return True
    return False


def primfacs(n):
   i = 2
   primfac = []
   while i * i <= n:
       while n % i == 0:
           primfac.append(int(i))
           n = n / i
       i = i + 1
   if n > 1:
       primfac.append(int(n))
   return primfac


def custom_function(x, n):
    return (math.pow(x, 2) + 1) % n


# functions

def mobius(n):
    num_of_multipliers = 0
    if n == 1:
        return 1
    for i in range(1, n + 1):
        if n % i == 0 and is_prime(i):
            if n % math.pow(i, 2) == 0:
                return 0
            else:
                num_of_multipliers += 1
    if num_of_multipliers % 2 == 0:
        return 1
    return -1


def euler(n):
    num_of_coprime = 0
    for i in range(1, n):
        if gcd(n, i) == 1 or gcd(n, i) == -1:
            num_of_coprime += 1
    return num_of_coprime


def system_solution(system):
    M = 1
    x = 0
    for i in range(0, len(system)):
        M *= system[i][1]
    for i in range(0, len(system)):
        x += (Decimal(M/system[i][1]))*find_N(Decimal(M/system[i][1]), system[i][1])*system[i][0]
    return x%M


def legendre_symbol(a, p):
    a = a % p
    if a == 0:
        return 0
    if find_X(a, p):
        return 1
    return -1


def jacobi_symbol(a, m):
    prime_multipliers = primfacs(m)
    number = 1
    for i in range(0, len(prime_multipliers)):
        number *= legendre_symbol(a, prime_multipliers[i])
    return number


def pollards_rho(n):
    x, y = 2, 2
    d = 1
    while d == 1:
        x = custom_function(x, n)
        y = custom_function(custom_function(y, n), n)
        d = gcd(math.fabs(x - y), n)

    if d == n:
        return -1
    return d


def euclid(a, b):
    if b == 0:
        return a, 1, 0
    else:
        d, xx, yy = euclid(b, a % b)
        x = yy
        y = xx - (a / b) * yy
        return d, x, y


def xab(N, n, alpha, beta, x, a, b):
    value = x % 3
    if value == 0:
        return Decimal(math.pow(x, 2))%N, (a*2)%n, (b*2)%n
    elif value == 1:
        return (x*alpha)%N, (a+1)%n, b
    else:
        return (x*beta)%N, a, (b+1)%n


def logarithm(g, p, beta):
    n = (p - 1)/2
    x, a, b = g*beta, 1, 1
    X, A, B = x, a, b
    for i in range(1, p):
        x, a, b = xab(p, n, g, beta, x, a, b)
        X, A, B = xab(p, n, g, beta, X, A, B)
        X, A, B = xab(p, n, g, beta, X, A, B)
        if x == X:
            res = (euclid(B - b, n)[1] * (a - A)) % n
            return res + n
    return -1


def cipolla_mult(a, b, power, p):
    i = 0
    a_begin = a
    n = 1
    while i < power-2:
        a_prev = a
        a = (a*a_begin+n*b) % p
        n = (a_prev+n*a_begin) % p
        i = i + 1
    return Decimal(a*a_begin+n*b+(a+n*a_begin)*int(math.sqrt(b))) % p


def cipolla(a, p):
    if legendre_symbol(a, p) != 1:
        return -1
    for i in range(1, p+1):
        value = math.pow(i, 2) - a
        if value < 0:
            value = p + value
        if value == 0:
            return i, -i
        if legendre_symbol(value, p) != 1:
            res = cipolla_mult(i, value, (p+1)/2, p)
            return res, p - res
    return -1


def solovay_strassen(n):
    k = N = n - 2
    for i in range(2, n):
        if gcd(n, i) > 1:
            return 0
        if math.pow(i, (n-1)/2) % n != legendre_symbol(i, n) and math.pow(i, (n-1)/2) % n != n+legendre_symbol(i, n):
            return 0
        k = k - 1
        if k == 0:
            return Decimal(1 - math.pow(2, -N))
    return 0


def extended_euklid(a, b):
    s = 0
    old_s = 1
    r = b
    old_r = a

    while r != 0:
        quotient = old_r / r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s

    if b != 0:
        bezout_t = (old_r - old_s * a)/b
    else:
        bezout_t = 0

    return old_s, bezout_t


def n_generation():
    primes = [i for i in range(500, 2500) if is_prime(i)]
    p = Fraction(random.choice(primes))
    q = Fraction(random.choice(primes))
    phi = int((p - 1) * (q - 1))
    e = 0
    for i in range(2, phi):
        if gcd(i, phi) == 1:
            e = i
            break
    if e == 0:
        return -1
    d = 0
    for i in range(2, phi):
        if np.mod((i*e), phi) == 1:
            d = i
            break
    if d == 0:
        return -1
    return int(p*q), e, d


def encryption(m, e, n):
    return np.mod(np.power(Fraction(m), e), n)


def decryption(c, d, n):
    return np.mod(np.power(c, d), n)


def key_generation(p, g):
    x = Fraction(random.randint(2, p))
    y = Fraction(np.mod(np.power(g, Fraction(x)), p))
    return x, y


def elgamal_encryption(m, p, g, y):
    k = Fraction(random.randint(2, p-1))
    a = Fraction(np.mod(np.power(Fraction(g), k), p))
    b = Fraction(np.mod(np.power(Fraction(y), k)*Fraction(m), p))
    return int(a), int(b)


def elgamal_decryption(a, b, x, p):
    return np.mod(np.power(a, p-1-x)*Fraction(b), p)


