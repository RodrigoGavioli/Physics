import numpy as np
import math
import matplotlib.pyplot as pl
from fontTools.varLib.plot import stops

G = 6.6738*10**(-11) #m^-3 kg^-1 s^-2
M = 1.9891*10**(30) #kg


v1 =30300 # m/s
l1 =147000000000 #m
def resolver_bhaskara(G,M,v1,l1):
    a = 1
    b = -((2*G*M)/(v1*l1))
    c =-(v1**2-(2*G*M)/(l1))
    raiz = b**2 - 4*a*c
    if raiz < 0:
        return("A raiz é negativa")
    v2 = (-b - np.sqrt(raiz)) / (2*a)
    return v2
#print("Este é o valor da velocidade no afelio:",resolver_bhaskara(G,M,v1,l1))
v2 = resolver_bhaskara(G, M, v1, l1)

def raio2(v1,l1,v2):
    l2 = (v1*l1)/v2
    return l2

l2 = raio2(v1,l1,v2)


a = (l1+l2)/2
b = np.sqrt(l1*l2)

def periodo_orbital(a,b,l1,v1):
    T = (2*np.pi*a*b)/(l1*v1)
    return T

#print("Em anos o periodo é:",periodo_orbital(a, b, l1, v1)/(60*60*24*365))

#Para fazer ao halley basta substituir os valores de v1 e l1, pode tambem deixar aberto
#como input sendo v1 = float("insira:"input()) e l1 = float("insira:"input())

##################################################################################################

#Exercicio 2

def catalan(n):
    C = ((4*n+2)/(n+2))
    return C

Serie = [1.0]
CatalanN = 1
for i in range(0,1000001):
    CatalanN = catalan(i)*CatalanN
    Serie.append(CatalanN)



#############################################################################################

#Exercicio 3


def madelung_constant(L):
    M = 0.0
    for i in range(-L, L + 1):
        for j in range(-L, L + 1):
            for k in range(-L, L + 1):
                if i == 0 and j == 0 and k == 0:
                    continue  # Ignorar a posição da origem

                distance = np.sqrt(i ** 2 + j ** 2 + k ** 2)
                sign = (-1) ** (i + j + k)
                M += sign / distance

    return M


# Definindo o tamanho da rede tridimensional
L = 10  # Pode aumentar para maior precisão, mas também aumenta o tempo de execução

M = madelung_constant(L)




########################################################################################################

#Exercicio 4
pascalbin = []


def binomio(valorn,valork):
    n = 1
    for i in range(1, valorn + 1):
        n = n * i
    k = 1
    for j in range(1, valork + 1):
        k = k * j
    dif = 1
    d = (valorn - valork)
    for l in range(1, d + 1):
        dif = dif * l
    if valork == 0:
        return 1
    else:
        bi = (n)/(k*dif)
        return bi
#Parte 1 feita

for i in range(21):
    valorn = i
    for j in range(valorn + 2):
        valork = j
        if binomio(valorn,valork) == int(binomio(valorn,valork)):
            pascalbin.append(binomio(valorn, valork))
#Parte 2 feita
valorn = 100
valork = 60

probabilidade = binomio(valorn,valork)/(2**(valorn))


#print("esta é a probabilidade:",probabilidade)
#Parte 3 feita



###########################################################################################

#Exercicio 5
ListaPrimos = [2]
n = 100000

for i in range(n): #i is the number we are going to check if it is prime
    for j in range(len(ListaPrimos)):
        if len(ListaPrimos) <= np.sqrt(i):
            if i % ListaPrimos[j] == 0:
                print(i,"is not prime")
                break
            else:
                print(i,"is prime")
                ListaPrimos.append(i)
                print("lista atualizada com:",i,ListaPrimos)
                break
        else:
            print(i,"is not prime")
            break
       
print(ListaPrimos)