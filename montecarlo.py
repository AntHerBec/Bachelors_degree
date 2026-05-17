# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 10:00:59 2025

@author: UO295695
"""

"Definimos condiciones del problema"

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as r
import scipy.optimize as opt

N0 = 100 # nº nucleos de Bi en el instante inicial
t = np.linspace(0, 50,101) # mallado de tiempo en dias
t_aleat = np.linspace(0, 50, 101) # mallado para generacion de numeros aleatorios
deltat = t[1] - t[0]

tau_bi = 7.8 # vida media del Bi en dias
tau_po = 196
omega_bi = 1 / tau_bi
omega_po = 1 / tau_po

N_bi = [N0] # nucleos de Bi a lo largo del tiempo
N_po = [0]

def prob_bi (t):
    return np.exp(-t * omega_bi) * omega_bi

def prob_po (t):
    return np.exp(-t * omega_po) * omega_po

def bi_real (t):
    return N0 * np.exp(-omega_bi * t)

def po_real (t):
    return N0 * omega_bi / (omega_po - omega_bi) * (-np.exp(-omega_po * t) + np.exp(-omega_bi * t))

def produccion_alphas (t):
    return -N0 * omega_bi * np.exp(-omega_bi * t) + N0 * omega_bi / (omega_po - omega_bi) * (omega_po * np.exp(-omega_po * t) - omega_bi * np.exp(-omega_bi * t))

def ajuste_exp(a, b, x):
    return b * np.exp(- a * x)

def ajuste_expdoble(a, b, c, x):
    return ((1 / a) / ((1 / b) - (1 / a))) * c * (np.exp(-x / a) - np.exp(-x / b))



bi_max = max(prob_bi(t))
po_max = max(prob_po(t))

'''
for i in range(1, len(t)):

    nbi = N_bi[-1]
    num_aleat = r.random(nbi)
    pbi = prob_bi[i]
    contador = 0

    for j in num_aleat:
        if j <= pbi: contador += 1

    N_bi.append(nbi - contador)

plt.scatter(t, N_bi, s = 1)

plt.show()
'''

total_aceptados = np.zeros(len(t_aleat))
total_aceptados2 = np.zeros(len(t_aleat))
total_aceptados3 = np.zeros(len(t_aleat))
num_muestras = 1000

for k in range(num_muestras):

    aceptados = np.zeros(len(t_aleat))
    aceptados2 = np.zeros(len(t_aleat))
    aceptados3 = np.zeros(len(t_aleat))


    u = r.random(len(t_aleat))
    x = u * (t_aleat[-1] - t_aleat[0])

    x.sort()

    u2 = r.random(len(t_aleat))
    x2 = u2 * (t_aleat[-1] - t_aleat[0])

    x2.sort()

    u3 = r.random(len(t_aleat))
    x3 = u3 * (t_aleat[-1] - t_aleat[0])

    x3.sort()

    w = prob_bi(x) / bi_max
    w2 = prob_po(x2) / po_max - prob_bi(x2) / bi_max
    w3 = prob_po(x3) / po_max

    for n in range(N0):

        y = r.random(len(u))
        y2 = r.random(len(u))
        y3 = r.random(len(u))

        for i in range(len(u)):

            if w[i] <= y[i]: aceptados[i] += 1
            if w2[i] <= y2[i]: aceptados2[i] += 1
            if w3[i] <= y3[i]: aceptados3[i] += 1


    total_aceptados += aceptados
    total_aceptados2 += aceptados2
    total_aceptados3 += aceptados3

aceptados_list = total_aceptados / num_muestras
aceptados2_list = total_aceptados2 / num_muestras
aceptados3_list = total_aceptados3 /num_muestras

numbi = N0 - aceptados_list
plt.scatter(t_aleat, numbi, s = 1)
plt.scatter(t_aleat, N0 - aceptados2_list, s = 1)
plt.scatter(t_aleat, aceptados3_list, s = 1)


plt.show()


plt.scatter(t_aleat, N0 - aceptados_list, s = 1)
plt.plot(t, bi_real(t), color = 'red')

plt.show()

plt.scatter(t_aleat, N0 - aceptados2_list, s = 1)
plt.plot(t, po_real(t), color = 'red')

plt.show()

plt.scatter(t_aleat, aceptados3_list, s = 1)
plt.plot(t, N0 - bi_real(t) - po_real(t), color = 'red')

plt.show()


"maxima produccion de alphas"

print(t[np.argmax(aceptados_list - aceptados2_list)])




"ajuste"

ajuste = opt.curve_fit(ajuste_exp, t, numbi)
print(ajuste[0])

func = ajuste_exp(ajuste[0][1], ajuste[0][0], t)

plt.plot(t, func)
plt.plot(t, bi_real(t))
plt.show()




'''
max_alfas = []

for j in range(100,110):

    p = np.poly1d(coefs)
    func = p(t)
    der_func = np.gradient(func, t)

    for i in range(len(t)):

        if der_func[i] == min(der_func) and t[i] < 2000:

            max_alfas.append(t[i])

print(min(max_alfas), max(max_alfas))
print(sum(max_alfas) / len(max_alfas))

der_real = np.gradient(po_real(t) - bi_real(t), t)

for i in range(len(t)):

    if der_real[i] == min(der_real):

        print(t[i])
'''
'''
plt.plot(t, - bi_real(t) + po_real(t))
plt.plot(t, (func - func2) * N0)
plt.show()
'''

plt.plot(t, - bi_real(t) + po_real(t), color = 'red')
plt.scatter(t_aleat, aceptados_list - aceptados2_list, s=1)
plt.show()

