# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 09:05:55 2025

@author: UO295695
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as r
import time as t
import scipy.special as spe

def generacion_poisson(nu, n):
    
    k = np.ones(n)
    A = np.ones(n)
    u = r.random(n)
    
    A = u * A
    exp = np.exp(- nu)
    
    for i in range(n):
    
        while A[i] >= exp:
            
            k[i] += 1
            u = r.random()
            A[i] = u * A[i]
    
    return k - np.ones(n)

def poisson_func(nu, x):
    return np.exp(-nu) * nu ** x / spe.gamma(x+1)

def promedio(list):
    return sum(list) / len(list)

def fotomultiplicador(nu, elec_inicial, numero_dinodos, reps):

    z = 0

    elec_final_simulado = []
    elec_final_esperado = []
    eficiencia_reps_simulado = []
    eficiencia_reps_esperado = []

    while z < reps:

        numero_elec_simulado = elec_inicial
        numero_elec_esperado = elec_inicial

        eficiencias_simulado = []
        eficiencias_esperado = []

        for j in range(numero_dinodos):

            anterior_simulado = numero_elec_simulado
            anterior_esperado = numero_elec_esperado
            numero_elec_simulado = int(sum(generacion_poisson(nu[j], numero_elec_simulado)))
            numero_elec_esperado = int(sum(r.poisson(nu[j], numero_elec_esperado)))

            if anterior_esperado != 0:

                if anterior_simulado != 0:

                    eficiencia_simulado = numero_elec_simulado / anterior_simulado * 100
                    eficiencias_simulado.append(eficiencia_simulado)
                    eficiencia_esperado = numero_elec_esperado / anterior_esperado * 100
                    eficiencias_esperado.append(eficiencia_esperado)

                elif anterior_simulado == 0:

                    eficiencias_simulado.append(0)
                    eficiencia_esperado = numero_elec_esperado / anterior_esperado * 100
                    eficiencias_esperado.append(eficiencia_esperado)

            elif anterior_esperado == 0:

                if anterior_simulado != 0:

                    eficiencia_simulado = numero_elec_simulado / anterior_simulado * 100
                    eficiencias_simulado.append(eficiencia_simulado)
                    eficiencias_esperado.append(0)

                elif anterior_simulado == 0:

                    eficiencias_simulado.append(0)
                    eficiencias_esperado.append(0)


        eficiencias_simulado = np.array(eficiencias_simulado)
        eficiencias_esperado = np.array(eficiencias_esperado)

        elec_final_simulado.append(numero_elec_simulado)
        eficiencia_reps_simulado.append(eficiencias_simulado)
        elec_final_esperado.append(numero_elec_esperado)
        eficiencia_reps_esperado.append(eficiencias_esperado)

        z += 1

    eficiencia_dinodos_simulado = promedio(eficiencia_reps_simulado)
    eficiencia_dinodos_esperado = promedio(eficiencia_reps_esperado)

    return [elec_final_simulado, elec_final_esperado, eficiencia_dinodos_simulado, eficiencia_dinodos_esperado]


"Comprobación del método"

nu = [1, 5, 10]
size = [100, 1000, 10000, 100000]

for i in nu:
    for j in size:
        valores_poisson = generacion_poisson(i, j)
        valores_esperados = r.poisson(i, j)

        x = np.linspace(0, max(valores_esperados)+2, 100)
        plt.hist(valores_poisson, label = "generador", bins = np.arange(0,max(valores_poisson)+2)-0.5, density=True, color = 'red', alpha = 0.5)
        plt.hist(valores_esperados, label = "random.poisson", bins = np.arange(0,max(valores_esperados)+2)-0.5, density = True, histtype='step', linestyle = '--', linewidth = 2, color = 'black')
        plt.plot(x, poisson_func(i, x), label = 'curva teórica', color = 'green')

        plt.title('Comparación entre generador y valores esperados ({} valores, media = {})'.format(j, i), fontsize='large')
        plt.xlabel('Probabilidad de que el valor generado sea x', fontsize='large')
        plt.ylabel('x, valores posibles', fontsize='large')
        plt.legend(fontsize='large')
        plt.show()


elec_0 = 1
n_din = 6
nreps = 1000

nu = 5 * np.ones(n_din)
rend_total = 1
for i in nu: rend_total *= i

elec_sim, elec_esp, efic_sim, efic_esp = fotomultiplicador(nu, elec_0, n_din, nreps)


plt.hist(elec_sim, bins = 50, alpha = 0.5, label = 'simulación (generador)')
plt.hist(elec_esp, bins = 50, histtype='step', linewidth = 1, color = 'black', label = 'esperado (random.poisson)')

plt.title(r'Número final de electrones (inicialmente {} $e^-$) para {} repeticiones y {} dínodos de $\nu$=5'.format(elec_0, nreps, n_din), fontsize='large')
plt.xlabel(r'nº $e^-$ final', fontsize='large')
plt.ylabel('Frecuencia', fontsize='large')
plt.legend(fontsize='large')
plt.show()

print(r'Rendimiento simulado en cada dínodo = {:.0f}%, {:.0f}%, {:.0f}%, {:.0f}%, {:.0f}% y {:.0f}%, respectivamente'.format(
    efic_sim[0], efic_sim[1], efic_sim[2], efic_sim[3], efic_sim[4], efic_sim[5]))
print(r'Rendimiento con random.poisson en cada dínodo = {:.0f}%, {:.0f}%, {:.0f}%, {:.0f}%, {:.0f}% y {:.0f}%, respectivamente'.format(
    efic_esp[0], efic_esp[1], efic_esp[2], efic_esp[3], efic_esp[4], efic_esp[5]))
print(r'Rendimiento total simulado = {:.0f}; Rendimiento total con random.poisson = {:.0f}. Rendimiento total esperado = {:.0f}.  Errores Relativos = {:.2f}% (simulado), {:.2f}% (random.poisson)'.format(
    promedio(elec_sim) / elec_0, promedio(elec_esp) / elec_0, rend_total, abs(promedio(elec_sim) - rend_total) / rend_total * 100, abs(promedio(elec_esp) - rend_total) / rend_total * 100), 2*'\n')


"optimizacion de la posicion del dinodo nu = 8"

elec_0 = 1
n_din = 6
nreps = 500

for l in range(n_din):

    nu = 5 * np.ones(n_din)
    nu[l] = 8

    rend_total = 1
    for i in nu: rend_total *= i

    elec_sim, elec_esp, efic_sim, efic_esp = fotomultiplicador(nu, elec_0, n_din, nreps)

    plt.hist(elec_sim, bins=50, alpha=0.5, label='simulación (generador)')
    plt.hist(elec_esp, bins=50, histtype='step', linewidth=1, color='black', label='esperado (random.poisson)')
    plt.show()

    print(
        r'Rendimiento simulado en cada dínodo = {:.0f}%, {:.0f}%, {:.0f}%, {:.0f}%, {:.0f}% y {:.0f}%, respectivamente'.format(
            efic_sim[0], efic_sim[1], efic_sim[2], efic_sim[3], efic_sim[4], efic_sim[5]))
    print(
        r'Rendimiento con random.poisson en cada dínodo = {:.0f}%, {:.0f}%, {:.0f}%, {:.0f}%, {:.0f}% y {:.0f}%, respectivamente'.format(
            efic_esp[0], efic_esp[1], efic_esp[2], efic_esp[3], efic_esp[4], efic_esp[5]))
    print(
        r'Rendimiento total simulado = {:.0f}; Rendimiento total con random.poisson = {:.0f}. Rendimiento total esperado = {:.0f}.  Errores Relativos = {:.2f}% (simulado), {:.2f}% (random.poisson)'.format(
            promedio(elec_sim) / elec_0, promedio(elec_esp) / elec_0, rend_total, abs(promedio(elec_sim) - rend_total) / rend_total * 100, abs(promedio(elec_esp) - rend_total) / rend_total * 100), 2*'\n')


"Eficiencia del fotomultiplicador"

elec_0 = 1
n_din = 6
nreps = [10, 32, 100, 316, 1000] # separaciones de aproximadamente 1/2 en log base 10
umbral = [1, 10000]

nu = 5 * np.ones(n_din)
nu2 = 5 * np.ones(n_din)
nu2[0] = 8

efi = []
efi2 = []

for i in nreps:

    elec_sim, elec_esp, efic_sim, efic_esp = fotomultiplicador(nu, elec_0, n_din, i)
    elec2_sim, elec2_esp, efic2_sim, efic2_esp = fotomultiplicador(nu2, elec_0, n_din, i)

    contador = i
    contador2 = i

    for j in umbral:

        for k in elec_sim:

            if k < j: contador = contador - 1

        for k in elec2_sim:

            if k < j: contador2 = contador2 - 1

        print(r'La eficiencia del fotomultiplicador (primer dínodo con nu = 5) para {} repeticiones y señal umbral = {} electrones es {:.1f}%'.format(i, j, contador / i * 100), '\n')
        print(r'La eficiencia del fotomultiplicador (primer dínodo con nu = 8) para {} repeticiones y señal umbral = {} electrones es {:.1f}%'.format(i, j, contador2 / i * 100), '\n')

        efi.append(contador / i * 100)
        efi2.append(contador2 / i * 100)

plt.semilogx(nreps, efi[0::len(umbral)], label = 'umbral = {}, nu[0] = {:.0f}'.format(umbral[0], nu[0]))
plt.semilogx(nreps, efi[1::len(umbral)], label = 'umbral = {}, nu[0] = {:.0f}'.format(umbral[1], nu[0]))
plt.semilogx(nreps, efi2[0::len(umbral)], label = 'umbral = {}, nu[0] = {:.0f}'.format(umbral[0], nu2[0]))
plt.semilogx(nreps, efi2[1::len(umbral)], label = 'umbral = {}, nu[0] = {:.0f}'.format(umbral[1], nu2[0]))

plt.title('Eficiencias del fotomultiplicador')
plt.xlabel('número de repeticiones')
plt.ylabel('eficiencia (%)')
plt.legend()
plt.show()















