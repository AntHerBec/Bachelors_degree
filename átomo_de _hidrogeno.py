# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 10:11:11 2025

@author: UO295695
"""

import numpy as np
import matplotlib.pyplot as plt
import time as t

"Definición de constantes"

c = 3e8
me = 0.511e6 / c**2
hbarra = 6.582e-16
a0 = 5.292e-11
V0 = 244

a = 10**(-10)

k = 2*me*V0*a**2/hbarra**2

L = [0, 1]
Z = 2

E0 = hbarra**2 / (2 * me * a0**2)
print(E0)

# Partimos desde el infinito hasta el origen

r = np.linspace(500, 1e-20, 3001) # todas las ecuaciones teniendo en cuenta que estamos usando una posición
                                 # en escala del radio de bohr, es decir, unidades de a0

deltar = abs((r[1]-r[0]))

def V(r, Z = 1, litio = False):

    '''
    El potencial se simplifica porque todas las constantes (excepto el 2) que aparecen son igual a
    a0, de manera que al expresar V en unidades de a0 tenemos 1 unidad de a0
    '''

    v = 2 * Z / r

    if litio == True:

        v = 2 / r

        for i in range(len(x)):

            if r[i] <= 1: v[i] = v[i] * 3 - 4

    return v

def psi_atomo(r, E, L = 0, Z = 1, litio = False):

    psi, phi = np.zeros(len(r)), np.zeros(len(r))

    phi[0] = 0
    phi[1] = 1e-5 # número pequeño cualquiera, parametro a modificar para ver si mejora

    func = L * (L+1) / r**2 - V(r, Z = Z, litio = litio) - E


    for i in range(2, len(r)):

        f_iter = func[i-1]
        phi[i] = phi[i-1] * (2 + deltar**2 * f_iter / (1 - deltar**2 * f_iter / 12)) - phi[i-2]

    psi = phi / (1 - deltar**2 * func / 12)

    return psi


def shooting_method_noprints (intervalos, x, L = 0, Z = 1, Etol = 1e-10, psitol = 1e-40, niter = 1000, paridad = 'par', litio = False):

    lista = []

    for i in range(1, len(intervalos)):

        inf, sup = intervalos[i-1], intervalos[i]

        psi_inf = psi_atomo(x, inf, L = L, Z = Z, litio = litio)[-1]
        psi_sup = psi_atomo(x, sup, L = L, Z = Z, litio = litio)[-1]

        prod_psi = psi_inf * psi_sup

        if prod_psi < 0:

            print('se puede aplicar bisección')

            for n in range(niter + 1):

                pmedio = (inf + sup)/2

                if abs(inf - sup) <= Etol:

                    print('se alcanzo la tolerancia de E')
                    print('diferencia entre inf y sup: {}'.format(abs(inf - sup)))
                    break

                psi_medio = psi_atomo(x, pmedio, L = L, Z = Z, litio = litio)[-1]
                psi_inf = psi_atomo(x, inf, L = L, Z = Z, litio = litio)[-1]

                if abs(psi_medio) <= psitol:

                    print('se alcanzo la tolerancia de psi')
                    break

                if psi_medio * psi_inf > 0: sup = pmedio

                else: inf = pmedio

            lista.append(pmedio)

    return lista
    
        
'''
E = - 13.6 * a0**2 * 2 * me / hbarra**2
    
psi, phi = psi_H(r, L_H, E)

plt.plot(r, psi)
plt.show()
'''
    
    
multiplos = np.linspace(E0*1.05, 0, 1001) * (-1)
    
multiplos1 = np.linspace(E0*1.05, E0*0.75, 1001) * (-1)
multiplos2 = np.linspace(E0*0.35, E0*0.2, 501) * (-1)
multiplos3 = np.linspace(E0*0.15, E0*0.1, 501) * (-1)

tlist = []
'''
for i in range(len(multiplos)):

    for j in range(i, len(multiplos)):

        t0 = t.time()
        E = shooting_method_noprints([multiplos[i], multiplos[j]], r)

        if E != None:

            for n in lista:

                if abs(n - E) <= 10 ** (-2): contador += 1

            if contador == 0:

                lista.append(float(f"{E:.13f}"))

            else:

                pass

        t1 = t.time()
        tlist.append(t1-t0)
        t0 = t1

        contador = 0
'''
'''
t0 = t.time()
E = shooting_method_noprints(multiplos1, r)
E.append(*shooting_method_noprints(multiplos2, r))
E.append(*shooting_method_noprints(multiplos3, r))
'''
for j in L:

    E = shooting_method_noprints(multiplos, r, L = j, Z = Z)

    lista_ev = [i for i in E] # para que los decimales no se redondeen

    print(lista_ev)

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.set_xlabel(r"radio (unidades de $a_0$)")
    ax1.set_ylabel(r"$\psi$")
    ax2.set_xlabel(r"radio (unidades de $a_0$)")
    ax2.set_ylabel(r"$\psi$")

    plt.tight_layout()

    for i in E:

        psi = psi_atomo(r, i, L = j, Z = Z)
        psi_norm = psi / max(psi)
        ax1.plot(r, psi_norm, label = "{:.2f} eV".format(i))
        ax2.plot(r, psi_norm, label = "{:.2f} eV".format(i))
        ax2.set_xlim(0,60)

    plt.legend()
    plt.show()

for j in L_li:

    E = shooting_method_noprints(multiplos, r, L = j, Z = Z, litio = False)

    lista_ev = [i for i in E] # para que los decimales no se redondeen

    print(lista_ev)

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.set_xlabel(r"radio (unidades de $a_0$)")
    ax1.set_ylabel(r"$\psi$")
    ax2.set_xlabel(r"radio (unidades de $a_0$)")
    ax2.set_ylabel(r"$\psi$")

    plt.tight_layout()

    for i in E:

        psi = psi_atomo(r, i, L = j, Z = Z, litio = True)
        psi_norm = psi / max(psi)
        ax1.plot(r, psi_norm, label = "{:.2f} eV".format(i))
        ax2.plot(r, psi_norm, label = "{:.2f} eV".format(i))
        ax2.set_xlim(0,60)

    plt.legend()
    plt.show()

'''
valores = np.array([-0.67971407, -1.35942813, -6.79714067, -7.47685473, -2.71885627, -3.39857033, -4.0782844, -4.75799847, -5.43771253, -6.1174266, -8.83628286, -9.51599693, -10.195711, -10.87542506, -11.55513913, -12.2348532, -12.91456726, -13.59428133])
valores2 = np.array([-1.0195711, -7.1369977, -3.0587133, -4.0782844, -5.0978555, -6.1174266, -9.1761399 -10.195711, -11.2152821, -12.2348532, -13.2544243])

for i in range(len(valores2)):

    print(psi_H(r, L_H, valores[i])[0][-1])
'''
'''
    plt.plot(r, psi_H(r, L_H, valores[i])[0])
    plt.show()
'''
