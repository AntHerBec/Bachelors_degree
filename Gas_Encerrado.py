###############################################
### 100 particulas con colision entre ellas ###
###############################################

import numpy as np
import numpy.random as ram
import matplotlib.pyplot as plt
from numpy.linalg import norm
import scipy.spatial as spt

"Datos del sistema"

l = 10 # lado de la caja cuadrada
m = 1 # masa del disco

"Datos de la simulación"
trep = 0.01       # tiempo de pausa entre representaciones
inter = 0.01      # intervalo temporal entre pasos
nump = 100  # nº particulas
n = 10000 # nº pasos
nrep = 10 # nº iteraciones hasta volver a pintar
d_choque = 0.25 # distancia de choque

"Condiciones iniciales"

rp = l * ram.rand(nump, 2)
angp, E0p = ram.random(nump) * 2 * np.pi, ram.exponential(0.5, nump)  # 0.5 = Kb*Ti = 0.01*50
direccion = np.array([np.cos(angp), np.sin(angp)])

vmodp = np.sqrt(E0p * 2 / m)
vp = direccion * vmodp
vp = vp.T

Ei = norm(vp, axis=1) ** 2 * m / 2
Emi = sum(Ei) / len(Ei)  # Energía media inicial
Ti = Emi / 0.01

k = nrep

"Funciones de la práctica"

def choque(r, va_inicial, vb_inicial):
    import numpy as np
    from numpy.linalg import norm

    ur = r / norm(r)  # vector unitario en la direccion ra - rb
    un = np.array([-ur[1], ur[0]])  # vector unitario en la direccion normal

    if np.dot(vb_inicial - va_inicial, ur) <= 0:  # comprobacion de que se van acercando

        # Definimos las velocidades de las particulas en las direcciones de ur y un
        var = np.dot(va_inicial, ur)
        van = np.dot(va_inicial, un)

        vbr = np.dot(vb_inicial, ur)
        vbn = np.dot(vb_inicial, un)

        va_final = vbr * ur + van * un  # nuevas velocidades tras choque, contrarias en direccion de ur
        vb_final = var * ur + vbn * un

        return [va_final, vb_final]

    else:

        return [va_inicial, vb_inicial]


def kdtree(rp):
    puntos = spt.cKDTree(rp)  # método que permite trabajar con los puntos en conjunto
    pares = puntos.query_pairs(d_choque)  # query_pairs encuentra las parejas de puntos a una distancia máxima d_choque

    return (pares)

"Simulación de las partículas"

fig1 = plt.figure(figsize=(8, 6), clear=True)

ax1 = fig1.add_subplot(111)  # Como si fuera 1 unico plt, facilmente modificable en el futuro
ax1.set_title('{} partículas colisionando en una caja'.format(nump))
ax1.set_xlabel('X')
ax1.set_ylabel('Y') # Títulos de ax1
ax1.set_xlim(0, l)
ax1.set_ylim(0, l)  # Limites de ax1

particula, = ax1.plot(rp[:, 0], rp[:, 1], 'bo', markersize=10) # Pintamos la partícula de manera que sea fácilmente animable como objeto 2D

pres_list = [] # Para el siguiente apartado, recoge los valores de presión en cada iteración

for i in range(n):

    presion = 0 # Para acumular la presión total en la iteración

    rp = rp + vp * inter # Actualizamos los valores de rp

    pares = kdtree(rp)

    for a in range(nump): # Rebotes en las paredes

        if (rp[a, 0] <= 0 and vp[a, 0] < 0) or (rp[a, 0] >= l and vp[a, 0] > 0):
            vp[a, 0] = - vp[a, 0]
            presion += m * abs(vp[a, 0]) / (2 * l * inter)  # ecuacion para P = F/(4*l) = m*dv/(4*l*dt) con dt
            # = intervalo entre iteraciones y dv = 2*|v_pared|

        if (rp[a, 1] <= 0 and vp[a, 1] < 0) or (rp[a, 1] >= l and vp[a, 1] > 0):
            vp[a, 1] = - vp[a, 1]
            presion += m * abs(vp[a, 1]) / (2 * l * inter)

    pres_list.append(presion)

    for par in pares: # Choques entre partículas

        rba = rp[par[1]] - rp[par[0]]
        vp[par[0]], vp[par[1]] = choque(rba, vp[par[0]], vp[par[1]])

    if k == nrep: # Volvemos a pintar las partículas cada nrep iteraciones

        plt.pause(trep)
        particula.set_data(rp[:, 0], rp[:, 1])

        k = 0

    k += 1

"Análisis físico"

T = 50  # temperatura inicial en K
Kb = 0.01  # cte de Boltzmann

E = norm(vp, axis=1) ** 2 * m / 2

Emf = sum(E) / len(E)  # Energía media final
Tf = Emf / Kb # Temperatura final

fig2 = plt.figure(figsize=(14, 6), clear=True) # Gráfica con histogramas de Ei y Ef

ax2_1 = fig2.add_subplot(121)
plt.xlabel('E inicial')
plt.ylabel('número de partículas')
plt.title('Energía cinética inicial de {} partículas'.format(nump))
ax2_1.hist(Ei, bins=25, ec='black')

ax2_2 = fig2.add_subplot(122)
plt.xlabel('E final')
plt.ylabel('número de partículas')
plt.title('Energía cinética final de {} partículas'.format(nump))
ax2_2.hist(E, bins=25, ec='black')

fig3 = plt.figure(figsize=(14, 6), clear=True) # Gráfica con valores de P
ax3 = fig3.add_subplot(111)

# para hacer una representación útil de los valores de la presión, tomamos el promedio de esta cada 100 valores:
num_promedio = 100
num_intervalos = int(n / num_promedio)
t_intervalos = np.linspace(0, int(n * inter / nrep), num_intervalos)  # tiempo que (idealmente, sin delay al hacer los cálculos) dura la simulación dividido entre el número de intervalos.

pres_promedio = []

for i in range(0, num_intervalos):  # para acumular los valores promedio en una lista

    pres_promedio.append(sum(pres_list[(i * num_promedio):((i + 1) * num_promedio)]) / len(
        pres_list[(i * num_promedio):((i + 1) * num_promedio)]))

ax3.plot(t_intervalos, pres_promedio)
plt.title('promedios (de {} valores) de presión a lo largo del tiempo'.format(num_promedio))
plt.xlabel('t')
plt.ylabel('P (promedios)')
plt.show()

P = sum(pres_list) / len(pres_list) # Presión media

V = l ** 2 - nump * 0.125 ** 2 * np.pi  # Consideramos como volumen (2-dimensional) el área de la caja menos el área conjunta de las partículas. Como chocan a 0.25 unidades, su radio es 0.125 unidades

relacion = P * V / (Kb * Tf * nump)  # Expresión que, si el gas es ideal, es igual a 1


print('\n','', 26 * '-', '\n', 'Valores numéricos de T y P', '\n', 26 * '-','\n') # Datos numéricos: T se conserva, y P es tal que PV = NkbT

print('+ La temperatura se conserva:')
print('    La temperatura inicial es:    {:.4f}'.format(Ti))
print('    La temperatura final es:      {:.4f}'.format(Tf),'\n')
print('+ La presión media total es:    {:.3f}'.format(P),'\n')
print('+ La relación de los gases ideales (PV = NKbT) se cumple:     PV/(NKbT) = {:.3f}'.format(relacion))
