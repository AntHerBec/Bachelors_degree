
import numpy as np
import matplotlib.pyplot as plt

"Datos para simulación"

n = 10000
nrep = 10
interv = 0.01

"Datos para mallado"

puntos_mallado = 1000

c = 3e8
dxp = 400e-9
dtp = dxp / c
dx = 10e-9 # m
dt = dx / (2 * c) # s

malla = np.arange(puntos_mallado) * dx * 1e6  # en micras

"Datos para pulso"

E0 = 1
t0p = 5 * dtp  # (cuando aparece el pulso)
xp = 250  # (donde aparece el pulso)



E = np.zeros(puntos_mallado)
H = np.zeros(puntos_mallado)

t = 0.0  # t inicial
k = nrep

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)  # Como si fuera 1 unico plt, facilmente modificable en el futuro

ax.set_xlim(0, puntos_mallado * dx * 1e6) # Limites de ax
ax.set_ylim(-3, 3)

ax.plot(0.5 * puntos_mallado * dx * 1e6 * np.ones(200), np.linspace(-12, 12, 200), 'black') # Franja negra para separar el medio según su epsilon relativa

campoelec, = ax.plot(malla, E, label = r'E ($\frac{V}{m}$)') # Campos iniciales
campomagn, = ax.plot(malla, H, label = r'H ($T$)')

ax.set_title("Propagación 1D de una onda EM") # Detalles de la gráfica
ax.set_xlabel(r'Posición ($\mu m$)')
ax.set_ylabel(r'Amplitud de la onda')
ax.legend()

"Conductividad y permitividad"

epsilon0, sigma1, sigma2 = 8.85 * 1e-12, 0, 4000
epsilon_r1, epsilon_r2 = 1, 4

epsilonr = [epsilon_r1 for i in range(int(puntos_mallado / 2))] # Determinamos qué zonas del medio tienen cada parámetro
sigma = [sigma1 for i in range(int(puntos_mallado / 2))]

for i in range(int(puntos_mallado / 2)): epsilonr.append(epsilon_r2)
for i in range(int(puntos_mallado / 2)): sigma.append(sigma2)

epsilonr = np.array(epsilonr) # Los pasamos a arrays para poder operar con ellos
sigma = np.array(sigma)

"Variables para evolución temporal de E"

cd = 1 / (2 * epsilonr)
ctc = sigma * dt / (2 * epsilonr * epsilon0)
ca = 1 - ctc / (1 + ctc)
cb = 1 / ((1 + ctc) * 2 * epsilonr)

tuberia1 = np.zeros(np.round(2 * np.sqrt(epsilon_r1)).astype(int)) # creación de tuberías según las epsilon relativas
tuberia2 = np.zeros(np.round(2 * np.sqrt(epsilon_r2)).astype(int))

"Bucle para animar la evolución temporal de la onda"

for j in range(n):

    t = t + dt # actualizamos el valor del tiempo

    E[1:] = ca[1:] * E[1:] - cb[1:] * (H[1:] - H[:-1]) # Ecuación de evol. temporal de E

    E_pulso = E0 * np.exp( -0.5 * ((t - t0p) / dtp) ** 2) # Pulso

    E[xp] = E[xp] + E_pulso

    E[0] = tuberia1[0] # Absorción en la pared izquierda
    tuberia1[:-1] = tuberia1[1:]
    tuberia1[-1] = E[1]

    H[:-1] = H[:-1] - 0.5 * (E[1:] - E[:-1]) # Ecuación de evol. temporal de H

    H[-1] = tuberia2[-1] # Absorción en la pared derecha
    tuberia2[1:] = tuberia2[:-1]
    tuberia2[0] = H[-2]

    if k == nrep: # Actualizamos gráfico cada nrep iteraciones

        plt.pause(interv)
        campoelec.set_data(malla, E)
        campomagn.set_data(malla, H)

        k = 0

    k += 1