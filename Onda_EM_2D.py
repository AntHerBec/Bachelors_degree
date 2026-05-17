########################################################
### Onda EM en 2D con medio dieléctrico y absorbente ###
########################################################

import numpy as np
import matplotlib.pyplot as plt

"Datos para simulación"

n = 10000
nrep = 10
interv = 0.01

"Datos para mallado"

puntos_mallado = 401

c = 3e8
dxp = 40e-9
dtp = dxp / c

"Datos para pulso"

E0 = 1
xp = (100, 200)
t0p = 5 * dtp
dx = 10e-9 # m
dt = dx / (2 * c) # s

t = 0.0  # t inicial
k = nrep

fig = plt.figure(figsize=(7, 8))
ax1 = fig.add_subplot(111)  # Como si fuera 1 unico plt, facilmente modificable en el futuro

ax1.set_title("Propagación en 2D de una onda EM en un medio dieléctrico, conductor y absorbente")
ax1.set_xlabel(r'Posición ($\mu m$)')
ax1.set_ylabel(r'Amplitud de la onda ($\frac{V}{m}$)')

malla = np.arange(puntos_mallado) * dx * 1e6 # en micras

Y, X = np.meshgrid(malla, malla) #Creamos el mallado bidimensional

Ez = np.zeros((401, 401)) # Definimos los campos eléctrico y magnético
Hx = np.zeros((401, 401))
Hy = np.zeros((401, 401))

# Hy = H[a, :], Hx = H[:, a]

nivelesdecalor = np.linspace(-0.1, 0.1, 21) # Número de tonos en los que se divide la barra de color

cs1 = ax1.contourf(X, Y, np.clip(Ez, -0.1, 0.1), nivelesdecalor, cmap='viridis') # Comando para pintar el mapa de calor

bar1 = plt.colorbar(cs1, label = r'$\frac{V}{m}$', orientation = 'horizontal') # Barra de color, indica el valor del campo en la zona con dicho color


"Conductividad y permitividad"

epsilon0, sigma1, sigma2 = 8.85 * 1e-12, 0, 4000
epsilon_r1, epsilon_r2 = 1, 4

epsilonr = [epsilon_r1 for i in range(int(puntos_mallado / 2))] # Determinamos qué zonas del medio tienen cada parámetro
sigma = [sigma1 for i in range(int(puntos_mallado / 2))]

for i in range(int(puntos_mallado / 2)): epsilonr.append(epsilon_r2)
for i in range(int(puntos_mallado / 2)): sigma.append(sigma2)

matriz = np.ones((400, 400)) # Al ser 2D asignamos los valores de epsilon y sigma ayudándonos de una matriz de unos

epsilonr = (np.array(epsilonr) * matriz).T # Los pasamos a arrays en forma de matriz para poder trabajar con ellos
sigma = (np.array(sigma) * matriz).T


"Variables para evolución temporal de E"

cd = 1 / (2 * epsilonr)
ctc = sigma * dt / (2 * epsilonr * epsilon0)

ca = 1 - ctc / (1 + ctc)
cb = 1 / ((1 + ctc) * 2 * epsilonr)

var1 = np.zeros(401)  # variables para poder simular la absorbancia del medio
var2 = np.zeros(401)
var3 = np.zeros(401)
var4 = np.zeros(401)
var5 = np.zeros(401)
var6 = np.zeros(401)
var7 = np.zeros(401)
var8 = np.zeros(401)

for j in range(n):  # Ez = E conserva su valor en los bordes de la simulación

    t = t + dt

    Ez[1:, 1:] = ca * Ez[1:, 1:] + cb * (Hy[1:, 1:] - Hy[:-1, 1:]) - cb * (Hx[1:, 1:] - Hx[1:, :-1]) # Ecuación de evol. temporal de Ez

    Ez[0, :] = var1[:] # Absorción en la pared izquierda
    var1[:] = var2[:]
    var2[:] = Ez[1, :]

    Ez[:, 0] = var3[:] #  Absorción en la pared inferior
    var3[:] = var4[:]
    var4[:] = Ez[:, 1]

    Ez[:, -1] = var5[:] # Absorción en la pared superior
    var5[:] = var6[:]
    var6[:] = Ez[:, -2]

    Ez[-1, :] = var7[:] #  Absorción en la pared derecha
    var7[:] = var8[:]
    var8[:] = Ez[-2, :]

    E_pulso = E0 * np.exp(-1 / 2 * ((t - t0p) / dtp) ** 2) # Pulso

    Ez[xp[0], xp[1]] = E_pulso

    Hx[:, :-1] = Hx[:, :-1] - 0.5 * (Ez[:, 1:] - Ez[:, :-1]) # Ecuación de evol. temporal de Hx


    Hy[:-1, :] = Hy[:-1, :] + 0.5 * (Ez[1:, :] - Ez[:-1, :]) # Ecuación de evol. temporal de Hy


    if k == nrep: # Actualizamos gráfico cada nrep iteraciones

        ax1.cla()

        ax1.set_title("Propagación 2D de una onda EM")
        ax1.set_xlabel(r'Posición ($\mu m$)')
        ax1.set_ylabel(r'Posición ($\mu m$)')


        ax1.contourf(X, Y, np.clip(Ez, -0.1, 0.1), nivelesdecalor, cmap='viridis')

        ax1.axvline(x= 0.5 * puntos_mallado * dx * 1e6 , linewidth=1, linestyle="--", color='black')

        plt.pause(interv)

        k = 0

    k += 1

'''
### Código pintando los 3 campos ###

fig = plt.figure(figsize = (18,8))
ax1 = fig.add_subplot(131)  # Como si fuera 1 unico plt, facilmente modificable en el futuro
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

ax1.set_xlabel(r'Posición ($\mu m$)')
ax1.set_ylabel(r'Posición ($\mu m$)')
ax2.set_xlabel(r'Posición ($\mu m$)')
ax2.set_ylabel(r'Posición ($\mu m$)')
ax3.set_xlabel(r'Posición ($\mu m$)')
ax3.set_ylabel(r'Posición ($\mu m$)')

puntos_mallado = 401
malla = np.arange(puntos_mallado) * dx

Y, X = np.meshgrid(malla, malla)

Ez = np.zeros((401, 401))
Hx = np.zeros((401, 401))
Hy = np.zeros((401, 401))

# Hy = H[a, :], Hx = H[:, a]

nivelesdecalor = np.linspace(-0.1, 0.1, 21)

cs1 = ax1.contourf(X, Y, np.clip(Ez, -0.1, 0.1), nivelesdecalor, cmap='viridis')
cs2 = ax2.contourf(X, Y, np.clip(Hy, -0.1, 0.1), nivelesdecalor, cmap='viridis')
cs3 = ax3.contourf(X, Y, np.clip(Hx, -0.1, 0.1), nivelesdecalor, cmap='viridis')

bar1 = plt.colorbar(cs1, label = r'$\frac{V}{m}$', orientation = 'horizontal')
bar2 = plt.colorbar(cs2, label = r'$T$', orientation = 'horizontal')
bar3 = plt.colorbar(cs3, label = r'$T$', orientation = 'horizontal')

"Conductividad y permitividad"

epsilon0, sigma1, sigma2 = 8.85 * 1e-12, 0, 4000
epsilon_r1, epsilon_r2 = 1, 4

epsilonr = [epsilon_r1 for i in range(int(puntos_mallado / 2))]
sigma = [sigma1 for i in range(int(puntos_mallado / 2))]

for i in range(int(puntos_mallado / 2)): epsilonr.append(epsilon_r2)
for i in range(int(puntos_mallado / 2)): sigma.append(sigma2)

matriz = np.ones((400, 400))

epsilonr = (np.array(epsilonr) * matriz).T
sigma = (np.array(sigma) * matriz).T

cd = 1 / (2 * epsilonr)
ctc = sigma * dt / (2 * epsilonr * epsilon0)

ca = 1 - ctc / (1 + ctc)
cb = 1 / ((1 + ctc) * 2 * epsilonr)

var1 = np.zeros(401)  # variables para poder simular la absorbancia del medio
var2 = np.zeros(401)
var3 = np.zeros(401)
var4 = np.zeros(401)
var5 = np.zeros(401)
var6 = np.zeros(401)
var7 = np.zeros(401)
var8 = np.zeros(401)

for j in range(n):  # Ez = E se mantiene a 0 en X = 0 e Y = 0, Hy = 0 en Y = 400, Hx = 0 en X = 400

    t = t + dt

    Ez[1:, 1:] = ca * Ez[1:, 1:] + cb * (Hy[1:, 1:] - Hy[:-1, 1:]) - cb * (Hx[1:, 1:] - Hx[1:, :-1])

    Ez[0, :] = var1[:]
    var1[:] = var2[:]
    var2[:] = Ez[1, :]

    Ez[:, 0] = var3[:]
    var3[:] = var4[:]
    var4[:] = Ez[:, 1]

    Ez[:, -1] = var5[:]
    var5[:] = var6[:]
    var6[:] = Ez[:, -2]

    Ez[-1, :] = var7[:]
    var7[:] = var8[:]
    var8[:] = Ez[-2, :]

    E_pulso = E0 * np.exp(-1 / 2 * ((t - t0p) / dtp) ** 2)

    Ez[xp[0], xp[1]] = E_pulso

    Hx[:, :-1] = Hx[:, :-1] - 0.5 * (Ez[:, 1:] - Ez[:, :-1])

    Hy[:-1, :] = Hy[:-1, :] + 0.5 * (Ez[1:, :] - Ez[:-1, :])

    if k == nrep:
        ax1.cla()
        ax2.cla()
        ax3.cla()

        ax1.contourf(X, Y, np.clip(Ez, -0.1, 0.1), nivelesdecalor, cmap='viridis')
        ax1.set_xlabel(r'Posición ($\mu m$)')
        ax1.set_ylabel(r'Posición ($\mu m$)')
        
        ax2.contourf(X, Y, np.clip(Hy, -0.1, 0.1), nivelesdecalor, cmap='viridis')
        ax2.set_xlabel(r'Posición ($\mu m$)')
        ax2.set_ylabel(r'Posición ($\mu m$)')
        
        ax3.contourf(X, Y, np.clip(Hx, -0.1, 0.1), nivelesdecalor, cmap='viridis')
        ax3.set_xlabel(r'Posición ($\mu m$)')
        ax3.set_ylabel(r'Posición ($\mu m$)')
        
        ax1.axvline(x=0.5 * puntos_mallado * 1e-8, linewidth=1, linestyle="--", color='black')
        ax2.axvline(x=0.5 * puntos_mallado * 1e-8, linewidth=1, linestyle="--", color='black')
        ax3.axvline(x=0.5 * puntos_mallado * 1e-8, linewidth=1, linestyle="--", color='black')

        plt.pause(interv)

        k = 0

    k += 1
'''