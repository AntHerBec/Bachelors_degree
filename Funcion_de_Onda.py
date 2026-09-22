import numpy as np
import matplotlib.pyplot as plt

"Parámetros Iniciales"

xmin = -5
xmax = 5
puntos_mallado = 1001

dx = 0.01 # m
dt = 5e-5 # s

malla = np.linspace(-5, 5, 1001) # mallado del sistema

"Datos de la simulación"

t = 0.0 # tiempo inicial
n = 5000  # nº pasos sistema
nrep = 50 # cada cuantas iteraciones queremos pintar
tpausa = 0.01 # cuanto tiempo esperara el programa para pintar

k = nrep

"funciones recurrentes"


def norma_func(phi_r, phi_i): # Calcula la norma de la funcion de onda
    norma_phi = sum(phi_r ** 2 + phi_i ** 2) * dx

    return norma_phi


def normalizar(phi_r, phi_i): # Normaliza la funcion de onda
    norma_phi = norma_func(phi_r, phi_i)

    cte = 1 / np.sqrt(norma_phi)

    return (phi_r * cte, phi_i * cte)


def seg_der(f):
    # Consideramos que mas alla de los extremos del intervalo la funcion se anula, calculamos la segunda derivada

    d2f = np.zeros(len(f))

    d2f[1:-1] = (f[2:] + f[0:-2] - 2 * f[1:-1]) / dx ** 2
    d2f[0] = (f[1] - 2 * f[0]) / dx ** 2
    d2f[-1] = (f[-2] - 2 * f[-1]) / dx ** 2

    return d2f


def seg_der_period(f): # Considerando que la función es periodica en los bordes, calculamos la segunda derivada
    d2f = np.zeros(len(f))

    d2f[1:-1] = (f[2:] + f[0:-2] - 2 * f[1:-1]) / dx ** 2
    d2f[0] = (f[1] + f[-1] - 2 * f[0]) / dx ** 2
    d2f[-1] = (f[0] + f[-2] - 2 * f[-1]) / dx ** 2

    return d2f


def prim_der(f): # Calculamos la primera derivada
    df = np.zeros(len(f))

    df[1:-1] = 0.5 * (f[2:len(f)] - f[0:-2]) / dx
    df[0] = 0.5 * (f[1] - f[-1]) / dx
    df[-1] = 0.5 * (f[2] - f[0]) / dx

    return df


def val_esp_x(phi_r, phi_i, malla): # calculamos el valor esperado de la posicion o potencias de esta
    valesp = sum((phi_r ** 2 + phi_i ** 2) * malla) * dx

    return valesp


def val_esp_p(phi_r, phi_i): # calculamos el valor esperado de p
    valesp = sum(phi_r * prim_der(phi_i) - phi_i * prim_der(phi_r)) * dx

    return valesp


def val_esp_p2(phi_r, phi_i): # calculamos el valor esperado de p^2
    valesp = - sum(phi_i * seg_der_period(phi_i) + phi_r * seg_der_period(phi_r)) * dx

    return valesp


def hermite(n, x): # calculamos el polinomio de hermite de grado dado n
    H0 = 1
    H1 = 2 * x

    if n >= 2:

        for i in range(n - 1):
            Hn = 2 * x * H1 - 2 * (i + 1) * H0
            H0 = H1
            H1 = Hn

        return Hn

    if n == 0: return H0

    if n == 1: return H1

"Parámetros del potencial"
x0 = 0
sigma = 0.5
k0 = 10
w = 4

V = 0.5 * w ** 2 * malla ** 2 # potencial V para la simulación

phi_r = np.exp(- 0.5 * ((malla - x0) / sigma)**2)*np.cos(k0 * malla) # Definimos las partes real e imaginaria de la función de onda

phi_i = np.exp(- 0.5 * ((malla - x0) / sigma)**2)*np.sin(k0 * malla)

norma_0 = norma_func(phi_r,phi_i) # Calculamos la norma inicial de la función de onda

phi_r, phi_i = normalizar(phi_r, phi_i) # Normalizamos la función para la representación

fig = plt.figure(figsize = (12,7))
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(234)
ax3 = fig.add_subplot(235)
ax4 = fig.add_subplot(236)



ax.plot(malla, V/max(V), "black", linestyle = "--") # Evolución temporal de la función de onda
parte_real, = ax.plot(malla, phi_r, label = 'parte real')
parte_imag, = ax.plot(malla, phi_i, label = 'parte imaginaria')
modulo, = ax.plot(malla, np.sqrt(phi_r**2 + phi_i**2),"black", label = 'módulo de la onda')
ax.set_title(r'Propagación de la función de onda ($\psi$) bajo un potencial  V = $\frac{\omega ^2 x^2}{2}$')
ax.legend(loc = 'upper right')


norm = norma_func(phi_r, phi_i) # Cambios en la norma de la función de onda
norma_list = [norm - 1]
t_list = [t]
norma, = ax2.plot(norma_list, 'black')
ax2.set_title(r'cambios en la norma de $\psi$ a lo largo del tiempo') # La norma 1 corresponde con y = 0 en el grafico, es decir, no cambia la norma

xesp_list = [val_esp_x(phi_r, phi_i, malla)]
pesp_list = [val_esp_p(phi_r, phi_i)]

xesp, = ax3.plot(t_list, xesp_list, label = "<x>") # Evolución temporal de los valores esperados de x y p
pesp, = ax3.plot(t_list, pesp_list, label = "<p>")
ax3.set_title('Valores esperados de x y p')
ax3.legend(loc = 'upper right')


deltax_list = [val_esp_x(phi_r, phi_i, malla**2) - val_esp_x(phi_r, phi_i, malla)**2]
deltap_list = [val_esp_p2(phi_r, phi_i) - val_esp_p(phi_r, phi_i)**2]
dxdp_list = [deltax_list[0] * deltap_list[0]]

deltax, = ax4.plot(t_list, deltax_list, label = r"$\Delta x$") # Evolución temporal de las dispersiones de x y p, y comprobación del principio de incertidumbre
deltap, = ax4.plot(t_list, deltap_list, label = r"$\Delta p$")
dxdp, = ax4.plot(t_list, dxdp_list, label = r"$\Delta x \Delta p$")
ax4.set_title('Dispersiones de x y p. Principio de incertidumbre')
ax4.legend(loc = 'upper right')


for i in range(n): # Bucle que va actualizando las partes de la funcion de onda cada iteración y las graficas cada nrep repeticiones

    t = t + dt

    phi_i_un_medio = phi_i + 0.25 * dt * seg_der_period(phi_r) - 0.5 * V * dt * phi_r
    phi_r = phi_r - 0.5 * dt * seg_der_period(phi_i_un_medio) + V * dt * phi_i_un_medio
    phi_i = phi_i_un_medio + 0.25 * dt * seg_der_period(phi_r) - 0.5 * V * dt * phi_r


    if k == nrep:

        plt.pause(tpausa)

        parte_real.set_data(malla, phi_r)
        parte_imag.set_data(malla, phi_i)
        modulo.set_data(malla, np.sqrt(phi_r**2 + phi_i**2))


        t_list.append(t)

        norm = norma_func(phi_r, phi_i) - 1
        norma_list.append(norm)

        norma.set_data(t_list, norma_list)

        ax2.set_ylim(min(norma_list), max(norma_list))
        ax2.set_xlim(0,t)

        xesp_list.append(val_esp_x(phi_r, phi_i, malla))
        pesp_list.append(val_esp_p(phi_r, phi_i))

        xesp.set_data(t_list, xesp_list)
        pesp.set_data(t_list, pesp_list)

        ax3.set_xlim(0,t)
        ax3.set_ylim(-10, 10)


        deltax_list.append(np.sqrt(val_esp_x(phi_r, phi_i, malla**2) - (val_esp_x(phi_r, phi_i, malla))**2))
        deltap_list.append(np.sqrt(val_esp_p2(phi_r, phi_i) - (val_esp_p(phi_r, phi_i))**2))
        dxdp_list.append(deltax_list[-1] * deltap_list[-1])

        deltax.set_data(t_list, deltax_list)
        deltap.set_data(t_list, deltap_list)
        dxdp.set_data(t_list, dxdp_list)
        ax4.set_xlim(0,t)
        ax4.set_ylim(0, 2)

        k = 0

    k += 1

phi_r, phi_i = normalizar(phi_r, phi_i)

norma_f = norma_func(phi_r,phi_i) # Calculamos la norma final de la función normalizada

plt.close(fig)

print(('La norma inicial de la función de onda es: {}').format(norma_0)) # Escribimos las normas inicial y final
print(('La norma final de la función de onda es: {}').format(norma_f))


"Hermite"


t = 0.0 # Reseteamos el tiempo inicial a 0

n_herm = int(input('Introduce el grado del polinomio de Hermite a utilizar: ')) # Selección del grado del polinomio de hermite

phi_n_r = np.exp(- malla ** 2 * w / 2) * hermite(n_herm, np.sqrt(w) * malla) # Definimos las partes real e imaginaria en el instante inicial
phi_n_i = np.zeros(len(malla))

## A partir de aquí, los pasos son análogos a los de la anterior parte ##
norma_0 = norma_func(phi_n_r, phi_n_i)

phi_n_r, phi_n_i = normalizar(phi_n_r, phi_n_i)

fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(234)
ax3 = fig.add_subplot(235)
ax4 = fig.add_subplot(236)

parte_real, = ax.plot(malla, phi_n_r, label='parte real')
parte_imag, = ax.plot(malla, phi_n_i, label='parte imaginaria')
modulo, = ax.plot(malla, np.sqrt(phi_n_r ** 2 + phi_n_i ** 2), "black", label=r'módulo de $\psi$')
ax.set_title('Función de onda con polinomio de hermite de grado {}'.format(n_herm))
ax.legend(loc='upper right')

norm = norma_func(phi_n_r, phi_n_i)
norma_list = [norm - 1]
t_list = [t]
norma, = ax2.plot(norma_list, 'black')
ax2.set_title(r'cambios en la norma de $\psi$')


xesp_list = [val_esp_x(phi_n_r, phi_n_i, malla)]
pesp_list = [val_esp_p(phi_n_r, phi_n_i)]

xesp, = ax3.plot(t_list, xesp_list, label="<x>")
pesp, = ax3.plot(t_list, pesp_list, label="<p>")
ax3.set_title('Valores esperados de x y p')
ax3.legend(loc='upper right')

deltax_list = [val_esp_x(phi_n_r, phi_n_i, malla ** 2) - val_esp_x(phi_n_r, phi_n_i, malla) ** 2]
deltap_list = [val_esp_p2(phi_n_r, phi_n_i) - val_esp_p(phi_n_r, phi_n_i) ** 2]
dxdp_list = [deltax_list[0] * deltap_list[0]]

deltax, = ax4.plot(t_list, deltax_list, label=r"$\Delta x$")
deltap, = ax4.plot(t_list, deltap_list, label=r"$\Delta p$")
dxdp, = ax4.plot(t_list, dxdp_list, label=r"$\Delta x \Delta p$")
ax4.set_title('Dispersiones de x y p. Principio de incertidumbre')
ax4.legend(loc='upper right')

for i in range(n):

    t = t + dt

    phi_n_i_un_medio = phi_n_i + 0.25 * dt * seg_der_period(phi_n_r) - 0.5 * V * dt * phi_n_r
    phi_n_r = phi_n_r - 0.5 * dt * seg_der_period(phi_n_i_un_medio) + V * dt * phi_n_i_un_medio
    phi_n_i = phi_n_i_un_medio + 0.25 * dt * seg_der_period(phi_n_r) - 0.5 * V * dt * phi_n_r

    if k == nrep:
        plt.pause(tpausa)

        parte_real.set_data(malla, phi_n_r)
        parte_imag.set_data(malla, phi_n_i)
        modulo.set_data(malla, np.sqrt(phi_n_r ** 2 + phi_n_i ** 2))

        t_list.append(t)

        norm = norma_func(phi_n_r, phi_n_i) - 1
        norma_list.append(norm)

        norma.set_data(t_list, norma_list)

        ax2.set_ylim(min(norma_list), max(norma_list))
        ax2.set_xlim(0, t)

        xesp_list.append(val_esp_x(phi_n_r, phi_n_i, malla))
        pesp_list.append(val_esp_p(phi_n_r, phi_n_i))

        xesp.set_data(t_list, xesp_list)
        pesp.set_data(t_list, pesp_list)

        ax3.set_xlim(0, t)
        ax3.set_ylim(min(min(pesp_list), min(xesp_list)), max(max(pesp_list), max(xesp_list)))

        deltax_list.append(np.sqrt(val_esp_x(phi_n_r, phi_n_i, malla ** 2) - (val_esp_x(phi_n_r, phi_n_i, malla)) ** 2))
        deltap_list.append(np.sqrt(val_esp_p2(phi_n_r, phi_n_i) - (val_esp_p(phi_n_r, phi_n_i)) ** 2))
        dxdp_list.append(deltax_list[-1] * deltap_list[-1])

        deltax.set_data(t_list, deltax_list)
        deltap.set_data(t_list, deltap_list)
        dxdp.set_data(t_list, dxdp_list)
        ax4.set_xlim(0, t)
        ax4.set_ylim(0, 20)

        k = 0

    k += 1

phi_r, phi_i = normalizar(phi_n_r, phi_n_i)

norma_f = norma_func(phi_n_r, phi_n_i)

print(('La norma inicial de la función de onda es: {}').format(norma_0))
print(('La norma final de la función de onda es: {}').format(norma_f))

plt.close(fig)
