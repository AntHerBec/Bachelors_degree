"Definimos condiciones del problema"

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as r
from scipy.optimize import curve_fit


N0 = 1000 # nº nucleos de Bi en el instante inicial
t = np.linspace(0, 2000, 201) # mallado de tiempo en dias
deltat = t[1] - t[0]

tau_bi = 7.8 # vida media del Bi en dias
tau_po = 196
omega_bi = 1 / tau_bi
omega_po = 1 / tau_po

def prob_bi (t):
    return np.exp(-t * omega_bi) * omega_bi

def prob_po (t):
    return np.exp(-t * omega_po) * omega_po

def bi_real (n, t):
    return n * np.exp(-omega_bi * t)

def po_real (n, t):
    return n * omega_bi / (omega_po - omega_bi) * (-np.exp(-omega_po * t) + np.exp(-omega_bi * t))

def numero_nucleos(N_0, t, num_reps, nucleo = 'bi'):

    long = len(t)
    total_aceptados = np.zeros(long)

    bi_max = max(prob_bi(t))
    po_max = max(prob_po(t))

    for k in range(num_reps):

        aceptados = np.zeros(long)

        u = r.random(long)
        x = u * (t[-1] - t[0]) + t[0]

        x.sort()

        if nucleo == 'bi': w = prob_bi(x) / bi_max
        elif nucleo == 'po': w = prob_po(x) / po_max - prob_bi(x) / bi_max
        else: raise ValueError('El valor introducido para nucleo no es válido.')

        for n in range(N_0):

            y = r.random(long)

            for i in range(long):

                if w[i] > y[i]: aceptados[i] += 1

        total_aceptados += aceptados

    return total_aceptados / num_reps

def bi_ajuste (a, b, t):
    return b * np.exp(-a * t)

def po_ajuste(t, a, b, c):
    return a / (b - a) * c * (np.exp(-a * t) - np.exp(-b * t))

def plot_y_errores(N0, t, nreps, errores):

    bi = numero_nucleos(N0, t, nreps, 'bi')
    po = numero_nucleos(N0, t, nreps, 'po')
    pb = N0 - bi - po

    fig1 = plt.figure(figsize = (10,6))
    ax1 = fig1.add_subplot(111)

    ax1.scatter(t, bi, s = 1, label = r'$^{{210}} Bi $')
    ax1.scatter(t, po, s = 1, label = r'$^{{210}} Po $')
    ax1.scatter(t, pb, s = 1, label = r'$^{{206}} Pb $')

    ax1.set_title(r'Desintegraciones para {} núcleos durante {} días (promedio de {} repeticiones)'.format(N0, int(t[-1]), nreps), fontsize = 'large')
    ax1.set_xlabel('t (s)', fontsize = 'large')
    ax1.set_ylabel('nº de núcleos', fontsize = 'large')
    plt.legend(fontsize = 'large', loc = 'center right')
    plt.show()

    fig2 = plt.figure(figsize=(10, 6))
    ax2 = fig2.add_subplot(111)

    ax2.scatter(t, bi , s = 1, label = 'simulación')
    ax2.set_title(r'nº $^{{210}} Bi$ para {} núcleos iniciales durante {} días (promedio de {} repeticiones)'.format(N0, int(t[-1]), nreps), fontsize = 'large')
    ax2.plot(t, bi_real(N0, t), color = 'red', label = 'teórica')
    ax2.set_xlabel('t (s)', fontsize = 'large')
    ax2.set_ylabel('nº de núcleos', fontsize = 'large')
    plt.show()

    fig3 = plt.figure(figsize=(10, 6))
    ax3 = fig3.add_subplot(111)

    ax3.scatter(t, po, s = 1, label = 'simulación')
    ax3.set_title(r'nº de $^{{210}} Po$ para {} núcleos iniciales durante {} días (promedio de {} repeticiones)'.format(N0, int(t[-1]), nreps), fontsize = 'large')
    ax3.plot(t, po_real(N0, t), color = 'red', label = 'teórica')
    ax3.set_xlabel('t (s)', fontsize = 'large')
    ax3.set_ylabel('nº de núcleos', fontsize = 'large')
    plt.show()

    fig4 = plt.figure(figsize=(10, 6))
    ax4 = fig4.add_subplot(111)

    ax4.scatter(t, pb, s = 1, label = 'simulación')
    ax4.set_title(r'nº de $^{{206}} Pb$ para {} núcleos iniciales durante {} días (promedio de {} repeticiones)'.format(N0, int(t[-1]), nreps), fontsize = 'large')
    ax4.plot(t, N0 - bi_real(N0, t) - po_real(N0, t), color = 'red', label = 'teórica')
    ax4.set_xlabel('t (s)', fontsize = 'large')
    ax4.set_ylabel('nº de núcleos', fontsize = 'large')
    plt.show()

    errores.append(abs(bi- bi_real(N0, t)))
    errores.append(abs(po - po_real(N0, t)))
    errores.append(abs(pb - (N0 - bi_real(N0, t) - po_real(N0, t))))




"Variando repeticiones"

reps = [1, 5, 20, 100]
N0 = 1000 # nº nucleos de Bi en el instante inicial
t = np.linspace(0, 1000, 2001) # mallado de tiempo en dias

errores_bi = []
errores_po = []
errores_pb = []

for i in reps:

    errores_reps = []
    plot_y_errores(N0, t, i, errores_reps)
    errores_bi.append(errores_reps[0])
    errores_po.append(errores_reps[1])
    errores_pb.append(errores_reps[2])

errores_reps = [errores_bi, errores_po, errores_pb]
nombres = [r'$^{{210}} Bi$', r'$^{{210}} Po$', r'$^{{206}} Pb$']

for i in range(len(errores_reps)):

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title('errores absolutos entre simulación y curva teórica del {}'.format(nombres[i]), fontsize = 'large')
    ax.set_xlabel('t (s)', fontsize = 'large')
    ax.set_ylabel('error absoluto', fontsize = 'large')

    for j in range(len(reps)):
        ax.semilogy(t, errores_reps[i][j], label = '{} repeticiones'.format(reps[j]))
        
    plt.legend()
    plt.show()


"Variando la muestra inicial"

reps = 30
N = [100, 500, 1000] # nº nucleos de Bi en el instante inicial
t = np.linspace(0, 1000, 1001) # mallado de tiempo en dias

errores_bi = []
errores_po = []
errores_pb = []

teoricas_bi = []
teoricas_po = []
teoricas_pb = []

for i in N:

    errores_N = []
    plot_y_errores(i, t, reps, errores_N)
    errores_bi.append(errores_N[0])
    errores_po.append(errores_N[1])
    errores_pb.append(errores_N[2])

    teoricas_bi.append(bi_real(i, t))
    teoricas_po.append(po_real(i, t))
    teoricas_pb.append(i - bi_real(i, t) - po_real(i, t))

teoricas = [teoricas_bi, teoricas_po, teoricas_pb]


errores_N = [errores_bi, errores_po, errores_pb]

for i in range(len(errores_N)):

    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot(121)
    ax1.set_title('errores absolutos entre simulación y curva teórica del {}'.format(nombres[i]), fontsize = 'large')
    ax1.set_xlabel('t (s)', fontsize = 'large')
    ax1.set_ylabel('error absoluto', fontsize = 'large')

    ax2 = fig.add_subplot(122)
    ax2.set_title('errores relativos entre simulación y curva teórica del {}'.format(nombres[i]), fontsize='large')
    ax2.set_xlabel('t (s)', fontsize='large')
    ax2.set_ylabel('error relativo (%)', fontsize='large')

    for j in range(len(N)):
        ax1.semilogy(t, errores_N[i][j], label = '{} núcleos'.format(N[j]))
        ax2.semilogy(t, errores_N[i][j] / max(teoricas[i][j]) * 100, label = '{} núcleos'.format(N[j]))

    plt.tight_layout()
    plt.legend()
    plt.show()


"Variando el instante de tiempo final"

reps = 30
N0 = 1000 # nº nucleos de Bi en el instante inicial
tiempo = [np.linspace(0, 200, 501), np.linspace(0, 500, 501), np.linspace(0, 2000, 501), np.linspace(0, 5000, 501)] # mallado de tiempo en dias

errores_bi = []
errores_po = []
errores_pb = []

for i in tiempo:

    errores_tiempo = []
    plot_y_errores(N0, i, reps, errores_tiempo)
    errores_bi.append(errores_tiempo[0])
    errores_po.append(errores_tiempo[1])
    errores_pb.append(errores_tiempo[2])


errores_tiempo = [errores_bi, errores_po, errores_pb]

for i in range(len(errores_tiempo)):

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title('errores absolutos entre simulación y curva teórica del {}'.format(nombres[i]), fontsize='large')
    ax.set_xlabel('t (s)', fontsize='large')
    ax.set_ylabel('error absoluto', fontsize='large')

    for j in range(len(tiempo)):
        ax.semilogy(tiempo[j], errores_tiempo[i][j], label='instante final en {:.0f} s'.format(tiempo[j][-1]))

    plt.legend()
    plt.show()

"Variando la densidad del mallado"

reps = 30
N0 = 1000 # nº nucleos de Bi en el instante inicial
mallado = [np.linspace(0, 1000, 101), np.linspace(0, 1000, 301), np.linspace(0, 1000, 1001), np.linspace(0, 1000, 5001)] # mallado de tiempo en dias

errores_bi = []
errores_po = []
errores_pb = []

for i in mallado:

    errores_mallado = []
    plot_y_errores(N0, i, reps, errores_mallado)
    errores_bi.append(errores_mallado[0])
    errores_po.append(errores_mallado[1])
    errores_pb.append(errores_mallado[2])


errores_mallado = [errores_bi, errores_po, errores_pb]

for i in range(len(errores_mallado)):

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    ax.set_title('errores absolutos entre simulación y curva teórica del {}'.format(nombres[i]), fontsize='large')
    ax.set_xlabel('t (s)', fontsize='large')
    ax.set_ylabel('error absoluto', fontsize='large')

    for j in range(len(mallado)):
        ax.semilogy(mallado[j], errores_mallado[i][j], label='mallado de {:.0f} puntos'.format(len(mallado[j])))

    plt.legend()
    plt.show()

"Máxima producción de alfas"

for i in [101, 501, 1001, 3001]:
    for j in [1, 10, 30]:
        t = np.linspace(0, 1000, i)
        po = numero_nucleos(1000, t, j, 'po')
        teor = po_real(1000, t)

        fig1 = plt.figure(figsize = (10,6))
        ax1 = fig1.add_subplot(111)

        ax1.scatter(t, po, s = 1, label = 'simulación')
        ax1.scatter(t[np.argmax(po)], max(po), s = 50, label = 'máximo simulado')
        ax1.plot(t, teor, color='red', label='teórica')
        ax1.scatter(t[np.argmax(teor)], max(teor), s = 50, label = 'máximo teórico')

        ax1.set_title(r'Máxima producción de partículas $\alpha$ para 1000 núcleos durante 1000 días (promedio de {} repeticiones, {} puntos de mallado)'.format(j, i), fontsize = 'large')
        ax1.set_xlabel('t (s)', fontsize = 'large')
        ax1.set_ylabel('nº de núcleos', fontsize = 'large')
        plt.legend(fontsize = 'large', loc = 'center right')
        plt.show()

        print(r'Simulación de 1000 núcleos durante 1000 días (promedio de {} repeticiones, {} puntos de mallado)'.format(j, i), '\n', 40*'-')
        print(r'Máximo simulado = {:.2f} núcleos, máximo teórico = {:.2f} núcleos. Error relativo = {:.1f}% '.format(max(po), max(teor), abs(max(po) - max(teor)) / max(teor) * 100))
        print(r't del máximo simulado = {:.2f} días, t del máximo teórico = {:.2f} días. Error relativo = {:.1f}% '.format(t[np.argmax(po)], t[np.argmax(teor)], abs(t[np.argmax(po)] - t[np.argmax(teor)]) / t[np.argmax(teor)] * 100), 2*'\n')

"Estimación de vida media"

for i in [101, 501, 1001]:
    for j in [50, 100]:
        t = np.linspace(0, j, i)
        bi = numero_nucleos(1000, t, 100, 'bi')
        po = numero_nucleos(1000, t, 100, 'po')
        biteor = bi_real(1000, t)
        poteor = po_real(1000, t)

        ajuste = curve_fit(bi_ajuste, t, bi)
        ajuste2 = curve_fit(po_ajuste, t, po, p0 = [1/10, 1/200, 1000])

        func = bi_ajuste(ajuste[0][1], ajuste[0][0], t)
        func2 = po_ajuste(t, ajuste2[0][0], ajuste2[0][1], ajuste2[0][2])

        fig1 = plt.figure(figsize=(10, 6))
        ax1 = fig1.add_subplot(111)

        ax1.plot(t, func, label= 'ajuste bi')
        ax1.plot(t, biteor, color='red', label='bi teórica')
        ax1.plot(t, func2, label = 'ajuste po')
        ax1.plot(t, poteor, color = 'green', label = 'po teórica')

        ax1.set_title(r'Ajuste de las desintegraciones para 1000 núcleos durante {} días (promedio de 100 repeticiones, {} puntos de mallado)'.format(j, i), fontsize='large')
        ax1.set_xlabel('t (s)', fontsize='large')
        ax1.set_ylabel('nº de núcleos', fontsize='large')
        plt.legend(fontsize='large', loc='center right')
        plt.show()

        print(r'Simulación de 1000 núcleos durante {} días (promedio de 100 repeticiones, {} puntos de mallado)'.format(j, i), '\n', 80 * '-')
        print(r'Bi: N0 simulado = {:.0f} núcleos, N0 teórico = {:.0f} núcleos. Error relativo = {:.1f}% '.format(ajuste[0][0], max(biteor), abs(ajuste[0][0] - max(biteor)) / max(biteor) * 100))
        print(r'Bi: Vida media simulada = {:.2f} días, vida media teórica = {:.1f} días. Error relativo = {:.1f}% '.format(1 / ajuste[0][1], tau_bi, abs(1 / ajuste[0][1] - tau_bi) / tau_bi * 100))
        print(r'Po: N0 simulado = {:.0f}, N0 teórico = {:.0f}. Error relativo = {:.1f}% '.format(ajuste2[0][2], max(biteor), abs(ajuste2[0][2] - max(biteor)) / max(biteor) * 100))
        print(r'Po: Vida media simulada = {:.2f} días, vida media teórica = {:.0f} días. Error relativo = {:.1f}% '.format(1 / ajuste2[0][1], tau_po, abs(1 / ajuste2[0][1] - tau_po) / tau_po * 100), 2 * '\n')




