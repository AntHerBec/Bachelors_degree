
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

E0 = hbarra**2 / (2 * me * a0**2)


def C(u, alpha):

    Cu = []

    for i in u:

        if abs(i) < 0.5:

            Cu.append(- k * alpha)

        else:

            Cu.append(k * (1 - alpha))

    return np.array(Cu)

def psi_pozo(alpha, x, metodo = 'dfinit', paridad = 'par'):

    if metodo == 'numerov':

        psi, phi = np.zeros(len(x)), np.zeros(len(x))
        deltax = x[1] - x[0]

        phi[0] = 0
        phi[1] = 1e-40 # número pequeño cualquiera, parametro a modificar para ver si mejora
        func = C(x, alpha)

        for i in range(2, len(x)):

            phi[i] = 2 * phi[i - 1] - phi[i - 2] + (deltax**2 * func[i - 1] / (1 - deltax ** 2 * func[i - 1] / 12)) * phi[i - 1]

        psi = phi * 12 / (12 - deltax ** 2 * func)

    elif metodo == 'dfinit':

        xpos = np.linspace(0, x[-1], int((len(x) + 1)/2))
        psipos, dpsi = np.zeros(len(xpos)), np.zeros(len(xpos))
        psineg = np.zeros(len(xpos) - 1)
        funcpos = C(xpos, alpha)

        if paridad == 'par':

            psipos[0] = 1
            dpsi[0] = 0
            simetria = 1

        elif paridad == 'impar':

            psipos[0] = 0
            dpsi[0] = 1
            simetria = -1

        else:
            raise ValueError("El valor para la paridad no es válido. Por favor, introduzca en su lugar 'par' o 'impar'")
            return None


        for i in range(1, len(xpos)):

            psipos[i] = psipos[i - 1] + dpsi[i - 1] * (xpos[i] - xpos[i - 1])
            dpsi[i] = dpsi[i - 1] + funcpos[i - 1] * psipos[i - 1] * (xpos[i] - xpos[i - 1])
            psineg[-i] = simetria * psipos[i]

        psi = [*psineg, *psipos]

    else:
        raise ValueError("El valor para el método no es válido. Por favor, introduzca en su lugar 'dfinit' o 'numerov'")
        return None

    return np.array(psi)

def V(r, Z = 1, litio = False):

    '''
    El potencial se simplifica porque todas las constantes (excepto el 2) que aparecen son igual a
    a0, de manera que al expresar V en unidades de a0 tenemos 1 unidad de a0
    '''

    v = 2 * Z / r

    if litio == True:

        v = 2 / r

        for i in range(len(r)):

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

def shooting_method (malla_alphas, x, alphatol = 1.5e-16, psitol = 1e-35 , niter = 1000, metodo = 'dfinit', paridad = 'par', L = 0, Z = 1, litio = False, atomo = True):

    lista = []

    if atomo == True:

        for i in range(1, len(malla_alphas)):

            inf, sup = malla_alphas[i - 1], malla_alphas[i]

            psi_inf = psi_atomo(x, inf, L=L, Z=Z, litio=litio)[-1]
            psi_sup = psi_atomo(x, sup, L=L, Z=Z, litio=litio)[-1]

            prod_psi = psi_inf * psi_sup

            if prod_psi < 0:

                for n in range(niter + 1):

                    pmedio = (inf + sup) / 2

                    if abs(inf - sup) <= alphatol:
                        break

                    psi_medio = psi_atomo(x, pmedio, L=L, Z=Z, litio=litio)[-1]
                    psi_inf = psi_atomo(x, inf, L=L, Z=Z, litio=litio)[-1]

                    if abs(psi_medio) <= psitol:
                        break

                    if psi_medio * psi_inf > 0:
                        sup = pmedio

                    else:
                        inf = pmedio

                lista.append(pmedio)

    else:

        for i in range(1, len(malla_alphas)):

            inf, sup = malla_alphas[i-1], malla_alphas[i]

            prod_psi = psi_pozo(inf, x, metodo = metodo, paridad = paridad)[-1] * psi_pozo(sup, x, metodo = metodo, paridad = paridad)[-1]

            if prod_psi < 0:

                for n in range(niter + 1):

                    pmedio = (inf + sup)/2

                    if abs(inf - sup) <= alphatol:
                        break

                    psi_medio = psi_pozo(pmedio, x, metodo = metodo, paridad = paridad)[-1]
                    psi_inf = psi_pozo(inf, x, metodo = metodo, paridad = paridad)[-1]

                    if abs(psi_medio) <= psitol:
                        break

                    elif psi_medio * psi_inf < 0: sup = pmedio

                    elif psi_medio * psi_inf > 0: inf = pmedio

                lista.append(pmedio)

    return lista

"Resultados del pozo"

u = np.linspace(-3, 3, 3001) # mallados siempre simétricos respecto al 0, para simplificar las funciones

multiplos = np.linspace(0, 1, 101)

fig = plt.figure()

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

valores_analiticos = np.array([0.098036, 0.383114, 0.808305])

alfas = shooting_method(multiplos, u, atomo = False)
alfas.append(*shooting_method(multiplos, u, paridad = 'impar', atomo = False))
alfas.append(alfas[1])
alfas.remove(alfas[1])

alfas_round = np.round(np.array(alfas), 6)
print("Valores de alpha para diferencias finitas:")
print(alfas_round)
errores_pozo_dif = (abs(valores_analiticos - alfas_round) / valores_analiticos * 100)

print("Errores relativos usando diferencias finitas:")
print(errores_pozo_dif)

ax1.plot(u, psi_pozo(alfas[0], u)/max(abs(psi_pozo(alfas[0], u))), label = r'$\alpha$ = {:.6f}'.format(alfas[0]))
ax1.plot(u, psi_pozo(alfas[1], u, paridad = 'impar')/max(abs(psi_pozo(alfas[1], u, paridad = 'impar'))), label = r'$\alpha$ = {:.6f}'.format(alfas[1]))
ax1.plot(u, psi_pozo(alfas[2], u)/max(abs(psi_pozo(alfas[2], u))), label = r'$\alpha$ = {:.6f}'.format(alfas[2]))

ax1.set_title("Diferencias Finitas")
ax1.set_xlabel("u (unidades de a)")
ax1.set_ylabel(r"$\psi$ normalizada")
ax1.set_xlim(-4,4)
ax1.set_ylim(-2,2)

alfas = shooting_method(multiplos, u, metodo = 'numerov', atomo = False)
alfas_round = np.round(np.array(alfas), 6)

print("\n", "Valores de alpha para Numerov:")
print(alfas_round)
errores_pozo_num = (abs(valores_analiticos - alfas_round) / valores_analiticos * 100)

print("Errores relativos usando Numerov:")
print(errores_pozo_num)

ax2.plot(u, psi_pozo(alfas[0], u, metodo = 'numerov')/max(abs(psi_pozo(alfas[0], u, metodo = 'numerov'))), label = r'$\alpha$ = {:.6f}'.format(alfas[0]))
ax2.plot(u, -psi_pozo(alfas[1], u, metodo = 'numerov')/max(abs(psi_pozo(alfas[1], u, metodo = 'numerov'))), label = r'$\alpha$ = {:.6f}'.format(alfas[1]))
ax2.plot(u, -psi_pozo(alfas[2], u, metodo = 'numerov')/max(abs(psi_pozo(alfas[2], u, metodo = 'numerov'))), label = r'$\alpha$ = {:.6f}'.format(alfas[2]))

ax2.set_title("Numerov")
ax2.set_xlabel("u (unidades de a)")
ax2.set_ylabel(r"$\psi$ normalizada")
ax2.set_xlim(-4,4)
ax2.set_ylim(-2,2)

ax1.axvline(-0.5, color = 'black', linestyle = 'dashed')
ax1.axvline(0.5, color = 'black', linestyle = 'dashed')
ax2.axvline(-0.5, color = 'black', linestyle = 'dashed')
ax2.axvline(0.5, color = 'black', linestyle = 'dashed')

plt.tight_layout()

ax1.legend()
ax2.legend()
plt.show()


segmentos = [101, 501, 1001, 3001, 5001]

err_dif_alpha1 = []
err_dif_alpha2 = []
err_dif_alpha3 = []
err_num_alpha1 = []
err_num_alpha2 = []
err_num_alpha3 = []
t_dif = []
t_num = []

for i in segmentos:

    u = np.linspace(-5, 5, i)

    multiplos = np.linspace(0, 1, 101)


    print("\n", "Mallado de posiciones de {} segmentos".format(i), "\n", "-"*40)

    t0 = t.time()

    alfas = shooting_method(multiplos, u, atomo=False)
    alfas.append(*shooting_method(multiplos, u, paridad='impar', atomo=False))
    alfas.append(alfas[1])
    alfas.remove(alfas[1])

    t1 = t.time()
    t_dif.append(t1-t0)

    alfas = np.array(alfas)
    alfas_round = np.round(alfas, 6)
    print("Valores de alpha para diferencias finitas:")
    print(alfas_round)
    errores_pozo_dif = np.round((abs(valores_analiticos - alfas) / valores_analiticos * 100), 3)

    err_dif_alpha1.append(errores_pozo_dif[0])
    err_dif_alpha2.append(errores_pozo_dif[1])
    err_dif_alpha3.append(errores_pozo_dif[2])


    print("Errores relativos usando diferencias finitas:")
    print(errores_pozo_dif)

    t0 = t.time()

    alfas = shooting_method(multiplos, u, metodo='numerov', atomo=False)

    t1 = t.time()
    t_num.append(t1 - t0)

    alfas = np.array(alfas)
    alfas_round = np.round(alfas, 6)

    print("\n", "Valores de alpha para Numerov:")
    print(alfas_round)
    errores_pozo_num = np.round((abs(valores_analiticos - alfas) / valores_analiticos * 100), 3)

    print("Errores relativos usando Numerov:")
    print(errores_pozo_num)

    err_num_alpha1.append(errores_pozo_num[0])
    err_num_alpha2.append(errores_pozo_num[1])
    err_num_alpha3.append(errores_pozo_num[2])


plt.plot(segmentos, err_dif_alpha1, label = "diferencias finitas")
plt.plot(segmentos, err_num_alpha1, label = "Numerov")
plt.title(r"Errores relativos con respecto a $\alpha$ = {}, pozo de anchura 10a".format(valores_analiticos[0]))
plt.xlabel("Número de segmentos")
plt.ylabel("Error relativo cometido (%)")
plt.legend()
plt.show()

plt.plot(segmentos, err_dif_alpha2, label = "diferencias finitas")
plt.plot(segmentos, err_num_alpha2, label = "Numerov")
plt.title(r"Errores relativos con respecto a $\alpha$ = {}, pozo de anchura 10a".format(valores_analiticos[1]))
plt.xlabel("Número de segmentos")
plt.ylabel("Error relativo cometido (%)")
plt.legend()
plt.show()

plt.plot(segmentos, err_dif_alpha3, label = "diferencias finitas")
plt.plot(segmentos, err_num_alpha3, label = "Numerov")
plt.title(r"Errores relativos con respecto a $\alpha$ = {}, pozo de anchura 10a".format(valores_analiticos[2]))
plt.xlabel("Número de segmentos")
plt.ylabel("Error relativo cometido (%)")
plt.legend()
plt.show()

plt.plot(segmentos, t_dif, label = "diferencias finitas")
plt.plot(segmentos, t_num, label = "Numerov")
plt.title(r"tiempo de ejecución de cada método")
plt.xlabel("Número de segmentos")
plt.ylabel("t(s)")
plt.legend()
plt.show()

extremos_mallado = [3, 4, 5]

err_dif_alpha1 = []
err_dif_alpha2 = []
err_dif_alpha3 = []
err_num_alpha1 = []
err_num_alpha2 = []
err_num_alpha3 = []

for i in extremos_mallado:

    u = np.linspace(-i, i, 5001)

    multiplos = np.linspace(0, 1, 101)

    print("\n", "pozo de anchura {}a".format(2*i), "\n", "-" * 20)
    alfas = shooting_method(multiplos, u, atomo=False)
    alfas.append(*shooting_method(multiplos, u, paridad='impar', atomo=False))
    alfas.append(alfas[1])
    alfas.remove(alfas[1])

    alfas = np.array(alfas)
    alfas_round = np.round(alfas, 6)
    print("Valores de alpha para diferencias finitas:")
    print(alfas_round)
    errores_pozo_dif = np.round((abs(valores_analiticos - alfas) / valores_analiticos * 100), 3)

    err_dif_alpha1.append(errores_pozo_dif[0])
    err_dif_alpha2.append(errores_pozo_dif[1])
    err_dif_alpha3.append(errores_pozo_dif[2])

    print("Errores relativos usando diferencias finitas:")
    print(errores_pozo_dif)

    alfas = shooting_method(multiplos, u, metodo='numerov', atomo=False)
    alfas = np.array(alfas)
    alfas_round = np.round(alfas, 6)

    print("\n", "Valores de alpha para Numerov:")
    print(alfas_round)
    errores_pozo_num = np.round((abs(valores_analiticos - alfas) / valores_analiticos * 100), 3)

    print("Errores relativos usando Numerov:")
    print(errores_pozo_num)

    err_num_alpha1.append(errores_pozo_num[0])
    err_num_alpha2.append(errores_pozo_num[1])
    err_num_alpha3.append(errores_pozo_num[2])

plt.plot(extremos_mallado, err_dif_alpha1, label="diferencias finitas")
plt.plot(extremos_mallado, err_num_alpha1, label="Numerov")
plt.title(r"Errores relativos con respecto a $\alpha$ = {}, 5001 segmentos".format(valores_analiticos[0]))
plt.xlabel("Anchura del pozo")
plt.ylabel("Error relativo cometido (%)")
plt.legend()
plt.show()

plt.plot(extremos_mallado, err_dif_alpha2, label="diferencias finitas")
plt.plot(extremos_mallado, err_num_alpha2, label="Numerov")
plt.title(r"Errores relativos con respecto a $\alpha$ = {}, 5001 segmentos".format(valores_analiticos[1]))
plt.xlabel("Anchura del pozo")
plt.ylabel("Error relativo cometido (%)")
plt.legend()
plt.show()

plt.plot(extremos_mallado, err_dif_alpha3, label="diferencias finitas")
plt.plot(extremos_mallado, err_num_alpha3, label="Numerov")
plt.title(r"Errores relativos con respecto a $\alpha$ = {}, 5001 segmentos".format(valores_analiticos[2]))
plt.xlabel("Anchura del pozo")
plt.ylabel("Error relativo cometido (%)")
plt.legend()
plt.show()




"Resultados de átomos"


r = np.linspace(300, 1e-10, 3001)
deltar = abs((r[1]-r[0]))


L = [0, 1, 2]
L_li = [0, 1]
Z = [1, 2]

multiplo_li1 = np.linspace(0.45, 0.35, 101) * (-1)
multiplo_li2 = np.linspace(0.35, 0.2, 101) * (-1)
multiplo_li3 = np.linspace(0.2, 0.1, 101) * (-1)



multiplo1 = np.linspace(1.05, 0.95, 101) * (-1)
multiplo2 = np.linspace(0.27, 0.23, 101) * (-1)
multiplo3 = np.linspace(0.13, 0.09, 101) * (-1)

errores = []
errores_litio = []

contador = 1

for j in L:
    for n in Z:

        if j == 0:
            multiplo = multiplo1
            esperado = -1
        if j == 1:
            multiplo = multiplo2
            esperado = -0.25
        if j == 2:
            multiplo = multiplo3
            esperado = -0.1111

        E = shooting_method(multiplo, r, L = j, Z = n)
        E.append(shooting_method(multiplo / 4, r, L=j, Z=n)[0])
        E.append(shooting_method(multiplo / 9, r, L=j, Z=n)[0])
        E.append(shooting_method(multiplo / 16, r, L=j, Z=n)[0])


        errores.append(abs(np.array(E) - np.array([esperado, esperado / 4, esperado / 9, esperado / 16])) / abs(np.array([esperado, esperado / 4, esperado / 9, esperado / 16])) * 100)

        lista_ev = [i for i in E[:4]] # para que los decimales no se redondeen



        plt.subplot(3, 2, contador)

        plt.title(r"$\psi$ para los 4 primeros estados ligados (Z = {}, L = {})".format(n,j))
        plt.xlabel(r"radio (unidades de $a_0$)")
        plt.ylabel(r"$\psi$")


        for i in E:

            psi = psi_atomo(r, i, L = j, Z = n)
            psi_norm = psi / max(psi)
            plt.plot(r, psi_norm, label = "{:.4f} eV".format(i*E0*n**2))

        plt.legend(loc="upper right")
        contador += 1


plt.legend()
plt.tight_layout()
plt.show()

for i in range(len(errores)):

    print("errores relativos para la gráfica número {}".format(i + 1))
    print(errores[i])

E = shooting_method(multiplo_li1, r, L = 0, litio = True)
E.append(shooting_method(multiplo_li2, r, L = 1, litio = True)[0])
E.append(shooting_method(multiplo_li3, r, L = 0, litio = True)[0])

valores_experimentales = np.array([-5.392, -3.544, -2.019])
errores_litio.append(abs(np.array(E) * E0 - valores_experimentales) / abs(valores_experimentales) * 100)

print("Errores relativos obtenidos para el litio:")
print(errores_litio)


fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title(r"$\psi$ para los 3 primeros estados ligados".format(n, j))
ax1.set_xlabel(r"radio (unidades de $a_0$)")
ax1.set_ylabel(r"$\psi$")

plt.tight_layout()

for i in E:

    if i == E[1]: j = 1
    else: j = 0
    psi = psi_atomo(r, i, L = j, litio = True)
    psi_norm = psi / max(psi)
    ax1.plot(r, psi_norm, label = "{:.4f} eV".format(i*E0))
    ax1.set_xlim(-5,50)


plt.legend()
plt.show()




