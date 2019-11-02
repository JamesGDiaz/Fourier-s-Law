""" Proyecto de la clase de Mecanica de Fluidos II
    Alumnos: Alejandro Camacho Ake
             Franco Giordani Diaz

    Ecuacion de conduccion de calor en una dimension
    mediante el metodo de diferencias finitas de 
    primer orden.
"""
import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm as cm
from progressbar import Counter, Timer, ProgressBar


copper = '#6b2e0d'

""" Propiedades del material """
# Del oro Au Au Au
# rho = 19300  # Densidad [kg m^-3]
# K = 318  # Conductividad Termica [W m^-1 K^-1]
# Cp = 130  # Capacidad Calorifica [J kg^-1 K^-1]

# Del cobre CUCUCUCU
rho = 8960  # Densidad [kg m^-3]
K = 410  # Conductividad Termica [W m^-1 K^-1]
Cp = 390  # Capacidad Calorifica [J kg^-1 K^-1]
alpha = K/(rho*Cp)  # Difusividad Termica [m^2 s^-1]
print(alpha)

""" Definicion del dominio """
# Espacial
L = 1  # Longitud unitaria [m]
nx = 50  # nodos
dx = L/nx
x = np.linspace(0, L, nx+1)
# Temporal
t_end = 3000  # [s]
dt = 1.666666666667  # 100ms
nt = int(t_end/dt)+1
dTdt = np.empty(nx+1)
t = np.linspace(0, t_end, nt)

Tdist = np.zeros((len(t), len(x)))

""" Condiciones Iniciales """
T0 = 20
T = np.exp(x*4.605170185)
""" Condiciones de Frontera """
TA = 1
TB = 100
T[0] = TA
T[-1] = TB

""" Simulacion """
# Condicion de estabilidad
courant = alpha*dt/(dx**2)
if(courant > 0.5):
    print(
        f"Courant: {courant:0.3f} (>= 0.5).\nLa solucion NO es estable. Abortando...")
    quit()
else:
    print(f"CFL: {courant:0.3f}. La solucion es estable.")

widgets = ['Solved: ', Counter(), ' steps (', Timer(), ')']
pbar = ProgressBar(widgets=widgets)
for i in pbar((i for i in range(nt))):  # tiempo
    plt.clf()
    for n in range(1, nx):  # espacio
        dTdt[n] = alpha*(-(T[n]-T[n-1])/(dx**2) + (T[n+1]-T[n])/(dx**2))
    dTdt[0] = 0  # en la frontera A
    dTdt[nx] = 0  # en la frontera B
    T = T + dTdt*dt  # incremento (T + dTdt*dt = T+dT)
    Tdist[i, :] = T

""" widgets = ['Renderizando 1/2: ', Counter(), f' de {nt}',
           ' [', Timer(format='%s'), ']']
pbar = ProgressBar(widgets=widgets)
for i in pbar(i for i in range(nt)):
    plt.clf()
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.35)
    plt.plot(x, Tdist[i, :], c=copper, linewidth='3.0')
    plt.text(0.49, 116, f"t: {i*dt:0.2f}s", size=15,
             ha="center", va="center",
             bbox=dict(boxstyle="round",
                       ec=(0.8, 0.8, 0.5),
                       fc=(.95, .95, .95),
                       )
             )
    plt.axis([0, L, min(T0, TA, TB)-10, max(T0, TA, TB)+10])
    plt.xlabel('Distancia [m]')
    plt.ylabel('Temperatura [°C]')
    plt.savefig(f"./copper/dist-{i}-{nt}.png")

widgets = ['Renderizando 2/2: ', Counter(), f' de {nt}',
           ' [', Timer(format='%s'), ']']
pbar = ProgressBar(widgets=widgets)
for i in pbar(i for i in range(nt)):
    grid = [Tdist[i, :], Tdist[i, :]]
    plt.clf()
    plt.figure(1)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(False)
    frame.axes.get_yaxis().set_visible(False)
    plt.imshow(grid, origin='upper', interpolation='gaussian',
               cmap=cm.afmhot, vmin=0, vmax=102)
    plt.savefig(f"./copper/wire-{i}-{nt}.png") """

fig, ax = plt.subplots(figsize=(8, 5))
plt.axis([0, L, min(T0, TA, TB)-10, max(T0, TA, TB)+10])
plt.xlabel('Distancia [m]')
plt.ylabel('Temperatura [°C]')
lines = [0, 50, 120, 260, 450, 750, 1750]
for i in range(len(lines)):
    plt.plot(x, Tdist[lines[i], :], c=copper, linewidth='1')
    ax.annotate(f't={lines[i]*dt:0.2f}s',
                xy=(x[15+(i*5)], Tdist[lines[i], 15+(i*5)]-0.6), xycoords='data',
                xytext=(0.05+(0.09*i), 60+(6*i)), textcoords='data',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc,angleA=0,armA=60,rad=10"))

plt.show()
