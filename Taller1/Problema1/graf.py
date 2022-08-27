import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True

epidemia = np.loadtxt(
    "c:/Users/Esteban/Desktop/Uni Esteban/Simulaciones/RepositorioGrupoSimulaciones/Taller1/Problema1//datosP1a.txt"
)

tiempo = epidemia[:, 0]
susceptibles = epidemia[:, 1]
infectados = epidemia[:, 2]
recuperados = epidemia[:, 3]

slim = np.loadtxt(
    "c:/Users/Esteban/Desktop/Uni Esteban/Simulaciones/RepositorioGrupoSimulaciones/Taller1/Problema1/datosP1c.txt"
)

R0 = slim[:, 0]
sinf = slim[:, 1]

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 4))

# Grafica de SIR
ax1.plot(tiempo, susceptibles, "b", label="Susceptibles")
ax1.plot(tiempo, infectados, "r", label="Infectados")
ax1.plot(tiempo, recuperados, "g", label="Recuperados")
ax1.set_title("Modelo SIR sin muertes o nacimientos")
ax1.set_xlabel("Tiempo")
ax1.set_ylabel("Susceptibles, Infectados, Recuperados")
ax1.legend(loc=(0.0185, 0.7))
ax1.set_xlim(xmin=0, xmax=70)
ax1.set_ylim(ymin=0, ymax=1)

# Grafica de S limite
ax2.plot(R0, sinf, "c", label=r"No. de Susceptibles si t $\rightarrow \infty$")
ax2.set_title(r"$s_{\infty}$ vs $R_{0}$")
ax2.legend()
ax2.set_xlabel(r"R$_{0}=\beta / \gamma$", fontsize=15)
ax2.set_ylabel(r"s$_{\infty}$", fontsize=15)
ax2.set_xlim(xmin=R0[0], xmax=R0[-1])
ax2.set_ylim(ymin=0, ymax=0.014)
