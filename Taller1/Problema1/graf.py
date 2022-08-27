import matplotlib.pyplot as plt
import numpy as np

epidemia = np.loadtxt(
    "c:/Users/Esteban/Desktop/Uni Esteban/Simulaciones/RepositorioGrupoSimulaciones/Taller1/Problema1//datosP1a.txt"
)
tiempo = epidemia[:, 0]
susceptibles = epidemia[:, 1]
infectados = epidemia[:, 2]
recuperados = epidemia[:, 3]

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(tiempo, susceptibles, "b", label="Susceptibles")
ax.plot(tiempo, infectados, "r", label="Infectados")
ax.plot(tiempo, recuperados, "g", label="Recuperados")
ax.set_title("Modelo SIR sin muertes o nacimientos")
ax.set_xlabel("Tiempo")
ax.set_ylabel("Susceptibles, Infectados, Recuperados")
ax.legend(loc=(0.05, 0.7))
ax.set_xlim(xmin=0, xmax=70)
ax.set_ylim(ymin=0, ymax=1)
ax.tight_layout()
