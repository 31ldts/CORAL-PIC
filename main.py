from analyze_interactions import *
'''import matplotlib.pyplot as plt
import numpy as np

# Datos de ejemplo
categorias = ['A', 'B', 'C', 'D']
grupos = [[20, 35, 30, 35], [25, 32, 34, 20]]

# Crear la figura y los ejes
fig, ax = plt.subplots()

# Posición de las barras
indices = np.arange(len(categorias))

# Ancho de las barras
barras = []
labels = ['Grupo 1', 'Grupo 2']
for index, grupo in enumerate(grupos):
    if index == 0:
        barras.append(ax.bar(indices, grupo, label=labels[index]))
    else:
        barras.append(ax.bar(indices, grupo, bottom=grupos[index-1], label=labels[index]))


# Personalización de la gráfica
ax.set_xlabel('Categorías')
ax.set_ylabel('Valores')
ax.set_title('Gráfico de Barras Apiladas')
ax.set_xticks(indices)
ax.set_xticklabels(categorias)
ax.legend()

# Mostrar la gráfica
plt.show()'''


change_directory("outputs")

matrix = analyze_files(directory="structures", protein=True, ligand=True)
#save_matrix(matrix=matrix, filename='cdu01 (prot-lig).csv')
matrix = analyze_files(directory="structures", protein=True, ligand=True, subunit=False)
#save_matrix(matrix=matrix, filename='cdu01 (prot-lig (without subunits)).csv')
matrix = analyze_files(directory="structures", protein=True, ligand=False)
#save_matrix(matrix=matrix, filename='cdu01 (prot).csv')
matrix = analyze_files(directory="structures", protein=False, ligand=True)
#save_matrix(matrix=matrix, filename='cdu01 (lig).csv')

selection = select_reactives(matrix=matrix, selectedItems=5, axis='columns')
#save_matrix(matrix=selection, filename='cdu02 (cols).csv')
selection = select_reactives(matrix=matrix, selectedItems=5, axis='rows')
#save_matrix(matrix=selection, filename='cdu02 (rows).csv')

transposed = transpose_matrix(matrix=selection)
#save_matrix(matrix=transposed, filename='cdu03.csv')

removed = select_reactives(matrix=matrix, threshold = 15, axis="columns")
#save_matrix(matrix=removed, filename='cdu04 (cols).csv')
removed = select_reactives(matrix=matrix, threshold = 15, axis="rows")
#save_matrix(matrix=removed, filename='cdu04 (rows).csv')

plot_matrix(matrix=matrix, plot_name="cdu05 (Original)", axis="columns", label_x="Drugs", label_y="Interactions", title="Interactions of the drugs", stacked=False, save=False)
#plot_matrix(matrix=removed, plotName="cdu05 (Removed Rows)", axis="rows", labelX="Residue", labelY="Interactions", title="Interactions of the residues that are over the threshold")

filtered = filter_by_interaction(matrix=matrix, interactions=[5])
#save_matrix(matrix=matrix, filename='cdu08 (5).csv')
