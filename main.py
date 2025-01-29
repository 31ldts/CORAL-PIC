from analyze_interactions import *

analizador = AnalyzeInteractions()
analizador.change_directory("outputs", mode='output')
analizador.change_directory("input", mode='input')

matrix = analizador.analyze_files(directory="ichem-ifp", mode='ichem', activity_file="ichem-ifp-activities.csv")
#matrix = analizador.filter_chain(matrix=matrix, chain=None, subpocket_path='subpockets.csv', subpockets=["S1'", "S2"], save="chain.csv")
analizador.heatmap(interaction_data=matrix, title="ichem-ifp-min", mode="min", save=False)
analizador.heatmap(interaction_data=matrix, title="ichem-ifp-max", mode="max", save=False)
analizador.heatmap(interaction_data=matrix, title="ichem-ifp-mean", mode="mean", save=False)
analizador.heatmap(interaction_data=matrix, title="ichem-ifp-count", mode="count", save=False)
analizador.heatmap(interaction_data=matrix, title="ichem-ifp-percentage", mode="percent", save=False)

#analizador.save_matrix(matrix=matrix, filename='_ichem.csv')

'''matrix = analizador.analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=True)
# save_matrix(matrix=matrix, filename='prueba2.csv')
matrix = analizador.analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=False, ligand=True, subunit=True)
# save_matrix(matrix=matrix, filename='prueba3.csv')
matrix = analizador.analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=False, subunit=True)
# save_matrix(matrix=matrix, filename='prueba4.csv')
matrix = analizador.analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=False)
# save_matrix(matrix=matrix, filename='prueba5.csv')
matrix = analizador.analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=False, ligand=True, subunit=False)
# save_matrix(matrix=matrix, filename='prueba6.csv')

matrix = analizador.analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=False)
selection = analizador.sort_matrix(matrix=matrix, axis='columns', residue_chain=True)
# save_matrix(matrix=selection, filename='prueba7.csv')
selection = analizador.sort_matrix(matrix=matrix, selected_items=5, axis='rows')
# save_matrix(matrix=selection, filename='prueba8.csv')

transposed = analizador.transpose_matrix(matrix=selection, save='3.1 transpose.csv')
analizador.save_matrix(matrix=transposed, filename='prueba9n.csv')

removed = analizador.sort_matrix(matrix=matrix, thr_activity= 6.45)
# save_matrix(matrix=removed, filename='prueba10.csv')
removed = analizador.sort_matrix(matrix=matrix, axis="columns") # Ordena según las columnas que  tengan más interacciones
# save_matrix(matrix=removed, filename='prueba11.csv')

filtered = analizador.filter_by_interaction(matrix=matrix, interactions=[5])
# save_matrix(matrix=filtered, filename='prueba12.csv')

filtered = analizador.filter_by_subunit(matrix=matrix, subunits=['A'])
# save_matrix(matrix=filtered, filename='prueba13.csv')
matrix = analizador.analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=False)
filtered = analizador.filter_by_subunit(matrix=matrix, subunits=['B'])
# save_matrix(matrix=filtered, filename='prueba14.csv')

analizador.plot_matrix(matrix=matrix, plot_name="prueba15", axis="columns", stacked=True)

matrix = analizador.analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=True)
analizador.plot_matrix(matrix=matrix, plot_name="cprueba16", axis="rows", stacked=True)
analizador.plot_matrix(matrix=filtered, plot_name="prueba17", axis="rows", stacked=True)
analizador.plot_matrix(matrix=analizador.remove_empty_axis(filtered), plot_name="prueba17", axis="rows", stacked=True)
'''