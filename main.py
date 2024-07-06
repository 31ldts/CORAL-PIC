from analyze_interactions import *

change_directory("outputs")

matrix = analyze_files(directory="ichem-ifp", protein=True, ligand=True, subunit=False)
#save_matrix(matrix=matrix, filename='cdu01 (prot-lig (without subunits)).csv')
matrix = analyze_files(directory="ichem-ifp", protein=True, ligand=False)
#save_matrix(matrix=matrix, filename='cdu01 (prot).csv')
matrix = analyze_files(directory="ichem-ifp", protein=False, ligand=True, subunit=False)
#save_matrix(matrix=matrix, filename='cdu01 (lig).csv')

selection = sort_matrix(matrix=matrix, axis='columns', residue_chain=True)
#save_matrix(matrix=selection, filename='cdu02 (cols).csv')
selection = sort_matrix(matrix=matrix, selected_items=5, axis='rows')
#save_matrix(matrix=selection, filename='cdu02 (rows).csv')

transposed = transpose_matrix(matrix=selection)
#save_matrix(matrix=transposed, filename='cdu03.csv')

removed = sort_matrix(matrix=matrix, thr_interactions = 15, axis="columns")
#save_matrix(matrix=removed, filename='cdu04 (cols).csv')
removed = sort_matrix(matrix=matrix, axis="columns")
#save_matrix(matrix=removed, filename='cdu04 (rows).csv')

plot_matrix(matrix=matrix, plot_name="cdu05 (Original)", axis="columns", stacked=False, save=False)
plot_matrix(matrix=removed, plot_name="cdu05 (Removed Rows)", axis="rows", stacked=True, save=False)

filtered = filter_by_interaction(matrix=matrix, interactions=[5])
#save_matrix(matrix=matrix, filename='cdu08 (5).csv')

aux = 0