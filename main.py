from analyze_interactions import *

change_directory("outputs")


matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv")
save_matrix(matrix=matrix, filename='0 aux.csv')

matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=True)
save_matrix(matrix=matrix, filename='1.1 prot-lig-subunits.csv')
matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=False, ligand=True, subunit=True)
save_matrix(matrix=matrix, filename='1.2 !prot-lig-subunits.csv')
matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=False, subunit=True)
save_matrix(matrix=matrix, filename='1.3 prot-!lig-subunits.csv')
matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=False)
save_matrix(matrix=matrix, filename='1.4 prot-lig-!subunits.csv')
matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=False, ligand=True, subunit=False)
save_matrix(matrix=matrix, filename='1.5 !prot-lig-!subunits.csv')

matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=True)
#selection = sort_matrix(matrix=matrix, axis='columns', residue_chain=True)
#save_matrix(matrix=selection, filename='2.1 sort by residue chain.csv')
selection = sort_matrix(matrix=matrix, selected_items=5, axis='rows')
save_matrix(matrix=selection, filename='2.2 sort and select bests 5.csv')

transposed = transpose_matrix(matrix=selection)
save_matrix(matrix=transposed, filename='3.1 transpose.csv')

removed = sort_matrix(matrix=matrix, thr_activity= 6.45, axis="columns")
save_matrix(matrix=removed, filename='4.1 sort and select by threshold.csv')
removed = sort_matrix(matrix=matrix, axis="columns")
save_matrix(matrix=removed, filename='4.2 sort an axis.csv')

filtered = filter_by_interaction(matrix=matrix, interactions=[5])
save_matrix(matrix=filtered, filename='5.1 filter by interaction.csv')

filtered = filter_by_subunit(matrix=matrix, subunits=['A'])
save_matrix(matrix=filtered, filename='6.1 filter by subunit (with subunits).csv')
matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=False)
filtered = filter_by_subunit(matrix=matrix, subunits=['B'])
save_matrix(matrix=filtered, filename='6.2 filter by subunit (with subunits).csv')

plot_matrix(matrix=matrix, plot_name="cdu05 (Original)", axis="columns", stacked=True, save=True)
plot_matrix(matrix=removed, plot_name="cdu05 (Removed Rows)", axis="rows", stacked=True, save=False)