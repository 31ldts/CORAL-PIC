from analyze_interactions import AnalyzeInteractions as new_api

analyzer = new_api()

analyzer.change_directory("output", mode=analyzer.OUTPUT)
analyzer.change_directory("input", mode=analyzer.INPUT)

analyzer.set_config(heat_max_cols=50)

data = analyzer.analyze_files(directory="386_arpeggio", mode=analyzer.ARPEGGIO, activity_file="386_Mpro_nc_February_2025.csv")
'''
for mode in [analyzer.COUNT, analyzer.MEAN]:
    analyzer.heatmap(interaction_data=data, title=mode, mode=mode, save=False)
'''
activity = analyzer.sort_matrix(interaction_data=data, thr_activity=7.1)
activity = analyzer.remove_empty_axis(interaction_data=activity)

for mode in [analyzer.COUNT, analyzer.MEAN]:
    analyzer.heatmap(interaction_data=activity, title=mode, mode=mode, save=False)
