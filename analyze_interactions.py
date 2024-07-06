import csv
import os
from matplotlib import pyplot as plt
from typing import Union
import numpy as np
from matplotlib.ticker import MaxNLocator
import mplcursors
import copy
from collections import Counter
import re

saving_directory = os.getcwd()

def change_directory(path: str) -> None:
    """
    Changes the saving directory to a subdirectory within the current directory.

    Args:
        path (str): The name of the subdirectory to change to.

    Returns:
        None

    Raises:
        ValueError: If the subdirectory does not exist.
    """
    global saving_directory

    # Get the full path of the new saving directory
    saving_directory = os.path.join(os.getcwd(), path)

    # Check if the new saving directory exists
    if not os.path.exists(saving_directory):
        raise ValueError("The saving directory must exist inside the project.")

def transpose_matrix(matrix: list) -> list:
    """
    Transposes the given matrix.

    Args:
        matrix (list of lists): The matrix to transpose.

    Returns:
        list of lists: The transposed matrix.
    """
    verify_dimensions(matrix=matrix)

    # Use list comprehension to transpose the matrix
    return [[row[i] for row in matrix] for i in range(len(matrix[0]))]

def save_matrix(matrix: list, filename: str) -> None:
    """
    Saves the matrix to a CSV file in the specified directory.

    Args:
        matrix (list of lists): The matrix to be saved.
        filename (str): The name of the file to save the matrix in.

    Returns:
        None
    """
    global saving_directory

    verify_dimensions(matrix=matrix)

    # Create the CSV file
    with open(os.path.join(saving_directory, filename), 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        
        # Write the rows of the matrix to the CSV file
        for row in matrix:
            csv_writer.writerow(row)

def analyze_files(directory: str, activity_file: str=None, protein: bool=True, ligand: bool=True, subunit: bool=True) -> list:
    """
    Analyzes interaction data files in a specified directory, categorizing interactions
    based on protein and ligand atoms involved.

    Args:
        directory (str): The directory path containing interaction data files.
        activity_file (str): Path to activity file.
        protein (bool): Flag indicating whether to include protein atoms in analysis.
        ligand (bool): Flag indicating whether to include ligand atoms in analysis.
        subunit (bool): Flag indicating whether to differentiate residues from different chains.

    Returns:
        list: A matrix representing categorized interactions between residues and files.

    Raises:
        FileNotFoundError: If the specified directory does not exist or cannot be accessed.
    """
    
    def label_matrix(matrix: list, rows: list, columns: list, activity_file: str) -> list:
        """
        Labels the matrix with residue names and file names for clarity.

        Args:
            matrix (list): The 2D list representing interaction data.
            rows (list): List of residue names.
            columns (list): List of file names (drugs).
            activity_file (str): Path to activity file.

        Returns:
            list: The labeled matrix with headers.
        """
        # Copy rows and columns to avoid modifying input parameters
        rows = [row.replace("\t", "") for row in rows]
        columns = [""] + columns[:]  # Include an empty string at the beginning

        if activity_file is not None:
            if not os.path.isfile(activity_file):
                raise FileNotFoundError(f"The file '{activity_file}' does not exist.")
            
            # Read data from CSV into a dictionary
            data_dict = {}
            with open(activity_file, newline='') as csvfile:
                csvreader = csv.reader(csvfile)
                next(csvreader)  # Skip header
                for key, value in csvreader:
                    data_dict[key] = value
            
            # Update column names with data from dictionary
            for i in range(1, len(columns)):
                drug_name = columns[i]
                if drug_name in data_dict:
                    columns[i] = f"{drug_name}_{data_dict[drug_name]}"
                else:
                    columns[i] = f"{drug_name}_0"

        # Insert columns as the first row in the matrix
        matrix.insert(0, columns)
        
        # Insert rows as the first element of each row in the matrix
        for i, row in enumerate(matrix[1:], start=1):
            row.insert(0, rows[i-1])
        
        return matrix
    
    def modify_cell(text: str, interaction: str, atoms: str) -> str:
        """
        Modifies the cell content based on interaction type and atoms involved.

        Args:
            text (str): The current content of the cell.
            interaction (str): The type of interaction.
            atoms (str): The atoms involved in the interaction.

        Returns:
            str: The modified cell content with updated interaction details.
        """
        # Map interaction types to specific codes
        interaction_map = {
            "Hydrophobic": '1',
            "Aromatic_Face/Face": '2',
            "Aromatic_Edge/Face": '3',
            "HBond_PROT": '4',
            "HBond_LIG": '5',
            "Ionic_PROT": '6',
            "Ionic_LIG": '7'
        }
        interaction = interaction_map.get(interaction, '8')

        if text == "-":
            return f"{interaction} |{atoms}|"

        content = text.replace(", ", "").split("|")[:-1]
        exists = False
        cell = ''

        # Check if interaction type already exists in the cell
        for index, segment in enumerate(content):
            if index % 2 == 0 and interaction == segment.strip():
                content[index + 1] += f" {atoms}"
                exists = True
                break

        # Reconstruct the cell content
        for index, segment in enumerate(content):
            if index % 2 == 0:
                cell += segment
            else:
                cell += f"|{segment}|, "
        cell = cell.rstrip(', ')

        # Add new interaction if it does not exist
        if not exists:
            cell += f", {interaction} |{atoms}|"
        return cell

    def read_txt_file(file_name: str) -> list:
        """
        Reads the content of a text file and returns it as a list of strings, with each string representing a line from the file.

        Args:
            file_name (str): The name of the text file to read.

        Returns:
            list: A list containing each line of the file as a string.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            Exception: For any other unexpected errors during file reading.
        """
        string_list = []  # Initialize an empty list to store lines of the file
        try:
            with open(file_name, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    string_list.append(line.strip())  # Remove newline character and add the line to the list
        except FileNotFoundError:
            print(f"Error: The file '{file_name}' does not exist.")
            raise  # Raise the error to notify the caller
        except Exception as e:
            print(f"Error: An unexpected error occurred while reading '{file_name}': {e}")
            raise  # Raise the error to notify the caller

        return string_list

    def adjust_subunits(matrix: list, subunits: list) -> list:
        def reduce_strings_in_even_indices(atoms, subunits):
            pattern = r'[() ]'
            result = list(filter(None, re.split(pattern, atoms)))
            # Filtrar las celdas con índice par
            even_index_strings = [result[i] for i in range(0, len(result), 2)]
            
            # Contar las apariciones de cada string
            string_count = Counter(even_index_strings)
            
            # Crear una nueva lista manteniendo solo una aparición de cada string
            # que aparece en índices pares si aparece más de 3 veces
            text = ''
            added_strings = set()
            
            for i in range(len(result)):
                string = result[i]
                if i % 2 == 0:  # Índice par
                    if string_count[string] >= subunits:
                        if string not in added_strings:
                            text = text + string + " "
                            added_strings.add(string)
                        # No agregamos el string si ya fue agregado
                    else:
                        text = text + string+'('+result[i+1]+') '
            
            return '|' + text[:-1] + '|'
        
        for row in range(len(matrix)):
            for column in range(len(matrix[row])):
                cell = matrix[row][column]
                text = ''
                if cell != '-':
                    sections = cell.split("|")
                    for i in range(1, len(sections), 2):
                        # Dividir segun '|' y quedarme con los impares
                        text = text + sections[i-1] + reduce_strings_in_even_indices(atoms=sections[i], subunits=len(subunits)) + ' '
                        # Diccionario
                    matrix[row][column] = text[:-1]
        return matrix

    files = os.listdir(directory)
    matrix = []
    aa = {}
    cont = 0
    subunits_set = set()
    for index, file in enumerate(files):
        file_path = os.path.join(directory, file)
        if os.path.isfile(file_path):
            content = read_txt_file(file_path)
            for line in content:
                elements = line.split("|")
                if len(elements) == 10:
                    interaction = elements[0].strip().replace("\t", "")
                    residue = elements[3].strip().replace("\t", "")
                    if not subunit:
                        sections = residue.split("-")
                        residue = sections[0]
                        subunits_set.add(sections[1])
                    
                    # Determine which atoms to include based on flags
                    if protein and ligand:
                        atoms = f"{elements[1].strip().replace(chr(9), '')}-{elements[4].strip().replace(chr(9), '')}"
                    elif protein:
                        atoms = elements[1].strip().replace("\t", "")
                    else:
                        atoms = elements[4].strip().replace("\t", "")

                    if not subunit:
                        atoms = atoms + "(" + sections[1] + ")"

                    # Add new residue to dictionary if it doesn't exist
                    if residue not in aa:
                        aa[residue] = cont
                        cont += 1
                    column = aa[residue]

                    # Ensure the matrix has enough columns
                    if len(matrix) <= column:
                        new_column = ["-"] * len(files)
                        matrix.append(new_column)

                    cell = matrix[column][index]
                    cell = modify_cell(text=cell, interaction=interaction, atoms=atoms)
                    matrix[column][index] = cell

        files[index] = file.replace(".txt", "")

    if not subunit:
        matrix = adjust_subunits(matrix=matrix, subunits=list(subunits_set))
    return label_matrix(matrix=matrix, rows=list(aa.keys()), columns=files, activity_file=activity_file)

def verify_dimensions(matrix: list):
    """
    Verifies the dimensions of the matrix to ensure it contains at least 2 rows and 2 columns.

    Args:
        matrix (list of lists): The matrix to be verified.

    Raises:
        ValueError: If the matrix has fewer than 2 rows or any row has fewer than 2 columns.

    Example:
        If matrix = [['-', '1'], ['2', '-']], verify_dimensions(matrix) will pass.
    """
    if len(matrix) < 2 or any(len(row) < 2 for row in matrix):
        raise ValueError("There are not interactions on the matrix.")

def sort_matrix(matrix: list, axis: str, thr_interactions: int=None, selected_items: int=None, count: bool=False, residue_chain: bool=False) -> list:
    """
    Sorts and selects reactive rows or columns from a matrix based on interactions.

    Args:
        matrix (list of lists): The matrix containing interaction data.
        axis (str): Specifies whether to select rows ('rows') or columns ('columns').
        thr_interactions (int, optional): Minimum number of interactions to select a row/column.
        selected_items (int, optional): Number of top rows/columns to select based on interactions.
        count (bool, optional): If True, returns counts of interactions instead of the matrix.
        residue_chain (bool, optional): If True, orders the resulting matrix based on the residues' order in the chain.

    Returns:
        list of lists: Selected rows or columns from the matrix based on the specified criteria.

    Raises:
        ValueError: If both thr_interactions and selected_items are provided simultaneously.
        ValueError: If the matrix dimensions are insufficient.
    """

    def get_interactions(cell: str) -> int:
        """
        Counts the number of interactions in a cell of formatted interaction data.

        Args:
            cell (str): The cell containing interaction data formatted with '|' separators.

        Returns:
            int: The total number of interactions found in the cell.
        """
        interactions = 0
        sections = cell.split("|")
        for index in range(1, len(sections), 2):
            interactions += len(sections[index].split(" "))
        return interactions
    
    def sort_by_residue(matrix: list) -> list:
        """
        Sorts the matrix based on the residue indices in the first column.

        Args:
            matrix (list of lists): The matrix to be sorted.

        Returns:
            list of lists: The matrix sorted by residue indices.
        """
        # Validate matrix dimensions
        verify_dimensions(matrix=matrix)

        # Determine if sorting is needed by rows or columns
        axis = get_residues_axis(matrix=matrix)
        
        # If sorting is needed by columns, transpose the matrix first
        if axis == 'columns':
            matrix = transpose_matrix(matrix=matrix)

        # Separate the header from the data rows
        header = matrix[0]
        data_rows = matrix[1:]

        # Sort the data rows based on the residue indices
        sorted_data_rows = sorted(data_rows, key=lambda row: int(row[0].replace(" ", "")[3:]))

        # Combine the header with the sorted data rows
        sorted_matrix = [header] + sorted_data_rows

        # If sorting was by columns, transpose the sorted matrix back
        if axis == 'columns':
            sorted_matrix = transpose_matrix(matrix=sorted_matrix)

        return sorted_matrix

    # Validate matrix dimensions
    verify_dimensions(matrix=matrix)

    # Raise an error if both thr_interactions and selected_items are provided simultaneously
    if thr_interactions is not None and selected_items is not None:
        raise ValueError("You cannot select by 'threshold' and by 'selected_items' at the same time.")

    # Transpose the matrix if axis is 'columns'
    if axis == 'columns':
        matrix = transpose_matrix(matrix=matrix)
    
    # Initialize a dictionary to store the number of interactions per row
    reactives = {}
    for row in range(1, len(matrix)):
        for column in range(1, len(matrix[row])):
            cell = matrix[row][column]
            interactions = get_interactions(cell)
            reactives[row] = reactives.get(row, 0) + interactions
    
    # If count is True, return the counts of interactions instead of the matrix
    if count:
        data = [list(reactives.keys()), list(reactives.values())]
        for index in reactives.keys():
            data[0][index-1] = matrix[index][0]
        return data
    # Select rows based on the specified criteria
    elif thr_interactions is None and selected_items is None:
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)]
    elif thr_interactions is not None:
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True) if value >= thr_interactions]
    else:
        selected_items = selected_items if selected_items < len(matrix) else len(matrix)
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)[:selected_items]]

    # Create the selection matrix with the selected rows
    selection = [matrix[0]] + [matrix[row] for row in reactives]

    # Transpose the selection matrix back if axis was 'columns'
    if axis == 'columns':
        selection = transpose_matrix(matrix=selection)
    
    # Sort the selection by residue chain if residue_chain is True
    if residue_chain:
        selection = sort_by_residue(matrix=selection)

    return selection

def plot_matrix(matrix: list, plot_name: str, axis: str, label_x: str = None, label_y: str = "Number of intermolecular interactions", title: str = "Protein-drug interactions", stacked: bool = False, save: bool = False) -> None:
    
    """
    Plots a bar chart based on selected rows or columns of a matrix and saves it as a PNG file.

    Args:
        matrix (list of lists): The matrix containing interaction data.
        plot_name (str): The name of the plot to be saved (without extension).
        axis (str): Specifies whether to select rows ('rows') or columns ('columns').
        label_x (str, optional): Label for the X-axis. Defaults to "PDB complexes".
        label_y (str, optional): Label for the Y-axis. Defaults to "Number of intermolecular interactions".
        title (str, optional): Title of the plot. Defaults to "Protein-drug interactions".
        stacked (bool, optional): If True, creates a stacked bar chart. Defaults to False.
        save (bool, optional): If True, saves the plot as a PNG file. Defaults to False.

    Returns:
        None

    Example:
        If matrix = [['', 'file1', 'file2'], ['residue1', '1 |atom1|, 2 |atom2|']],
        plot_matrix(matrix, 'interaction_plot', 'rows', label_x='Residues', label_y='Interactions',
                    title='Interactions per Residue')
    """

    def get_interactions(cell: str) -> list:
        """
        Extracts and counts interactions from a cell string.

        Args:
            cell (str): The cell string containing interaction data.

        Returns:
            list: A list of interaction counts for each interaction type.
        """
        interactions = [0] * 8
        sections = cell.split("|")
        for index in range(1, len(sections), 2):
            interaction = int(sections[index - 1].replace(" ", "").replace(",", ""))
            interactions[interaction - 1] += len(sections[index].split(" "))
        return interactions

    def stack_reactives(matrix: list, axis: str) -> tuple[list, list]:
        """
        Accumulates interaction counts for rows or columns and returns the stacked data.

        Args:
            matrix (list of lists): The matrix containing interaction data.
            axis (str): Specifies whether to select rows ('rows') or columns ('columns').

        Returns:
            tuple: A tuple containing the stacked data and indices.
        """
        verify_dimensions(matrix=matrix)
        if axis == 'columns':
            matrix = transpose_matrix(matrix)
        
        reactives = {row: [0] * 8 for row in range(1, len(matrix))}
        indices = [matrix[row][0] for row in range(1, len(matrix))]
        
        for row in range(1, len(matrix)):
            for column in range(1, len(matrix[row])):
                cell = matrix[row][column]
                interactions = get_interactions(cell)
                for i in range(8):
                    reactives[row][i] += interactions[i]

        result_list = list(reactives.values())
        return result_list, indices

    # Create a new figure
    fig, ax = plt.subplots(num=plot_name, figsize=(12, 6))

    if stacked:
        bars = []
        data, indices = stack_reactives(matrix=matrix, axis=axis)
        labels = [
            "Hydrophobic", "Aromatic_Face/Face", "Aromatic_Edge/Face", "HBond_PROT", "HBond_LIG", "Ionic_PROT", "Ionic_LIG", "Other_Interactions"
        ]

        transposed_data = transpose_matrix(data)
        bottoms = [0] * len(indices)
        for index, group in enumerate(transposed_data):
            bars.append(ax.bar(indices, group, bottom=bottoms, label=labels[index]))
            bottoms = [i + j for i, j in zip(bottoms, group)]

        ax.set_xticks(range(len(indices)))
        ax.set_xticklabels(indices)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)

        max_y = max([sum(col) for col in data])
    else:
        data = sort_matrix(matrix, axis, count=True)
        ax.bar(data[0], data[1])
        
        max_y = max(data[1])
    
    ax.set_ylim(0, max_y * 1.1)

    ax.set_ylabel(label_y)
    if label_x is  None:
        residues_axis = get_residues_axis(matrix=matrix)
        label_x = "Interacting protein residues" if residues_axis == axis else "PDB complexes"
    ax.set_xlabel(label_x)
    ax.set_title(title)

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.xticks(rotation=90, ha='center')
    plt.tight_layout()

    # Add interactive cursors to display percentages for stacked bars
    if stacked:
        cursor = mplcursors.cursor(bars, hover=True)

        @cursor.connect("add")
        def on_add(sel):
            index = sel.index
            total = sum(transposed_data[i][index] for i in range(len(transposed_data)))
            percentages = [transposed_data[i][index] / total * 100 if total != 0 else 0 for i in range(len(transposed_data))]
            annotation_text = "\n".join([f"{labels[i]}: {percentages[i]:.1f}%" for i in range(len(labels))])
            sel.annotation.set_text(annotation_text)
            sel.annotation.get_bbox_patch().set(fc="white", alpha=0.9)  # Set background to white with 90% opacity

    if not save:
        plt.show()
    else:
        plt.savefig(os.path.join(saving_directory, plot_name + '.png'))
        plt.close(fig)  # Close the figure after saving to avoid display overlap

def filter_by_interaction(matrix: list, interactions: list) -> list:
    """
    Filters the matrix based on specified interactions.

    Args:
        matrix (list): The matrix to filter, represented as a list of lists.
        interactions (list): List of valid interactions (numbers 1 to 7).

    Returns:
        list: Filtered matrix with updated interaction information.

    Raises:
        ValueError: If matrix dimensions are invalid or if no desired interactions are found.
    """

    def validate_list(interactions: list) -> None:
        """
        Validates the interactions list to ensure it contains valid numbers (1 to 7) without duplicates.

        Args:
            interactions (list): List of integers representing interactions.

        Returns:
            None

        Raises:
            ValueError: If any number is outside the range 1 to 7 or if there are duplicates.
        """
        
        # Valid numbers are between 1 and 7
        valid_numbers = set(range(1, 8))

        # Check if all numbers in the list are valid
        for num in interactions:
            if num not in valid_numbers:
                raise ValueError(f"The number {num} is not valid. It should be a number from 1 to 7.")

        # Check for duplicates in the list
        if len(set(interactions)) != len(interactions):
            raise ValueError("The list should not contain duplicate numbers.")

    # Validate matrix dimensions
    verify_dimensions(matrix=matrix)

    # Ensure the interactions list is valid
    validate_list(interactions=interactions)

    filtered = copy.deepcopy(matrix)

    # Flag to track if any changes were made in the matrix
    changes = False

    # Iterate through each cell in the matrix
    for i in range(1, len(filtered)):
        for j in range(1, len(filtered[i])):
            cell = filtered[i][j]
            
            # Process non-empty cells
            if cell != '-':
                sections = cell.split(", ")
                cell = ""
                
                # Iterate through each section in the cell
                for section in sections:
                    # Check if the first number in the section is in the valid interactions
                    if int(section.split(" ")[0]) in interactions:
                        if cell == '':
                            cell = section
                        else:
                            cell += ', ' + section
                        changes = True
                
                # Update the matrix cell based on filtered sections
                if cell == '':
                    filtered[i][j] = '-'
                else:
                    filtered[i][j] = cell
    
    # If no changes were made, raise an error
    if not changes:
        raise ValueError("The matrix does not have any of the desired interactions.")

    return filtered

def get_residues_axis(matrix: list) -> str:
    """
    Determines whether the residues' axis is in the rows or columns of the matrix.

    Args:
        matrix (list of lists): The matrix containing interaction data.

    Returns:
        str: 'columns' if residues' axis is in the columns, 'rows' if residues' axis is in the rows.

    Raises:
        ValueError: If the residues' axis cannot be determined.
    """
    # Validate matrix dimensions
    verify_dimensions(matrix=matrix)
    
    # Count the spaces in specific positions to determine the axis
    count_0_1 = matrix[0][1].count(' ')
    count_1_0 = matrix[1][0].count(' ')

    # If the counts are equal, the axis cannot be determined
    if count_0_1 == count_1_0:
        raise ValueError("Cannot determine the residues' axis.")
    # If the count in the first row, second column is 1, the axis is 'columns'
    elif count_0_1 == 1:
        return 'columns'
    # If the count in the second row, first column is 1, the axis is 'rows'
    elif count_1_0 == 1:
        return 'rows'
    # If neither condition is met, the axis cannot be determined
    else:
        raise ValueError("Cannot determine the residues' axis.")
