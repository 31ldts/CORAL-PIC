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

##############
# Exceptions #
##############

class TypeMismatchException(Exception):
    def __init__(self, variable_name, expected_types, actual_type):
        expected_types_str = ", ".join([t.__name__ for t in expected_types])
        self.message = f"Variable '{variable_name}' has type {actual_type.__name__}, expected one of ({expected_types_str})."
        super().__init__(self.message)

class FileOrDirectoryNotFoundException(Exception):
    def __init__(self, path, message="File or directory not found"):
        self.path = path
        self.message = f"{message}: '{path}'"
        super().__init__(self.message)

class EmptyDirectoryException(Exception):
    def __init__(self, path):
        self.path = path
        self.message = f"Directory '{path}' is empty"
        super().__init__(self.message)

###################
# Private Methods #
###################

def _check_variable_types(variables, expected_types, variable_names):
    # Check that all lists have the same length
    if len(variables) != len(expected_types) or len(variables) != len(variable_names):
        raise ValueError("The lists of variables, expected types, and variable names must all have the same length.")
    
    for i, variable in enumerate(variables):
        expected_type = expected_types[i]
        variable_name = variable_names[i]
        
        # Check if the variable's type matches one of the expected types
        if not isinstance(variable, expected_type if isinstance(expected_type, tuple) else (expected_type,)):
            actual_type = type(variable)
            raise TypeMismatchException(variable_name, expected_type if isinstance(expected_type, tuple) else (expected_type,), actual_type)

def _get_residues_axis(matrix: list) -> str:
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
    _verify_dimensions(matrix=matrix)
    
    # Count the spaces in specific positions to determine the axis
    count_0_1 = matrix[0][1].count('(')
    count_1_0 = matrix[1][0].count('(')

    # If the counts are equal, the axis cannot be determined
    if count_0_1 == count_1_0:
        raise ValueError("Cannot determine the residues' axis.")
    # If the count in the first row, second column is 1, the axis is 'columns'
    elif count_0_1 == 1:
        return 'rows'
    # If the count in the second row, first column is 1, the axis is 'rows'
    elif count_1_0 == 1:
        return 'columns'
    # If neither condition is met, the axis cannot be determined
    else:
        raise ValueError("Cannot determine the residues' axis.")

def _verify_dimensions(matrix: list):
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

def _remove_void(matrix: list):
    def _remove_void_rows(matrix: list) -> list:
        changes = 0
        for row in range(1, len(matrix)):
            if all(column == '-' for column in matrix[row - changes][1:]):
                matrix.pop(row - changes)
                changes += 1
        return matrix

    # Remove empty rows
    matrix = _remove_void_rows(matrix)
    
    # Transpose the matrix, remove empty columns (which are now rows)
    matrix = transpose_matrix(matrix=matrix)
    matrix = _remove_void_rows(matrix)
    
    # Transpose back to restore original format
    matrix = transpose_matrix(matrix=matrix)
    
    return matrix

##################
# Public Methods #
##################

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

def transpose_matrix(matrix: list, save: str = None) -> list:
    """
    Transposes the given matrix (list of lists).

    Args:
        matrix (list): A 2D list (matrix) to transpose.
        save (str, optional): If provided, saves the transposed matrix to the specified file path.

    Returns:
        list: The transposed matrix.
    """
    # Check types of matrix and save parameters
    _check_variable_types(
        variables=[matrix, save], 
        expected_types=[list, (str, None.__class__)], 
        variable_names=['matrix', 'save']
    )

    # Ensure matrix has consistent dimensions
    _verify_dimensions(matrix=matrix)

    # Transpose the matrix using list comprehension
    transposed = [[row[i] for row in matrix] for i in range(len(matrix[0]))]

    # Save the transposed matrix if a file path is provided
    if save:
        save_matrix(matrix=transposed, filename=save)

    return transposed

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

    _verify_dimensions(matrix=matrix)

    # Create the CSV file
    with open(os.path.join(saving_directory, filename), 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        
        # Write the rows of the matrix to the CSV file
        for row in matrix:
            csv_writer.writerow(row)

def analyze_files(directory: str, activity_file: str = None, protein: bool = True, ligand: bool = True, subunit: bool = False, save: str = None) -> list:
    """
    Analyzes interaction data files in a specified directory, categorizing interactions
    based on protein and ligand atoms involved.

    Args:
        directory (str): Path to the directory containing interaction data files.
        activity_file (str, optional): Path to the activity file (CSV) for labeling data.
        protein (bool, optional): Include protein atoms in the analysis if True.
        ligand (bool, optional): Include ligand atoms in the analysis if True.
        subunit (bool, optional): Differentiate between subunits if True.
        save (str, optional): Path to save the resulting matrix (optional).

    Returns:
        list: A matrix categorizing interactions between residues and files.

    Raises:
        FileNotFoundError: If the specified directory or activity file does not exist.
        EmptyDirectoryException: If the specified directory is empty.
    """
    
    def label_matrix(matrix: list, rows: list, columns: list, activity_file: str) -> list:
        """
        Adds headers to the interaction matrix with residue names and file names.

        Args:
            matrix (list): 2D list representing interaction data.
            rows (list): List of residue names for row labeling.
            columns (list): List of file names for column labeling.
            activity_file (str): Path to the activity file for activity-based labeling.

        Returns:
            list: The labeled matrix.
        """
        rows = [row.replace("\t", "") for row in rows]
        columns = [""] + columns[:]  # Add an empty string at the start for residue names
        
        if activity_file:
            if not os.path.isfile(activity_file):
                raise FileNotFoundError(f"The file '{activity_file}' does not exist.")
            
            # Read activity data into a dictionary
            data_dict = {}
            with open(activity_file, newline='') as csvfile:
                csvreader = csv.reader(csvfile)
                try:
                    next(csvreader)  # Skip header
                except:
                    raise ValueError(f"The CSV file '{activity_file}' is missing a header.")
                for key, value in csvreader:
                    data_dict[key] = str(round(float(value), 3))

            if not data_dict:
                raise ValueError(f"The CSV file '{activity_file}' must contain at least one row of data.")

            # Update column names with activity data
            for i in range(1, len(columns)):
                drug_name = columns[i]
                columns[i] = f"{drug_name} ({data_dict.get(drug_name, '0')})"

        # Insert headers into the matrix
        matrix.insert(0, columns)
        for i, row in enumerate(matrix[1:], start=1):
            row.insert(0, rows[i-1])

        return matrix

    def modify_cell(text: str, interaction: str, atoms: str) -> str:
        """
        Updates cell content by adding interaction type and involved atoms.

        Args:
            text (str): Current cell content.
            interaction (str): Interaction type.
            atoms (str): Atoms involved in the interaction.

        Returns:
            str: Updated cell content.
        """
        interaction_map = {
            "Hydrophobic": '1', "Aromatic_Face/Face": '2', "Aromatic_Edge/Face": '3',
            "HBond_PROT": '4', "HBond_LIG": '5', "Ionic_PROT": '6', "Ionic_LIG": '7'
        }
        interaction_code = interaction_map.get(interaction, '8')

        if text == "-":
            return f"{interaction_code} |{atoms}|"
        
        content = text.replace(", ", "").split("|")[:-1]
        exists = False
        cell = ''

        for index, segment in enumerate(content):
            if index % 2 == 0 and interaction_code == segment.strip():
                content[index + 1] += f" {atoms}"
                exists = True
                break

        cell = ', '.join(f"{content[i]}|{content[i+1]}|" for i in range(0, len(content), 2))
        
        if not exists:
            cell += f", {interaction_code} |{atoms}|"
        
        return cell

    def read_txt_file(file_name: str) -> list:
        """
        Reads a text file and returns its content as a list of lines.

        Args:
            file_name (str): The name of the text file to read.

        Returns:
            list: List of lines from the file.

        Raises:
            FileNotFoundError: If the file does not exist.
            Exception: For any unexpected errors during reading.
        """
        try:
            with open(file_name, 'r') as file:
                return [line.strip() for line in file.readlines()]
        except FileNotFoundError:
            print(f"Error: The file '{file_name}' does not exist.")
            raise
        except Exception as e:
            print(f"Error: An unexpected error occurred while reading '{file_name}': {e}")
            raise

    def adjust_subunits(matrix: list, subunits: list) -> list:
        """
        Adjusts matrix data to reflect residue counts for subunits.

        Args:
            matrix (list): The interaction matrix.
            subunits (list): List of subunit identifiers.

        Returns:
            list: Matrix with adjusted subunit information.
        """
        def reduce_strings_in_even_indices(atoms: str, subunits: int) -> str:
            pattern = r'[() ]'
            result = list(filter(None, re.split(pattern, atoms)))
            even_index_strings = [result[i] for i in range(0, len(result), 2)]
            string_count = Counter(even_index_strings)

            text = ''
            added_strings = set()
            for i in range(len(result)):
                string = result[i]
                if i % 2 == 0:
                    if string_count[string] >= subunits:
                        if string not in added_strings:
                            text += f"{string} "
                            added_strings.add(string)
                    else:
                        text += f"{string}({result[i+1]}) "
            return '|' + text.strip() + '|'
        
        for row in range(len(matrix)):
            for column in range(len(matrix[row])):
                cell = matrix[row][column]
                if cell != '-':
                    sections = cell.split("|")
                    text = ''.join(
                        sections[i-1] + reduce_strings_in_even_indices(sections[i], len(subunits))
                        for i in range(1, len(sections), 2)
                    )
                    matrix[row][column] = text
        return matrix

    # Validate input types
    _check_variable_types(
        variables=[directory, activity_file, protein, ligand, subunit, save],
        expected_types=[str, (str, None.__class__), bool, bool, bool, (str, None.__class__)],
        variable_names=['directory', 'activity_file', 'protein', 'ligand', 'subunit', 'save']
    )

    # Check if the directory exists
    if not os.path.exists(directory):
        raise FileOrDirectoryNotFoundException(path=directory)
    
    files = os.listdir(directory)
    if not files:
        raise EmptyDirectoryException(path=directory)

    matrix = []
    aa = {}
    cont = 0
    subunits_set = set()
    
    # Analyze each file in the directory
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
                    
                    atoms = f"{elements[1].strip()}-{elements[4].strip()}" if protein and ligand else elements[1].strip() if protein else elements[4].strip()
                    if not subunit:
                        atoms += f"({sections[1]})"

                    if residue not in aa:
                        aa[residue] = cont
                        cont += 1
                    column = aa[residue]

                    # Ensure matrix size and modify cell
                    if len(matrix) <= column:
                        matrix.append(["-"] * len(files))

                    matrix[column][index] = modify_cell(matrix[column][index], interaction, atoms)

        files[index] = file.replace(".txt", "")

    if not subunit:
        matrix = adjust_subunits(matrix, list(subunits_set))
    matrix = label_matrix(matrix, rows=list(aa.keys()), columns=files, activity_file=activity_file)

    # Save the matrix if specified
    if save:
        save_matrix(matrix=matrix, filename=save)

    return matrix

def sort_matrix(matrix: list, axis: str = 'rows', thr_interactions: int = None, thr_activity: float = None, selected_items: int = None, count: bool = False, residue_chain: bool = False, save: str = None) -> list:
    """
    Sorts and selects reactive rows or columns from a matrix based on interactions.

    Args:
        matrix (list of lists): The matrix containing interaction data.
        axis (str, optional): Specifies whether to select rows ('rows') or columns ('columns'). Defaults to 'rows'.
        thr_interactions (int, optional): Minimum number of interactions to select a row/column.
        thr_activity (float, optional): Activity threshold to select rows/columns if activity values are given.
        selected_items (int, optional): Number of top rows/columns to select based on interactions.
        count (bool, optional): If True, returns the count of interactions instead of the matrix.
        residue_chain (bool, optional): If True, sorts the resulting matrix based on residue order in the chain.
        save (str, optional): File path to save the resulting matrix.

    Returns:
        list of lists: The selected rows or columns based on the specified criteria.

    Raises:
        ValueError: If both `thr_interactions` and `selected_items` are provided simultaneously.
        ValueError: If the matrix dimensions are insufficient.
    """

    def get_interactions(cell: str) -> int:
        """
        Counts the number of interactions in a cell formatted with interaction data.

        Args:
            cell (str): Cell containing interaction data, formatted with '|' separating values.

        Returns:
            int: Total number of interactions in the cell.
        """
        interactions = 0
        sections = cell.split("|")
        for index in range(1, len(sections), 2):
            interactions += len(sections[index].split(" "))
        return interactions
    
    def sort_by_residue(matrix: list) -> list:
        """
        Sorts the matrix based on residue indices in the first column.

        Args:
            matrix (list of lists): The matrix to be sorted.

        Returns:
            list of lists: The matrix sorted by residue indices.
        """
        # Validate matrix dimensions
        _verify_dimensions(matrix=matrix)

        # Determine whether to sort by rows or columns
        axis = _get_residues_axis(matrix=matrix)
        
        # If sorting by columns, transpose the matrix first
        if axis == 'columns':
            matrix = transpose_matrix(matrix=matrix)

        # Separate the header from the data rows
        header = matrix[0]
        data_rows = matrix[1:]

        # Sort the data rows based on residue indices
        sorted_data_rows = sorted(data_rows, key=lambda row: int(row[0].replace(" ", "")[3:]))

        # Combine the header with the sorted data rows
        sorted_matrix = [header] + sorted_data_rows

        # If sorting was by columns, transpose the sorted matrix back
        if axis == 'columns':
            sorted_matrix = transpose_matrix(matrix=sorted_matrix)

        return sorted_matrix

    # Check variable types to ensure correct input
    _check_variable_types(
        variables=[matrix, axis, thr_interactions, thr_activity, selected_items, count, residue_chain, save], 
        expected_types=[list, str, (int, None.__class__), (float, None.__class__), (int, None.__class__), bool, bool, (str, None.__class__)], 
        variable_names=['matrix', 'axis', 'thr_interactions', 'thr_activity', 'selected_items', 'count', 'residue_chain', 'save']
    )

    # Validate matrix dimensions
    _verify_dimensions(matrix=matrix)

    # Raise an error if `thr_interactions` and `selected_items` are provided simultaneously
    if thr_interactions is not None and selected_items is not None:
        raise ValueError("Cannot select by both 'thr_interactions' and 'selected_items' at the same time.")
    if thr_interactions is not None and thr_activity is not None:
        raise ValueError("Cannot select by both 'thr_interactions' and 'thr_activity' at the same time.")
    if thr_activity is not None and selected_items is not None:
        raise ValueError("Cannot select by both 'thr_activity' and 'selected_items' at the same time.")
    
    # Transpose the matrix if operating on columns
    if axis == 'columns':
        matrix = transpose_matrix(matrix=matrix)
    
    # Initialize a dictionary to store interaction counts per row/column
    reactives = {}
    for row in range(1, len(matrix)):
        for column in range(1, len(matrix[row])):
            cell = matrix[row][column]
            interactions = get_interactions(cell)
            reactives[row] = reactives.get(row, 0) + interactions
    
    # If `count` is True, return the interaction counts
    if count:
        data = [list(reactives.keys()), list(reactives.values())]
        for index in reactives.keys():
            original_value = matrix[index][0]
            split_value = original_value.split('_')[0].strip()  # Splits and removes spaces
            data[0][index-1] = split_value
        return data
    
    # Select rows/columns based on the provided criteria
    elif thr_interactions is not None:
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True) if value >= thr_interactions]
    elif _get_residues_axis(matrix) == "columns" and thr_activity is not None:
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: float(matrix[item[0]][0].split(" (")[1].replace(")", "")), reverse=True) if float(matrix[key][0].split(" (")[1].replace(")", "")) >= thr_activity]
    elif selected_items:
        selected_items = min(selected_items, len(matrix))
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)[:selected_items]]
    else:
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)]
    
    # Create the selection matrix with the chosen rows/columns
    selection = [matrix[0]] + [matrix[row] for row in reactives]

    # Transpose the selection back if it was initially transposed for columns
    if axis == 'columns':
        selection = transpose_matrix(matrix=selection)
    
    # Sort the selection by residue chain if specified
    if residue_chain:
        selection = sort_by_residue(matrix=selection)

    # Save the matrix if a save path is provided
    if save:
        save_matrix(matrix=selection, filename=save)

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
        _verify_dimensions(matrix=matrix)
        if axis == 'columns':
            matrix = transpose_matrix(matrix)
        
        reactives = {row: [0] * 8 for row in range(1, len(matrix))}
        indices = [matrix[row][0].split('_')[0].strip() for row in range(1, len(matrix))]
        
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
        residues_axis = _get_residues_axis(matrix=matrix)
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

def filter_by_interaction(matrix: list, interactions: list, save: str = None) -> list:
    """
    Filters a matrix based on specified interaction types.

    Args:
        matrix (list): The matrix to filter, represented as a list of lists.
        interactions (list): List of valid interaction types (numbers 1 to 7) to retain in the matrix.
        save (str, optional): Path to save the filtered matrix. Defaults to None.

    Returns:
        list: The filtered matrix with only the specified interactions retained.

    Raises:
        ValueError: If the matrix dimensions are invalid, or if no desired interactions are found.
    """

    def validate_list(interactions: list) -> None:
        """
        Validates the interaction list to ensure it contains unique numbers between 1 and 7.

        Args:
            interactions (list): List of integers representing interaction types.

        Raises:
            ValueError: If any number is outside the range of 1 to 7 or if there are duplicates.
        """
        # Valid interactions are numbers from 1 to 7
        valid_numbers = set(range(1, 8))

        # Check if all numbers in the list are within the valid range
        for num in interactions:
            if num not in valid_numbers:
                raise ValueError(f"Invalid interaction: {num}. Must be a number between 1 and 7.")

        # Ensure the list contains no duplicate values
        if len(set(interactions)) != len(interactions):
            raise ValueError("The interaction list must not contain duplicates.")

    # Validate types of the matrix, interactions, and save parameters
    _check_variable_types(
        variables=[matrix, interactions, save], 
        expected_types=[list, list, (str, None.__class__)], 
        variable_names=['matrix', 'interactions', 'save']
    )

    # Validate that the matrix has appropriate dimensions
    _verify_dimensions(matrix=matrix)

    # Validate that the interaction list contains valid values
    validate_list(interactions=interactions)

    # Create a deep copy of the matrix to avoid modifying the original
    filtered = copy.deepcopy(matrix)

    # Track whether any interactions were filtered
    changes = False

    # Iterate through each cell in the matrix (skipping the header row/column)
    for i in range(1, len(filtered)):
        for j in range(1, len(filtered[i])):
            cell = filtered[i][j]
            
            # If the cell is not empty ('-'), process it
            if cell != '-':
                sections = cell.split(", ")
                cell = ""
                
                # Iterate through the sections in the cell to filter valid interactions
                for section in sections:
                    # Check if the first number in the section is in the valid interaction list
                    if int(section.split(" ")[0]) in interactions:
                        # Add the section to the cell if it contains a valid interaction
                        if cell == '':
                            cell = section
                        else:
                            cell += ', ' + section
                        changes = True
                
                # Update the cell with the filtered sections or set it to '-' if empty
                filtered[i][j] = cell if cell != '' else '-'
    
    # If no changes were made, raise an error indicating no matching interactions were found
    if not changes:
        raise ValueError("No matching interactions were found in the matrix.")

    # Save the filtered matrix to a file if a save path is provided
    if save:
        save_matrix(matrix=filtered, filename=save)

    return filtered

def filter_by_subunit(matrix: list, subunits: list, save: str = None) -> list:
    """
    Filters a matrix based on specified subunits, removing rows or columns
    that do not contain the desired subunits.

    Args:
        matrix (list): The matrix to filter, represented as a list of lists.
        subunits (list): List of valid subunits to retain in the matrix.
        save (str, optional): Path to save the filtered matrix. Defaults to None.

    Returns:
        list: The filtered matrix with only the specified subunits retained.

    Raises:
        ValueError: If the matrix dimensions are invalid or if no desired subunits are found.
    """
    
    def get_subunits_location(matrix: list) -> str:
        """
        Determines the location of subunits in the matrix (residues or interactions).

        Args:
            matrix (list): The matrix being analyzed.

        Returns:
            str: 'residues' if the first column indicates residue data, otherwise 'interactions'.
        """
        return 'residues' if len(matrix[1][0].split('-')) == 2 else 'interactions'
    
    # Check types of the matrix, subunits, and save parameters
    _check_variable_types(
        variables=[matrix, subunits, save], 
        expected_types=[list, list, (str, None.__class__)], 
        variable_names=['matrix', 'subunits', 'save']
    )

    # Validate the dimensions of the matrix
    _verify_dimensions(matrix=matrix)

    # Create a deep copy of the matrix to avoid modifying the original
    filtered = copy.deepcopy(matrix)

    # Determine the axis of residues in the matrix
    axis = _get_residues_axis(matrix=filtered)

    # Transpose the matrix if the axis is columns
    if axis == "columns":
        filtered = transpose_matrix(matrix=filtered)
    
    # Determine whether the matrix contains residues or interactions
    subunitsLocation = get_subunits_location(matrix=filtered)

    # Initialize change tracking variables
    changes = 0

    # Filter based on residue locations
    if subunitsLocation == 'residues':
        for index in range(1, len(filtered)):
            sections = filtered[index - changes][0].split("-")
            # Remove rows without valid sections or subunits
            if len(sections) != 2 or sections[1] not in subunits:
                filtered.pop(index - changes)
                changes += 1
    else:
        # Iterate through each cell in the matrix to filter interactions
        for i in range(1, len(filtered)):
            for j in range(1, len(filtered[i])):
                cell = filtered[i][j]
                
                # Process non-empty cells
                if cell != '-':
                    sections = cell.split(", ")
                    cell = ""
                    
                    # Iterate through each section in the cell
                    for section in sections:
                        separators = section.split('|')[:-1]

                        # Filter out unwanted interactions
                        for index in range(1, len(separators)):
                            if index % 2 != 0:
                                interactions = separators[index].split(' ')
                                subchanges = 0
                                
                                # Remove interactions not in the valid subunits
                                for interaction in range(len(interactions)):
                                    if interactions[interaction - subchanges][-2] not in subunits:
                                        changes += 1
                                        interactions.pop(interaction - subchanges)
                                        subchanges += 1

                                # Rebuild the cell if there are valid interactions
                                if interactions:
                                    cell += separators[index - 1] + '|'
                                    cell += ' '.join(interactions) + ' |, '
                                    
                    # Update the cell with filtered interactions
                    cell = cell[:-2]  # Remove trailing comma and space
                    filtered[i][j] = cell if cell else '-'
    
    # Raise an error if no changes were made (no desired subunits found)
    if changes == 0:
        raise ValueError("The matrix does not contain any of the desired subunits.")

    # Transpose the matrix back if it was originally in columns
    if axis == "columns":
        filtered = transpose_matrix(matrix=filtered)

    # Remove any empty rows or columns
    filtered = _remove_void(matrix=filtered)

    # Save the filtered matrix to a file if a save path is provided
    if save:
        save_matrix(matrix=filtered, filename=save)

    return filtered