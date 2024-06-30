import csv
import os
from matplotlib import pyplot as plt
from files_mediator import *
from typing import Union

saving_directory = os.getcwd()

"""
Changes the saving directory to a subdirectory within the current directory.

Args:
    path (str): The name of the subdirectory to change to.

Returns:
    None

Raises:
    ValueError: If the subdirectory does not exist.
"""
def change_directory(path: str) -> None:

    global saving_directory

    # Get the full path of the new saving directory
    saving_directory = os.path.join(os.getcwd(), path)

    # Check if the new saving directory exists
    if not os.path.exists(saving_directory):
        raise ValueError("The saving directory must exist inside the project.")

"""
Transposes the given matrix.

Args:
    matrix (list of lists): The matrix to transpose.

Returns:
    list of lists: The transposed matrix.
"""
def transpose_matrix(matrix: list) -> list:
    verify_dimensions(matrix=matrix)

    # Use list comprehension to transpose the matrix
    return [[row[i] for row in matrix] for i in range(len(matrix[0]))]

"""
Saves the matrix to a CSV file in the specified directory.

Args:
    matrix (list of lists): The matrix to be saved.
    filename (str): The name of the file to save the matrix in.

Returns:
    None
"""
def save_matrix(matrix: list, filename: str) -> None:
    global saving_directory

    verify_dimensions(matrix=matrix)

    # Create the CSV file
    with open(os.path.join(saving_directory, filename), 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        
        # Write the rows of the matrix to the CSV file
        for row in matrix:
            csv_writer.writerow(row)

"""
Analyzes interaction data files in a specified directory, categorizing interactions
based on protein and ligand atoms involved.

Args:
    directory (str): The directory path containing interaction data files.
    protein (bool): Flag indicating whether to include protein atoms in analysis.
    ligand (bool): Flag indicating whether to include ligand atoms in analysis.

Returns:
    list: A matrix representing categorized interactions between residues and files.

Raises:
    FileNotFoundError: If the specified directory does not exist or cannot be accessed.
"""
def analyze_files(directory: str, protein: bool=True, ligand: bool=True, subunit: bool=True) -> list:

    """
    Labels the matrix with residue names and file names for clarity.

    Args:
        matrix (list): The 2D list representing interaction data.
        rows (list): List of residue names.
        columns (list): List of file names (drugs).

    Returns:
        list: The labeled matrix with headers.
    """
    def label_matrix(matrix: list, rows: list, columns: list) -> list:
        for index, row in enumerate(rows):
            matrix[index].insert(0, row.replace("\t", ""))
        columns.insert(0, "")
        matrix.insert(0, columns)
        return matrix
    
    """
    Modifies the cell content based on interaction type and atoms involved.

    Args:
        text (str): The current content of the cell.
        interaction (str): The type of interaction.
        atoms (str): The atoms involved in the interaction.

    Returns:
        str: The modified cell content with updated interaction details.
    """
    def modify_cell(text: str, interaction: str, atoms: str) -> str:
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
        interaction = interaction_map.get(interaction, interaction)

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

    files = os.listdir(directory)
    matrix = []
    aa = {}
    cont = 0
    for index, file in enumerate(files):
        file_path = os.path.join(directory, file)
        if os.path.isfile(file_path):
            content = read_txt_file(file_path)  # Assume read_txt_file is defined elsewhere
            for line in content:
                elements = line.split("|")
                if len(elements) == 10:
                    interaction = elements[0].strip().replace("\t", "")
                    residue = elements[3].strip().replace("\t", "")
                    if not subunit:
                        residue = residue.split("-")[0]
                    
                    # Determine which atoms to include based on flags
                    if protein and ligand:
                        atoms = f"{elements[1].strip().replace(chr(9), '')}-{elements[4].strip().replace(chr(9), '')}"
                    elif protein:
                        atoms = elements[1].strip().replace("\t", "")
                    elif ligand:
                        atoms = elements[4].strip().replace("\t", "")
                    else:
                        atoms = ""

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

    return label_matrix(matrix=matrix, rows=list(aa.keys()), columns=files)

"""
Counts the number of interactions in a cell of formatted interaction data.

Args:
    cell (str): The cell containing interaction data formatted with '|' separators.

Returns:
    int: The total number of interactions found in the cell.

Example:
    If cell = "1 |atom1|, 2 |atom2|, 3 |atom3|", get_interactions(cell) returns 3.
"""
def get_interactions(cell: str) -> int:
    interactions = 0
    sections = cell.split("|")
    
    # Iterate over every other section starting from the second one
    for index in range(1, len(sections), 2):
        # Split the section by spaces and count the number of elements
        interactions += len(sections[index].split(" "))
    
    return interactions

"""
Verifies the dimensions of the matrix to ensure it contains at least 2 rows and 2 columns.

Args:
    matrix (list of lists): The matrix to be verified.

Raises:
    ValueError: If the matrix has fewer than 2 rows or any row has fewer than 2 columns.

Example:
    If matrix = [['-', '1'], ['2', '-']], verify_dimensions(matrix) will pass.
"""
def verify_dimensions(matrix: list):
    if len(matrix) < 2 or any(len(row) < 2 for row in matrix):
        raise ValueError("There are not interactions on the matrix.")

"""
Transposes a given matrix.

Args:
    matrix (list of lists): The matrix to be transposed.

Returns:
    list of lists: The transposed matrix.

Example:
    If matrix = [[1, 2, 3], [4, 5, 6]], transpose_matrix(matrix) returns [[1, 4], [2, 5], [3, 6]].
"""
def transpose_matrix(matrix: list):
    return [[row[i] for row in matrix] for i in range(len(matrix[0]))]

"""
Selects reactive rows or columns from a matrix based on interactions.

Args:
    matrix (list of lists): The matrix containing interaction data.
    axis (str): Specifies whether to select rows ('rows') or columns ('columns').
    threshold (int, optional): Minimum number of interactions to select a row/column.
    selectedItems (int, optional): Number of top rows/columns to select based on interactions.

Returns:
    list of lists: Selected rows or columns from the matrix based on the specified criteria.

Raises:
    ValueError: If both threshold and selectedItems are provided simultaneously.
    ValueError: If the matrix dimensions are insufficient.

Example:
    If matrix = [['', 'file1', 'file2'], ['residue1', '1 |atom1|, 2 |atom2|', '3 |atom3|, 4 |atom4|']], 
    select_reactives(matrix, 'rows', threshold=2) returns [['', 'file1', 'file2'], ['residue1', '1 |atom1|, 2 |atom2|']].

"""
def select_reactives(matrix: list, axis: str, threshold: int=None, selectedItems:int =None) -> list:
    verify_dimensions(matrix=matrix)

    if threshold is not None and selectedItems is not None:
        raise ValueError("You cannot select by 'threshold' and by 'selectedItems' at the same time.")

    if axis == 'columns':
        matrix = transpose_matrix(matrix=matrix)
    
    reactives = {}
    for row in range(1, len(matrix)):
        for column in range(1, len(matrix[row])):
            cell = matrix[row][column]
            interactions = get_interactions(cell)
            reactives[row] = reactives.get(row, 0) + interactions
    
    if threshold is None and selectedItems is None:
        data = [list(reactives.keys()), list(reactives.values())]
        for index in reactives.keys():
            data[0][index-1] = matrix[index][0]
        return data
    elif threshold is not None:
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True) if value >= threshold]
    else:
        selectedItems = selectedItems if selectedItems < len(matrix) else len(matrix)
        reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)[:selectedItems]]

    selection = [matrix[0]] + [matrix[row] for row in reactives]

    if axis == 'columns':
        selection = transpose_matrix(matrix=selection)

    return selection

"""
Plots a bar chart based on selected rows or columns of a matrix and saves it as a PNG file.

Args:
    matrix (list of lists): The matrix containing interaction data.
    plotName (str): The name of the plot to be saved (without extension).
    axis (str): Specifies whether to select rows ('rows') or columns ('columns').
    labelX (str, optional): Label for the X-axis.
    labelY (str, optional): Label for the Y-axis.
    title (str, optional): Title of the plot.

Returns:
    None

Example:
    If matrix = [['', 'file1', 'file2'], ['residue1', '1 |atom1|, 2 |atom2|']], 
    plot_matrix(matrix, 'interaction_plot', 'rows', labelX='Residues', labelY='Interactions', 
                title='Interactions per Residue')

"""
def plot_matrix(matrix: list, plotName: str, axis:str, labelX: str=None, labelY: str=None, title: str=None) -> None:
    global saving_directory
    data = select_reactives(matrix, axis)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_ylabel(labelY)
    ax.set_xlabel(labelX)
    ax.set_title(title)

    ax.bar(data[0], data[1])
    
    # Rotate X-axis labels to prevent overlap
    plt.xticks(rotation=90, ha='center')
    
    # Adjust layout to prevent labels from being cut off
    plt.tight_layout()
    
    # Save the plot as a PNG file
    plt.savefig(os.path.join(saving_directory, plotName + '.png'))
