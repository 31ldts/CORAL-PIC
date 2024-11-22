import csv
import os
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import mplcursors
import copy
import re
import json

###########
# Globals #
###########

# Labels for interaction types
INTERACTION_LABELS = [
    "Hydrophobic", "Aromatic_Face/Face", "Aromatic_Edge/Face", "HBond_PROT", "HBond_LIG", 
    "Ionic_PROT", "Ionic_LIG", "Metal Acceptor", "Pi/Cation", "Other_Interactions"
]

# List of colors
COLORS = [
    "#ff6384", "#36a2eb", "#ffce56", "#4bc0c0", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b"
]

# List of predictors
PREDICTORS = [
    "ichem", "arpeggio"
]

# Arpeggio interaction entities
ARPEGGIO_INT_ENT = [
    "INTER"
]

ARPEGGIO_CONT = [
    "covalent", "hbond", "aromatic", "hydrophobic", "polar", "ionic", "xbond", "metal", 
    "carbonyl", "CARBONPI", "CATIONPI", "DONORPI", "HALOGENPI", "METSULPHURPI", 
    "AMIDEAMIDE", "AMIDERING"
]

ARPEGGIO_TYPE = [
    "plane-plane"
]

ARPEGGIO_COLORS = [
    "#ff6384", "#36a2eb", "#ffce56", "#4bc0c0", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ff9da7", "#1a9ceb", "#fdca52"
]


# Constants for delimiters
SAME_DELIM = ', '       # Delimiter used to separate interactions of the same type.
DIFF_DELIM = '; '       # Delimiter used to separate interactions of different types.
GROUP_DELIM = '|'        # Delimiter used to group interactions of the same type.

# Constants for cell values
EMPTY_CELL = ''         # Represents an empty cell value from the beginning.
EMPTY_DASH_CELL = '-'   # Represents a cell that is considered empty due to filtering.

###########
# Lambdas #
###########

is_not_empty_or_dash = lambda cell: not (cell == EMPTY_DASH_CELL or cell == EMPTY_CELL)


saving_directory = os.getcwd()

class AnalyzeInteractions:

    def __init__(self):
        """Initializes the class with the current working directory and default settings."""
        self.saving_directory = os.getcwd()  # Set the default saving directory
        self.interaction_labels = INTERACTION_LABELS  # Default interaction labels
        self.colors = COLORS  # Default color configuration

    ###################
    # Private Methods #
    ###################

    def _check_variable_types(self, variables: list, expected_types: list, variable_names: list[str]):
        """
        Check if variables match their expected types.

        Args:
            variables (list): List of variables to check.
            expected_types (list): List of expected types (can include tuples).
            variable_names (list[str]): List of variable names for error messages.

        Raises:
            ValueError: If the lengths of input lists don't match.
            TypeMismatchException: If a variable type doesn't match the expected one.
        """
        if len(variables) != len(expected_types) or len(variables) != len(variable_names):
            raise ValueError("The lists of variables, expected types, and variable names must all have the same length.")

        for i, variable in enumerate(variables):
            expected_type = expected_types[i]
            variable_name = variable_names[i]

            if not isinstance(variable, expected_type if isinstance(expected_type, tuple) else (expected_type,)):
                actual_type = type(variable)
                raise TypeMismatchException(variable_name, expected_type if isinstance(expected_type, tuple) else (expected_type,), actual_type)

    def _get_residues_axis(self, matrix: list[list[str]]) -> str:
        """
        Determine whether residues' axis is in rows or columns.

        Args:
            matrix (list[list[str]]): Matrix with interaction data.

        Returns:
            str: 'columns' if residues are in columns, 'rows' otherwise.

        Raises:
            ValueError: If the axis cannot be determined.
        """
        self._verify_dimensions(matrix=matrix)

        # Check '(' occurrence in specific cells to determine axis
        count_0_1 = matrix[0][1].count('(')
        count_1_0 = matrix[1][0].count('(')

        if count_0_1 == count_1_0:
            raise ValueError("Cannot determine the residues' axis.")
        elif count_0_1 == 1:
            return 'rows'
        elif count_1_0 == 1:
            return 'columns'
        else:
            raise ValueError("Cannot determine the residues' axis.")

    def _verify_dimensions(self, matrix: list[list[str]]) -> None:
        """
        Verify matrix dimensions to ensure at least 2 rows and 2 columns.

        Args:
            matrix (list[list[str]]): Matrix to verify.
        
        Returns:
            None

        Raises:
            ValueError: If the matrix is too small or any row is too short.
        """
        if len(matrix) < 2 or any(len(row) < 2 for row in matrix):
            raise ValueError("The matrix must have at least 2 rows and 2 columns.")

    #################################
    # Public Methods: Configuration #
    #################################

    def change_directory(
            self, 
            path: str
            ) -> None:
        """
        Changes the saving directory to a specified subdirectory within the project.

        Args:
            path (str): Name of the subdirectory to switch to.

        Returns:
            None

        Raises:
            ValueError: If the subdirectory does not exist within the project.
        """
        
        self._check_variable_types(
            variables=[path],
            expected_types=[str],
            variable_names=['path']
        )

        # Construct the new path
        new_path = os.path.join(os.getcwd(), path)

        # Verify the new path exists
        if not os.path.exists(new_path):
            raise ValueError("The specified directory must exist inside the project.")

        # Update the saving directory
        self.saving_directory = new_path

    def set_plot_config(
            self, 
            interactions: list[str] = None, 
            colors: list[str] = None, 
            reset: bool = False,
            mode: str = None
            ) -> None:
        """
        Updates interaction labels and colors or resets them to default values.

        Args:
            interactions (list[str], optional): List of interaction labels to update.
            colors (list[str], optional): List of colors in hexadecimal format.
            reset (bool, optional): If True, resets the configuration to default values.

        Returns:
            None

        Raises:
            InvalidColorException: If any color in the provided list is not a valid hexadecimal value.
        """

        def is_valid_hex_color(color: str) -> bool:
            """
            Validates if a string is a valid hexadecimal color.

            Args:
                color (str): Color string to validate.

            Returns:
                bool: True if the color is a valid hex value, False otherwise.
            """
            return bool(re.match(r'^#([0-9a-fA-F]{6}|[0-9a-fA-F]{3})$', color))

        def reset_configuration() -> None:
            """
            Resets interaction labels and colors to their default values.

            Returns:
                None
            """
            self.interaction_labels = INTERACTION_LABELS
            self.colors = COLORS

        self._check_variable_types(
            variables=[interactions, colors, reset, mode],
            expected_types=[(list, None.__class__), (list, None.__class__), bool, (str, None.__class__)],
            variable_names=['interactions', 'colors', 'reset', 'mode']
        )

        # Perform reset if requested
        if reset:
            reset_configuration()
            return

        # Update interaction labels if provided
        if interactions:
            self.interaction_labels = interactions

        # Update colors if provided and valid
        if colors:
            invalid_colors = [color for color in colors if not is_valid_hex_color(color)]
            if invalid_colors:
                raise InvalidColorException(invalid_colors)
            self.colors = colors

        if mode:
            if mode == 'ichem':
                self.interaction_labels = INTERACTION_LABELS  # Default interaction labels
                self.colors = COLORS  # Default color configuration
            elif mode == 'arpeggio':
                self.interaction_labels = ARPEGGIO_CONT + ARPEGGIO_TYPE  # Default interaction labels
                self.colors = ARPEGGIO_COLORS  # Default color c
            else:
                raise InvalidPredictorException(predictor=mode)

    #################################
    # Public Methods: Functionality #
    #################################

    def analyze_files(self, 
        directory: str, 
        predictor: str, 
        activity_file: str = None, 
        protein: bool = True, 
        ligand: bool = True, 
        subunit: bool = False, 
        save: str = None
    ) -> list[list[str]]:
        """
        Analyzes interaction data files in a specified directory, categorizing interactions
        based on protein and ligand atoms involved.

        Args:
            directory (str): Path to the directory containing interaction data files.
            predictor (str): Label to differentiate predictors' files.
            activity_file (str, optional): Path to the activity file (CSV) for labeling data.
            protein (bool, optional): Include protein atoms in the analysis if True.
            ligand (bool, optional): Include ligand atoms in the analysis if True.
            subunit (bool, optional): Differentiate between subunits if True.
            save (str, optional): Path to save the resulting matrix (optional).

        Returns:
            list[list[str]]: A matrix categorizing interactions between residues and files.

        Raises:
            FileNotFoundError: If the specified directory or activity file does not exist.
            EmptyDirectoryException: If the specified directory is empty.
        """

        def label_matrix(
            matrix: list[list[str]], 
            rows: list[str], 
            columns: list[str], 
            activity_file: str,
            correction: list[str] = None
        ) -> list[list[str]]:
            """
            Adds headers to the interaction matrix with residue names and file names.

            Args:
                matrix (list[list[str]]): 2D list representing interaction data.
                rows (list[str]): List of residue names for row labeling.
                columns (list[str]): List of file names for column labeling.
                activity_file (str): Path to the activity file for activity-based labeling.

            Returns:
                list[list[str]]: The labeled matrix.
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
                if correction:
                    for i in range(1, len(columns)):
                        drug_name = columns[i]
                        columns[i] = f"{drug_name} ({data_dict.get(correction[i-1], '0')})"
                else:
                    for i in range(1, len(columns)):
                        drug_name = columns[i]
                        columns[i] = f"{drug_name} ({data_dict.get(drug_name, '0')})"

            # Insert headers into the matrix
            matrix.insert(0, columns)
            for i, row in enumerate(matrix[1:], start=1):
                row.insert(0, rows[i-1])

            return matrix

        def modify_cell(
                text: str, 
                interaction: str, 
                atoms: str,
                interaction_labels: list,
                ) -> str:
            """
            Updates the cell content by adding the interaction type and the involved atoms.

            Args:
                text (str): Current cell content.
                interaction (str): The interaction type.
                atoms (str): Atoms involved in the interaction.

            Returns:
                str: The updated cell content.
            """
            # Create an interaction_map based on the global INTERACTION_LABELS list
            interaction_map = {label: str(index + 1) for index, label in enumerate(interaction_labels)}

            # Assign the interaction code based on the interaction_map, or default to the last value
            interaction_code = interaction_map.get(interaction, str(len(interaction_labels)))

            # If the text is empty, add the interaction with the provided atoms
            if text == "":
                return f"{interaction_code} {GROUP_DELIM}{atoms}{GROUP_DELIM}"

            # Split the cell content and remove empty parts, keeping existing interactions
            content = text.replace(DIFF_DELIM, "").split(GROUP_DELIM)[:-1]
            exists = False
            cell = ''

            # Check if the interaction already exists and add atoms to the corresponding interaction
            for index, segment in enumerate(content):
                if index % 2 == 0 and interaction_code == segment.strip():
                    content[index + 1] += f"{SAME_DELIM}{atoms}"
                    exists = True
                    break

            # Rebuild the cell content, preserving existing interactions
            cell = DIFF_DELIM.join(f"{content[i]}{GROUP_DELIM}{content[i+1]}{GROUP_DELIM}" for i in range(0, len(content), 2))
            
            # If the interaction didn't exist, append it at the end
            if not exists:
                cell += f"{DIFF_DELIM}{interaction_code} {GROUP_DELIM}{atoms}{GROUP_DELIM}"
            
            return cell

        def read_file(
                file_name: str
                ) -> list[str]:
            """
            Reads a file and returns its content as a list of lines.

            Args:
                file_name (str): The name of the file to read.

            Returns:
                list[str]: List of lines from the file.

            Raises:
                FileNotFoundError: If the file does not exist.
                Exception: For any unexpected errors during reading.
            """
            try:
                with open(file_name, 'r') as file:
                    if file_name.split('.')[-1] == "json":
                        return json.load(file)
                    else:
                        return [line.strip() for line in file.readlines()]
            except FileNotFoundError:
                print(f"Error: The file '{file_name}' does not exist.")
                raise
            except Exception as e:
                print(f"Error: An unexpected error occurred while reading '{file_name}': {e}")
                raise

        def adjust_subunits(
                matrix: list[list[str]], 
                ) -> list[list[str]]:
            """
            Adjusts matrix data to reflect residue counts for subunits.

            Args:
                matrix (list[list[str]]): The interaction matrix.

            Returns:
                list[list[str]]: Matrix with adjusted subunit information.
            """
            
            def remove_duplicate_atoms(
                    atoms: str
                    ) -> str:
                """
                Removes duplicate atoms from the input string.

                Args:
                    atoms (str): A string containing atom interactions separated by SAME_DELIM.

                Returns:
                    str: A formatted string with unique atoms, enclosed in '|' delimiters.
                """
                # Split the input string by SAME_DELIM to get individual atoms.
                sections = atoms.split(SAME_DELIM)

                # Use a set to keep only unique atoms in their original order.
                unique_atoms = list(dict.fromkeys(sections))

                # Join the unique atoms back into a single string separated by SAME_DELIM.
                text = SAME_DELIM.join(unique_atoms)

                # Return the formatted string enclosed in '|' delimiters.
                return f"{GROUP_DELIM}{text}{GROUP_DELIM}"

            for row in range(len(matrix)):
                for column in range(len(matrix[row])):
                    cell = matrix[row][column]
                    if cell != '':
                        sections = cell.split(GROUP_DELIM)
                        text = ''.join(
                            sections[i-1] + remove_duplicate_atoms(sections[i])
                            for i in range(1, len(sections), 2)
                        )
                        matrix[row][column] = text
            return matrix

        def sort_interactions(
                matrix: list[list[str]]
                ) -> list[list[str]]:
            """
            Sorts the interactions within each cell of the matrix in ascending order based on the number that follows the initial space.

            Args:
                matrix (list[list[str]]): The matrix to be sorted, represented as a list of lists.

            Returns:
                list[list[str]]: The sorted matrix with interactions in each cell ordered in ascending order.
            """
            for row_index, row in enumerate(matrix):
                for cell_index, cell in enumerate(row):
                    if cell != '':
                        # Split the cell into individual interactions
                        interactions = cell.split(DIFF_DELIM)

                        # Sort the interactions based on the number that follows the initial space
                        if len(interactions) > 1:
                            interactions = sorted(interactions, key=lambda x: int(x.split(" ")[0]))

                        # Join the sorted interactions back and update the cell
                        matrix[row_index][cell_index] = DIFF_DELIM.join(interactions)

            return matrix

        def validate_string(
                input_string: str
                ) -> bool:
            """
            Validates the given string to ensure it follows a specific format:
            - The first three characters must represent an amino acid abbreviation.
            - A space follows the amino acid.
            - After the space(s), there must be a sequence of digits followed by a dash and additional characters.

            Args:
                input_string (str): The string to be validated.

            Returns:
                bool: True if the string matches the required format, False otherwise.
            """
            # Regular expression to validate the required format:
            # - Three uppercase letters (amino acid abbreviation)
            # - Followed by one or more spaces
            # - One or more digits, a dash, and then anything after it
            pattern = r'^[A-Z]{3} +\d+-.+$'
            
            # Match the input string with the regular expression pattern
            return bool(re.match(pattern, input_string))

        def get_protein_ligand(begin: dict, end: dict) -> tuple[dict, dict]:
            if begin["label_comp_type"] == "P":
                return begin, end
            else:
                return end, begin

        # Validate input types
        self._check_variable_types(
            variables=[directory, predictor, activity_file, protein, ligand, subunit, save],
            expected_types=[str, str, (str, None.__class__), bool, bool, bool, (str, None.__class__)],
            variable_names=['directory', 'predictor', 'activity_file', 'protein', 'ligand', 'subunit', 'save']
        )

        # Check if the directory exists
        if not os.path.exists(directory):
            raise FileOrDirectoryNotFoundException(path=directory)
        
        # Check if the predictor is registered
        if predictor not in PREDICTORS:
            raise InvalidPredictorException(predictor=predictor)
        else:
            self.set_plot_config(mode=predictor)

        files = os.listdir(directory)
        if not files:
            raise EmptyDirectoryException(path=directory)

        ligands = [None] * len(files)
        matrix = []
        aa = {}
        cont = 0
        subunits_set = set()
        
        if predictor == 'ichem':
            # Analyze each file in the directory
            for index, file in enumerate(files):
                file_path = os.path.join(directory, file)
                if os.path.isfile(file_path):
                    content = read_file(file_path)
                    for line in content:
                        elements = line.split(GROUP_DELIM)
                        if len(elements) == 10:
                            interaction = elements[0].strip().replace("\t", "")
                            residue = elements[3].strip().replace("\t", "")
                            if validate_string(residue):
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
                                    matrix.append([""] * len(files))

                                matrix[column][index] = modify_cell(text=matrix[column][index], interaction=interaction, atoms=atoms, interaction_labels=INTERACTION_LABELS)
                        
                ligands[index] = file.replace(".txt", "")
            files = None
        elif predictor == 'arpeggio':
            # Analyze each file in the directory
            for index, file in enumerate(files):
                file_path = os.path.join(directory, file)
                if os.path.isfile(file_path):
                    content = read_file(file_path)
                    # Filter to obtain entries with interacting_entities == INTER
                    inter_set = [elem for elem in content if elem["interacting_entities"] in ARPEGGIO_INT_ENT]
                    # Filter to obtain entries with the desired contact or type
                    inter_set = [
                        elem for elem in inter_set
                        for contact in elem["contact"]
                        if contact in ARPEGGIO_CONT or elem["type"] in ARPEGGIO_TYPE
                    ]
                    for inter in inter_set:
                        if inter["type"] in ARPEGGIO_TYPE:
                            contact = inter["type"]
                        else:
                            contact = [cont for cont in inter["contact"] if cont in ARPEGGIO_CONT]
                        protein, ligand = get_protein_ligand(begin=inter["bgn"], end=inter["end"])
                        residue = protein["label_comp_id"] + " " + str(protein["auth_seq_id"])
                        prot_atom = protein["auth_atom_id"]
                        prot_subunit = protein["auth_asym_id"]
                        ligand_code = ligand["label_comp_id"]
                        lig_atom = ligand["auth_atom_id"]
                        
                        subunits_set.add(prot_subunit)
                        atoms = f"{prot_atom}-{lig_atom}" if protein and ligand else prot_atom if protein else lig_atom

                        if subunit:
                            residue += "-" + prot_subunit
                        else:
                            atoms += f"({prot_subunit})"

                        if residue not in aa:
                            aa[residue] = cont
                            cont += 1
                        column = aa[residue]

                        # Ensure matrix size and modify cell
                        if len(matrix) <= column:
                            matrix.append([""] * len(files))

                        for interaction in contact:
                            matrix[column][index] = modify_cell(text=matrix[column][index], interaction=interaction, atoms=atoms, interaction_labels=ARPEGGIO_CONT+ARPEGGIO_TYPE)
 
                # En este caso no vale, se debe detectar el combre del ligando    
                files[index] = file.replace(".json", "")
                ligands[index] = ligand_code

        if not subunit:
            matrix = adjust_subunits(matrix=matrix)
        matrix = sort_interactions(matrix=matrix)
        matrix = label_matrix(matrix=matrix, rows=list(aa.keys()), columns=ligands, activity_file=activity_file, correction=files)

        # Save the matrix if specified
        if save:
            self.save_matrix(matrix=matrix, filename=save)

        return matrix
    
    def filter_by_interaction(self, 
            matrix: list[list[str]], 
            interactions: list[int], 
            save: str = None
            ) -> list[list[str]]:
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

        def validate_list(interactions: list[int]) -> None:
            """
            Validates the interaction list to ensure it contains unique numbers between 1 and 7.

            Args:
                interactions (list): List of integers representing interaction types.

            Raises:
                ValueError: If any number is outside the range of 1 to 7 or if there are duplicates.
            """
            # Valid interactions are numbers from 1 to 7
            valid_numbers = set(range(1, len(self.interaction_labels) + 1))

            # Check if all numbers in the list are within the valid range
            for num in interactions:
                if num not in valid_numbers:
                    raise ValueError(f"Invalid interaction: {num}. Must be a number between 1 and 7.")

            # Ensure the list contains no duplicate values
            if len(set(interactions)) != len(interactions):
                raise ValueError("The interaction list must not contain duplicates.")

        # Validate types of the matrix, interactions, and save parameters
        self._check_variable_types(
            variables=[matrix, interactions, save], 
            expected_types=[list, list, (str, None.__class__)], 
            variable_names=['matrix', 'interactions', 'save']
        )

        # Validate that the matrix has appropriate dimensions
        self._verify_dimensions(matrix=matrix)

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
                if is_not_empty_or_dash(cell):
                    sections = cell.split(DIFF_DELIM)
                    cell = ""
                    
                    # Iterate through the sections in the cell to filter valid interactions
                    for section in sections:
                        # Check if the first number in the section is in the valid interaction list
                        if int(section.split(' ')[0]) in interactions:
                            # Add the section to the cell if it contains a valid interaction
                            if cell == '':
                                cell = section
                            else:
                                cell += DIFF_DELIM + section
                            changes = True
                    
                    # Update the cell with the filtered sections or set it to '-' if empty
                    filtered[i][j] = cell if cell != '' else EMPTY_DASH_CELL
        
        # If no changes were made, raise an error indicating no matching interactions were found
        if not changes:
            raise ValueError("No matching interactions were found in the matrix.")

        # Save the filtered matrix to a file if a save path is provided
        if save:
            self.save_matrix(matrix=filtered, filename=save)

        return filtered

    def filter_by_subunit(self,
            matrix: list[list[str]], 
            subunits: list[str], 
            save: str = None
            ) -> list[list[str]]:
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
        
        def get_subunits_location(matrix: list[list[str]]) -> str:
            """
            Determines the location of subunits in the matrix (residues or interactions).

            Args:
                matrix (list): The matrix being analyzed.

            Returns:
                str: 'residues' if the first column indicates residue data, otherwise 'interactions'.
            """
            return 'residues' if len(matrix[1][0].split('-')) == 2 else 'interactions'
        
        # Check types of the matrix, subunits, and save parameters
        self._check_variable_types(
            variables=[matrix, subunits, save], 
            expected_types=[list, list, (str, None.__class__)], 
            variable_names=['matrix', 'subunits', 'save']
        )

        # Validate the dimensions of the matrix
        self._verify_dimensions(matrix=matrix)

        # Create a deep copy of the matrix to avoid modifying the original
        filtered = copy.deepcopy(matrix)

        # Determine the axis of residues in the matrix
        axis = self._get_residues_axis(matrix=filtered)

        # Transpose the matrix if the axis is columns
        if axis == "columns":
            filtered = self.transpose_matrix(matrix=filtered)
        
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
                    if is_not_empty_or_dash(cell=cell):
                        sections = cell.split(DIFF_DELIM)
                        cell = ""
                        
                        # Iterate through each section in the cell
                        for section in sections:
                            separators = section.split(GROUP_DELIM)[:-1]

                            # Filter out unwanted interactions
                            for index in range(1, len(separators)):
                                if index % 2 != 0:
                                    interactions = separators[index].split(SAME_DELIM)
                                    subchanges = 0
                                    
                                    # Remove interactions not in the valid subunits
                                    for interaction in range(len(interactions)):
                                        if len(interactions[interaction-subchanges].split('(')) > 1 and interactions[interaction - subchanges][-2] not in subunits:
                                            changes += 1
                                            interactions.pop(interaction - subchanges)
                                            subchanges += 1

                                    # Rebuild the cell if there are valid interactions
                                    if interactions:
                                        cell += separators[index - 1] + GROUP_DELIM
                                        cell += SAME_DELIM.join(interactions) + ' ' + GROUP_DELIM + DIFF_DELIM
                                        
                        # Update the cell with filtered interactions
                        cell = cell[:-2]  # Remove trailing comma and space
                        filtered[i][j] = cell if cell else EMPTY_DASH_CELL
        
        # Raise an error if no changes were made (no desired subunits found)
        if changes == 0:
            raise ValueError("The matrix does not contain any of the desired subunits.")

        # Transpose the matrix back if it was originally in columns
        if axis == "columns":
            filtered = self.transpose_matrix(matrix=filtered)

        # Save the filtered matrix to a file if a save path is provided
        if save:
            self.save_matrix(matrix=filtered, filename=save)

        return filtered

    def plot_matrix(self,
        matrix: list[list[str]],
        plot_name: str,
        axis: str,
        label_x: str = None,
        label_y: str = "Number of intermolecular interactions",
        title: str = "Protein-drug interactions",
        stacked: bool = False,
        save: bool = False,
        show_pie_chart: bool = False,
        colors: list[str] = None
    ) -> None:
        """
        Plots a bar chart or pie chart based on selected rows or columns of a matrix and saves it as a PNG file.

        Args:
            matrix (list of lists): The matrix containing interaction data.
            plot_name (str): The name of the plot to be saved (without extension).
            axis (str): Specifies whether to select rows ('rows') or columns ('columns').
            label_x (str, optional): Label for the X-axis. Defaults to "PDB complexes".
            label_y (str, optional): Label for the Y-axis. Defaults to "Number of intermolecular interactions".
            title (str, optional): Title of the plot. Defaults to "Protein-drug interactions".
            stacked (bool, optional): If True, creates a stacked bar chart. Defaults to False.
            save (bool, optional): If True, saves the plot as a PNG file. Defaults to False.
            show_pie_chart (bool, optional): If True, shows a pie chart instead of a bar chart. Defaults to False.
            colors (list, optional): List of colors to use for each interaction type. Defaults to None.

        Returns:
            None
        """

        if colors is None:
            colors = self.colors

        # Ensure the number of colors matches the number of interaction labels
        if len(colors) < len(self.interaction_labels):
            raise ValueError(f"Not enough colors provided. Expected at least {len(self.interaction_labels)} colors, but got {len(colors)}.")

        def get_interactions(cell: str) -> list:
            """
            Extracts and counts interactions from a cell string.

            Args:
                cell (str): The cell string containing interaction data.

            Returns:
                list: A list of interaction counts for each interaction type.
            """
            interactions = [0] * len(self.interaction_labels)
            sections = cell.split(GROUP_DELIM)
            for index in range(1, len(sections), 2):
                interaction = int(sections[index - 1].replace(DIFF_DELIM, "").replace(" ", ""))
                interactions[interaction - 1] += len(sections[index].split(SAME_DELIM))
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
            self._verify_dimensions(matrix=matrix)
            if axis == 'columns':
                matrix = self.transpose_matrix(matrix)
            
            reactives = {row: [0] * len(self.interaction_labels) for row in range(1, len(matrix))}
            indices = [matrix[row][0].split('_')[0].strip() for row in range(1, len(matrix))]
            
            for row in range(1, len(matrix)):
                for column in range(1, len(matrix[row])):
                    cell = matrix[row][column]
                    interactions = get_interactions(cell)
                    for i in range(len(self.interaction_labels)):
                        reactives[row][i] += interactions[i]

            result_list = list(reactives.values())
            return result_list, indices

        # Calculate stacked data if necessary
        if stacked or show_pie_chart:
            data, indices = stack_reactives(matrix=matrix, axis=axis)
            transposed_data = self.transpose_matrix(data)

        # Plot pie chart if requested
        if show_pie_chart:
            # Sum up the total interactions for each type
            total_interactions = [sum(transposed_data[i]) for i in range(len(transposed_data))]

            # Filter out interactions with zero counts
            non_zero_interactions = [(label, total, colors[i]) for i, (label, total) in enumerate(zip(self.interaction_labels, total_interactions)) if total > 0]

            # Prepare pie chart data
            if non_zero_interactions:
                labels_pie, sizes, pie_colors = zip(*non_zero_interactions)
            else:
                labels_pie, sizes, pie_colors = [], [], []

            fig, ax_pie = plt.subplots(figsize=(10, 6))

            # Plotting the pie chart without labels around it
            wedges, texts, autotexts = ax_pie.pie(sizes, labels=None, colors=pie_colors, autopct='', startangle=140)

            # Set the title of the pie chart
            ax_pie.set_title('Interaction Percentages')

            total = 0
            for interaction in sizes:
                total += interaction

            # Adding a legend with all possible interaction labels, regardless of their values
            labels = [label + " (" + str(round(count/total*100, 2)) +"%)" for label, count in zip(self.interaction_labels, total_interactions) if count != 0]
            ax_pie.legend(labels, title="Interaction Types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

        else:
            # Create a new figure for the bar chart
            fig, ax = plt.subplots(num=plot_name, figsize=(12, 6))

            if stacked:
                bars = []
                bottoms = [0] * len(indices)
                for index, group in enumerate(transposed_data):
                    bars.append(ax.bar(indices, group, bottom=bottoms, label=self.interaction_labels[index], color=colors[index]))
                    bottoms = [i + j for i, j in zip(bottoms, group)]

                ax.set_xticks(range(len(indices)))
                if self._get_residues_axis(matrix=matrix) == axis:
                    ax.set_xticklabels(indices)
                else:
                    ax.set_xticklabels(item.split()[0] for item in indices)
                # Adding legend for all labels used in the bar chart
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)

                max_y = max([sum(col) for col in data])
            else:
                data = self.sort_matrix(matrix, axis, count=True)
                ax.bar(data[0], data[1], color=colors[0] if len(colors) > 0 else None)

                max_y = max(data[1])

            ax.set_ylim(0, max_y * 1.1)
            ax.set_ylabel(label_y)
            if label_x is None:
                residues_axis = self._get_residues_axis(matrix=matrix)
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
                    annotation_text = "\n".join([f"{self.interaction_labels[i]}: {percentages[i]:.1f}%" for i in range(len(self.interaction_labels))])
                    sel.annotation.set_text(annotation_text)
                    sel.annotation.get_bbox_patch().set(fc="white", alpha=0.9)  # Set background to white with 90% opacity

        # Show or save the plot
        if not save:
            plt.show()
        else:
            plt.savefig(os.path.join(saving_directory, plot_name + '.png'))
            plt.close(fig)  # Close the figure after saving to avoid display overlap

    def remove_empty_axis(
            self, 
            matrix: list[list[str]],
            save: str = None
            ) -> list[list[SyntaxError]]:
        """
        Remove empty rows and columns from a given matrix.

        Args:
            matrix (list[list[str]]): The input matrix from which empty rows and columns will be removed.
            save (str, optional): The filename to save the cleaned matrix. If None, the matrix will not be saved.

        Returns:
            list[list[str]]: The cleaned matrix with empty rows and columns removed.

        Raises:
            TypeMismatchException: If the types of input variables do not match the expected types.
            ValueError: If the dimensions of the matrix are not valid or if the matrix is too small or any row is too short.
        """
        
        def _remove_empty_rows(matrix: list[list[str]]) -> list[list[str]]:
            """
            Helper function to remove empty rows from the matrix.

            Args:
                matrix (list[list[str]]): The matrix from which empty rows will be removed.

            Returns:
                list[list[str]]: The matrix with empty rows removed.
            """
            changes = 0
            for row in range(1, len(matrix)):
                if all(not is_not_empty_or_dash(cell=column) for column in matrix[row - changes][1:]):
                    matrix.pop(row - changes)
                    changes += 1
            return matrix

        self._check_variable_types(
            variables=[matrix, save], 
            expected_types=[list, (str, None.__class__)], 
            variable_names=['matrix', 'save']
        )

        self._verify_dimensions(matrix=matrix)

        matrix = _remove_empty_rows(matrix=matrix)
        
        # Transpose the matrix, remove empty columns (which are now rows)
        matrix = self.transpose_matrix(matrix=matrix)
        matrix = _remove_empty_rows(matrix=matrix)
        
        # Transpose back to restore original format
        matrix = self.transpose_matrix(matrix=matrix)

        if save:
            self.save_matrix(matrix=matrix, filename=save)
        
        return matrix

    def save_matrix(
            self, 
            matrix: list[list[str]], 
            filename: str
            ) -> None:
        """
        Saves the matrix to a CSV file in the specified directory.

        Args:
            matrix (list[list[str]]): The matrix to be saved.
            filename (str): The name of the file to save the matrix in.

        Returns:
            None

        Raises:
            ValueError: If the dimensions of the matrix are not valid or if the lengths of input lists don't match.
            InvalidFileExtensionException: If the filename does not end with '.csv'.
            TypeMismatchException: If a variable type doesn't match the expected one.
        """

        self._check_variable_types(
            variables=[matrix, filename], 
            expected_types=[list, str], 
            variable_names=['matrix', 'filename']
        )

        # Verify that the filename ends with .csv
        if not filename.endswith('.csv'):
            raise InvalidFileExtensionException(filename)

        self._verify_dimensions(matrix=matrix)

        # Create the CSV file
        file_path = os.path.join(self.saving_directory, filename)
        with open(file_path, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)

            # Write the rows of the matrix to the CSV file
            for row in matrix:
                csv_writer.writerow(row)

    def sort_matrix(self,
        matrix: list[list[str]],
        axis: str = 'rows',
        thr_interactions: int = None,
        thr_activity: float = None,
        selected_items: int = None,
        count: bool = False,
        residue_chain: bool = False,
        save: str = None
    ) -> list[list[str]]:

        """
        Sorts and selects reactive rows or columns from a matrix based on interactions.

        Args:
            matrix (list[list[str]]): The matrix containing interaction data.
            axis (str, optional): Specifies whether to select rows ('rows') or columns ('columns'). Defaults to 'rows'.
            thr_interactions (int, optional): Minimum number of interactions to select a row/column.
            thr_activity (float, optional): Activity threshold to select rows/columns if activity values are given.
            selected_items (int, optional): Number of top rows/columns to select based on interactions.
            count (bool, optional): If True, returns the count of interactions instead of the matrix.
            residue_chain (bool, optional): If True, sorts the resulting matrix based on residue order in the chain.
            save (str, optional): File path to save the resulting matrix.

        Returns:
            list[list[str]]: The selected rows or columns based on the specified criteria.

        Raises:
            ValueError: If both `thr_interactions` and `selected_items` are provided simultaneously.
            ValueError: If the matrix dimensions are insufficient.
        """
        def get_interactions(
                cell: str
                ) -> int:
            """
            Counts the number of interactions in a cell formatted with interaction data.

            Args:
                cell (str): Cell containing interaction data, formatted with '|' separating values.

            Returns:
                int: Total number of interactions in the cell.
            """
            interactions = 0
            sections = cell.split(GROUP_DELIM)
            for index in range(1, len(sections), 2):
                interactions += len(sections[index].split(SAME_DELIM))
            return interactions

        def sort_by_residue(
                matrix: list[list[str]]
                ) -> list[list[str]]:
            """
            Sorts the matrix based on residue indices in the first column.

            Args:
                matrix (list[list[str]]): The matrix to be sorted.

            Returns:
                list[list[str]]: The matrix sorted by residue indices.
            """
            # Validate matrix dimensions
            self._verify_dimensions(matrix=matrix)

            # Determine whether to sort by rows or columns
            axis = self._get_residues_axis(matrix=matrix)
            
            # If sorting by columns, transpose the matrix first
            if axis == 'columns':
                matrix = self.transpose_matrix(matrix=matrix)

            # Separate the header from the data rows
            header = matrix[0]
            data_rows = matrix[1:]

            # Sort the data rows based on residue indices
            sorted_data_rows = sorted(data_rows, key=lambda row: int(row[0].replace(" ", "")[3:].split("-")[0]))

            # Combine the header with the sorted data rows
            sorted_matrix = [header] + sorted_data_rows

            # If sorting was by columns, transpose the sorted matrix back
            if axis == 'columns':
                sorted_matrix = self.transpose_matrix(matrix=sorted_matrix)

            return sorted_matrix

        self._check_variable_types(
            variables=[matrix, axis, thr_interactions, thr_activity, selected_items, count, residue_chain, save], 
            expected_types=[list, str, (int, None.__class__), (float, None.__class__), (int, None.__class__), bool, bool, (str, None.__class__)], 
            variable_names=['matrix', 'axis', 'thr_interactions', 'thr_activity', 'selected_items', 'count', 'residue_chain', 'save']
        )

        self._verify_dimensions(matrix=matrix)

        if axis not in ['rows', 'columns']:
            raise InvalidAxisException(axis)
    
        # Raise an error if `thr_interactions` and `selected_items` are provided simultaneously
        if thr_interactions is not None and selected_items is not None:
            raise ValueError("Cannot select by both 'thr_interactions' and 'selected_items' at the same time.")
        if thr_interactions is not None and thr_activity is not None:
            raise ValueError("Cannot select by both 'thr_interactions' and 'thr_activity' at the same time.")
        if thr_activity is not None and selected_items is not None:
            raise ValueError("Cannot select by both 'thr_activity' and 'selected_items' at the same time.")
        
        # Ensure the correct axis is selected for activity-based selection
        if thr_activity is not None:
            axis = 'rows' if self._get_residues_axis(matrix=matrix) == 'columns' else 'columns'

        # Transpose the matrix if operating on columns
        if axis == 'columns':
            matrix = self.transpose_matrix(matrix=matrix)
        
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
        elif self._get_residues_axis(matrix) == "columns" and thr_activity is not None:
            reactives = [key for key, value in sorted(reactives.items(), key=lambda item: float(matrix[item[0]][0].split(" (")[1].replace(")", "")), reverse=True) if float(matrix[key][0].split(" (")[1].replace(")", "")) >= thr_activity]
        elif selected_items:
            selected_items = min(selected_items, len(matrix))
            reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)[:selected_items]]
        else:
            reactives = [key for key, value in sorted(reactives.items(), key=lambda item: item[1], reverse=True)]
        
        # Create the selection matrix with the chosen rows/columns
        selection = [matrix[0]] + [matrix[row] for row in reactives]

        # Sort the selection by residue chain if specified
        if residue_chain:
            selection = sort_by_residue(matrix=matrix)

        # Transpose the selection back if it was initially transposed for columns
        if axis == 'columns':
            selection = self.transpose_matrix(matrix=selection)

        if save:
            self.save_matrix(matrix=selection, filename=save)

        return selection
    
    def transpose_matrix(
            self, 
            matrix: list[list[str]], 
            save: str = None
            ) -> list[list[str]]:
        """
        Transposes the given matrix.

        Args:
            matrix (list[list[str]]): A matrix to transpose.
            save (str, optional): If provided, saves the transposed matrix to the specified file path.

        Returns:
            list[list[str]]: The transposed matrix.
        
        Raises:
            TypeMismatchException: If the types of the provided arguments are incorrect.
            ValueError: If the dimensions of the matrix are not valid.
        """
        
        self._check_variable_types(
            variables=[matrix, save], 
            expected_types=[list, (str, None.__class__)], 
            variable_names=['matrix', 'save']
        )

        self._verify_dimensions(matrix=matrix)

        # Transpose the matrix using list comprehension
        transposed = [[row[i] for row in matrix] for i in range(len(matrix[0]))]

        if save:
            self.save_matrix(matrix=transposed, filename=save)

        return transposed

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

class InvalidColorException(Exception):
    """
    Exception raised when an invalid hexadecimal color is provided.
    """
    def __init__(self, invalid_colors: list[str]):
        self.invalid_colors = invalid_colors
        self.message = f"Invalid hexadecimal color(s) detected: {', '.join(invalid_colors)}"
        super().__init__(self.message)

class InvalidFileExtensionException(Exception):
    """Exception raised when the file extension is not .csv."""
    def __init__(self, filename, message="Invalid file extension. Only '.csv' files are allowed"):
        self.filename = filename
        self.message = f"{message}: '{filename}'"
        super().__init__(self.message)

class InvalidAxisException(Exception):
    """Exception raised for invalid axis values."""
    def __init__(self, axis_value):
        self.axis_value = axis_value
        self.message = f"Invalid axis value: '{axis_value}'. Expected 'rows' or 'columns'."
        super().__init__(self.message)

class InvalidPredictorException(Exception):
    """Exception raised for invalid predictor values."""
    def __init__(self, predictor):
        self.predictor_value = predictor
        expected_values = ""
        for value in PREDICTORS:
            expected_values += "\n\t- " + value
        self.message = f"Invalid predictor value: '{predictor}'. Expected:{expected_values}"
        super().__init__(self.message)
