from analyze_interactions import *
from colorama import Fore, init
import os
import sys

# Directory for test support files
test_support_directory = "tests_support_files"
total_tests = 0
total_passed_tests = 0

def load_and_compare_csv(directory: str, filename: str, result: list[list[str]]) -> bool:
    """
    Load a CSV file from the specified directory into a matrix and compare it with the provided result matrix.

    Args:
        directory (str): Path to the directory containing the CSV file.
        filename (str): Name of the CSV file to load.
        result (List[List[str]]): The matrix to compare with the loaded content.

    Returns:
        bool: True if the loaded matrix matches the result, False otherwise.

    Raises:
        FileOrDirectoryNotFoundException: If the directory or file does not exist.
        TypeMismatchException: If the result is not of type List[List[str]].
    """
    # Validate 'result' type
    if not isinstance(result, list) or not all(isinstance(row, list) for row in result):
        raise TypeMismatchException(
            variable_name='result',
            expected_types=[list],
            actual_type=type(result)
        )

    # Construct the full path to the file
    file_path = os.path.join(directory, filename)

    # Check if the file exists
    if not os.path.isfile(file_path):
        raise FileOrDirectoryNotFoundException(file_path)

    # Load the CSV content into a matrix
    csv_matrix: list[list[str]] = []
    with open(file_path, 'r', newline='') as file:
        reader = csv.reader(file)
        csv_matrix = [row for row in reader]

    # Compare the loaded matrix with the provided result matrix
    return csv_matrix == result

def create_directory(directory_name: str) -> None:
    """
    Create a directory if it doesn't exist. Stops the program if creation fails.

    Args:
        directory_name (str): Name of the directory to create.

    Raises:
        SystemExit: If the directory cannot be created due to an error.
    """
    try:
        # Create the directory if it does not exist
        os.makedirs(directory_name, exist_ok=True)
        print(f"Directory '{directory_name}' created or already exists.")
    except OSError as e:
        print(f"Critical Error: Could not create directory '{directory_name}'. {e}")
        sys.exit(1)  # Stop the program immediately with a non-zero exit code

def check_test_result(result: bool, expected: bool, test_case_index: int, message: str) -> int:
    """
    Check if the test result matches the expected output and print the corresponding message.

    Args:
        result (bool): The actual result of the test case.
        expected (bool): The expected result of the test case.
        test_case_index (int): The index of the test case.

    Returns:
        int: 1 if the test passed, 0 if it failed.
    """
    if result == expected:
        print(f"{Fore.LIGHTGREEN_EX}\tTest case {test_case_index} passed correctly{message}")
        return 1
    else:
        print(f"{Fore.LIGHTRED_EX}\tTest case {test_case_index} failed. Expected {'pass' if expected else 'fail'}, but got {'pass' if result else 'fail'}.")
        return 0

def print_final_result(passed_tests: int, total_tests: int) -> None:
    """
    Print the final result of the test cases.

    Args:
        passed_tests (int): The number of tests that passed.
        total_tests (int): The total number of tests executed.
    """
    result_message = f"Passed tests ({passed_tests}/{total_tests}).\n"
    print(f"{Fore.LIGHTGREEN_EX if passed_tests == total_tests else Fore.LIGHTRED_EX}{result_message}")

def run_tests(
    analyzer: AnalyzeInteractions, 
    method_name: str, 
    test_cases: list[dict], 
    expected_output: list[bool], 
    exceptions: tuple, 
    validate_result
) -> None:
    """
    Execute test cases and validate the results.

    Args:
        analyzer (AnalyzeInteractions): An instance of the AnalyzeInteractions class.
        method_name (str): The name of the method being tested.
        test_cases (list[dict]): List of dictionaries containing test parameters.
        expected_output (list[bool]): Expected results for each test case.
        exceptions (tuple): Exception types expected during tests.
        validate_result (callable): Function to validate the result for each test case.
    """
    print(f"Testing {method_name}:")
    passed_tests = 0

    for i, (case, expected) in enumerate(zip(test_cases, expected_output), start=1):
        message = "."
        try:
            # Call the method dynamically and store the result if any
            case['result'] = getattr(analyzer, method_name)(**case)

            # Validate the result using custom logic
            result = validate_result(analyzer, case)
            if result and 'save' in case:
                if isinstance(case['save'], str):
                    exceptions = exceptions + (FileOrDirectoryNotFoundException,)
                    result = load_and_compare_csv(directory=analyzer.saving_directory,
                                                  filename=case['save'],
                                                  result=case['result']);
                
        except exceptions as e:
            result = False
            message = f": {e.args[0]}"

        # Check if the test passed and accumulate the count of passed tests
        passed_tests += check_test_result(
            result=result, expected=expected, test_case_index=i, message=message
        )

    # Print the final summary of tests
    print_final_result(passed_tests=passed_tests, total_tests=len(test_cases))
    global total_tests, total_passed_tests
    total_tests += len(test_cases)
    total_passed_tests += passed_tests

def validate_change_directory_result(analyzer: AnalyzeInteractions, case: dict) -> bool:
    """Validation logic for change_directory."""
    return os.path.isdir(case.get("path", ""))

def test_change_directory(analyzer: AnalyzeInteractions) -> None:
    """
    Test the change_directory method of the AnalyzeInteractions class.
    
    Args:
        analyzer (AnalyzeInteractions): An instance of the AnalyzeInteractions class.

    Returns:
        None
    """

    cases = [
        {"path": test_support_directory},
        {"path": "wrong_dir"},
        {"path": 1}
    ]
    expected = [True, False, False]
    exceptions = (ValueError, TypeMismatchException)

    run_tests(
        analyzer=analyzer,
        method_name="change_directory",
        test_cases=cases,
        expected_output=expected,
        exceptions=exceptions,
        validate_result=validate_change_directory_result
    )

def validate_set_plot_config_result(analyzer: AnalyzeInteractions, case: dict) -> bool:
    """Validation logic for set_plot_config."""
    if 'reset' in case:
        return (case['reset'] and analyzer.interaction_labels == INTERACTION_LABELS and analyzer.colors == COLORS) or not case['reset']
    if 'colors' in case and 'interactions' in case:
        return case['colors'] == analyzer.colors and case['interactions'] == analyzer.interaction_labels
    if 'colors' in case:
        return case['colors'] == analyzer.colors
    if 'interactions' in case:
        return case['interactions'] == analyzer.interaction_labels
    return False

def test_set_plot_config(analyzer: AnalyzeInteractions) -> None:
    """
    Test the set_plot_config method of the AnalyzeInteractions class.
    
    Args:
        analyzer (AnalyzeInteractions): An instance of the AnalyzeInteractions class.

    Returns:
        None
    """

    interactions = ["Hydrophobic", "Aromatic_Face/Face"]
    colors_valids = ["#ff6384", "#36a2eb"]
    colors_no_valids = ["ff6384", "#36x2eb"]
    interactions_bigger = ["Hydrophobic", "Aromatic_Face/Face", "Other"]
    colors_bigger = ["#ff6384", "#36a2eb", "#36a45b"]

    cases = [
        {"interactions": interactions, "colors": colors_valids, "reset": False},  # Valid configuration
        {"interactions": 123, "colors": colors_valids, "reset": False},  # Invalid types
        {"interactions": interactions, "colors": True, "reset": False},
        {"interactions": interactions, "colors": colors_valids, "reset": "False"},
        {"interactions": interactions, "colors": colors_no_valids, "reset": False},
        {"colors": colors_valids, "reset": False},  # Partial updates
        {"interactions": interactions, "reset": False},
        {"interactions": interactions, "colors": colors_valids},
        {"reset": True},  # Reset cases
        {"reset": False},
        {},  # Empty input
        {"interactions": interactions_bigger, "colors": colors_valids},  # Size mismatches
        {"interactions": interactions, "colors": colors_bigger}
    ]

    expected = [
        True,  # Valid configuration
        False, False, False, False,  # Invalid types or invalid colors
        True, True, True,  # Partial updates
        True, True,  # Reset cases
        False,  # Empty input (no valid config)
        True, True  # Size mismatches allowed
    ]

    exceptions = (ValueError, TypeMismatchException, InvalidColorException)

    run_tests(
        analyzer=analyzer,
        method_name="set_plot_config",
        test_cases=cases,
        expected_output=expected,
        exceptions=exceptions,
        validate_result=validate_set_plot_config_result
    )

def validate_transpose_matrix_result(analyzer: AnalyzeInteractions, case: dict) -> bool:
    """Validation logic for transpose_matrix."""

    matrix = case['matrix']
    transpose = case['result']

    # Check if dimensions are compatible
    if len(matrix) != len(transpose[0]) or len(matrix[0]) != len(transpose):
        return False

    # Check if elements match for the transpose condition
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] != transpose[j][i]:
                return False

    return True

def test_transpose_matrix(analyzer: AnalyzeInteractions) -> None:
    """
    Test the transpose_matrix method of the AnalyzeInteractions class.
    
    Args:
        analyzer (AnalyzeInteractions): An instance of the AnalyzeInteractions class.

    Returns:
        None
    """

    matrix = [['abc', 'def', 'hij'], ['klm', 'nop', 'qrs']]
    cases = [
        {"matrix": matrix, "save": 'test.csv'},  # Valid configuration
        {"matrix": matrix, "save": False},  # Invalid types
        {"matrix": False, "save": 'test.csv'},
        {"matrix": matrix}, # Partial updates
        {"save": 'test.csv'},
        {}  # Empty input
    ]

    expected = [
        True,  # Valid configuration
        False, False, # Invalid types or invalid colors
        True, False, # Partial updates
        False  # Empty input (no valid config)
    ]

    exceptions = (ValueError, TypeMismatchException, TypeError)

    run_tests(
        analyzer=analyzer,
        method_name="transpose_matrix",
        test_cases=cases,
        expected_output=expected,
        exceptions=exceptions,
        validate_result=validate_transpose_matrix_result
    )

def validate_save_matrix_result(analyzer: AnalyzeInteractions, case: dict) -> bool:
    """Validation logic for transpose_matrix."""

    return load_and_compare_csv(directory=analyzer.saving_directory,
                                filename=case['filename'],
                                result=case['matrix'])

def test_save_matrix(analyzer: AnalyzeInteractions) -> None:
    """
    Test the save_matrix method of the AnalyzeInteractions class.
    
    Args:
        analyzer (AnalyzeInteractions): An instance of the AnalyzeInteractions class.

    Returns:
        None
    """

    matrix = [['abc', 'def', 'hij'], ['klm', 'nop', 'qrs']]
    cases = [
        {"matrix": matrix, "filename": 'test.csv'},  # Valid configuration
        {"matrix": matrix, "filename": False},  # Invalid types
        {"matrix": False, "filename": 'test.csv'},
        {"matrix": matrix}, # Partial updates
        {"filename": 'test.csv'},
        {},  # Empty input
        {"matrix": matrix, "filename": 'test.txt'}   # Invalid file extension
    ]

    expected = [
        True,  # Valid configuration
        False, False, # Invalid types or invalid colors
        False, False, # Partial updates
        False,  # Empty input (no valid config)
        False   # Invalid file extension
    ]

    exceptions = (ValueError, TypeMismatchException, TypeError, InvalidFileExtensionException)

    run_tests(
        analyzer=analyzer,
        method_name="save_matrix",
        test_cases=cases,
        expected_output=expected,
        exceptions=exceptions,
        validate_result=validate_save_matrix_result
    )


def validate_remove_empty_axis_result(analyzer: AnalyzeInteractions, case: dict) -> bool:
    """Validation logic for transpose_matrix."""

    empty_chars = [EMPTY_CELL, EMPTY_DASH_CELL]

    def is_empty_row(row: list[str]) -> bool:
        return all(cell in empty_chars for cell in row[1:])

    def is_empty_column(column: tuple[str]) -> bool:
        return all(cell in empty_chars for cell in column[1:])

    matrix = case['matrix']
    
    # Remove empty rows
    matrix = [row for row in matrix if not is_empty_row(row)]

    # Transpose the filtered matrix to easily remove empty columns
    matrix = list(zip(*matrix))

    # Remove empty columns
    matrix = [column for column in matrix if not is_empty_column(column)]

    # Transpose back to the original orientation and convert tuples to lists
    matrix = [list(row) for row in zip(*matrix)]

    # Compare the final_matrix with the expected result
    return matrix == case['result']

def test_remove_empty_axis(analyzer: AnalyzeInteractions) -> None:
    """
    Test the remove_empty_axis method of the AnalyzeInteractions class.
    
    Args:
        analyzer (AnalyzeInteractions): An instance of the AnalyzeInteractions class.

    Returns:
        None
    """

    matrix_1 = [['abc', '123', 'hij', '456'], 
                ['klm', EMPTY_DASH_CELL, 'qrs', EMPTY_CELL],
                ['klm', EMPTY_DASH_CELL, 'qrs', EMPTY_CELL]]
    matrix_2 = [['123','456','789'], 
                ['abc', 'def', 'hij'], 
                ['123', EMPTY_CELL, EMPTY_DASH_CELL], 
                ['klm', 'nop', 'qrs']]
    cases = [
        {"matrix": matrix_1, "save": 'test.csv'},  # Valid configuration
        {"matrix": matrix_1, "save": False},  # Invalid types
        {"matrix": False, "save": 'test.csv'},
        {"matrix": matrix_2}, # Partial updates
        {"save": 'test.csv'},
        {},  # Empty input
    ]

    expected = [
        True,  # Valid configuration
        False, False, # Invalid types or invalid colors
        True, False, # Partial updates
        False,  # Empty input (no valid config)
        False   # Invalid file extension
    ]

    exceptions = (ValueError, TypeMismatchException, TypeError, InvalidFileExtensionException)

    run_tests(
        analyzer=analyzer,
        method_name="remove_empty_axis",
        test_cases=cases,
        expected_output=expected,
        exceptions=exceptions,
        validate_result=validate_remove_empty_axis_result
    )


# Initialize colorama
init(autoreset=True)
create_directory(test_support_directory)
analyzer = AnalyzeInteractions()

# Run the tests
test_change_directory(analyzer)
test_set_plot_config(analyzer)
test_transpose_matrix(analyzer)
test_save_matrix(analyzer)
test_remove_empty_axis(analyzer)

result_message = f"Total passed tests ({total_passed_tests}/{total_tests}).\n"
print(f"{Fore.LIGHTGREEN_EX if total_passed_tests == total_tests else Fore.LIGHTRED_EX}{result_message}")