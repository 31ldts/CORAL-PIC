from analyze_interactions import *
from colorama import Fore, init

# Directory for test support files
test_support_directory = "tests_support_files"

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
    result_message = f"Tests passed ({passed_tests}/{total_tests}).\n"
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
            # Dynamically call the method using its name
            getattr(analyzer, method_name)(**case)
            result = validate_result(analyzer, case)  # Call custom validation logic
        except exceptions as e:
            result = False
            message = f": {e.args[0]}"

        passed_tests += check_test_result(result=result, expected=expected, test_case_index=i, message=message)

    print_final_result(passed_tests=passed_tests, total_tests=len(test_cases))

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
        {"path": "outputs"},
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

# Initialize colorama
init(autoreset=True)
analyzer = AnalyzeInteractions()

# Run the tests
test_change_directory(analyzer)
test_set_plot_config(analyzer)
