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

def test_change_directory(analyzer: AnalyzeInteractions) -> None:
    """
    Test the change_directory method of the AnalyzeInteractions class.
    
    Args:
        analyzer (AnalyzeInteractions): An instance of the AnalyzeInteractions class.
    """
    def run_tests(test_cases, expected_output, exceptions, test_name) -> None:
        """
        Execute test cases and validate results.

        Args:
            test_cases (list): List of dictionaries containing test parameters.
            expected_output (list): Expected results for each test case.
            exceptions (tuple): Exception types expected during tests.
            test_name (str): Name of the test being executed.
        """
        print(f"Testing {test_name}:")
        passed_tests = 0

        for i, (case, expected) in enumerate(zip(test_cases, expected_output), start=1):
            message = "."
            try:
                analyzer.change_directory(**case)
                result = True
            except exceptions as e:
                result = False
                message = ": " + e.args[0] + "."

            passed_tests += check_test_result(result=result, expected=expected, test_case_index=i, message=message)
        
        print_final_result(passed_tests=passed_tests, total_tests=len(test_cases))

    # Test cases for change_directory
    cases = [
        {"path": "outputs"},
        {"path": "wrong_dir"},
        {"path": 1}
    ]
    expected = [True, False, False]
    exceptions = (ValueError, TypeMismatchException)

    run_tests(cases, expected, exceptions, "change_directory")

def test_set_plot_config(analyzer: AnalyzeInteractions) -> None:
    """
    Test the set_plot_config method of the AnalyzeInteractions class.
    
    Args:
        analyzer (AnalyzeInteractions): An instance of the AnalyzeInteractions class.
    """
    def run_tests(test_cases, expected_output, exceptions, test_name) -> None:
        """
        Execute test cases for set_plot_config and validate results.

        Args:
            test_cases (list): List of dictionaries containing test parameters.
            expected_output (list): Expected results for each test case.
            exceptions (tuple): Exception types expected during tests.
            test_name (str): Name of the test being executed.
        """
        print(f"Testing {test_name}:")
        passed_tests = 0

        for i, (case, expected) in enumerate(zip(test_cases, expected_output), start=1):
            message = "."
            try:
                analyzer.set_plot_config(**case)

                # Initialize result as False
                result = False  

                # Check conditions to set result to True
                if 'reset' in case:
                    result = (case['reset'] and analyzer.interaction_labels == INTERACTION_LABELS and analyzer.colors == COLORS) or not case['reset']
                elif 'colors' in case and 'interactions' in case:
                    result = (case['colors'] == analyzer.colors and case['interactions'] == analyzer.interaction_labels)
                elif 'colors' in case:
                    result = case['colors'] == analyzer.colors
                elif 'interactions' in case:
                    result = case['interactions'] == analyzer.interaction_labels

            except exceptions as e:
                result = False
                message = f": {e.args[0]}"

            passed_tests += check_test_result(result=result, expected=expected, test_case_index=i, message=message)
        
        print_final_result(passed_tests=passed_tests, total_tests=len(test_cases))

    # Test data
    interactions = ["Hydrophobic", "Aromatic_Face/Face"]
    colors_valids = ["#ff6384", "#36a2eb"]
    colors_no_valids = ["ff6384", "#36x2eb"]
    interactions_bigger = ["Hydrophobic", "Aromatic_Face/Face", "Other"]
    colors_bigger = ["#ff6384", "#36a2eb", "#36a45b"]

    # Test cases for set_plot_config
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
    
    # Execute the tests
    run_tests(cases, expected, exceptions, "set_plot_config")

# Initialize colorama for colored output
init(autoreset=True)

# Create an instance of AnalyzeInteractions
analyzer = AnalyzeInteractions()

# Execute all test functions
test_change_directory(analyzer=analyzer)
test_set_plot_config(analyzer=analyzer)
