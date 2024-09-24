from analyze_interactions import *
from colorama import Fore, init
import os

test_support_directory = "tests_support_files"
change_directory("outputs")

def tests_cdu_01():
    def run_tests(test_cases, expected_output, exceptions, test_name, save=False):
        print(f"Running {test_name}:")
        passed_tests = 0

        for i, (case, expected) in enumerate(zip(test_cases, expected_output), start=1):
            
            try:
                output = analyze_files(**case)
                result = True
                if save:
                    save_matrix(matrix=output, filename=f'{test_name}-{i}.csv')
            except exceptions as e:
                result = False
                print(f"{Fore.LIGHTYELLOW_EX}\tTest case {i} failed with error: {e}")            

            # Check if the test result matches the expected output
            if result == expected:
                print(f"{Fore.LIGHTGREEN_EX}\tTest case {i} passed correctly.")
                passed_tests += 1
            else:
                print(f"{Fore.LIGHTRED_EX}\tTest case {i} failed. Expected {'pass' if expected else 'fail'}, but got {'pass' if result else 'fail'}.")

        # Print the final result
        result_message = f"Tests passed ({passed_tests}/{len(test_cases)}).\n"
        if passed_tests == len(test_cases):
            print(f"{Fore.LIGHTGREEN_EX}{result_message}")
        else:
            print(f"{Fore.LIGHTRED_EX}{result_message}")
    
    # Cdu 01.T1
    test_1_cases = [
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"directory": 1, "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": True, "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": "True", "subunit": True},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "protein": True, "ligand": True, "subunit": True}
    ]
    expected_1 = [True, False, False, False, True, True]
    exceptions_1 = TypeMismatchException

    run_tests(test_1_cases, expected_1, exceptions_1, "Cdu 01.T1")

    # Cdu 01.T2
    test_2_cases = [
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"directory": os.path.join(test_support_directory, "cdu 01. T2: Case 3"), "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"directory": os.path.join(test_support_directory, "cdu 01. T2: Case 4"), "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"directory": os.path.join(test_support_directory, "cdu 01. T2: Case 5"), "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True}
    ]
    expected_2 = [True, False, False, True, True, False]
    exceptions_2 = (FileOrDirectoryNotFoundException, EmptyDirectoryException, TypeError)

    run_tests(test_2_cases, expected_2, exceptions_2, "Cdu 01.T2")

    # Cdu 01.T3
    test_3_cases = [
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": "activities.csv", "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": os.path.join(test_support_directory, "cdu 01. T3: Case 4.csv"), "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": os.path.join(test_support_directory, "cdu 01. T3: Case 5.csv"), "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": os.path.join(test_support_directory, "cdu 01. T3: Case 6.csv"), "protein": True, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": os.path.join(test_support_directory, "cdu 01. T3: Case 7.csv"), "protein": True, "ligand": True, "subunit": True}
    ]
    expected_3 = [True, False, True, False, False, True, True]
    exceptions_3 = (FileOrDirectoryNotFoundException, EmptyDirectoryException, TypeError, FileNotFoundError, ValueError)

    run_tests(test_3_cases, expected_3, exceptions_3, "Cdu 01.T3", True)

    # Cdu 01.T4
    test_4_cases = [
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": False, "ligand": False, "subunit": False},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": False, "subunit": False},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": False},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": False, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": False, "ligand": True, "subunit": False},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": False, "ligand": True, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": False, "ligand": False, "subunit": True},
        {"directory": "ichem-ifp", "activity_file": "ichem-ifp-activities.csv", "protein": True, "ligand": True, "subunit": True},
    ]
    expected_4 = [True, True, True, True, True, True, True, True]
    exceptions_4 = (FileOrDirectoryNotFoundException, EmptyDirectoryException, TypeError, FileNotFoundError, ValueError)

    run_tests(test_4_cases, expected_4, exceptions_4, "Cdu 01.T4", True)

def tests_cdu_02():
    def run_tests(test_cases, expected_output, exceptions, test_name, save=False):
        print(f"Running {test_name}:")
        passed_tests = 0

        for i, (case, expected) in enumerate(zip(test_cases, expected_output), start=1):
            
            try:
                output = sort_matrix(**case)
                result = True
                if save:
                    save_matrix(matrix=output, filename=f'{test_name}-{i}.csv')
            except exceptions as e:
                result = False
                print(f"{Fore.LIGHTYELLOW_EX}\tTest case {i} failed with error: {e}")            

            # Check if the test result matches the expected output
            if result == expected:
                print(f"{Fore.LIGHTGREEN_EX}\tTest case {i} passed correctly.")
                passed_tests += 1
            else:
                print(f"{Fore.LIGHTRED_EX}\tTest case {i} failed. Expected {'pass' if expected else 'fail'}, but got {'pass' if result else 'fail'}.")

        # Print the final result
        result_message = f"Tests passed ({passed_tests}/{len(test_cases)}).\n"
        if passed_tests == len(test_cases):
            print(f"{Fore.LIGHTGREEN_EX}{result_message}")
        else:
            print(f"{Fore.LIGHTRED_EX}{result_message}")
    
    matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=False)

    # Cdu 02.T1
    test_1_cases = [
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_activity": 6.1, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "selected_items": 4, "count": True, "residue_chain": True},
        {"matrix": 'matrix', "axis": "columns", "thr_interactions": 7, "thr_activity": 6.1, "selected_items": 4, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": 23, "thr_interactions": 7, "thr_activity": 6.1, "selected_items": 4, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "fila", "selected_items": 4, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 2.3, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "thr_activity": 6, "selected_items": 4, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "thr_activity": True, "selected_items": 4, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "thr_activity": 6.1, "selected_items": 4.8, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "thr_activity": 6.1, "selected_items": 'False', "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "thr_activity": 6.1, "selected_items": 4, "count": 4, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "thr_activity": 6.1, "selected_items": 4, "count": True, "residue_chain": 'True'},
        {"matrix": matrix, "thr_interactions": 7, "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "count": True, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "residue_chain": True},
        {"matrix": matrix, "axis": "columns", "thr_interactions": 7, "count": True}
    ]
    expected_1 = [True, True, True, False, False, True, False, False, False, False, False, False, False, True, True, True, True]
    exceptions_1 = TypeMismatchException
    
    run_tests(test_1_cases, expected_1, exceptions_1, "Cdu 02.T1")

def tests_cdu_03():
    def run_tests(test_cases, expected_output, exceptions, test_name, save=False):
        print(f"Running {test_name}:")
        passed_tests = 0

        for i, (case, expected) in enumerate(zip(test_cases, expected_output), start=1):
            
            try:
                output = transpose_matrix(**case)
                result = True
                if save:
                    save_matrix(matrix=output, filename=f'{test_name}-{i}.csv')
            except exceptions as e:
                result = False
                print(f"{Fore.LIGHTYELLOW_EX}\tTest case {i} failed with error: {e}")            

            # Check if the test result matches the expected output
            if result == expected:
                print(f"{Fore.LIGHTGREEN_EX}\tTest case {i} passed correctly.")
                passed_tests += 1
            else:
                print(f"{Fore.LIGHTRED_EX}\tTest case {i} failed. Expected {'pass' if expected else 'fail'}, but got {'pass' if result else 'fail'}.")

        # Print the final result
        result_message = f"Tests passed ({passed_tests}/{len(test_cases)}).\n"
        if passed_tests == len(test_cases):
            print(f"{Fore.LIGHTGREEN_EX}{result_message}")
        else:
            print(f"{Fore.LIGHTRED_EX}{result_message}")
    
    matrix = analyze_files(directory="ichem-ifp", activity_file="ichem-ifp-activities.csv", protein=True, ligand=True, subunit=False)

    # Cdu 03.T1
    test_1_cases = [
        {"matrix": matrix},
        {"matrix": 12}
    ]
    expected_1 = [True, False, False, False, True, True]
    exceptions_1 = TypeMismatchException

    run_tests(test_1_cases, expected_1, exceptions_1, "Cdu 03.T1")

# Initialize colorama
init(autoreset=True)

# Run all tests
#tests_cdu_01()
#tests_cdu_02()
tests_cdu_03()