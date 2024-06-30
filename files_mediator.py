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
def read_txt_file(file_name: str) -> list:
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