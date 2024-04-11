import os
import json
import sys

def extract_functions_from_py(file_path):
    """
    Extracts functions from a .py file.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    return extract_functions_from_lines(lines)

def extract_functions_from_ipynb(file_path):
    """
    Extracts functions from a .ipynb file.
    """
    with open(file_path, 'r') as file:
        notebook = json.load(file)

    function_defs = []
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") == "code":
            lines = cell.get("source", [])
            function_defs.extend(extract_functions_from_lines(lines))

    return function_defs

def extract_functions_from_lines(lines):
    """
    Extracts function definitions from a list of lines.
    """
    function_defs = []
    in_function = False
    for line in lines:
        if line.startswith('def '):
            in_function = True
            function_defs.append(line)
        elif in_function and (line.startswith(' ') or line.startswith('\t')):
            function_defs.append(line)
        else:
            if in_function:
                function_defs.append('\n\n')
            in_function = False

    return function_defs

def extract_functions_from_files(directory,outfile='../extracted_functions.py'):
    """
    Extracts functions from .py and .ipynb files in the given directory and writes them to a new Python file.
    """
    all_function_defs = []

    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if file.endswith('.py'):
                all_function_defs.extend(extract_functions_from_py(file_path))
            elif file.endswith('.ipynb'):
                all_function_defs.extend(extract_functions_from_ipynb(file_path))

    with open(outfile, 'w') as output_file:
        output_file.writelines(all_function_defs)

# Example usage
#extract_functions_from_files('/home/matthew.schmitz/Matthew/code/v1-analysis/','v1_extracted_functions.py')
extract_functions_from_files(*sys.argv[1:])
