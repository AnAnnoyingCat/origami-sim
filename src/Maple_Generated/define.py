import os
import re
import argparse

def add_double_declarations(file_content):
    lines = file_content.splitlines()
    declared = set()
    updated_lines = []

    for line in lines:
        match = re.match(r"\s*(t\d+)\s*=", line)
        if match:
            var = match.group(1)
            if var not in declared:
                line = f"double {line.strip()}"
                declared.add(var)
            else:
                line = line.strip()
        updated_lines.append(line)

    return '\n'.join(updated_lines)

def process_and_save_file(input_path):
    base, ext = os.path.splitext(os.path.basename(input_path))
    output_path = os.path.join(os.path.dirname(input_path), f"{base}-defined.cpp")

    with open(input_path, "r") as infile:
        content = infile.read()
    updated_content = add_double_declarations(content)

    with open(output_path, "w") as outfile:
        outfile.write(updated_content)

    print(f"Processed file saved as: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Add double declarations to variables in a file.")
    parser.add_argument("filepath", help="Path to the input file")
    args = parser.parse_args()

    process_and_save_file(args.filepath)

if __name__ == "__main__":
    main()
