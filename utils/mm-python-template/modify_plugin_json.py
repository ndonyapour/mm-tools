"""Need to reformat plugin.json to have short lines and indents,
also replace single quotes with double since cookiecutter seems 
to add single quotes but json needs double quotes
"""

from pathlib import Path
import json
import sys
import re

# read in path to json from command line
json_path = Path(sys.argv[1])
input_json = json_path.read_text()

# Regex pattern to find single quotes with characters on left and right
pattern = re.compile(r"(?<=\S)'(?=\S)")

# Find all matches in the input JSON
matches = pattern.finditer(input_json)

indices_to_remove = []
skip_chars = [':', '}', '{', ',']
for match in matches:
    # Get the start and end positions of the match
    start_pos = match.start()
    end_pos = match.end()

    # Determine the left and right characters
    left_char = input_json[start_pos - 1] if start_pos > 0 else ''
    right_char = input_json[end_pos] if end_pos < len(input_json) else ''

    # Check if left or right character contains ':', '}', or '{'
    if left_char not in skip_chars and right_char not in skip_chars:
        # Replace the single quote with an empty character in the input_json
        indices_to_remove.append(start_pos)

# Remove the single quotes from the input_json
for index in reversed(indices_to_remove):
    input_json = input_json[:index] + input_json[index + 1:]

# Replace single quotes that are surrounding keys or values with double quotes
cleaned_json = json.loads(input_json.replace("'", '"'))

json_path.write_text(json.dumps(cleaned_json, indent=2))
