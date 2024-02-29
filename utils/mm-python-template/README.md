# Initial Setup
'''
sudo apt update
sudo apt install pipx
pipx ensurepath
pipx install poetry
pipx install bump2version
pipx install cookiecutter
pipx install pre-commit # for code sanitation
cd ~/mm-tools
pre-commit install
'''

# Set Mount Directory Environment Variables
'''
cd ~/mm-tools/utils/mm-python-template
export tool_mount_dir=extract_pdbbind_refined
chmod +x mount_dir.sh
./mount_dir.sh
'''

# Create New Tool
'''
python read_cwl_inputs_outputs.py ~/mm-workflows/cwl_adapters/extract_pdbbind_refined.cwl
'''

# Reformat Plugin.json
'''
python modify_plugin_json.py pdbbind-refined-v2020-plugin/plugin.json
'''

# Modify cookiecutter.json as needed
'''
python modify_base_template.py
pipx run cookiecutter . --no-input
cd new_tool_created
'''
# Modify pyproject.toml
# Modify Dockerfile
# Modify python source file, remove argparse, copy source code etc…
# Build container
'''
./build-docker.sh
'''

# Create Non-Default Environment Variables
'''
export query='(Kd_Ki == "Kd") and (value < 0.001)'
export min_row=0
export max_row=1
export convert_Kd_dG=True
'''

# Run container & Extract Default Parameters to Environment Variables
'''
python run-plugin.py
./run-plugin.sh
'''

# Add tests
'''
pytest
'''

# Post Development
'''
mv new_tool_created ../.
pipx run bump2version dev --dry-run —verbose
'''