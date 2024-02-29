import json
import os
import subprocess

def load_json_from_file(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

def extract_keys_and_defaults(data, defaults={}):
    keys = []
    for key, value in data.items():
        keys.append(key)
        if isinstance(value, dict):
            if "default" in value:
                defaults[key] = value["default"]

    return keys, defaults

def create_env_key(key):
    return f'env_{key}'


input_file_path = 'plugin.json'
json_data = load_json_from_file(input_file_path)

input_keys, default_input_keys = extract_keys_and_defaults(json_data["inputs"])
output_keys, default_output_keys = extract_keys_and_defaults(json_data["outputs"])
all_keys = input_keys + output_keys
all_defaults = {**default_input_keys, **default_output_keys}

# only use defaults if they are not already set in the environment
defaults = {key: all_defaults[key] for key in all_defaults if create_env_key(key) not in os.environ}

# now add the defaults to the environment
for key, value in defaults.items():
    os.environ[create_env_key(key)] = str(value)

docker_cmd = "docker run -v $tool_mount_path_input:/inpdir \
           -v $tool_mount_path_output:/outdir \
            --user $(id -u):$(id -g) \
            {{cookiecutter.container_id}}:{{cookiecutter.plugin_version}} \
            {{cookiecutter.base_command}}"

keys_to_env_keys = {key: create_env_key(key) for key in all_keys}
prefixes = {key: json_data["inputs"][key].get("prefix", "") for key in all_keys}
prefix_to_env_key = zip(prefixes.values(), keys_to_env_keys.values())

for prefix, env_key in prefix_to_env_key:
    if prefix != "":
        docker_cmd += f' --{key} ${env_key}'
    else:
        docker_cmd += f' ${env_key}'

# now run in subprocess with the environment set, capture any output and errors
process = subprocess.Popen(docker_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=os.environ)
stdout, stderr = process.communicate()
# print the output and errors
print(stdout.decode('utf-8'))
print(stderr.decode('utf-8'))
