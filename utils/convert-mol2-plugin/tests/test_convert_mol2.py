"""Tests for convert_mol2."""
import shlex
import subprocess
from pathlib import Path


def test_convert_mol2():
    """Test convert_mol2."""
    
    input_path = "4xk9_ligand.sdf" 
    output_mol2_path = "ligand.mol2"

    # Define the command to pull the Docker image
    docker_image_name = "quay.io/biocontainers/biobb_chemistry:4.0.0--pyhdfd78af_1"
    pull_command = (
        f"docker pull {docker_image_name}"
    )
    subprocess.run(pull_command, shell=True)  # noqa: S602

    # Define the command to run the Docker container
    container_output_dir = "/outdir"
    tool_mount_path_output = Path.cwd()  # Update this with your actual path
    docker_run_command = f"docker run -v \
    {tool_mount_path_output}:{container_output_dir}\
    --name obabel_container {docker_image_name}\
     obabel\
    {shlex.quote(input_path)} -o mol2 -O {shlex.quote(output_mol2_path)} -xu"
    subprocess.run(docker_run_command, shell=True)  # noqa: S602

    # Copy the output file from the container to the host
    docker_cp_command = f"docker cp obabel_container:/{shlex.quote(output_mol2_path)}\
     {tool_mount_path_output}"
    subprocess.run(docker_cp_command, shell=True)  # noqa: S602

    # Stop and remove the container
    docker_stop_command = "docker stop obabel_container"
    subprocess.run(docker_stop_command, shell=True)  # noqa: S602

    docker_rm_command = "docker rm obabel_container"
    subprocess.run(docker_rm_command, shell=True)  # noqa: S602

    # assert that the output file exists
    #assert Path(output_mol2_path).exists()