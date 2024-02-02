#!/bin/bash

version=$(<VERSION)

container_input_dir="/inddir"
container_output_dir="/outdir"

docker run -v $tool_mount_path_input:${container_input_dir} \
           -v $tool_mount_path_output:${container_output_dir} \
            mm/sanitize-ligand-plugin:${version} \
            --pattern ${pattern} \
            --indir ${container_input_dir} \
            --outdir ${container_output_dir}
