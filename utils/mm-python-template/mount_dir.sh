#!/bin/bash

# Set environment mount directory variables
mkdir -p ~/mm-tools-mount 
export root_mount_dir=~/mm-tools-mount
export tool_mount_path=$root_mount_dir/$tool_mount_dir # will use this in run-plugin.sh
mkdir -p $tool_mount_path
export tool_mount_path_input=$tool_mount_path/inpdir
mkdir -p $tool_mount_path_input
export tool_mount_path_output=$tool_mount_path/outdir
mkdir -p $tool_mount_path_output

echo "root_mount_dir: $root_mount_dir"
echo "tool_mount_path: $tool_mount_path"
echo "tool_mount_path_input: $tool_mount_path_input"
echo "tool_mount_path_output: $tool_mount_path_output"
