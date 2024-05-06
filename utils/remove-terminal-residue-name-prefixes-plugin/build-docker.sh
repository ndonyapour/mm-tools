#!/bin/bash

version=$(<VERSION)
docker build . -t polusai/remove-terminal-residue-name-prefixes-plugin:${version}
