#!/bin/bash

version=$(<VERSION)
docker build . -t mm/sanitize-ligand-plugin:${version}
