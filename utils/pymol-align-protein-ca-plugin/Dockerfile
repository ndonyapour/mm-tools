FROM condaforge/mambaforge

RUN apt-get update && apt-get install libgl1-mesa-glx -y

ENV EXEC_DIR="/opt/executables"
ENV POLUS_LOG="INFO"
RUN mkdir -p ${EXEC_DIR}

# Work directory defined in the base container
# WORKDIR ${EXEC_DIR}

COPY pyproject.toml ${EXEC_DIR}
COPY VERSION ${EXEC_DIR}
COPY README.md ${EXEC_DIR}
COPY CHANGELOG.md ${EXEC_DIR}

# Install needed packages here
RUN mamba install -c conda-forge mdanalysis
RUN mamba install -c conda-forge pymol-open-source

COPY src ${EXEC_DIR}/src

RUN pip3 install ${EXEC_DIR} --no-cache-dir

CMD ["--help"]
