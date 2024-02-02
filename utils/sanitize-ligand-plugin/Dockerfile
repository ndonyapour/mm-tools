FROM condaforge/mambaforge

# environment variables defined
ENV EXEC_DIR="/opt/executables"

RUN mkdir -p ${EXEC_DIR}
COPY pyproject.toml ${EXEC_DIR}
COPY VERSION ${EXEC_DIR}
COPY README.md ${EXEC_DIR}
COPY CHANGELOG.md ${EXEC_DIR}

RUN conda install -c conda-forge rdkit --yes

RUN pip install filepattern

RUN conda init bash

RUN mamba clean --all --yes

ADD Dockerfile ${EXEC_DIR}

COPY src ${EXEC_DIR}/src

RUN pip3 install ${EXEC_DIR} --no-cache-dir
# Default command. Additional arguments are provided through the command line
ENTRYPOINT ["python3", "-m", "polus.mm.utils"]
CMD ["--help"]
