#!/bin/bash

PYTHON_VERSIONS=("cp38-cp38" "cp39-cp39" "cp310-cp310" "cp311-cp311" "cp312-cp312" "cp313-cp313")

cd /io
for PYVER in "${PYTHON_VERSIONS[@]}"
do
    PYTHON_BIN="/opt/python/${PYVER}/bin"
    ${PYTHON_BIN}/pip install --upgrade pip setuptools wheel auditwheel pybind11
    ${PYTHON_BIN}/python setup.py bdist_wheel

    for whl in dist/*.whl
    do
        ${PYTHON_BIN}/auditwheel repair "$whl" --plat manylinux2014_x86_64 -w dist/
    done
done
