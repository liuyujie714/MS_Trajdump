from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os, sys

__version__ = '0.1'

src=[
    'MS_dump/MSTrjParser.cpp',
    'MS_dump/PyModule.cpp',
]

# tell compiler cpp must need g++
LINKERS = []
DEFINES = [('VERSION_INFO', __version__), ('_CRT_SECURE_NO_WARNINGS', 1)]
if os.name == 'posix':
    os.environ['CC'] = 'g++' 
    LINKERS=['-O3', '-std=c++17']
elif os.name == 'nt':
    LINKERS=['/O2', '-std:c++17']
    pyver=sys.version_info
    if pyver.major == 3 and pyver.minor == 13 and pyver.micro == 4:
        DEFINES.append(('Py_GIL_DISABLED', 1))

ext_modules = [
    Pybind11Extension(
        "PyMSDump",
        src, 
        include_dirs=["MS_dump", "MS_dump/eigen-3.3.9"], 
        define_macros=DEFINES,
        extra_compile_args=LINKERS,
    ),
]

description = ''
try:
    description = open("README.md", 'r', encoding='utf-8').read()
except:
    pass

setup(
    name="PyMSDump",                # 模块名称
    version=__version__,            # 版本号
    author="YujieLiu",
    author_email="",
    description="A parser for Materials Studio .trj format",
    long_description=description,
    long_description_content_type="text/markdown",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,                  # 非压缩文件包
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
    ],
)
