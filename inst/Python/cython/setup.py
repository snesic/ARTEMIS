from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

setup(
    name="TSW_scoreMat",
    ext_modules=cythonize([
        Extension(
            name="TSW_scoreMat",
            sources=["TSW_scoreMat.pyx"],
            include_dirs=[numpy.get_include()],
            extra_compile_args=["-O3"]
        )
    ],
    compiler_directives={
        "boundscheck": False,
        "wraparound": False,
        "nonecheck": True,
        "cdivision": True,
        "language_level": 3
    })
)