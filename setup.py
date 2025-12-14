import glob
import os
import sys

from Cython.Build import cythonize
from setuptools import Extension, setup

# Fix for macOS build with PyPy or specific environments
if sys.platform == "darwin":
    os.environ["LDFLAGS"] = os.environ.get("LDFLAGS", "").replace(
        "-Wl,-Bsymbolic-functions", ""
    )

extensions = [
    Extension(
        "scphylo.external._mltd",
        ["src/scphylo/external/_mltd.pyx", "src/scphylo/external/mltd/mltd.cpp"],
        include_dirs=["src/scphylo/external/mltd"],
        language="c++",
        extra_compile_args=["-std=c++11"],
    ),
    Extension(
        "scphylo.external._scistree",
        ["src/scphylo/external/_scistree.pyx"]
        + glob.glob("src/scphylo/external/scistree/*.cpp"),
        include_dirs=["src/scphylo/external/scistree"],
        language="c++",
        extra_compile_args=["-std=c++11"],
    ),
    Extension(
        "scphylo.external._scite",
        ["src/scphylo/external/_scite.pyx"]
        + glob.glob("src/scphylo/external/scite/*.cpp"),
        include_dirs=["src/scphylo/external/scite"],
        language="c++",
        extra_compile_args=["-std=c++11"],
    ),
    Extension(
        "scphylo.external._scprob",
        ["src/scphylo/external/_scprob.pyx", "src/scphylo/external/scprob/main.cpp"],
        include_dirs=["src/scphylo/external/scprob"],
        language="c++",
        extra_compile_args=["-std=c++11"],
    ),
]

setup(
    ext_modules=cythonize(extensions),
)
