"""Setup Module."""

import os
import sys
from pathlib import Path
from sys import platform

from setuptools import find_packages, setup
from setuptools.extension import Extension

if platform == "linux" or platform == "linux2":
    os.environ["CC"] = "g++"
elif platform == "darwin":
    os.environ["CC"] = "clang++"

extensions = [
    Extension(
        "scphylo.external._mltd",
        sources=[
            "src/scphylo/external/_mltd.pyx",
            "src/scphylo/external/mltd/mltd.cpp",
        ],
        include_dirs=["src/scphylo/external/mltd"],
        extra_compile_args=["-std=c++11"],
        language="c++",
    ),
    Extension(
        "scphylo.external._scprob",
        sources=[
            "src/scphylo/external/_scprob.pyx",
            "src/scphylo/external/scprob/main.cpp",
        ],
        include_dirs=["src/scphylo/external/scprob"],
        extra_compile_args=["-std=c++11"],
        language="c++",
    ),
    Extension(
        "scphylo.external._scistree",
        sources=[
            "src/scphylo/external/_scistree.pyx",
            "src/scphylo/external/scistree/Utils.cpp",
            "src/scphylo/external/scistree/Utils2.cpp",
            "src/scphylo/external/scistree/Utils3.cpp",
            "src/scphylo/external/scistree/Utils4.cpp",
            "src/scphylo/external/scistree/UtilsNumerical.cpp",
            "src/scphylo/external/scistree/RerootTreeUtils.cpp",
            "src/scphylo/external/scistree/TreeBuilder.cpp",
            "src/scphylo/external/scistree/UnWeightedGraph.cpp",
            "src/scphylo/external/scistree/MarginalTree.cpp",
            "src/scphylo/external/scistree/RBT.cpp",
            "src/scphylo/external/scistree/PhylogenyTreeBasic.cpp",
            "src/scphylo/external/scistree/PhylogenyTree.cpp",
            "src/scphylo/external/scistree/BioSequenceMatrix.cpp",
            "src/scphylo/external/scistree/BinaryMatrix.cpp",
            "src/scphylo/external/scistree/GenotypeMatrix.cpp",
            "src/scphylo/external/scistree/ScistGenotype.cpp",
            "src/scphylo/external/scistree/ScistPerfPhyUtils.cpp",
            "src/scphylo/external/scistree/ScistPerfPhyImp.cpp",
            "src/scphylo/external/scistree/ScistDoublet.cpp",
            "src/scphylo/external/scistree/ScistErrRateInf.cpp",
            "src/scphylo/external/scistree/main.cpp",
        ],
        include_dirs=["src/scphylo/external/scistree"],
        extra_compile_args=["-O3", "-std=c++11", "-c"],
        language="c++",
    ),
    Extension(
        "scphylo.external._scite",
        sources=[
            "src/scphylo/external/_scite.pyx",
            "src/scphylo/external/scite/matrices.cpp",
            "src/scphylo/external/scite/mcmcBinTreeMove.cpp",
            "src/scphylo/external/scite/mcmc.cpp",
            "src/scphylo/external/scite/mcmcTreeMove.cpp",
            "src/scphylo/external/scite/output.cpp",
            "src/scphylo/external/scite/rand.cpp",
            "src/scphylo/external/scite/scoreBinTree.cpp",
            "src/scphylo/external/scite/scoreTree.cpp",
            "src/scphylo/external/scite/treelist.cpp",
            "src/scphylo/external/scite/trees.cpp",
            "src/scphylo/external/scite/findBestTrees.cpp",
        ],
        include_dirs=["src/scphylo/external/scite"],
        extra_compile_args=["-O3", "-std=c++11", "-c"],
        language="c++",
    ),
]


def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                sfile = path + ".cpp"
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


try:
    from Cython.Build import cythonize

    HAS_CYTHON = True
except ImportError:
    HAS_CYTHON = False

CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 1 if HAS_CYTHON else 0)))
if CYTHONIZE:
    if not HAS_CYTHON:
        sys.stderr.write(
            "Cannot find Cython. Have you installed all the requirements?\n"
            "Try pip install Cython\n"
        )
        sys.exit(1)
    compiler_directives = {"language_level": 2, "embedsignature": True}
    extensions = cythonize(extensions, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(extensions)

if __name__ == "__main__":
    setup(
        name="scphylo-tools",
        version="0.0.5",
        ext_modules=extensions,
        description=Path("README.rst").read_text("utf-8").split("\n")[3],
        package_dir={"": "src"},
        packages=find_packages(where="src"),
    )
