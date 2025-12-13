"""Setup Module."""

import os
import sys
from pathlib import Path
from sys import platform

from setuptools import setup
from setuptools.extension import Extension

if platform == "linux" or platform == "linux2":
    os.environ["CC"] = "g++"
elif platform == "darwin":
    os.environ["CC"] = "clang++"
extensions = [
    Extension(
        "scphylo.external._mltd",
        sources=["scphylo/external/_mltd.pyx", "scphylo/external/mltd/mltd.cpp"],
        include_dirs=["scphylo/external/mltd"],
        extra_compile_args=["-std=c++11"],
        language="c++",
    ),
    Extension(
        "scphylo.external._scprob",
        sources=[
            "scphylo/external/_scprob.pyx",
            "scphylo/external/scprob/main.cpp",
        ],
        include_dirs=["scphylo/external/scprob"],
        extra_compile_args=["-std=c++11"],
        language="c++",
    ),
    Extension(
        "scphylo.external._scistree",
        sources=[
            "scphylo/external/_scistree.pyx",
            "scphylo/external/scistree/Utils.cpp",
            "scphylo/external/scistree/Utils2.cpp",
            "scphylo/external/scistree/Utils3.cpp",
            "scphylo/external/scistree/Utils4.cpp",
            "scphylo/external/scistree/UtilsNumerical.cpp",
            "scphylo/external/scistree/RerootTreeUtils.cpp",
            "scphylo/external/scistree/TreeBuilder.cpp",
            "scphylo/external/scistree/UnWeightedGraph.cpp",
            "scphylo/external/scistree/MarginalTree.cpp",
            "scphylo/external/scistree/RBT.cpp",
            "scphylo/external/scistree/PhylogenyTreeBasic.cpp",
            "scphylo/external/scistree/PhylogenyTree.cpp",
            "scphylo/external/scistree/BioSequenceMatrix.cpp",
            "scphylo/external/scistree/BinaryMatrix.cpp",
            "scphylo/external/scistree/GenotypeMatrix.cpp",
            "scphylo/external/scistree/ScistGenotype.cpp",
            "scphylo/external/scistree/ScistPerfPhyUtils.cpp",
            "scphylo/external/scistree/ScistPerfPhyImp.cpp",
            "scphylo/external/scistree/ScistDoublet.cpp",
            "scphylo/external/scistree/ScistErrRateInf.cpp",
            "scphylo/external/scistree/main.cpp",
        ],
        include_dirs=["scphylo/external/scistree"],
        extra_compile_args=["-O3", "-std=c++11", "-c"],
        language="c++",
    ),
    Extension(
        "scphylo.external._scite",
        sources=[
            "scphylo/external/_scite.pyx",
            "scphylo/external/scite/matrices.cpp",
            "scphylo/external/scite/mcmcBinTreeMove.cpp",
            "scphylo/external/scite/mcmc.cpp",
            "scphylo/external/scite/mcmcTreeMove.cpp",
            "scphylo/external/scite/output.cpp",
            "scphylo/external/scite/rand.cpp",
            "scphylo/external/scite/scoreBinTree.cpp",
            "scphylo/external/scite/scoreTree.cpp",
            "scphylo/external/scite/treelist.cpp",
            "scphylo/external/scite/trees.cpp",
            "scphylo/external/scite/findBestTrees.cpp",
        ],
        include_dirs=["scphylo/external/scite"],
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
        ext_modules=extensions,
        description=Path("README.rst").read_text("utf-8").split("\n")[3],
    )
