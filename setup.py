"""Setup Module."""

import os
import sys
from pathlib import Path
from sys import platform

from setuptools import find_packages, setup
from setuptools.extension import Extension

try:
    from scphylo import __author__, __email__, __maintainer__, __version__
except ImportError:
    __author__ = ", ".join(["Farid Rashidi"])
    __maintainer__ = ", ".join(["Farid Rashidi"])
    __email__ = ", ".join(["farid.rsh@gmail.com"])
    __version__ = "0.0.2"

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


CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0)))
if CYTHONIZE:
    try:
        from Cython.Build import cythonize
    except ImportError:
        sys.stderr.write(
            "Cannot find Cython. Have you installed all the requirements?\n"
            "Try pip install -r requirements.txt\n"
        )
        sys.exit(1)
    compiler_directives = {"language_level": 2, "embedsignature": True}
    extensions = cythonize(extensions, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(extensions)

if __name__ == "__main__":
    setup(
        name="scphylo-tools",
        entry_points="""
            [console_scripts]
            scphylo=scphylo.commands.scphylo:cli
        """,
        use_scm_version=True,
        setup_requires=["setuptools_scm"],
        python_requires=">=3.6",
        install_requires=[
            r.strip() for r in Path("requirements.txt").read_text("utf-8").splitlines()
        ],
        ext_modules=extensions,
        extras_require={
            "dev": [
                "pre-commit",
                "pytest-cov",
                "pytest-xdist",
                "bumpversion",
            ],
            "docs": [
                r.strip()
                for r in (Path("docs") / "requirements.txt")
                .read_text("utf-8")
                .splitlines()
                if not r.startswith("-r")
            ],
        },
        platforms=["Linux", "MacOSX"],
        packages=find_packages(),
        include_package_data=True,
        author=__author__,
        author_email=__email__,
        email=__email__,
        maintainer=__maintainer__,
        maintainer_email=__email__,
        version=__version__,
        description=Path("README.rst").read_text("utf-8").split("\n")[3],
        long_description=Path("README.rst").read_text("utf-8"),
        long_description_content_type="text/x-rst; charset=UTF-8",
        license="BSD",
        url="https://github.com/faridrashidi/scphylo-tools",
        project_urls={
            "Documentation": "https://scphylo-tools.readthedocs.io/en/latest",
            "Source Code": "https://github.com/faridrashidi/scphylo-tools",
        },
        download_url="https://github.com/faridrashidi/scphylo-tools",
        keywords=[
            "tumor phylogeny",
            "single cell",
            "tools",
        ],
        classifiers=[
            "License :: OSI Approved :: BSD License",
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Framework :: Jupyter",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Visualization",
        ],
    )
