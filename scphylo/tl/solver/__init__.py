"""Solver Module."""

from scphylo.tl.solver._bnb import bnb
from scphylo.tl.solver._cardelino import cardelino
from scphylo.tl.solver._dendro import dendro
from scphylo.tl.solver._gpps import gpps
from scphylo.tl.solver._grmt import grmt
from scphylo.tl.solver._onconem import onconem
from scphylo.tl.solver._phiscs import (
    phiscs_readcount,
    phiscsb,
    phiscsb_bulk,
    phiscsi,
    phiscsi_bulk,
)
from scphylo.tl.solver._sbm import sbm
from scphylo.tl.solver._sciphi import sciphi
from scphylo.tl.solver._scistree import iscistree, rscistree, scistree
from scphylo.tl.solver._scite import infscite, scite
from scphylo.tl.solver._siclonefit import siclonefit
from scphylo.tl.solver._sphyr import sphyr
from scphylo.tl.solver.booster import booster
from scphylo.tl.solver.huntress import huntress

__all__ = (
    bnb,
    cardelino,
    dendro,
    gpps,
    grmt,
    onconem,
    phiscs_readcount,
    phiscsb,
    phiscsb_bulk,
    phiscsi,
    phiscsi_bulk,
    sbm,
    sciphi,
    iscistree,
    rscistree,
    scistree,
    infscite,
    scite,
    siclonefit,
    sphyr,
    booster,
    huntress,
)
