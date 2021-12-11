"""Input/Output Module."""

from scphylo.io._genotype import read, write
from scphylo.io._tree import to_png

__all__ = (
    read,
    write,
    to_png,
)
