"""quadratic_sos — exact integral lengths of sums of squares in real quadratic rings.

This package implements the algorithms described in the paper
*Computing Integral Length of Sums of Squares and Pythagoras Elements
in Real Quadratic Rings*.
"""

from quadratic_sos.ring import QuadraticInteger, RealQuadraticOrder
from quadratic_sos.length import exact_length, is_representable, decomposition
from quadratic_sos.pythagoras import pythagoras_number, pythagoras_element
from quadratic_sos.experiments import search_profile, validation_sweep

__all__ = [
    "QuadraticInteger",
    "RealQuadraticOrder",
    "exact_length",
    "is_representable",
    "decomposition",
    "pythagoras_number",
    "pythagoras_element",
    "search_profile",
    "validation_sweep",
]
