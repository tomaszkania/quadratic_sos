"""Exact integral lengths of sums of squares in real quadratic rings.

The package implements the algorithms described in the accompanying TOMS
submission, *Computing Integral Length of Sums of Squares and Pythagoras
Elements in Real Quadratic Rings*.
"""

from quadratic_sos.enumeration import is_representable
from quadratic_sos.experiments import search_profile, validation_sweep
from quadratic_sos.length import decomposition, exact_length
from quadratic_sos.pythagoras import pythagoras_element, pythagoras_number
from quadratic_sos.ring import QuadraticInteger, RealQuadraticOrder

__version__ = "0.2.0"

__all__ = [
    "QuadraticInteger",
    "RealQuadraticOrder",
    "decomposition",
    "exact_length",
    "is_representable",
    "pythagoras_element",
    "pythagoras_number",
    "search_profile",
    "validation_sweep",
]
