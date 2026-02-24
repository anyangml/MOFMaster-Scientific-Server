"""
Base module for MOF tools.

Provides common imports, constants, and base classes used across all tools.
"""

import os
from io import StringIO
from typing import List, Optional

from pydantic import BaseModel, Field, field_validator, ValidationError, ConfigDict

import ase.io
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS, LBFGS, FIRE
from ase.filters import FrechetCellFilter
from ase.constraints import FixSymmetry
from deepmd.calculator import DP
from deepmd.pt.infer.deep_eval import DeepProperty

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

__all__ = ['Atoms', 'EMT', 'BFGS', 'LBFGS', 'FIRE', 'FrechetCellFilter', 'FixSymmetry',
           'BaseModel', 'Field', 'field_validator', 'ValidationError', 
           'ConfigDict', 'Optional', 'List', 'StringIO', 'ase', 'DP', 'DeepProperty', 'os', 'DATA_DIR']
