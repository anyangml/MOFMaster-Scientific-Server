"""
Structure Fetcher Tool

Retrieve MOF structural information from the QMOF database using the MOF ID.
"""

import json
from pydantic import ValidationError
from .base import BaseModel, Field, Optional, os, DATA_DIR

class FetchStructureInput(BaseModel):
    """Input model for structure fetching."""
    mof_id: str = Field(..., description="The unique identifier for the material (e.g., qmof-8b5bb88)")

class FetchStructureOutput(BaseModel):
    """Output model for structure fetching."""
    success: bool = Field(..., description="Whether the fetch was successful")
    atoms_dict: Optional[dict] = Field(None, description="ASE Atoms object as dictionary")
    metadata: Optional[dict] = Field(None, description="Metadata dictionary for the MOF")
    error: Optional[str] = Field(None, description="Error message if fetch failed")
    message: str = Field(..., description="Human-readable result message")

# Global caches for large JSON files to speed up subsequent queries
_qmof_metadata_cache = None
_qmof_structs_cache = None

def fetch_structure(mof_id: str) -> dict:
    """
    Fetch MOF structural information and metadata from the QMOF database.
    
    Args:
        mof_id: The unique identifier for the material (e.g., qmof-8b5bb88)
        
    Returns:
        Dictionary containing the structure and metadata
    """
    global _qmof_metadata_cache, _qmof_structs_cache
    
    try:
        validated_input = FetchStructureInput(mof_id=mof_id)
        mof_id = validated_input.mof_id
        
        # Load caches if not loaded
        if _qmof_metadata_cache is None:
            qmof_path = os.path.join(DATA_DIR, "qmof.json")
            if not os.path.exists(qmof_path):
                raise FileNotFoundError(f"{qmof_path} not found.")
            with open(qmof_path, "r") as f:
                _qmof_metadata_cache = {entry['qmof_id']: entry for entry in json.load(f)}
                
        if _qmof_structs_cache is None:
            qmof_structs_path = os.path.join(DATA_DIR, "qmof_structure_data.json")
            if not os.path.exists(qmof_structs_path):
                raise FileNotFoundError(f"{qmof_structs_path} not found.")
            with open(qmof_structs_path, "r") as f:
                _qmof_structs_cache = {entry['qmof_id']: entry['structure'] for entry in json.load(f)}
                
        # Fetch data
        if mof_id not in _qmof_metadata_cache or mof_id not in _qmof_structs_cache:
            return FetchStructureOutput(
                success=False,
                error="MOF ID not found",
                message=f"MOF ID '{mof_id}' not found in the database."
            ).model_dump()
            
        metadata = _qmof_metadata_cache[mof_id]
        struct_dict = _qmof_structs_cache[mof_id]
        
        # We lazily import pymatgen to avoid initial cost or failure if uninstalled
        from pymatgen.core import Structure
        from pymatgen.io.ase import AseAtomsAdaptor
        
        pmg_struct = Structure.from_dict(struct_dict)
        ase_atoms = AseAtomsAdaptor.get_atoms(pmg_struct)
        
        # Convert ASE Atoms to dictionary format expected by other tools
        atoms_dict = {
            "numbers": ase_atoms.get_atomic_numbers().tolist(),
            "positions": ase_atoms.get_positions().tolist(),
            "cell": ase_atoms.get_cell().tolist(),
            "pbc": ase_atoms.get_pbc().tolist()
        }
        
        return FetchStructureOutput(
            success=True,
            atoms_dict=atoms_dict,
            metadata=metadata,
            message=f"Successfully fetched structure and metadata for {mof_id}."
        ).model_dump()
        
    except ValidationError as e:
        return FetchStructureOutput(
            success=False,
            error="Input validation error",
            message=f"Input validation error: {str(e)}"
        ).model_dump()
    except Exception as e:
        return FetchStructureOutput(
            success=False,
            error=str(e),
            message=f"Error fetching structure: {str(e)}"
        ).model_dump()
