import os
import sys
from pathlib import Path
from typing import Dict, Any, List
import json

# Try to import scientific libs
try:
    import ase.io
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS
except ImportError:
    print("Warning: ASE not found. Some tools will fail.", file=sys.stderr)

from mcp.server.fastmcp import FastMCP

# Initialize FastMCP server
mcp = FastMCP("MOFMaster-Scientific-Server")

# Configuration - Shared data directory
# In a real cloud setup, this would be a cloud bucket or shared volume
DATA_DIR = Path(os.getenv("DATA_DIR", "./data"))
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Sample MOF database
SAMPLE_MOF_DB = [
    {
        "mof_name": "HKUST-1",
        "formula": "Cu3(BTC)2",
        "description": "Copper-based MOF with high surface area",
        "tags": ["copper", "high surface area", "paddle-wheel"],
        "cif_filename": "HKUST-1.cif",
        "properties": {"surface_area_m2g": 1850, "pore_volume_cm3g": 0.75},
    },
    {
        "mof_name": "MOF-5",
        "formula": "Zn4O(BDC)3",
        "description": "Zinc-based MOF, one of the first MOFs discovered",
        "tags": ["zinc", "BDC", "cubic"],
        "cif_filename": "MOF-5.cif",
        "properties": {"surface_area_m2g": 3800, "pore_volume_cm3g": 1.55},
    },
    {
        "mof_name": "UiO-66",
        "formula": "Zr6O4(OH)4(BDC)6",
        "description": "Zirconium-based MOF with exceptional stability",
        "tags": ["zirconium", "stable", "water-stable"],
        "cif_filename": "UiO-66.cif",
        "properties": {"surface_area_m2g": 1187, "pore_volume_cm3g": 0.44},
    },
    {
        "mof_name": "MIL-101",
        "formula": "Cr3F(H2O)2O[(O2C)-C6H4-(CO2)]3",
        "description": "Chromium-based MOF with very high surface area",
        "tags": ["chromium", "very high surface area", "mesoporous"],
        "cif_filename": "MIL-101.cif",
        "properties": {"surface_area_m2g": 4100, "pore_volume_cm3g": 2.15},
    },
]

@mcp.tool()
def search_mof_db(query_string: str) -> str:
    """
    Search for MOF structures in the database.
    Returns: JSON string containing MOF metadata.
    """
    query_lower = query_string.lower()
    for mof in SAMPLE_MOF_DB:
        if (
            query_lower in mof["mof_name"].lower()
            or query_lower in mof["formula"].lower()
            or query_lower in mof["description"].lower()
            or any(query_lower in tag.lower() for tag in mof["tags"])
        ):
            cif_path = DATA_DIR / mof["cif_filename"]
            if not cif_path.exists():
                cif_content = f"""data_{mof['mof_name']}
_cell_length_a    26.343
_cell_length_b    26.343
_cell_length_c    26.343
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'F m -3 m'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cu1 Cu 0.250 0.250 0.250
O1  O  0.200 0.200 0.200
C1  C  0.150 0.150 0.150
"""
                with open(cif_path, "w") as f:
                    f.write(cif_content)

            result = {
                "mof_name": mof["mof_name"],
                "formula": mof["formula"],
                "description": mof["description"],
                "cif_filepath": str(cif_path.absolute()),
                "properties": mof["properties"],
            }
            return json.dumps(result)
    
    return json.dumps({"error": f"No MOF found matching query: {query_string}"})

@mcp.tool()
def optimize_structure_ase(cif_filepath: str) -> str:
    """
    Optimize the geometry of a MOF structure using ASE.
    """
    try:
        atoms = ase.io.read(cif_filepath)
        atoms.calc = EMT()
        initial_energy = atoms.get_potential_energy()
        optimizer = BFGS(atoms, logfile=None)
        optimizer.run(fmax=0.05)
        final_energy = atoms.get_potential_energy()
        
        input_path = Path(cif_filepath)
        output_filename = input_path.stem + "_optimized.cif"
        output_path = DATA_DIR / output_filename
        ase.io.write(str(output_path), atoms, format="cif")

        return json.dumps({
            "optimized_cif_filepath": str(output_path.absolute()),
            "initial_energy_ev": float(initial_energy),
            "final_energy_ev": float(final_energy),
            "energy_change_ev": float(final_energy - initial_energy),
            "n_steps": optimizer.get_number_of_steps(),
            "converged": True,
        })
    except Exception as e:
        return json.dumps({"error": f"Optimization failed: {str(e)}", "cif_filepath": cif_filepath})

@mcp.tool()
def calculate_energy_force(cif_filepath: str) -> str:
    """
    Calculate the energy and forces of a structure using ASE.
    """
    try:
        atoms = ase.io.read(cif_filepath)
        atoms.calc = EMT()
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        max_force = float(abs(forces).max())

        return json.dumps({
            "energy_ev": float(energy),
            "max_force_ev_ang": max_force,
            "cif_filepath": str(Path(cif_filepath).absolute()),
            "n_atoms": len(atoms),
        })
    except Exception as e:
        return json.dumps({"error": f"Energy calculation failed: {str(e)}", "cif_filepath": cif_filepath})

if __name__ == "__main__":
    mcp.run()
