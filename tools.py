from io import StringIO
from pathlib import Path

# Optional ASE integration
try:
    import ase.io
    from ase.calculators.emt import EMT
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

# Internal database
SAMPLE_MOF_DB = [
    {"name": "HKUST-1", "formula": "Cu3(BTC)2", "surface_area": 1850},
    {"name": "MOF-5", "formula": "Zn4O(BDC)3", "surface_area": 3800},
    {"name": "UiO-66", "formula": "Zr6O4(OH)4(BDC)6", "surface_area": 1187},
]

def search_mofs(query: str) -> str:
    """Search for MOFs by name or formula."""
    results = [m for m in SAMPLE_MOF_DB if query.lower() in m["name"].lower() or query.lower() in m["formula"].lower()]
    return str(results) if results else f"No MOFs found for '{query}'"

def calculate_energy(data: str) -> str:
    """Calculate the potential energy of a structure (CIF string or path)."""
    if not ASE_AVAILABLE:
        return "Error: ASE library not installed."
    try:
        fileobj = StringIO(data) if "\n" in data else data
        atoms = ase.io.read(fileobj, format="cif" if "\n" in data else None)
        atoms.calc = EMT()
        return f"Energy: {atoms.get_potential_energy():.4f} eV"
    except Exception as e:
        return f"Error: {str(e)}"

def optimize_structure(name: str) -> str:
    """Perform structure optimization (standard placeholder)."""
    return f"Successfully initiated optimization for {name}"
