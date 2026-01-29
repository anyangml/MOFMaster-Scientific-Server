"""
Unit tests for MOF Tools with Pydantic validation.

Tests cover the refactored modular tool structure.
"""

import json
import pytest
from typing import Dict, Any

from tools import parse_structure, static_calculation, optimize_geometry
from tool_registry import (
    ToolRegistry,
    ToolCategory,
    ToolMetadata,
    get_registry,
)


class TestStructureParser:
    """Test parse_structure tool."""
    
    def test_parse_structure_valid_cif(self):
        """Test parsing valid CIF data."""
        cif_data = """data_test
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cu1 Cu 0.0 0.0 0.0
Cu2 Cu 0.5 0.5 0.5
"""
        result = parse_structure(cif_data, format="cif")
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["num_atoms"] == 2
        assert "Cu" in parsed["formula"]
        assert parsed["atoms_dict"] is not None
    
    def test_parse_structure_invalid_data(self):
        """Test parsing invalid data."""
        result = parse_structure("invalid data", format="cif")
        parsed = json.loads(result)
        
        assert parsed["success"] is False
        assert parsed["error"] is not None
    
    def test_parse_structure_empty_data(self):
        """Test parsing empty data."""
        result = parse_structure("   ", format="cif")
        parsed = json.loads(result)
        
        assert parsed["success"] is False
        assert "validation error" in parsed["message"].lower()


class TestStaticCalculation:
    """Test static_calculation tool."""
    
    def get_test_atoms_dict(self):
        """Helper to get test atoms dictionary."""
        return {
            "numbers": [29, 29],
            "positions": [[0, 0, 0], [2, 2, 2]],
            "cell": [[4, 0, 0], [0, 4, 0], [0, 0, 4]],
            "pbc": [True, True, True]
        }
    
    def test_static_calculation_basic(self):
        """Test basic static calculation."""
        atoms_dict = self.get_test_atoms_dict()
        result = static_calculation(atoms_dict)
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["total_energy"] is not None
        assert isinstance(parsed["total_energy"], (int, float))
    
    def test_static_calculation_with_forces(self):
        """Test static calculation with forces."""
        atoms_dict = self.get_test_atoms_dict()
        result = static_calculation(atoms_dict, compute_forces=True)
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["forces"] is not None
        assert len(parsed["forces"]) == 2  # 2 atoms
    
    def test_static_calculation_normalize_per_atom(self):
        """Test static calculation with per-atom normalization."""
        atoms_dict = self.get_test_atoms_dict()
        result = static_calculation(atoms_dict, normalize_per_atom=True)
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["energy_per_atom"] is not None
        assert isinstance(parsed["energy_per_atom"], (int, float))
    
    def test_static_calculation_with_virial(self):
        """Test static calculation with virial."""
        atoms_dict = self.get_test_atoms_dict()
        result = static_calculation(atoms_dict, compute_virial=True)
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        # Virial may or may not be computed depending on calculator
        assert "virial" in parsed
    
    def test_static_calculation_invalid_atoms_dict(self):
        """Test static calculation with invalid atoms dict."""
        result = static_calculation({"invalid": "data"})
        parsed = json.loads(result)
        
        assert parsed["success"] is False
        assert parsed["error"] is not None


class TestGeometryOptimization:
    """Test optimize_geometry tool."""
    
    def get_test_atoms_dict(self):
        """Helper to get test atoms dictionary."""
        return {
            "numbers": [29, 29],
            "positions": [[0, 0, 0], [2, 2, 2]],
            "cell": [[4, 0, 0], [0, 4, 0], [0, 0, 4]],
            "pbc": [True, True, True]
        }
    
    def test_optimize_geometry_basic(self):
        """Test basic geometry optimization."""
        atoms_dict = self.get_test_atoms_dict()
        result = optimize_geometry(atoms_dict)
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["optimized_atoms_dict"] is not None
        assert parsed["metadata"] is not None
    
    def test_optimize_geometry_metadata(self):
        """Test optimization metadata."""
        atoms_dict = self.get_test_atoms_dict()
        result = optimize_geometry(atoms_dict, fmax=0.05, max_steps=10)
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        metadata = parsed["metadata"]
        assert "converged" in metadata
        assert "final_fmax" in metadata
        assert "steps" in metadata
        assert "initial_energy" in metadata
        assert "final_energy" in metadata
    
    def test_optimize_geometry_different_optimizers(self):
        """Test different optimizer types."""
        atoms_dict = self.get_test_atoms_dict()
        
        for optimizer in ["BFGS", "LBFGS", "FIRE"]:
            result = optimize_geometry(atoms_dict, optimizer=optimizer, max_steps=5)
            parsed = json.loads(result)
            assert parsed["success"] is True
    
    def test_optimize_geometry_invalid_optimizer(self):
        """Test with invalid optimizer."""
        atoms_dict = self.get_test_atoms_dict()
        result = optimize_geometry(atoms_dict, optimizer="INVALID")
        parsed = json.loads(result)
        
        assert parsed["success"] is False
        assert "validation error" in parsed["message"].lower()
    
    def test_optimize_geometry_invalid_atoms_dict(self):
        """Test optimization with invalid atoms dict."""
        result = optimize_geometry({"invalid": "data"})
        parsed = json.loads(result)
        
        assert parsed["success"] is False
        assert parsed["error"] is not None

    def test_optimize_geometry_relax_cell(self):
        """Test geometry optimization with cell relaxation."""
        atoms_dict = self.get_test_atoms_dict()
        result = optimize_geometry(atoms_dict, relax_cell=True, max_steps=2)
        parsed = json.loads(result)
        
        # This may fail if FrechetCellFilter is not available, but logic is correct
        assert parsed["success"] is True
        assert "cell" in parsed["optimized_atoms_dict"]

    def test_optimize_geometry_fix_symmetry(self):
        """Test geometry optimization with fix_symmetry."""
        atoms_dict = self.get_test_atoms_dict()
        result = optimize_geometry(atoms_dict, fix_symmetry=True, max_steps=2)
        parsed = json.loads(result)
        
        assert parsed["success"] is True


class TestToolRegistry:
    """Test tool registry functionality."""
    
    def test_registry_initialization(self):
        """Test creating a new registry."""
        registry = ToolRegistry()
        assert len(registry) == 0
    
    def test_register_tool(self):
        """Test registering a tool."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        metadata = registry.register(
            name="test_tool",
            description="A test tool",
            category=ToolCategory.CALCULATION,
            function=dummy_func
        )
        
        assert len(registry) == 1
        assert "test_tool" in registry
        assert metadata.name == "test_tool"
    
    def test_register_duplicate_tool(self):
        """Test registering a duplicate tool name."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="test_tool",
            description="A test tool",
            category=ToolCategory.CALCULATION,
            function=dummy_func
        )
        
        with pytest.raises(ValueError, match="already registered"):
            registry.register(
                name="test_tool",
                description="Another tool",
                category=ToolCategory.UTILS,
                function=dummy_func
            )
    
    def test_get_by_category(self):
        """Test getting tools by category."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="calc_tool",
            description="Calculate",
            category=ToolCategory.CALCULATION,
            function=dummy_func
        )
        
        registry.register(
            name="utils_tool",
            description="Utility",
            category=ToolCategory.UTILS,
            function=dummy_func
        )
        
        calc_tools = registry.get_by_category(ToolCategory.CALCULATION)
        assert len(calc_tools) == 1
        assert calc_tools[0].name == "calc_tool"
    
    def test_unregister_tool(self):
        """Test unregistering a tool."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="test_tool",
            description="A test tool",
            category=ToolCategory.CALCULATION,
            function=dummy_func
        )
        
        assert len(registry) == 1
        result = registry.unregister("test_tool")
        assert result is True
        assert len(registry) == 0


class TestYAMLConfiguration:
    """Test YAML configuration loading."""
    
    def test_load_tool_definitions(self):
        """Test loading tool definitions from YAML."""
        from main import load_tool_definitions
        
        definitions = load_tool_definitions()
        
        # Check that we loaded the correct number of tools
        assert len(definitions) == 3
        
        # Check that all expected tools are present
        tool_names = [d['name'] for d in definitions]
        assert 'Static Calculation' in tool_names
        assert 'Structure Parser' in tool_names
        assert 'Geometry Optimization' in tool_names
    
    def test_tool_definition_structure(self):
        """Test that tool definitions have the correct structure."""
        from main import load_tool_definitions
        
        definitions = load_tool_definitions()
        
        for tool_def in definitions:
            # Check required fields
            assert 'name' in tool_def
            assert 'description' in tool_def
            assert 'category' in tool_def
            assert 'function_name' in tool_def
            assert 'tags' in tool_def
            assert 'version' in tool_def
            
            # Check types
            assert isinstance(tool_def['name'], str)
            assert isinstance(tool_def['description'], str)
            assert isinstance(tool_def['category'], str)
            assert isinstance(tool_def['function_name'], str)
            assert isinstance(tool_def['tags'], list)
            assert isinstance(tool_def['version'], str)
    
    def test_tool_registration_from_yaml(self):
        """Test that tools can be registered from YAML definitions."""
        from main import register_tools_in_registry
        
        # Clear registry first
        registry = get_registry()
        registry.clear()
        
        # Register tools from YAML
        register_tools_in_registry()
        
        # Verify all tools were registered
        assert len(registry) == 3


class TestIntegration:
    """Integration tests for complete workflows."""
    
    def test_parse_and_calculate_workflow(self):
        """Test parsing a structure and calculating energy."""
        # Parse structure
        cif_data = """data_test
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cu1 Cu 0.0 0.0 0.0
Cu2 Cu 0.5 0.5 0.5
"""
        parse_result = parse_structure(cif_data, format="cif")
        parsed = json.loads(parse_result)
        assert parsed["success"] is True
        
        # Calculate energy
        atoms_dict = parsed["atoms_dict"]
        calc_result = static_calculation(atoms_dict, compute_forces=True)
        calc = json.loads(calc_result)
        assert calc["success"] is True
        assert calc["total_energy"] is not None
    
    def test_full_workflow(self):
        """Test complete workflow: parse -> calculate -> optimize."""
        # Parse structure
        cif_data = """data_test
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cu1 Cu 0.0 0.0 0.0
Cu2 Cu 0.5 0.5 0.5
"""
        parse_result = parse_structure(cif_data, format="cif")
        parsed = json.loads(parse_result)
        assert parsed["success"] is True
        
        atoms_dict = parsed["atoms_dict"]
        
        # Static calculation
        calc_result = static_calculation(atoms_dict)
        calc = json.loads(calc_result)
        assert calc["success"] is True
        
        # Geometry optimization
        opt_result = optimize_geometry(atoms_dict, max_steps=5)
        opt = json.loads(opt_result)
        assert opt["success"] is True
        assert opt["metadata"]["converged"] is not None


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
