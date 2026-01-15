"""
Unit tests for MOF Tools with Pydantic validation.

Tests cover input validation, output validation, error handling, and tool registry.
"""

import json
import pytest
from typing import Dict, Any

import tools
from tools import (
    SearchMOFsInput,
    SearchMOFsOutput,
    CalculateEnergyInput,
    CalculateEnergyOutput,
    OptimizeStructureInput,
    OptimizeStructureOutput,
    MOFRecord,
)
from tool_registry import (
    ToolRegistry,
    ToolCategory,
    ToolMetadata,
    get_registry,
)
from pydantic import ValidationError


class TestPydanticModels:
    """Test Pydantic model validation."""
    
    def test_mof_record_valid(self):
        """Test valid MOF record creation."""
        record = MOFRecord(
            name="Test-MOF",
            formula="Cu3(BTC)2",
            surface_area=1850.0
        )
        assert record.name == "Test-MOF"
        assert record.formula == "Cu3(BTC)2"
        assert record.surface_area == 1850.0
    
    def test_mof_record_invalid_surface_area(self):
        """Test MOF record with invalid surface area."""
        with pytest.raises(ValidationError):
            MOFRecord(
                name="Test-MOF",
                formula="Cu3(BTC)2",
                surface_area=-100.0  # Negative surface area
            )
    
    def test_search_mofs_input_valid(self):
        """Test valid search input."""
        input_data = SearchMOFsInput(query="MOF-5")
        assert input_data.query == "MOF-5"
    
    def test_search_mofs_input_strips_whitespace(self):
        """Test that search input strips whitespace."""
        input_data = SearchMOFsInput(query="  MOF-5  ")
        assert input_data.query == "MOF-5"
    
    def test_search_mofs_input_empty_query(self):
        """Test search input with empty query."""
        with pytest.raises(ValidationError):
            SearchMOFsInput(query="   ")
    
    def test_search_mofs_input_too_long(self):
        """Test search input with query too long."""
        with pytest.raises(ValidationError):
            SearchMOFsInput(query="x" * 201)
    
    def test_calculate_energy_input_valid(self):
        """Test valid energy calculation input."""
        input_data = CalculateEnergyInput(data="some cif content")
        assert input_data.data == "some cif content"
    
    def test_calculate_energy_input_empty(self):
        """Test energy calculation input with empty data."""
        with pytest.raises(ValidationError):
            CalculateEnergyInput(data="   ")
    
    def test_optimize_structure_input_valid(self):
        """Test valid structure optimization input."""
        input_data = OptimizeStructureInput(name="test-structure")
        assert input_data.name == "test-structure"
    
    def test_optimize_structure_input_strips_whitespace(self):
        """Test that optimization input strips whitespace."""
        input_data = OptimizeStructureInput(name="  test-structure  ")
        assert input_data.name == "test-structure"
    
    def test_optimize_structure_input_empty(self):
        """Test optimization input with empty name."""
        with pytest.raises(ValidationError):
            OptimizeStructureInput(name="   ")


class TestSearchMOFs:
    """Test search_mofs tool."""
    
    def test_search_mofs_found(self):
        """Test searching for existing MOF."""
        result = tools.search_mofs("MOF-5")
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["count"] == 1
        assert len(parsed["results"]) == 1
        assert parsed["results"][0]["name"] == "MOF-5"
    
    def test_search_mofs_not_found(self):
        """Test searching for non-existing MOF."""
        result = tools.search_mofs("NonExistent")
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["count"] == 0
        assert len(parsed["results"]) == 0
        assert "No MOFs found" in parsed["message"]
    
    def test_search_mofs_multiple_results(self):
        """Test searching with multiple results."""
        result = tools.search_mofs("Cu")
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["count"] >= 1
    
    def test_search_mofs_case_insensitive(self):
        """Test case-insensitive search."""
        result = tools.search_mofs("mof-5")
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["count"] == 1
    
    def test_search_mofs_invalid_input(self):
        """Test search with invalid input."""
        result = tools.search_mofs("   ")
        parsed = json.loads(result)
        
        assert parsed["success"] is False
        assert "validation error" in parsed["message"].lower()


class TestCalculateEnergy:
    """Test calculate_energy tool."""
    
    def test_calculate_energy_ase_not_available(self):
        """Test energy calculation when ASE is not available."""
        # This test will vary depending on whether ASE is installed
        result = tools.calculate_energy("test data")
        parsed = json.loads(result)
        
        assert "success" in parsed
        assert "message" in parsed
    
    def test_calculate_energy_invalid_input(self):
        """Test energy calculation with invalid input."""
        result = tools.calculate_energy("   ")
        parsed = json.loads(result)
        
        assert parsed["success"] is False
        assert "validation error" in parsed["message"].lower()
    
    def test_calculate_energy_error_handling(self):
        """Test energy calculation error handling."""
        result = tools.calculate_energy("invalid cif data")
        parsed = json.loads(result)
        
        # Should handle gracefully
        assert "success" in parsed
        assert "message" in parsed


class TestOptimizeStructure:
    """Test optimize_structure tool."""
    
    def test_optimize_structure_valid(self):
        """Test structure optimization with valid input."""
        result = tools.optimize_structure("test-structure")
        parsed = json.loads(result)
        
        assert parsed["success"] is True
        assert parsed["structure_name"] == "test-structure"
        assert "Successfully initiated optimization" in parsed["message"]
    
    def test_optimize_structure_invalid_input(self):
        """Test optimization with invalid input."""
        result = tools.optimize_structure("   ")
        parsed = json.loads(result)
        
        assert parsed["success"] is False
        assert "validation error" in parsed["message"].lower()


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
            category=ToolCategory.SEARCH,
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
            category=ToolCategory.SEARCH,
            function=dummy_func
        )
        
        with pytest.raises(ValueError, match="already registered"):
            registry.register(
                name="test_tool",
                description="Another tool",
                category=ToolCategory.CALCULATION,
                function=dummy_func
            )
    
    def test_get_tool(self):
        """Test retrieving a tool."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="test_tool",
            description="A test tool",
            category=ToolCategory.SEARCH,
            function=dummy_func
        )
        
        metadata = registry.get("test_tool")
        assert metadata is not None
        assert metadata.name == "test_tool"
    
    def test_get_nonexistent_tool(self):
        """Test retrieving a non-existent tool."""
        registry = ToolRegistry()
        metadata = registry.get("nonexistent")
        assert metadata is None
    
    def test_get_by_category(self):
        """Test getting tools by category."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="search_tool",
            description="Search",
            category=ToolCategory.SEARCH,
            function=dummy_func
        )
        
        registry.register(
            name="calc_tool",
            description="Calculate",
            category=ToolCategory.CALCULATION,
            function=dummy_func
        )
        
        search_tools = registry.get_by_category(ToolCategory.SEARCH)
        assert len(search_tools) == 1
        assert search_tools[0].name == "search_tool"
    
    def test_get_by_tag(self):
        """Test getting tools by tag."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="tool1",
            description="Tool 1",
            category=ToolCategory.SEARCH,
            function=dummy_func,
            tags=["mof", "database"]
        )
        
        registry.register(
            name="tool2",
            description="Tool 2",
            category=ToolCategory.CALCULATION,
            function=dummy_func,
            tags=["energy", "mof"]
        )
        
        mof_tools = registry.get_by_tag("mof")
        assert len(mof_tools) == 2
    
    def test_unregister_tool(self):
        """Test unregistering a tool."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="test_tool",
            description="A test tool",
            category=ToolCategory.SEARCH,
            function=dummy_func
        )
        
        assert len(registry) == 1
        result = registry.unregister("test_tool")
        assert result is True
        assert len(registry) == 0
    
    def test_unregister_nonexistent_tool(self):
        """Test unregistering a non-existent tool."""
        registry = ToolRegistry()
        result = registry.unregister("nonexistent")
        assert result is False
    
    def test_clear_registry(self):
        """Test clearing the registry."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="tool1",
            description="Tool 1",
            category=ToolCategory.SEARCH,
            function=dummy_func
        )
        
        registry.register(
            name="tool2",
            description="Tool 2",
            category=ToolCategory.CALCULATION,
            function=dummy_func
        )
        
        assert len(registry) == 2
        registry.clear()
        assert len(registry) == 0
    
    def test_list_names(self):
        """Test listing tool names."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="tool1",
            description="Tool 1",
            category=ToolCategory.SEARCH,
            function=dummy_func
        )
        
        registry.register(
            name="tool2",
            description="Tool 2",
            category=ToolCategory.CALCULATION,
            function=dummy_func
        )
        
        names = registry.list_names()
        assert "tool1" in names
        assert "tool2" in names
        assert len(names) == 2
    
    def test_list_categories(self):
        """Test listing categories."""
        registry = ToolRegistry()
        
        def dummy_func():
            return "test"
        
        registry.register(
            name="search_tool",
            description="Search",
            category=ToolCategory.SEARCH,
            function=dummy_func
        )
        
        categories = registry.list_categories()
        assert categories[ToolCategory.SEARCH] == 1
        assert categories[ToolCategory.CALCULATION] == 0




if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
