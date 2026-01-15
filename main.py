"""
MOF Tools MCP Server - Main entry point.

Production-ready MCP server with formal tool registration and Pydantic validation.
"""

import sys
from typing import List

from mcp.server.fastmcp import FastMCP

from tool_registry import ToolRegistry, ToolCategory, get_registry
import tools


def initialize_server() -> FastMCP:
    """
    Initialize the FastMCP server with production settings.
    
    Returns:
        Configured FastMCP server instance
    """
    # Initialize server with host/port defaults for HTTP transport
    mcp = FastMCP(
        "mof-tools",
        host="0.0.0.0",
        port=8080,
        log_level="ERROR"  # Suppress info logs that might corrupt the protocol stream
    )
    
    return mcp


def register_tools_in_registry():
    """
    Register all available tools in the formal tool registry.
    
    This provides a centralized, organized way to manage tools with metadata.
    """
    registry = get_registry()
    
    # Define tools to register with their metadata
    tool_definitions = [
        {
            "name": "search_mofs",
            "description": "Search for Metal-Organic Frameworks by name or formula in the database",
            "category": ToolCategory.SEARCH,
            "function": tools.search_mofs,
            "requires_ase": False,
            "is_experimental": False,
            "tags": ["mof", "database", "search", "query"],
            "version": "1.0.0"
        },
        {
            "name": "calculate_energy",
            "description": "Calculate the potential energy of a structure from CIF file content or path using ASE",
            "category": ToolCategory.CALCULATION,
            "function": tools.calculate_energy,
            "requires_ase": True,
            "is_experimental": False,
            "tags": ["energy", "calculation", "ase", "cif"],
            "version": "1.0.0"
        },
        {
            "name": "optimize_structure",
            "description": "Perform structure optimization for a named MOF structure (placeholder implementation)",
            "category": ToolCategory.OPTIMIZATION,
            "function": tools.optimize_structure,
            "requires_ase": False,
            "is_experimental": True,
            "tags": ["optimization", "structure", "geometry"],
            "version": "1.0.0"
        }
    ]
    
    # Register all tools using a loop
    for tool_def in tool_definitions:
        registry.register(**tool_def)


def register_tools_with_mcp(mcp: FastMCP, registry: ToolRegistry):
    """
    Register tools from the registry with the FastMCP server.
    
    Args:
        mcp: FastMCP server instance
        registry: Tool registry containing registered tools
    """
    for tool_metadata in registry.get_all():
        # Register each tool function with FastMCP
        mcp.tool()(tool_metadata.function)


def print_registered_tools(registry: ToolRegistry):
    """
    Print information about registered tools for debugging/info purposes.
    
    Args:
        registry: Tool registry to display
    """
    print(f"\n=== MOF Tools Server ===", file=sys.stderr)
    print(f"Registered {len(registry)} tools:", file=sys.stderr)
    
    for tool_metadata in registry.get_all():
        status = "[EXPERIMENTAL]" if tool_metadata.is_experimental else ""
        ase_req = "[ASE Required]" if tool_metadata.requires_ase else ""
        print(
            f"  - {tool_metadata.name} ({tool_metadata.category.value}) "
            f"{status} {ase_req}",
            file=sys.stderr
        )
        print(f"    {tool_metadata.description}", file=sys.stderr)
    
    print(f"\nTools by category:", file=sys.stderr)
    for category, count in registry.list_categories().items():
        print(f"  {category.value}: {count} tool(s)", file=sys.stderr)
    print("", file=sys.stderr)


# Initialize server and registry
mcp = initialize_server()
register_tools_in_registry()
register_tools_with_mcp(mcp, get_registry())


if __name__ == "__main__":
    # Print registered tools info
    print_registered_tools(get_registry())
    
    # The 'run()' method automatically handles CLI arguments like --transport
    # To run as HTTP by default in code, we can pass transport="streamable-http"
    # or just let the user use 'python main.py --transport streamable-http'
    if len(sys.argv) == 1:
        mcp.run(transport="streamable-http")
    else:
        mcp.run()
