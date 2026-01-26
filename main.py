"""
MOF Tools MCP Server - Main entry point.

Production-ready MCP server with formal tool registration and Pydantic validation.
"""

import sys
from typing import List
from pathlib import Path
import yaml

from mcp.server.fastmcp import FastMCP
from dp.agent.server.calculation_mcp_server import CalculationMCPServer

from tool_registry import ToolRegistry, ToolCategory, get_registry
import tools


def load_tool_definitions(yaml_path: str = "tool_definitions.yaml") -> List[dict]:
    """
    Load tool definitions from a YAML configuration file.
    
    Args:
        yaml_path: Path to the YAML file containing tool definitions
        
    Returns:
        List of tool definition dictionaries
        
    Raises:
        FileNotFoundError: If the YAML file doesn't exist
        yaml.YAMLError: If the YAML file is malformed
    """
    config_path = Path(__file__).parent / yaml_path
    
    if not config_path.exists():
        raise FileNotFoundError(f"Tool definitions file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config.get('tools', [])


def initialize_server() -> CalculationMCPServer:
    """
    Initialize the CalculationMCPServer with production settings.
    
    Returns:
        Configured CalculationMCPServer instance
    """
    # CalculationMCPServer will initialize its own FastMCP instance internally
    calc_server = CalculationMCPServer(
        "mof-tools",
        host="0.0.0.0",
        port=50001,
        log_level="ERROR"
    )
    
    return calc_server


def register_tools_in_registry():
    """
    Register all available tools in the formal tool registry.
    
    Tool definitions are loaded from tool_definitions.yaml, which provides
    a centralized, organized way to manage tools with metadata.
    """
    registry = get_registry()
    
    # Load tool definitions from YAML configuration
    tool_definitions = load_tool_definitions()
    
    # Register all tools using a loop
    for tool_def in tool_definitions:
        # Get the function from the tools module using the function_name
        function_name = tool_def.pop('function_name')
        tool_function = getattr(tools, function_name)
        
        # Convert category string to ToolCategory enum
        category_str = tool_def.pop('category')
        category = ToolCategory[category_str]
        
        # Register the tool with the converted values
        registry.register(
            name=tool_def['name'],
            description=tool_def['description'],
            category=category,
            function=tool_function,
            tags=tool_def['tags'],
            version=tool_def['version']
        )


def register_tools_with_mcp(calc_server: CalculationMCPServer, registry: ToolRegistry):
    """
    Register tools from the registry with the CalculationMCPServer.
    
    This wrapper automatically adds "submit_", "query_job_status", and 
    "get_job_results" tools required by the Bohr Agent SDK.
    """
    for tool_metadata in registry.get_all():
        # Register each tool function with CalculationMCPServer
        calc_server.tool()(tool_metadata.function)
    return calc_server


def print_registered_tools(registry: ToolRegistry):
    """
    Print information about registered tools for debugging/info purposes.
    
    Args:
        registry: Tool registry to display
    """
    print(f"\n=== MOF Tools Server ===", file=sys.stderr)
    print(f"Registered {len(registry)} tools:", file=sys.stderr)
    
    for tool_metadata in registry.get_all():
        print(
            f"  - {tool_metadata.name} ({tool_metadata.category.value})",
            file=sys.stderr
        )
        print(f"    {tool_metadata.description}", file=sys.stderr)
    
    print(f"\nTools by category:", file=sys.stderr)
    for category, count in registry.list_categories().items():
        print(f"  {category.value}: {count} tool(s)", file=sys.stderr)
    print("", file=sys.stderr)


# Initialize server and registry
calc_server = initialize_server()
register_tools_in_registry()
register_tools_with_mcp(calc_server, get_registry())


def run_server():
    """Main entry point for running the server."""
    # Print registered tools info
    print_registered_tools(get_registry())
    
    # Use calc_server.run to include health check endpoint and correct tool patches
    # CalculationMCPServer.run handles host and port settings internally if passed correctly
    if len(sys.argv) == 1:
        calc_server.run(transport="streamable-http", host="0.0.0.0", port=50001)
    else:
        calc_server.run(host="0.0.0.0", port=50001)


if __name__ == "__main__":
    run_server()
