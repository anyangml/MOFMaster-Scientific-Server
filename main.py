"""
MOF Tools MCP Server - Main entry point.

Production-ready MCP server with formal tool registration and Pydantic validation.
"""

import sys
from typing import List
from pathlib import Path
import yaml

from mcp.server.fastmcp import FastMCP

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
            requires_ase=tool_def['requires_ase'],
            is_experimental=tool_def['is_experimental'],
            tags=tool_def['tags'],
            version=tool_def['version']
        )


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
