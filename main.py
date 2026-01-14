from mcp.server.fastmcp import FastMCP
import tools
import os

# Initialize server
# We set host/port in constructor to provide defaults for HTTP transport
mcp = FastMCP(
    "mof-tools", 
    host="0.0.0.0", 
    port=8080,
    log_level="ERROR" # Suppress info logs that might corrupt the protocol stream
)

# Register tools from tools module
mcp.tool()(tools.search_mofs)
mcp.tool()(tools.calculate_energy)
mcp.tool()(tools.optimize_structure)

if __name__ == "__main__":
    # The 'run()' method automatically handles CLI arguments like --transport
    # To run as HTTP by default in code, we can pass transport="streamable-http"
    # or just let the user use 'python main.py --transport streamable-http'
    import sys
    
    # If no transport is specified in CLI, default to streamable-http for your remote needs
    if len(sys.argv) == 1:
        mcp.run(transport="streamable-http")
    else:
        mcp.run()
