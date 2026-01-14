from mcp.server.fastmcp import FastMCP
import tools

# Initialize server with host and port configuration
mcp = FastMCP("mof-tools", host="0.0.0.0", port=8080)

# Register tools from tools module
mcp.tool()(tools.search_mofs)
mcp.tool()(tools.calculate_energy)
mcp.tool()(tools.optimize_structure)

if __name__ == "__main__":
    # Run using SSE transport
    mcp.run(transport="sse")
