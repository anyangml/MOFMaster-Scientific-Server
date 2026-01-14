# MOF Tools MCP Server

A professional, lightweight MCP server for Metal-Organic Framework (MOF) research. This server uses the **SSE (Server-Sent Events)** transport to provide scientific tools over HTTP.

## 📁 Repository Structure
- `main.py`: Server entrypoint and tool registration.
- `tools.py`: Core logic for scientific calculations and database queries.
- `pyproject.toml`: Dependency and package management.

## 🚀 Installation

```bash
# Install dependencies
pip install mcp[server] ase
```

## 🏃 Running the Server

To run as the modern HTTP server (Streamable HTTP):
```bash
python main.py
```

To run in the standard local mode (stdio):
```bash
python main.py --transport stdio
```

## 🧪 Testing the Server

### Method 1: MCP Inspector
```bash
# Connect inspector to your running HTTP server
npx @modelcontextprotocol/inspector http://localhost:8080/mcp
```


### Method 2: Manual HTTP Check
Since the server runs on HTTP, you can verify it's up with a simple `curl`:

```bash
curl http://localhost:8080/mcp
```

## 🛠️ Available Tools

| Tool | Description | Input |
| :--- | :--- | :--- |
| `search_mofs` | Search internal MOF database | `query` (name or formula) |
| `calculate_energy` | Calculate potential energy via ASE | `data` (CIF content or path) |
| `optimize_structure` | Structure optimization placeholder | `name` (MOF name) |

## 🔌 Integration
To connect your agent to this server, use the Streamable HTTP endpoint:
`URL: http://localhost:8080/mcp`
