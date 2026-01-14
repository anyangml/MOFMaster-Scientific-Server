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

Start the server using Python:

```bash
python main.py
```

The server will be live at `http://localhost:8080`.

## 🧪 Testing the Server

### Method 1: MCP Inspector (Recommended)
The official way to test MCP tools interactively. 

1. Start your server in one terminal:
   ```bash
   python main.py
   ```

2. In a second terminal, run the inspector pointing to the SSE URL:
   ```bash
   npx @modelcontextprotocol/inspector http://localhost:8080/sse
   ```


### Method 2: Manual HTTP Check
Since the server runs on HTTP/SSE, you can verify it's up with a simple `curl`:

```bash
curl http://localhost:8080/health
```

## 🛠️ Available Tools

| Tool | Description | Input |
| :--- | :--- | :--- |
| `search_mofs` | Search internal MOF database | `query` (name or formula) |
| `calculate_energy` | Calculate potential energy via ASE | `data` (CIF content or path) |
| `optimize_structure` | Structure optimization placeholder | `name` (MOF name) |

## 🔌 Integration
To connect your agent to this server, use the SSE endpoint:
`URL: http://localhost:8080/sse`
