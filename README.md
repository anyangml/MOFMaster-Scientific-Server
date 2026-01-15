# MOF Tools MCP Server

A **production-ready**, professional MCP server for Metal-Organic Framework (MOF) research. This server uses the **SSE (Server-Sent Events)** transport to provide scientific tools over HTTP with comprehensive **Pydantic validation** and a **formal tool registry system**.

## ✨ Key Features

- **🔒 Input/Output Validation**: All tool inputs and outputs are validated using Pydantic models
- **📋 Formal Tool Registry**: Centralized tool management with metadata, categories, and tags
- **🛡️ Production-Ready**: Industry-standard code with comprehensive error handling
- **📊 Type Safety**: Full type hints throughout the codebase
- **✅ Tested**: Comprehensive test suite with 33 tests

## 📁 Repository Structure

- `main.py`: Server entrypoint with tool registration and initialization
- `tools.py`: Core scientific tools with Pydantic validation models
- `tool_registry.py`: Formal tool registration system with metadata management
- `tests/test_tools.py`: Comprehensive test suite for tools and validation
- `pyproject.toml`: Dependency and package management

## 🚀 Installation

```bash
# Install dependencies
pip install mcp[server] ase pydantic
```

## 🏃 Running the Server

To run as the modern HTTP server (Streamable HTTP):
```bash
python main.py
```

The server will start on `http://0.0.0.0:8080` and display registered tools information:

```
=== MOF Tools Server ===
Registered 3 tools:
  - search_mofs (search)  
  - calculate_energy (calculation) [ASE Required]
  - optimize_structure (optimization) [EXPERIMENTAL]

Tools by category:
  search: 1 tool(s)
  calculation: 1 tool(s)
  optimization: 1 tool(s)
  analysis: 0 tool(s)
```

## 🧪 Testing

Run the comprehensive test suite:

```bash
# Run all tests
pytest tests/test_tools.py -v

# Run specific test class
pytest tests/test_tools.py::TestPydanticModels -v

# Run with coverage
pytest tests/test_tools.py --cov=tools --cov=tool_registry
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

All tools return validated JSON responses with comprehensive error handling.

### search_mofs

**Category**: Search  
**Description**: Search for Metal-Organic Frameworks by name or formula in the database

**Input**:
```json
{
  "query": "MOF-5"  // String, 1-200 characters
}
```

**Output**:
```json
{
  "success": true,
  "results": [
    {
      "name": "MOF-5",
      "formula": "Zn4O(BDC)3",
      "surface_area": 3800.0
    }
  ],
  "count": 1,
  "message": "Found 1 MOF(s)"
}
```

### calculate_energy

**Category**: Calculation  
**Requirements**: ASE library  
**Description**: Calculate the potential energy of a structure from CIF file content or path using ASE

**Input**:
```json
{
  "data": "<CIF file content or path>"  // String, non-empty
}
```

**Output**:
```json
{
  "success": true,
  "energy": -123.4567,
  "error": null,
  "message": "Energy: -123.4567 eV"
}
```

### optimize_structure

**Category**: Optimization  
**Status**: Experimental  
**Description**: Perform structure optimization for a named MOF structure (placeholder implementation)

**Input**:
```json
{
  "name": "HKUST-1"  // String, 1-200 characters
}
```

**Output**:
```json
{
  "success": true,
  "structure_name": "HKUST-1",
  "message": "Successfully initiated optimization for HKUST-1"
}
```

## 🔌 Integration

To connect your agent to this server, use the Streamable HTTP endpoint:

**URL**: `http://localhost:8080/mcp`

## 🏗️ Architecture

### Pydantic Validation

All tools use Pydantic models for comprehensive input/output validation:

- **Input Models**: Validate and sanitize user inputs
- **Output Models**: Ensure consistent response structure
- **Field Validators**: Custom validation logic (e.g., trimming whitespace, checking ranges)
- **Error Handling**: Graceful error messages returned as validated JSON

### Tool Registry System

The formal tool registry provides:

- **Metadata Management**: Description, category, version, tags for each tool
- **Categorization**: Tools organized by category (search, calculation, optimization, analysis)
- **Tag-based Discovery**: Find tools by tags (e.g., "mof", "energy", "database")
- **Dependency Tracking**: Mark tools that require specific libraries (e.g., ASE)
- **Status Indicators**: Mark experimental or production-ready tools

### Production-Ready Features

- ✅ Full type hints throughout
- ✅ Comprehensive error handling
- ✅ Input sanitization and validation
- ✅ Consistent JSON output format
- ✅ Detailed documentation
- ✅ Test coverage (33 tests)
- ✅ Modular, maintainable code structure

## 📚 Development

### Adding a New Tool

1. **Define Pydantic models** in `tools.py`:
```python
class MyToolInput(BaseModel):
    param: str = Field(..., description="Parameter description")

class MyToolOutput(BaseModel):
    success: bool
    result: str
```

2. **Implement the tool function** in `tools.py`:
```python
def my_tool(param: str) -> str:
    try:
        validated_input = MyToolInput(param=param)
        # ... tool logic ...
        output = MyToolOutput(success=True, result="...")
        return output.model_dump_json(indent=2)
    except ValidationError as e:
        # ... error handling ...
```

3. **Register the tool** in `main.py`:
```python
registry.register(
    name="my_tool",
    description="Tool description",
    category=ToolCategory.CALCULATION,
    function=tools.my_tool,
    tags=["tag1", "tag2"]
)
```

4. **Add tests** in `test_tools.py`:
```python
def test_my_tool():
    result = tools.my_tool("test")
    parsed = json.loads(result)
    assert parsed["success"] is True
```

## 📝 License

See LICENSE file for details.
