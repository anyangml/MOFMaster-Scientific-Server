"""
Tool Registry - Formal registration system for MCP tools.

This module provides a production-ready tool registration system with metadata,
validation, and automatic discovery.
"""

from typing import Callable, Dict, List, Optional
from dataclasses import dataclass, field
from enum import Enum


class ToolCategory(str, Enum):
    """Categories for organizing tools."""
    SEARCH = "search"
    CALCULATION = "calculation"
    OPTIMIZATION = "optimization"
    ANALYSIS = "analysis"


@dataclass
class ToolMetadata:
    """Metadata for a registered tool."""
    name: str
    description: str
    category: ToolCategory
    function: Callable
    requires_ase: bool = False
    is_experimental: bool = False
    tags: List[str] = field(default_factory=list)
    version: str = "1.0.0"
    
    def __post_init__(self):
        """Validate metadata after initialization."""
        if not self.name:
            raise ValueError("Tool name cannot be empty")
        if not self.description:
            raise ValueError("Tool description cannot be empty")
        if not callable(self.function):
            raise ValueError("Tool function must be callable")


class ToolRegistry:
    """
    Central registry for managing MCP tools.
    
    This class provides a formal system for registering, discovering, and
    retrieving tools with associated metadata.
    """
    
    def __init__(self):
        """Initialize the tool registry."""
        self._tools: Dict[str, ToolMetadata] = {}
        self._categories: Dict[ToolCategory, List[str]] = {
            category: [] for category in ToolCategory
        }
    
    def register(
        self,
        name: str,
        description: str,
        category: ToolCategory,
        function: Callable,
        requires_ase: bool = False,
        is_experimental: bool = False,
        tags: Optional[List[str]] = None,
        version: str = "1.0.0"
    ) -> ToolMetadata:
        """
        Register a new tool with metadata.
        
        Args:
            name: Unique name for the tool
            description: Human-readable description
            category: Tool category
            function: The callable function implementing the tool
            requires_ase: Whether the tool requires ASE library
            is_experimental: Whether the tool is experimental
            tags: Optional list of tags for categorization
            version: Tool version
            
        Returns:
            ToolMetadata: The registered tool metadata
            
        Raises:
            ValueError: If tool name already exists or validation fails
        """
        if name in self._tools:
            raise ValueError(f"Tool '{name}' is already registered")
        
        metadata = ToolMetadata(
            name=name,
            description=description,
            category=category,
            function=function,
            requires_ase=requires_ase,
            is_experimental=is_experimental,
            tags=tags or [],
            version=version
        )
        
        self._tools[name] = metadata
        self._categories[category].append(name)
        
        return metadata
    
    def get(self, name: str) -> Optional[ToolMetadata]:
        """
        Get tool metadata by name.
        
        Args:
            name: Tool name
            
        Returns:
            ToolMetadata if found, None otherwise
        """
        return self._tools.get(name)
    
    def get_all(self) -> List[ToolMetadata]:
        """
        Get all registered tools.
        
        Returns:
            List of all tool metadata
        """
        return list(self._tools.values())
    
    def get_by_category(self, category: ToolCategory) -> List[ToolMetadata]:
        """
        Get all tools in a specific category.
        
        Args:
            category: Tool category
            
        Returns:
            List of tool metadata in the category
        """
        tool_names = self._categories.get(category, [])
        return [self._tools[name] for name in tool_names]
    
    def get_by_tag(self, tag: str) -> List[ToolMetadata]:
        """
        Get all tools with a specific tag.
        
        Args:
            tag: Tag to search for
            
        Returns:
            List of tool metadata with the tag
        """
        return [
            metadata for metadata in self._tools.values()
            if tag in metadata.tags
        ]
    
    def list_names(self) -> List[str]:
        """
        Get list of all registered tool names.
        
        Returns:
            List of tool names
        """
        return list(self._tools.keys())
    
    def list_categories(self) -> Dict[ToolCategory, int]:
        """
        Get count of tools per category.
        
        Returns:
            Dictionary mapping categories to tool counts
        """
        return {
            category: len(tools)
            for category, tools in self._categories.items()
        }
    
    def unregister(self, name: str) -> bool:
        """
        Unregister a tool.
        
        Args:
            name: Tool name to unregister
            
        Returns:
            True if tool was unregistered, False if not found
        """
        if name not in self._tools:
            return False
        
        metadata = self._tools[name]
        self._categories[metadata.category].remove(name)
        del self._tools[name]
        
        return True
    
    def clear(self):
        """Clear all registered tools."""
        self._tools.clear()
        self._categories = {category: [] for category in ToolCategory}
    
    def __len__(self) -> int:
        """Return the number of registered tools."""
        return len(self._tools)
    
    def __contains__(self, name: str) -> bool:
        """Check if a tool is registered."""
        return name in self._tools
    
    def __repr__(self) -> str:
        """String representation of the registry."""
        return f"ToolRegistry(tools={len(self._tools)})"


# Global registry instance
_global_registry = ToolRegistry()


def get_registry() -> ToolRegistry:
    """
    Get the global tool registry instance.
    
    Returns:
        The global ToolRegistry instance
    """
    return _global_registry
