"""MCP Bridge: proxies Datagrok tool calls from Agent SDK to browser via WebSocket.

Tool definitions map to browser-side context class methods:
- TableViewContext (tableview-tools.ts)
- SQLGenerationContext (sql-tools.ts)
- ScriptGenerationContext (script-tools.ts)
"""
import asyncio
import json
import logging
from uuid import uuid4

from claude_agent_sdk import tool, create_sdk_mcp_server

logger = logging.getLogger(__name__)


# (name, description, {param: python_type})  — one entry per Datagrok tool
TOOL_DEFS = [
    # Table View tools (-> TableViewContext)
    ('add_viewer',
     'Add a new viewer to the table view with specified properties. '
     'Use describe_viewer first to understand available properties.',
     {'viewerType': str, 'viewerProperties': dict}),
    ('describe_viewer',
     'Get detailed info about a viewer type including all its properties.',
     {'viewerType': str}),
    ('list_current_viewers',
     'List all currently added viewers with IDs, types, and properties.',
     {}),
    ('adjust_viewer',
     'Adjust properties of an existing viewer by its ID.',
     {'viewerId': str, 'viewerProperties': dict}),
    ('find_viewers_by_type',
     'Find all viewers of a specific type and get their IDs and properties.',
     {'viewerType': str}),
    ('get_available_viewers',
     'Get list of all available viewer types.',
     {}),
    ('highlight_element',
     'Highlight a UI element to guide the user.',
     {'elementId': str, 'description': str}),
    ('list_ui_elements',
     'List all available UI elements with IDs and descriptions.',
     {}),

    # SQL tools (-> SQLGenerationContext)
    ('list_tables_in_schema',
     'List all table names in a schema with descriptions.',
     {'schemaName': str, 'catalogName': str}),
    ('describe_tables',
     'Get detailed column info for tables. Use schema.table format.',
     {'tables': list}),
    ('list_joins',
     'List foreign key relationships for specified tables.',
     {'tables': list}),
    ('try_sql',
     'Execute SQL query to test it. Returns row count and columns (limited to 10 rows).',
     {'sql': str, 'description': str}),
    ('list_all_schemas',
     'List all schemas in a catalog.',
     {'catalogName': str}),
    ('list_catalogs',
     'List all database catalogs available on the connection.',
     {}),
    ('find_similar_queries',
     'Find similar previously executed queries based on a prompt.',
     {'prompt': str}),

    # Script tools (-> ScriptGenerationContext)
    ('find_similar_script_samples',
     'Find similar script code samples from Datagrok codebase. '
     'Optionally filter by language: javascript, python, r, julia, octave.',
     {'description': str, 'language': str}),
    ('search_documentation',
     'Search Datagrok documentation and API references.',
     {'query': str, 'maxResults': int}),
    ('list_js_members',
     'List members of a JS namespace to verify API existence.',
     {'expression': str}),
    ('search_js_api_sources',
     'Search Datagrok JS API source code to verify a method exists.',
     {'query': str, 'maxResults': int}),
]

# Tool name sets for context-based filtering
TABLE_VIEW_TOOLS = {
    'add_viewer', 'describe_viewer', 'list_current_viewers', 'adjust_viewer',
    'find_viewers_by_type', 'get_available_viewers', 'highlight_element', 'list_ui_elements',
}
SQL_TOOLS = {
    'list_tables_in_schema', 'describe_tables', 'list_joins', 'try_sql',
    'list_all_schemas', 'list_catalogs', 'find_similar_queries',
}
SCRIPT_TOOLS = {
    'find_similar_script_samples', 'search_documentation', 'list_js_members', 'search_js_api_sources',
}


def mcp_tool_names(context: dict) -> list[str]:
    """Return mcp__datagrok__* tool names enabled for the given view context."""
    names = []
    if context.get('viewType') == 'TABLE_VIEW':
        names.extend(f'mcp__datagrok__{t}' for t in TABLE_VIEW_TOOLS)
    if context.get('connectionId'):
        names.extend(f'mcp__datagrok__{t}' for t in SQL_TOOLS)
    names.extend(f'mcp__datagrok__{t}' for t in SCRIPT_TOOLS)
    return names


class DatagrokMcpBridge:
    """Proxies MCP tool calls from the Agent SDK to the browser over WebSocket.

    The Agent SDK invokes Datagrok tools through the MCP server returned by
    ``mcp_server``.  Each call is forwarded to the browser as a ``tool_execute``
    WebSocket message; the browser executes it and replies with ``tool_result``.
    """

    def __init__(self, ws):
        self.ws = ws
        self.pending: dict[str, asyncio.Future] = {}
        self.mcp_server = self._build_mcp_server()

    def _build_mcp_server(self):
        """Create MCP server with all Datagrok tool definitions."""
        tools = []
        for name, description, params in TOOL_DEFS:
            tool_name = name
            input_schema = _build_schema(params)

            @tool(tool_name, description, input_schema)
            async def handler(args, _n=tool_name):
                return await self._proxy(_n, args)

            tools.append(handler)
        return create_sdk_mcp_server(name='datagrok', tools=tools)

    async def _proxy(self, tool_name: str, args: dict) -> dict:
        """Send tool call to browser over WebSocket, await result."""
        call_id = str(uuid4())
        loop = asyncio.get_running_loop()
        future = loop.create_future()
        self.pending[call_id] = future

        await self.ws.send(json.dumps({
            'type': 'tool_execute',
            'callId': call_id,
            'tool': tool_name,
            'args': args,
        }))

        try:
            result = await asyncio.wait_for(future, timeout=30.0)
        except asyncio.TimeoutError:
            result = f'Error: browser did not respond to {tool_name} within 30s'
        finally:
            self.pending.pop(call_id, None)

        return {'content': [{'type': 'text', 'text': str(result)}]}

    def resolve(self, call_id: str, result: str):
        """Called when browser sends back a tool_result message."""
        future = self.pending.get(call_id)
        if future and not future.done():
            future.set_result(result)

    def cancel_all(self):
        """Cancel all pending tool calls (called on WebSocket disconnect)."""
        for future in self.pending.values():
            if not future.done():
                future.set_result('Error: WebSocket disconnected')
        self.pending.clear()


_PY_TO_JSON = {str: 'string', int: 'integer', float: 'number', bool: 'boolean', dict: 'object', list: 'array'}


def _build_schema(params: dict) -> dict:
    """Convert {name: python_type} to JSON Schema dict for the @tool decorator."""
    if not params:
        return {'type': 'object', 'properties': {}}
    props = {k: {'type': _PY_TO_JSON.get(v, 'string')} for k, v in params.items()}
    return {'type': 'object', 'properties': props, 'required': list(params.keys())}
