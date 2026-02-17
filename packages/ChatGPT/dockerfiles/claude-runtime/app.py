"""Quart WebSocket server bridging Datagrok browser clients and Claude Agent SDK."""
import asyncio
import json
import logging
import os
import traceback
from collections import OrderedDict

from quart import Quart, jsonify, websocket, request
from hypercorn.asyncio import serve
from hypercorn.config import Config

from claude_agent_sdk import (
    query, ClaudeAgentOptions,
    SystemMessage, AssistantMessage, ResultMessage,
    TextBlock, ToolUseBlock,
)
from datagrok_mcp import DatagrokMcpBridge, mcp_tool_names

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Quart(__name__)

# Client sessionId -> SDK session_id (LRU, max 200 entries)
session_map: OrderedDict[str, str] = OrderedDict()
MAX_SESSIONS = 200

WORKSPACE = os.environ.get('CLAUDE_WORKSPACE', '/workspace')

SYSTEM_PROMPT = """\
You are a coding assistant for the Datagrok platform — a data analytics and visualization tool.
You help users explore data, create visualizations, write SQL queries, and build scripts.

## Workspace

The Datagrok public repository is mounted at /workspace (the CLAUDE_WORKSPACE directory).
Key paths inside /workspace:
- js-api/src/           — Core TypeScript API (datagrok-api): dataframe.ts, viewer.ts, grid.ts, etc.
- libraries/            — Shared TypeScript libraries (@datagrok-libraries/*)
- packages/             — 70+ extension packages (viewers, connectors, scientific tools)
- packages/ChatGPT/     — This package (AI integration)
- packages/ApiSamples/  — API usage examples
- packages/ApiTests/    — API test coverage
- help/                 — Documentation
- tools/                — CLI tool "grok" (datagrok-tools)

Use the built-in Read, Glob, and Grep tools to explore the codebase when answering
questions about Datagrok internals, API usage, or package implementation patterns.

## Tool Usage

When working with table views:
- Use describe_viewer first to understand available properties before adding viewers
- Use list_current_viewers to see what's already on screen
- Use adjust_viewer to modify existing viewers instead of adding duplicates

When working with SQL:
- Use list_tables_in_schema and describe_tables to understand the database structure
- Use try_sql to test queries before presenting them (limited to 10 rows)
- Always validate column names against the schema

When writing scripts:
- Use find_similar_script_samples to find relevant examples
- Use list_js_members to verify API methods exist before using them
- Use search_documentation for platform-specific guidance

Be concise. Show results, not process.\
"""

SDK_BUILTIN_TOOLS = ['Read', 'Glob', 'Grep', 'WebSearch', 'WebFetch']


async def _prompt_stream(message: str):
    """Wrap a string prompt as an async iterable.

    The SDK has a bug where string prompts call end_input() immediately,
    closing stdin before SDK MCP servers finish their handshake.  Passing
    the prompt as an async iterable forces the SDK to use stream_input(),
    which properly waits for the first result before closing.
    """
    yield {
        'type': 'user',
        'session_id': '',
        'message': {'role': 'user', 'content': message},
        'parent_tool_use_id': None,
    }


def store_session(client_id: str, sdk_id: str):
    if client_id in session_map:
        session_map.move_to_end(client_id)
    session_map[client_id] = sdk_id
    while len(session_map) > MAX_SESSIONS:
        session_map.popitem(last=False)


def get_sdk_session(client_id: str) -> str | None:
    sid = session_map.get(client_id)
    if sid is not None:
        session_map.move_to_end(client_id)
    return sid


# -- HTTP endpoints -----------------------------------------------------------

@app.route('/health', methods=['GET'])
async def health():
    return jsonify({'status': 'ok'}), 200


@app.errorhandler(404)
async def not_found(e):
    logger.warning(f'404 Not Found: {request.path}')
    return jsonify({'error': 'Not found', 'path': request.path}), 404


@app.errorhandler(500)
async def internal_error(e):
    logger.error(f'500 Internal Server Error: {e}')
    return jsonify({'error': 'Internal server error', 'message': str(e)}), 500


# -- WebSocket handler --------------------------------------------------------

@app.websocket('/ws')
async def ws_handler():
    ws = websocket._get_current_object()
    bridge = DatagrokMcpBridge(ws)
    incoming: asyncio.Queue = asyncio.Queue()

    async def receive_loop():
        """Read messages from client; route tool_results to bridge, queue the rest."""
        try:
            while True:
                raw = await ws.receive()
                try:
                    data = json.loads(raw)
                except json.JSONDecodeError:
                    await ws.send(json.dumps({
                        'type': 'error', 'sessionId': '', 'message': 'Invalid JSON',
                    }))
                    continue

                if data.get('type') == 'tool_result':
                    bridge.resolve(data.get('callId', ''), data.get('result', ''))
                else:
                    await incoming.put(data)
        except asyncio.CancelledError:
            pass
        except Exception:
            logger.debug('receive_loop ended')

    recv_task = asyncio.create_task(receive_loop())
    try:
        while True:
            data = await incoming.get()
            if data.get('type') == 'user_message':
                await _handle_user_message(ws, bridge, data)
            else:
                await ws.send(json.dumps({
                    'type': 'error',
                    'sessionId': data.get('sessionId', ''),
                    'message': f'Unknown message type: {data.get("type")}',
                }))
    except asyncio.CancelledError:
        pass
    except Exception:
        logger.error(f'ws_handler error:\n{traceback.format_exc()}')
    finally:
        recv_task.cancel()
        bridge.cancel_all()


async def _handle_user_message(ws, bridge: DatagrokMcpBridge, data: dict):
    session_id = data.get('sessionId', '')
    message = data.get('message', '')
    context = data.get('context', {})

    if not message:
        await ws.send(json.dumps({
            'type': 'error', 'sessionId': session_id, 'message': 'Empty message',
        }))
        return

    allowed_tools = SDK_BUILTIN_TOOLS + mcp_tool_names(context)
    sdk_session = get_sdk_session(session_id)

    options = ClaudeAgentOptions(
        system_prompt=SYSTEM_PROMPT,
        allowed_tools=allowed_tools,
        mcp_servers={'datagrok': bridge.mcp_server},
        permission_mode='bypassPermissions',
        max_turns=15,
        model='sonnet',
        include_partial_messages=True,
        cwd=WORKSPACE,
    )
    if sdk_session:
        options.resume = sdk_session

    try:
        async for event in query(prompt=_prompt_stream(message), options=options):
            await _forward_event(ws, session_id, event)
    except Exception as e:
        logger.error(f'Query error: {traceback.format_exc()}')
        await ws.send(json.dumps({
            'type': 'error', 'sessionId': session_id, 'message': str(e),
        }))


async def _forward_event(ws, session_id: str, event):
    """Forward an Agent SDK event to the browser as a WebSocket message.

    Event types from the SDK:
    - SystemMessage  (subtype='init')  -> store session mapping
    - AssistantMessage                 -> text chunks + tool use notifications
    - ResultMessage  (subtype='success'/'error_*') -> final/error
    """
    if isinstance(event, SystemMessage):
        if event.subtype == 'init':
            sdk_sid = event.data.get('session_id')
            if sdk_sid:
                store_session(session_id, sdk_sid)

    elif isinstance(event, AssistantMessage):
        for block in event.content:
            if isinstance(block, TextBlock):
                await ws.send(json.dumps({
                    'type': 'chunk', 'sessionId': session_id, 'content': block.text,
                }))
            elif isinstance(block, ToolUseBlock):
                await ws.send(json.dumps({
                    'type': 'tool_use',
                    'sessionId': session_id,
                    'tool': block.name,
                    'input': block.input,
                    'status': 'running',
                }))

    elif isinstance(event, ResultMessage):
        if not event.is_error:
            await ws.send(json.dumps({
                'type': 'final',
                'sessionId': session_id,
                'content': event.result or '',
                'usage': event.usage,
                'cost_usd': event.total_cost_usd,
                'num_turns': event.num_turns,
            }))
        else:
            await ws.send(json.dumps({
                'type': 'error',
                'sessionId': session_id,
                'message': event.result or event.subtype,
            }))


# -- Entry point --------------------------------------------------------------

async def run_server():
    config = Config()
    config.bind = ['0.0.0.0:5353']
    await serve(app, config)


if __name__ == '__main__':
    logger.info('Starting claude-runtime on port 5353')
    asyncio.run(run_server())
