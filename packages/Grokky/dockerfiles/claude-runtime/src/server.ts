import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {createNodeWebSocket} from '@hono/node-ws';
import {query} from '@anthropic-ai/claude-agent-sdk';
import type {SDKMessage} from '@anthropic-ai/claude-agent-sdk';
import type {UserMessage, OutgoingMessage, ToolInputs, McpInputs, ToolName, McpName} from './types';

const WORKSPACE = process.env['CLAUDE_WORKSPACE'] || '/workspace';
const MCP_SERVER_URL = process.env['MCP_SERVER_URL'] || '';
const PORT = 5355;
const MAX_SESSIONS = 200;

const DATAGROK_PROMPT = `\
You are Datagrok Code Assistant — an AI coding agent embedded in the Datagrok data analytics platform.
You have full access to the Datagrok public repository at ${WORKSPACE}, including the JS API source,
70+ extension packages, API samples, and documentation.

NEVER guess API methods, property names, or function signatures — ALWAYS search the codebase first.
Key lookup paths: js-api/src/ (core API), packages/ApiSamples/ (usage examples), help/ (docs).

## Code Execution

When the user asks you to **do** something (add a viewer, modify data, run a script, etc.):
1. FIRST search the codebase to find the correct API (method names, property names, option keys).
2. THEN emit the code in a \`\`\`datagrok-exec fenced block — this is the ONLY way code gets executed.

NEVER emit a \`\`\`datagrok-exec block with guessed or unverified API calls.
Regular \`\`\`javascript blocks are for explanations only and will NOT run.

### Globals in \`\`\`datagrok-exec blocks

| Variable | Type | When available |
|----------|------|----------------|
| \`grok\` | module | Always |
| \`ui\` | module | Always |
| \`DG\` | module | Always |
| \`view\` | \`DG.ViewBase\` | Always — check \`view.type\` for specific view type |
| \`t\` | \`DG.DataFrame\` | Only in TableView (\`view.type === 'TableView'\`) |

### Rules
- \`await\` works directly — the block runs in an async context.
- Keep blocks short and focused: one action per block.
- Do NOT import \`grok\`, \`ui\`, or \`DG\` — they are already in scope.
- Prefer the JS API over console commands or scripts.
- Keep responses concise.`;

const sessions = new Map<string, string>();

function storeSession(clientId: string, sdkId: string): void {
  sessions.delete(clientId);
  sessions.set(clientId, sdkId);
  if (sessions.size > MAX_SESSIONS)
    sessions.delete(sessions.keys().next().value!);
}

function getSession(clientId: string): string | undefined {
  const sid = sessions.get(clientId);
  if (sid !== undefined) {
    sessions.delete(clientId);
    sessions.set(clientId, sid);
  }
  return sid;
}

async function* promptStream(message: string) {
  yield {
    type: 'user' as const,
    session_id: '',
    message: {role: 'user' as const, content: message},
    parent_tool_use_id: null,
  };
}

function buildOptions(resume?: string) {
  return {
    systemPrompt: DATAGROK_PROMPT,
    allowedTools: ['Read', 'Glob', 'Grep', 'Edit', 'Write', 'Bash', 'WebSearch', 'WebFetch'],
    ...(MCP_SERVER_URL ? {
      mcpServers: {
        datagrok: {
          type: 'http' as const,
          url: MCP_SERVER_URL,
        },
      },
    } : {}),
    permissionMode: 'bypassPermissions' as const,
    allowDangerouslySkipPermissions: true,
    model: 'sonnet' as const,
    includePartialMessages: true,
    cwd: WORKSPACE,
    ...(resume ? {resume} : {}),
  };
}

const toolFormatters: {[K in ToolName]: (i: ToolInputs[K]) => string} = {
  Read: (i) => `Read ${i.file_path ?? ''}`,
  Write: (i) => `Write ${i.file_path ?? ''}`,
  Edit: (i) => `Edit ${i.file_path ?? ''}`,
  Bash: (i) => `Bash: ${(i.command ?? '').slice(0, 120)}`,
  Glob: (i) => `Glob ${i.pattern ?? ''} ${i.path ? 'in ' + i.path : ''}`.trim(),
  Grep: (i) => `Grep '${i.pattern ?? ''}' ${i.path ? 'in ' + i.path : ''}`.trim(),
  WebSearch: (i) => `WebSearch: ${i.query ?? ''}`,
  WebFetch: (i) => `WebFetch: ${i.url ?? ''}`,
};

const mcpFormatters: {[K in McpName]: (i: McpInputs[K]) => string} = {
  call_function: (i) => `Call ${i.name ?? 'function'}`,
  list_functions: (i) => `List functions${i.filter ? ': ' + i.filter : ''}`,
  get_function: (i) => `Get function ${i.id ?? ''}`,
  list_files: (i) => `list files ${i.path ?? ''}`.trim(),
  download_file: (i) => `download file ${i.path ?? ''}`.trim(),
  upload_file: (i) => `upload file ${i.path ?? ''}`.trim(),
};

function toolSummary(name: string, input: Record<string, unknown>): string {
  if (name in toolFormatters)
    return toolFormatters[name as ToolName](input as ToolInputs[ToolName]);
  const mcp = name.replace(/^mcp__[^_]+__/, '');
  if (mcp !== name) {
    if (mcp in mcpFormatters)
      return mcpFormatters[mcp as McpName](input as McpInputs[McpName]);
    return mcp.replace(/_/g, ' ');
  }
  return name;
}

function extractResult(evt: any): string | null {
  const r = evt.tool_use_result ?? evt.tool_result;
  if (!r)
    return null;
  if (typeof r === 'string')
    return r || null;
  if (Array.isArray(r))
    return r.filter((b: any) => b.type === 'text').map((b: any) => b.text).join('\n') || null;
  return null;
}

interface WsSender {
  send(data: string): void;
}

function emit(ws: WsSender, msg: OutgoingMessage): void {
  ws.send(JSON.stringify(msg));
}

function forwardEvent(ws: WsSender, sid: string, event: SDKMessage): void {
  const e = event as any;
  switch (event.type) {
  case 'system':
    if (event.subtype === 'init' && e.session_id)
      storeSession(sid, e.session_id);
    break;
  case 'assistant':
    for (const block of e.message?.content ?? []) {
      if (block.type === 'tool_use')
        emit(ws, {type: 'tool_activity', sessionId: sid, summary: toolSummary(block.name, block.input ?? {})});
    }
    break;
  case 'user': {
    const content = extractResult(e);
    if (content)
      emit(ws, {type: 'tool_result', sessionId: sid, content});
    break;
  }
  case 'stream_event':
    if (e.event?.delta?.type === 'text_delta' && e.event.delta.text)
      emit(ws, {type: 'chunk', sessionId: sid, content: e.event.delta.text});
    break;
  case 'tool_progress':
    emit(ws, {type: 'tool_activity', sessionId: sid, summary: `Running ${e.tool_name ?? ''}…`});
    break;
  case 'tool_use_summary':
    emit(ws, {type: 'tool_activity', sessionId: sid, summary: e.summary ?? ''});
    break;
  case 'result':
    if (e.subtype === 'success')
      emit(ws, {type: 'final', sessionId: sid, content: e.result || ''});
    else
      emit(ws, {type: 'error', sessionId: sid, message: (e.errors ?? []).join(', ') || e.subtype || 'unknown'});
    break;
  }
}

async function handleMessage(ws: WsSender, data: UserMessage): Promise<void> {
  const sid = data.sessionId ?? '';
  const message = data.message ?? '';
  if (!message)
    return emit(ws, {type: 'error', sessionId: sid, message: 'Empty message'});

  let gotResult = false;
  try {
    const opts = buildOptions(getSession(sid));
    for await (const event of query({prompt: promptStream(message), options: opts})) {
      if (event.type === 'result')
        gotResult = true;
      forwardEvent(ws, sid, event);
    }
  } catch (e: any) {
    if (!gotResult || !/exited with code/i.test(String(e.message)))
      emit(ws, {type: 'error', sessionId: sid, message: String(e.message || e)});
  }
}

const app = new Hono();
const {injectWebSocket, upgradeWebSocket} = createNodeWebSocket({app});

app.get('/health', (c) => c.json({status: 'ok'}));

app.get('/ws', upgradeWebSocket(() => {
  let busy = false;
  return {
    onMessage(evt, ws) {
      const sender = ws as unknown as WsSender;
      let data: any;
      try {
        data = JSON.parse(String(evt.data));
      } catch {
        return emit(sender, {type: 'error', sessionId: '', message: 'Invalid JSON'});
      }

      if (data.type !== 'user_message') {
        return emit(sender, {
          type: 'error', sessionId: data.sessionId ?? '',
          message: `Unknown type: ${data.type}`,
        });
      }

      if (busy) {
        return emit(sender, {
          type: 'error', sessionId: data.sessionId ?? '',
          message: 'Query already in progress',
        });
      }

      busy = true;
      handleMessage(sender, data).finally(() => { busy = false; });
    },
  };
}));

app.notFound((c) => c.json({error: 'Not found'}, 404));
app.onError((err, c) => c.json({error: String(err)}, 500));

if (!process.env['ANTHROPIC_API_KEY'])
  console.warn('ANTHROPIC_API_KEY is not set — Claude API calls will fail');

const server = serve({fetch: app.fetch, port: PORT});
injectWebSocket(server);
console.log(`claude-runtime listening on :${PORT}`);
