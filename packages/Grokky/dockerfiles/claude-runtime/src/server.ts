import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {createNodeWebSocket} from '@hono/node-ws';
import {query} from '@anthropic-ai/claude-agent-sdk';
import type {SDKMessage} from '@anthropic-ai/claude-agent-sdk';
import type {UserMessage, AbortMessage, InputResponseMessage, OutgoingMessage, ToolInputs, McpInputs, ToolName, McpName} from './types';

const WORKSPACE = process.env['CLAUDE_WORKSPACE'] || '/workspace';
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
- Keep responses concise.

## Clarifying ambiguous requests

You have the AskUserQuestion tool available. When the user's request could reasonably be interpreted in multiple ways, you MUST use AskUserQuestion to clarify before acting. Examples of ambiguous requests:
- "add a plot" → ask which plot type (scatter, histogram, bar chart, etc.)
- "filter the data" → ask which column and criteria
- "clean up the data" → ask what kind of cleanup
- "correlate columns" → ask which columns

NEVER guess when there are multiple valid options. Always ask first.

## Entity References

When your response mentions Datagrok entities (files, scripts, queries, connections, projects, spaces),
render them as interactive cards by emitting a \`\`\`datagrok-entities fenced block with a JSON array.
Each entry MUST have \`type\` and \`name\`. Server entities MUST include \`id\`. Files use \`connector\` + \`path\`.

Supported types and required/optional fields:

| type         | required                     | optional                        |
|--------------|------------------------------|---------------------------------|
| file         | connector, path, name        | isDirectory, size               |
| script       | id, name                     | language                        |
| query        | id, name                     | connectionName                  |
| connection   | id, name                     | dataSource                      |
| project      | id, name                     |                                 |
| space        | id, name                     |                                 |

Example — after listing files or creating a script, emit:

\`\`\`datagrok-entities
[{"type":"file","connector":"System:DemoFiles","path":"datasets/demog.csv","name":"demog.csv","size":12345}]
\`\`\`

Rules:
- ALWAYS use this block when referencing entities the user may want to open or inspect.
- You get \`id\` values from MCP tool responses — pass them through exactly.
- For files, \`connector\` is the connection name (e.g. "System:DemoFiles") and \`path\` is relative.
- Add surrounding text before/after the block as needed for context.
- Do NOT put entity info as plain text when a datagrok-entities block is appropriate.`;

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

function rewriteForDocker(url: string): string {
  return url.replace('//localhost', '//host.docker.internal');
}

function buildMcpHeaders(apiKey?: string, apiUrl?: string): Record<string, string> {
  const headers: Record<string, string> = {};
  if (apiKey) {
    headers['Authorization'] = apiKey;
    headers['x-user-api-key'] = apiKey;
  }
  if (apiUrl)
    headers['x-datagrok-api-url'] = apiUrl;
  return headers;
}

function apiUrlFromMcpUrl(mcpUrl: string): string | undefined {
  const idx = mcpUrl.indexOf('/docker/containers/proxy/');
  return idx > 0 ? mcpUrl.substring(0, idx) : undefined;
}

function buildOptions(resume?: string, apiKey?: string, mcpServerUrl?: string) {
  const mcpUrl = mcpServerUrl || '';
  const apiUrl = mcpUrl ? apiUrlFromMcpUrl(mcpUrl) : undefined;
  return {
    systemPrompt: DATAGROK_PROMPT,
    allowedTools: ['Read', 'Glob', 'Grep', 'Edit', 'Write', 'Bash', 'WebSearch', 'WebFetch', 'AskUserQuestion'],
    ...(mcpUrl ? {
      mcpServers: {
        datagrok: {
          type: 'http' as const,
          url: mcpUrl,
          headers: buildMcpHeaders(apiKey, apiUrl),
        },
      },
    } : {}),
    permissionMode: 'acceptEdits' as const,
    model: 'opus' as const,
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
  AskUserQuestion: (i) => `Asking: ${i.questions?.[0]?.question ?? 'a question'}`,
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

interface ActiveQuery {
  abortController: AbortController;
  queryHandle: ReturnType<typeof query> | null;
  pendingInputResolve: ((value: any) => void) | null;
}

const activeQueries = new Map<WsSender, ActiveQuery>();

function handleInputResponse(ws: WsSender, data: InputResponseMessage): void {
  const active = activeQueries.get(ws);
  if (active?.pendingInputResolve) {
    const resolve = active.pendingInputResolve;
    active.pendingInputResolve = null;
    resolve(data.value);
  }
}

async function handleMessage(ws: WsSender, data: UserMessage): Promise<void> {
  const sid = data.sessionId ?? '';
  const message = data.message ?? '';
  if (!message)
    return emit(ws, {type: 'error', sessionId: sid, message: 'Empty message'});

  const mcpUrl = rewriteForDocker(data.mcpServerUrl || '');
  const abortController = new AbortController();
  const active: ActiveQuery = {abortController, queryHandle: null, pendingInputResolve: null};
  activeQueries.set(ws, active);

  let gotResult = false;
  try {
    const opts = buildOptions(getSession(sid), data.apiKey, mcpUrl);
    const canUseTool = async (toolName: string, input: any) => {
      if (toolName === 'AskUserQuestion') {
        emit(ws, {type: 'input_request', sessionId: sid, toolName, input});
        const updatedInput = await new Promise<any>((resolve, reject) => {
          active.pendingInputResolve = resolve;
          abortController.signal.addEventListener('abort', () => {
            active.pendingInputResolve = null;
            reject(new Error('aborted'));
          }, {once: true});
        });
        return {behavior: 'allow' as const, updatedInput};
      }
      return {behavior: 'allow' as const, updatedInput: input};
    };
    const q = query({prompt: promptStream(message), options: {...opts, canUseTool, abortController}});
    active.queryHandle = q;
    for await (const event of q) {
      if (abortController.signal.aborted)
        break;
      if (event.type === 'result')
        gotResult = true;
      forwardEvent(ws, sid, event);
    }
  } catch (e: any) {
    if (!abortController.signal.aborted && (!gotResult || !/exited with code/i.test(String(e.message))))
      emit(ws, {type: 'error', sessionId: sid, message: String(e.message || e)});
  } finally {
    activeQueries.delete(ws);
  }
}

function handleAbort(ws: WsSender, data: AbortMessage): void {
  const active = activeQueries.get(ws);
  if (!active)
    return;
  active.abortController.abort();
  try {
    if (active.queryHandle)
      active.queryHandle.close();
  } catch { /* query may have already finished */ }
  emit(ws, {type: 'aborted', sessionId: data.sessionId});
}

const app = new Hono();
const {injectWebSocket, upgradeWebSocket} = createNodeWebSocket({app});

app.get('/health', (c) => c.json({status: 'ok'}));

app.get('/ws', upgradeWebSocket(() => {
  let busy = false;
  return {
    onMessage(evt: any, ws: any) {
      const sender = ws as unknown as WsSender;
      let data: any;
      try {
        data = JSON.parse(String(evt.data));
      } catch {
        return emit(sender, {type: 'error', sessionId: '', message: 'Invalid JSON'});
      }

      if (data.type === 'abort') {
        handleAbort(sender, data);
        return;
      }

      if (data.type === 'input_response') {
        handleInputResponse(sender, data);
        return;
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
