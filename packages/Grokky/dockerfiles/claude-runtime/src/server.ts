import * as fs from 'node:fs';
import * as path from 'node:path';
import {spawn} from 'node:child_process';
import type {ChildProcess} from 'node:child_process';
import {randomUUID} from 'node:crypto';
import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {createNodeWebSocket} from '@hono/node-ws';
import {query, createSdkMcpServer, tool as sdkTool} from '@anthropic-ai/claude-agent-sdk';
import type {SDKMessage, HookCallback} from '@anthropic-ai/claude-agent-sdk';
import {z} from 'zod/v4';
import type {UserMessage, AbortMessage, InputResponseMessage, OutgoingMessage, ToolInputs, McpInputs, ToolName, McpName} from './types';
import {ClaudeModel} from './types';
import {WORKSPACE} from './constants';
import {syncUserFiles} from './sync/orchestrator';
import {ensureUserDir} from './user/user-dir';
import {createPackageKnowledgeServer} from './package-knowledge-tool';
import {awaitWorkspaceSync, markQueryStart, markQueryEnd, startWorkspaceSync} from './sync/workspace';

const PORT = 5355;
const MAX_SESSIONS = 200;

type FenceMode = 'prose' | 'other';
interface FenceState { mode: FenceMode; carry: string; lineInProgress: boolean; }
const fenceStates = new Map<string, FenceState>();
const FENCE_RE = /^```([\w-]*)\s*$/;

const ENTITY_TYPE_MAP: Record<string, string> = {
  project: 'project', dataquery: 'query', dataconnection: 'connection',
  script: 'script', space: 'space', usergroup: 'group', group: 'group',
  user: 'user', fileinfo: 'file', file: 'file',
};

const entityRefSchema = z.object({
  '#type': z.string().optional(),
  type: z.string().optional(),
  name: z.string(),
  id: z.string().optional(),
  connector: z.string().optional(),
  path: z.string().optional(),
  isDirectory: z.boolean().optional(),
}).transform((e) => {
  const raw = (e['#type'] ?? e.type ?? '') as string;
  return {
    type: ENTITY_TYPE_MAP[raw.toLowerCase()] ?? raw.toLowerCase(),
    name: e.name,
    ...(e.id && {id: e.id}),
    ...(e.connector && {connector: e.connector}),
    ...(e.path && {path: e.path}),
    ...(e.isDirectory !== undefined && {isDirectory: e.isDirectory}),
  };
});

const BASH_EXEC_PROMPT = `\
Execute the given shell command using the Bash tool. \
Output ONLY the exact stdout of the command — no preamble, no explanation, nothing else.`;

const DATAGROK_PROMPT = `\
You are Datagrok Code Assistant — an AI coding agent embedded in the Datagrok data analytics platform.
You have access to the Datagrok public repository at \`workspace/\` (relative to your current working
directory), including the JS API source, API samples, and documentation. This view contains only
packages installed on this instance — uninstalled packages are not present.

ALWAYS use paths under \`workspace/\` (e.g. \`workspace/packages/Chem/...\`, \`workspace/js-api/src/...\`).
NEVER use \`/workspace/...\` as an absolute path — access is blocked at the hook level.

Personal user knowledge files (if any) live in the \`agents/\` directory in your current working
directory. Use Glob or \`ls agents/\` to discover them when relevant.

## Don't invent names

NEVER guess function names, parameter names, signatures, or JS API methods. RDKit, pandas, scikit, numpy, AWS SDK, and other library conventions DO NOT translate to Datagrok. If you cannot point to an exact name in the inlined skills below, in an MCP discovery result, or in \`workspace/js-api/src/\`, STOP and look it up before emitting code. Inventing names is the #1 cause of silent failures.

For a Datagrok function in a package not covered by an inlined skill, call \`list_functions(keyword)\` (MCP) to discover it; use \`get_function(name)\` only when you need full parameter details.

## Knowledge graph — query it before you grep

\`workspace/\` ships a queryable knowledge graph of the whole platform (packages, functions, scripts, queries, libraries, classes/methods, docs, changelog — and the edges between them). For ANY structural question ("what implements X?", "which package owns / tests / imports Y?", "where is Z documented?"), your FIRST action is a KG query, not Grep:

\`\`\`bash
python3 workspace/.kg/scripts/qq.py "MATCH (p:Package {name:'Chem'})-[:HAS_FEATURE]->(f:Feature) RETURN f.name LIMIT 10"
\`\`\`

First call self-installs the engine (~30s, one-time); after that, milliseconds. Read \`workspace/.kg/CLAUDE.md\` for the node/edge schema and query cookbook before composing non-trivial queries. Fall back only when the KG can't answer.

## Clarifying ambiguous requests

You have the AskUserQuestion tool available. Use it only when the **intent** itself is unclear — e.g. a request names no action, or two fundamentally different actions are equally plausible.`;

// Inlined into the system prompt. datagrok-exec is universal — it defines the contract for the
// datagrok_exec tool, which nearly every action-taking response uses. Everything else is loaded
// on demand via the Skill tool — skills' description triggers handle routing.
const INLINED_SKILL_NAMES = [
  'datagrok-exec',
  'datagrok-entities',
];

function loadInlinedSkills(): string {
  const sections: string[] = [];
  for (const name of INLINED_SKILL_NAMES) {
    const skillPath = `/app/plugin/skills/${name}/SKILL.md`;
    try {
      const raw = fs.readFileSync(skillPath, 'utf8');
      const body = raw.replace(/^---\n[\s\S]*?\n---\n+/, '').trim();
      sections.push(`### ${name}\n\n${body}`);
    } catch {
    }
  }
  return sections.join('\n\n---\n\n');
}

const INLINED_SKILLS = loadInlinedSkills();

// Build once at module load so the system prompt prefix is byte-stable across every turn and
// every user — required for Anthropic prompt-cache hits on the ~20-30 KB prefix.
const DATAGROK_SYSTEM_PROMPT = INLINED_SKILLS
  ? `${DATAGROK_PROMPT}\n\n## Inlined Skills\n\nThese skills are available in this context — invoke them directly without loading via the Skill tool. Each section gives the canonical signatures, conventions, and examples for one capability.\n\n${INLINED_SKILLS}`
  : DATAGROK_PROMPT;

const USER_WORKSPACE_PATTERN = /\/users\/[\w.-]+\/workspace/g;
const WORKSPACE_ACCESS_PATTERN = /\/workspace(?:[/"'\s\\]|$)/;

const blockWorkspaceAccess: HookCallback = async (input) => {
  if (input.hook_event_name !== 'PreToolUse')
    return {continue: true};
  const inputStr = JSON.stringify(input.tool_input ?? '');
  const stripped = inputStr.replace(USER_WORKSPACE_PATTERN, '');
  if (WORKSPACE_ACCESS_PATTERN.test(stripped)) {
    console.log(`PreToolUse: blocked /workspace reference in ${input.tool_name}`);
    return {
      decision: 'block',
      reason: 'Direct access to /workspace is blocked. Use cwd-relative `workspace/...` paths (e.g. `workspace/packages/Chem/...`) instead.',
    };
  }
  return {continue: true};
};

function buildSystemPrompt(mode?: string): string {
  if (mode === 'bash') return BASH_EXEC_PROMPT;
  if (mode === 'none') return '';
  return DATAGROK_SYSTEM_PROMPT;
}

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

// In-process MCP server whose datagrok_exec tool round-trips JS to the browser tab (via
// awaitBrowserInput) and returns the result. Synchronous, so Claude reports only what ran — fixes B2.
// Created per handleMessage so the ws/sid/active closure stays current.
function createBrowserExecServer(ws: WsSender, sid: string, active: ActiveQuery) {
  const asResult = (o: unknown) => ({content: [{type: 'text' as const, text: JSON.stringify(o)}]});
  return createSdkMcpServer({
    name: 'datagrok-browser',
    version: '1.0.0',
    tools: [
      sdkTool(
        'datagrok_exec',
        'Run JavaScript in the Datagrok tab to perform an action (add viewer, filter, open file, ' +
        'upload data, …). Returns {success, returnValue?, error?}. For informational questions ' +
        '("how do I…", "what is…") answer in plain text — do NOT call this tool.',
        {code: z.string().describe(
          'Async JS (await works). Globals: grok, ui, DG, view, t. Return a plain object confirming the ' +
          'action (shape per the datagrok-exec skill), or an HTMLElement only to render output in chat.',
        )},
        async ({code}) => {
          try {
            return asResult(await awaitBrowserInput(ws, sid, active, 'datagrok_exec', {code}));
          } catch (e: any) {
            return asResult({success: false, error: e.message});
          }
        },
      ),
      sdkTool(
        'datagrok_show_entities',
        'Display Datagrok entities (files, scripts, queries, connections, projects, spaces, groups, users) ' +
        'as interactive cards in the chat. ALWAYS call this tool when you have entity data from an MCP result — ' +
        'never list entities as plain text, markdown links, or a fenced block.',
        {entities: z.array(entityRefSchema).describe(
          'Array of entity refs. Each entry: ' +
          'type (file|script|query|connection|project|space|group|user) or the raw #type string from the MCP result — the tool normalizes it; ' +
          'name (for users, use friendlyName); ' +
          'id (required for all types except file); ' +
          'connector + path (required for files; connector is the connection name, e.g. "System:DemoFiles"); ' +
          'isDirectory (optional, files only). Pass all values from the MCP result exactly — never invent.',
        )},
        async ({entities}) => {
          try {
            return asResult(await awaitBrowserInput(ws, sid, active, 'datagrok_show_entities', {entities}));
          } catch (e: any) {
            return asResult({success: false, error: e.message});
          }
        },
      ),
    ],
  });
}

function buildMcpServers(browserExecServer: ReturnType<typeof createBrowserExecServer>, apiKey?: string, mcpServerUrl?: string, userId?: string): Record<string, any> | undefined {
  const servers: Record<string, any> = {};

  servers['datagrok-knowledge'] = createPackageKnowledgeServer(userId);
  servers['datagrok-browser'] = browserExecServer;

  const mcpUrl = mcpServerUrl || '';
  if (mcpUrl) {
    const apiUrl = apiUrlFromMcpUrl(mcpUrl);
    servers['datagrok'] = {
      type: 'http' as const,
      url: mcpUrl,
      headers: buildMcpHeaders(apiKey, apiUrl),
    };
  }

  return Object.keys(servers).length > 0 ? servers : undefined;
}

function buildOptions(
  browserExecServer: ReturnType<typeof createBrowserExecServer>,
  resume?: string, apiKey?: string, mcpServerUrl?: string,
  systemPromptMode?: string, userDir?: string,
  userId?: string,
  model?: ClaudeModel,
) {
  const systemPrompt = buildSystemPrompt(systemPromptMode);
  const mcpServers = buildMcpServers(browserExecServer, apiKey, mcpServerUrl, userId);
  // Bash and 'none' modes are minimal — don't pull in output-format skills.
  const loadPlugin = !systemPromptMode || (systemPromptMode !== 'bash' && systemPromptMode !== 'none');
  return {
    systemPrompt,
    allowedTools: ['Read', 'Glob', 'Grep', 'Edit', 'Write', 'Bash', 'WebSearch', 'WebFetch', 'AskUserQuestion'],
    // TODO: temporarily disabled — DB tools were getting picked up for unrelated prompts
    // (e.g. "convert GROKPEP-000002 sequence to molecule"). Re-enable once routing is tighter.
    disallowedTools: [
      'mcp__datagrok__db_list_catalogs', 'mcp__datagrok__db_list_schemas',
      'mcp__datagrok__db_list_tables', 'mcp__datagrok__db_describe_tables',
      'mcp__datagrok__db_list_joins', 'mcp__datagrok__db_try_sql',
    ],
    ...(loadPlugin ? {plugins: [{type: 'local' as const, path: '/app/plugin'}]} : {}),
    ...(mcpServers ? {mcpServers} : {}),
    strictMcpConfig: true,
    permissionMode: 'acceptEdits' as const,
    model: model ?? ClaudeModel.Sonnet,
    effort: 'low' as const,
    thinking: {type: 'disabled' as const},
    includePartialMessages: true,
    cwd: userDir || WORKSPACE,
    hooks: {
      PreToolUse: [{hooks: [blockWorkspaceAccess]}],
    },
    ...(resume ? {resume} : {}),
  };
}

const toolFormatters: {[K in ToolName]: (i: ToolInputs[K]) => string} = {
  Read: (i) => `Read ${i.file_path ?? ''}`,
  Write: (i) => `Write ${i.file_path ?? ''}`,
  Edit: (i) => `Edit ${i.file_path ?? ''}`,
  Bash: (_i) => 'Bash',
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
  index_codebase: (i) => `Index codebase ${i.path ?? ''}`.trim(),
  search_code: (i) => `Search code: ${i.query ?? ''}`,
  get_indexing_status: (i) => `Indexing status ${i.path ?? ''}`.trim(),
  clear_index: (i) => `Clear index ${i.path ?? ''}`.trim(),
  datagrok_exec: (_i) => 'Execute in browser',
  datagrok_show_entities: (i) => `Show ${(i.entities ?? []).length} entit${(i.entities ?? []).length === 1 ? 'y' : 'ies'}`,
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

interface WsSender {
  send(data: string): void;
}

function emit(ws: WsSender, msg: OutgoingMessage): void {
  ws.send(JSON.stringify(msg));
}

function emitChunk(ws: WsSender, sid: string, content: string): void {
  emit(ws, {type: 'chunk', sessionId: sid, content});
}

// Streams text from a Claude text_delta. Holds only the partial trailing line when it could
// still become a fence marker (starts with `` ` `` at line start); everything else is emitted
// immediately under the current mode. Complete lines are batched per same-mode group so a
// 50-line delta yields ~3 emits, not 50.
function emitFiltered(ws: WsSender, sid: string, text: string): void {
  const st = fenceStates.get(sid) ?? {mode: 'prose' as FenceMode, carry: '', lineInProgress: false};
  fenceStates.set(sid, st);

  const buf = st.carry + text;
  st.carry = '';
  const lastNl = buf.lastIndexOf('\n');

  if (lastNl >= 0) {
    const lines = buf.slice(0, lastNl).split('\n');
    let groupStart = 0;
    let groupMode = st.mode;
    for (let i = 0; i < lines.length; i++) {
      const couldBeFence = i > 0 || !st.lineInProgress;
      const fence = couldBeFence ? FENCE_RE.exec(lines[i]) : null;
      if (!fence) continue;

      if (i > groupStart)
        emitChunk(ws, sid, lines.slice(groupStart, i).join('\n') + '\n');

      emitChunk(ws, sid, lines[i] + '\n');
      st.mode = st.mode === 'prose' ? 'other' : 'prose';
      groupStart = i + 1;
      groupMode = st.mode;
    }
    if (groupStart < lines.length)
      emitChunk(ws, sid, lines.slice(groupStart).join('\n') + '\n');
    st.lineInProgress = false;
  }

  const partial = lastNl < 0 ? buf : buf.slice(lastNl + 1);
  if (partial.length === 0) return;

  if (!st.lineInProgress && partial.startsWith('`'))
    st.carry = partial;
  else {
    emitChunk(ws, sid, partial);
    st.lineInProgress = true;
  }
}

function flushFenceState(ws: WsSender, sid: string): void {
  const st = fenceStates.get(sid);
  if (st?.carry)
    emitChunk(ws, sid, st.carry);
  fenceStates.delete(sid);
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
  case 'stream_event':
    if (e.event?.delta?.type === 'text_delta' && e.event.delta.text)
      emitFiltered(ws, sid, e.event.delta.text);
    break;
  case 'tool_progress':
    emit(ws, {type: 'tool_activity', sessionId: sid, summary: `Running ${e.tool_name ?? ''}…`});
    break;
  case 'tool_use_summary':
    emit(ws, {type: 'tool_activity', sessionId: sid, summary: e.summary ?? ''});
    break;
  case 'result':
    flushFenceState(ws, sid);
    if (e.subtype === 'success')
      emit(ws, {type: 'final', sessionId: sid, content: e.result || '', ...(e.structured_output ? {structured_output: e.structured_output} : {})});
    else
      emit(ws, {type: 'error', sessionId: sid, message: (e.errors ?? []).join(', ') || e.subtype || 'unknown'});
    break;
  }
}

interface ActiveQuery {
  abortController: AbortController;
  queryHandle: ReturnType<typeof query> | null;
  // requestId → resolver, so parallel tool calls each await their own browser reply.
  pendingInputs: Map<string, (value: any) => void>;
}

// Round-trips a tool call to the browser: emits input_request, resolves on the matching
// input_response. The requestId keeps parallel calls from crossing wires; rejects on abort.
function awaitBrowserInput(ws: WsSender, sid: string, active: ActiveQuery, toolName: string, input: any): Promise<any> {
  const requestId = randomUUID();
  emit(ws, {type: 'input_request', sessionId: sid, requestId, toolName, input});
  return new Promise<any>((resolve, reject) => {
    active.pendingInputs.set(requestId, resolve);
    active.abortController.signal.addEventListener('abort', () => {
      if (active.pendingInputs.delete(requestId))
        reject(new Error('aborted'));
    }, {once: true});
  });
}

const activeQueries = new Map<string, ActiveQuery>();

function registerActiveQuery(sid: string, q: ActiveQuery): void {
  activeQueries.set(sid, q);
  markQueryStart();
}

function unregisterActiveQuery(sid: string): void {
  activeQueries.delete(sid);
  markQueryEnd();
}

function handleInputResponse(ws: WsSender, data: InputResponseMessage): void {
  const active = activeQueries.get(data.sessionId);
  if (!active)
    return;
  // Correlate by requestId; fall back to the sole pending request for older clients that omit it.
  const id = data.requestId ?? (active.pendingInputs.size === 1 ? active.pendingInputs.keys().next().value : undefined);
  const resolve = id !== undefined ? active.pendingInputs.get(id) : undefined;
  if (!resolve)
    return;
  active.pendingInputs.delete(id!);
  resolve(data.value);
}

async function handleMessage(ws: WsSender, data: UserMessage): Promise<void> {
  const sid = data.sessionId ?? '';
  let message = data.message ?? '';
  if (!message)
    return emit(ws, {type: 'error', sessionId: sid, message: 'Empty message'});

  // Don't block this turn on workspace git pull — it runs every 30 min in the background and
  // a stale read for one turn is fine.
  void awaitWorkspaceSync();

  const mcpUrl = rewriteForDocker(data.mcpServerUrl || '');
  const abortController = new AbortController();
  const active: ActiveQuery = {abortController, queryHandle: null, pendingInputs: new Map()};
  registerActiveQuery(sid, active);

  const userDir = data.apiKey ? await ensureUserDir(data.apiKey) : undefined;
  const userId = userDir ? path.basename(userDir) : undefined;

  const apiUrl = apiUrlFromMcpUrl(mcpUrl);
  if (apiUrl && data.apiKey) {
    // Fire-and-forget: file sync writes to disk in userDir; the model reads from disk on demand.
    syncUserFiles(apiUrl, data.apiKey).catch((e: any) =>
      console.warn('handleMessage: failed to sync user files:', e.message));
  }

  const DB_CLIENT_TOOLS = new Set([
    'mcp__datagrok__db_list_catalogs', 'mcp__datagrok__db_list_schemas',
    'mcp__datagrok__db_list_tables', 'mcp__datagrok__db_describe_tables',
    'mcp__datagrok__db_list_joins', 'mcp__datagrok__db_try_sql',
  ]);

  const browserExecServer = createBrowserExecServer(ws, sid, active);

  let gotResult = false;
  try {
    const existingSession = getSession(sid);
    const opts = buildOptions(browserExecServer, existingSession, data.apiKey, mcpUrl, data.systemPromptMode, userDir, userId, data.model);
    const canUseTool = async (toolName: string, input: any) => {
      if (toolName === 'AskUserQuestion' || DB_CLIENT_TOOLS.has(toolName)) {
        const updatedInput = await awaitBrowserInput(ws, sid, active, toolName, input);
        return {behavior: 'allow' as const, updatedInput};
      }
      return {behavior: 'allow' as const, updatedInput: input};
    };
    const outputFormat = data.outputSchema ? {type: 'json_schema' as const, schema: data.outputSchema as Record<string, unknown>} : undefined;
    const q = query({prompt: promptStream(message), options: {...opts, canUseTool, abortController, ...(outputFormat ? {outputFormat} : {})}});
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
    fenceStates.delete(sid);
    unregisterActiveQuery(sid);
  }
}

let authProc: ChildProcess | null = null;
let authWs: WsSender | null = null;

function handleAuthStart(ws: WsSender): void {
  if (authProc) {
    emit(ws, {type: 'auth_error', message: 'Authentication already in progress'});
    return;
  }
  authWs = ws;
  authProc = spawn('claude', ['auth', 'login'], {
    env: {...process.env, TERM: 'dumb', FORCE_COLOR: '0'},
    stdio: ['pipe', 'pipe', 'pipe'],
  });

  let urlSent = false;
  const onData = (data: Buffer) => {
    if (urlSent)
      return;
    const text = data.toString();
    const m = text.match(/https:\/\/claude\.com\/cai\/oauth[^\s]+/);
    if (m) {
      urlSent = true;
      emit(ws, {type: 'auth_url', url: m[0]});
    }
  };
  authProc.stdout?.on('data', onData);
  authProc.stderr?.on('data', onData);

  authProc.on('exit', (code) => {
    authProc = authWs = null;
    if (code === 0)
      emit(ws, {type: 'auth_done'});
    else
      emit(ws, {type: 'auth_error', message: `claude auth login exited with code ${code}`});
  });

  authProc.on('error', (err: Error) => {
    authProc = authWs = null;
    emit(ws, {type: 'auth_error', message: err.message});
  });
}

function handleAuthCode(ws: WsSender, code: string): void {
  if (!authProc?.stdin) {
    emit(ws, {type: 'auth_error', message: 'No authentication in progress'});
    return;
  }
  authProc.stdin.write(code + '\n');
}

function handleAbort(ws: WsSender, data: AbortMessage): void {
  const active = activeQueries.get(data.sessionId);
  if (!active)
    return;
  active.abortController.abort();
  try {
    if (active.queryHandle)
      active.queryHandle.close();
  } catch { /* query may have already finished */ }
  // Drop the SDK session: an aborted session can't be resumed. Leaving it here makes the
  // next message resume a dead session, so the following turn would fail to continue.
  sessions.delete(data.sessionId);
  flushFenceState(ws, data.sessionId);
  emit(ws, {type: 'aborted', sessionId: data.sessionId});
}

const app = new Hono();
const {injectWebSocket, upgradeWebSocket} = createNodeWebSocket({app});

app.get('/health', (c) => c.json({status: 'ok'}));

app.get('/ws', upgradeWebSocket(() => {
  return {
    onMessage(evt: any, ws: any) {
      const sender = ws as unknown as WsSender;
      let data: any;
      try {
        data = JSON.parse(String(evt.data));
      } catch {
        return emit(sender, {type: 'error', sessionId: '', message: 'Invalid JSON'});
      }

      if (data.apiKey)
        ensureUserDir(data.apiKey).catch((e: any) => console.warn('user-dir pre-hook:', e.message));

      if (data.type === 'abort') {
        handleAbort(sender, data);
        return;
      }

      if (data.type === 'input_response') {
        handleInputResponse(sender, data);
        return;
      }

      if (data.type === 'auth_start') {
        handleAuthStart(sender);
        return;
      }

      if (data.type === 'auth_code') {
        handleAuthCode(sender, data.code ?? '');
        return;
      }

      if (data.type === 'sync_user_files') {
        const mcpUrl = rewriteForDocker(data.mcpServerUrl || '');
        const apiUrl = apiUrlFromMcpUrl(mcpUrl);
        if (apiUrl && data.apiKey) {
          const scope = data.scope || 'all';
          const packageName = data.packageName;
          console.log(`sync_user_files: scope=${scope}, packageName=${packageName ?? '<none>'}`);
          (async () => {
            try {
              const result = await syncUserFiles(apiUrl, data.apiKey, scope, packageName);
              console.log(`sync_user_files: synced ${result.files.length} file(s) (scope=${scope})`);
              emit(sender, {type: 'sync_status', status: 'done', files: result.files});
            } catch (e: any) {
              console.warn('sync_user_files failed:', e.message);
              emit(sender, {type: 'sync_status', status: 'error', message: e.message});
            }
          })();
        }
        return;
      }

      if (data.type !== 'user_message') {
        return emit(sender, {
          type: 'error', sessionId: data.sessionId ?? '',
          message: `Unknown type: ${data.type}`,
        });
      }

      handleMessage(sender, data);
    },
  };
}));

app.notFound((c) => c.json({error: 'Not found'}, 404));
app.onError((err, c) => c.json({error: String(err)}, 500));

// Provider config arrives as container env, forwarded from the Grokky package credentials.
// Here we translate those into the env vars the Claude Agent SDK (which wraps Claude Code) reads
// at spawn. Field-name -> SDK-env mappings below mirror Claude Code's documented provider setup:
//   Bedrock  -> CLAUDE_CODE_USE_BEDROCK + AWS_REGION + (AWS_BEARER_TOKEN_BEDROCK | AWS_* IAM creds)
//              https://code.claude.com/docs/en/amazon-bedrock
//   Foundry  -> CLAUDE_CODE_USE_FOUNDRY + ANTHROPIC_FOUNDRY_RESOURCE + (ANTHROPIC_FOUNDRY_API_KEY | Entra ID)
//              https://code.claude.com/docs/en/microsoft-foundry
//   Anthropic-> ANTHROPIC_API_KEY
// The model aliases buildOptions() passes (sonnet/opus/haiku) resolve per provider via
// ANTHROPIC_DEFAULT_{OPUS,SONNET,HAIKU}_MODEL — Bedrock needs inference-profile ids, Foundry deployment names.
// Translate the injected credential fields into the SDK provider env, collecting any
// missing-required-credential problems so they surface in the container logs at startup.
function applyProviderConfig(): void {
  const e = process.env;
  const problems: string[] = [];

  const provider = e['provider'] || 'Anthropic';
  if (provider === 'Bedrock') {
    e['CLAUDE_CODE_USE_BEDROCK'] = '1';
    if (e['region'])
      e['AWS_REGION'] = e['region'];
    if (e['awsBearerToken'])
      e['AWS_BEARER_TOKEN_BEDROCK'] = e['awsBearerToken'];
    if (e['awsAccessKeyId'])
      e['AWS_ACCESS_KEY_ID'] = e['awsAccessKeyId'];
    if (e['awsSecretAccessKey'])
      e['AWS_SECRET_ACCESS_KEY'] = e['awsSecretAccessKey'];
    if (e['awsSessionToken'])
      e['AWS_SESSION_TOKEN'] = e['awsSessionToken'];
    if (!e['awsBearerToken'] && !(e['awsAccessKeyId'] && e['awsSecretAccessKey']))
      problems.push('Bedrock selected but no credentials — set awsBearerToken, or awsAccessKeyId + awsSecretAccessKey');
  }
  else if (provider === 'Microsoft Foundry') {
    e['CLAUDE_CODE_USE_FOUNDRY'] = '1';
    if (e['foundryResource'])
      e['ANTHROPIC_FOUNDRY_RESOURCE'] = e['foundryResource'];
    if (e['foundryApiKey'])
      e['ANTHROPIC_FOUNDRY_API_KEY'] = e['foundryApiKey'];
    if (!e['foundryResource'])
      problems.push('Microsoft Foundry selected but foundryResource is missing — required to reach the endpoint');
    if (!e['foundryApiKey'])
      problems.push('Microsoft Foundry selected without foundryApiKey — falls back to Entra ID, which is not configured in this container');
  }
  else {
    if (e['apiKey'])
      e['ANTHROPIC_API_KEY'] = e['apiKey'];
  }

  if (e['opusModel'])
    e['ANTHROPIC_DEFAULT_OPUS_MODEL'] = e['opusModel'];
  if (e['sonnetModel'])
    e['ANTHROPIC_DEFAULT_SONNET_MODEL'] = e['sonnetModel'];
  if (e['haikuModel'])
    e['ANTHROPIC_DEFAULT_HAIKU_MODEL'] = e['haikuModel'];

  for (var p of problems)
    console.warn(`[provider-config] ${p}`);
}

applyProviderConfig();

const usingBedrock = process.env['CLAUDE_CODE_USE_BEDROCK'] === '1';
const usingFoundry = process.env['CLAUDE_CODE_USE_FOUNDRY'] === '1';
const hasApiKey = !!process.env['ANTHROPIC_API_KEY'];
// Subscription auth requires the host's ~/.claude/.credentials.json to be mounted into the container at this path.
const hasSubscription = fs.existsSync('/home/grok/.claude/.credentials.json');
if (usingBedrock)
  console.log('Claude auth: using Amazon Bedrock');
else if (usingFoundry)
  console.log('Claude auth: using Microsoft Foundry');
else if (hasApiKey)
  console.log('Claude auth: using ANTHROPIC_API_KEY');
else if (hasSubscription)
  console.log('Claude auth: using subscription credentials at ~/.claude/.credentials.json');
else
  console.warn('Claude auth: no provider configured (no Bedrock/Foundry/ANTHROPIC_API_KEY and no ~/.claude/.credentials.json) — API calls will fail');

const server = serve({fetch: app.fetch, port: PORT});
injectWebSocket(server);
startWorkspaceSync();

console.log(`claude-runtime listening on :${PORT}`);
