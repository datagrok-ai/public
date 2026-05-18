import * as fs from 'node:fs';
import * as path from 'node:path';
import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {createNodeWebSocket} from '@hono/node-ws';
import {query} from '@anthropic-ai/claude-agent-sdk';
import type {SDKMessage, HookCallback} from '@anthropic-ai/claude-agent-sdk';
import type {UserMessage, AbortMessage, InputResponseMessage, OutgoingMessage, ToolInputs, McpInputs, ToolName, McpName} from './types';
import {ClaudeModel} from './types';
import {WORKSPACE} from './constants';
import {syncUserFiles} from './sync/orchestrator';
import {ensureUserDir} from './user/user-dir';
import {createPackageKnowledgeServer} from './package-knowledge-tool';
import {awaitWorkspaceSync, markQueryStart, markQueryEnd, startWorkspaceSync} from './sync/workspace';
const PORT = 5355;
const MAX_SESSIONS = 200;

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

## Clarifying ambiguous requests

You have the AskUserQuestion tool available. When the user's request could reasonably be interpreted in multiple ways (which plot type, which column, which kind of cleanup), you MUST use AskUserQuestion to clarify before acting. NEVER guess when there are multiple valid options.`;

// Inlined into the system prompt. Only datagrok-exec is universal (the fenced-block
// contract is needed for every response that emits code). Everything else is loaded
// on demand via the Skill tool — skills' description triggers handle routing.
const INLINED_SKILL_NAMES = [
  'datagrok-exec',
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

function buildMcpServers(apiKey?: string, mcpServerUrl?: string, userId?: string): Record<string, any> | undefined {
  const servers: Record<string, any> = {};

  servers['datagrok-knowledge'] = createPackageKnowledgeServer(userId);

  const mcpUrl = mcpServerUrl || '';
  if (mcpUrl) {
    const apiUrl = apiUrlFromMcpUrl(mcpUrl);
    servers['datagrok'] = {
      type: 'http' as const,
      url: mcpUrl,
      headers: buildMcpHeaders(apiKey, apiUrl),
    };
  }

  if (process.env['MILVUS_TOKEN']) {
    const env: Record<string, string> = {MILVUS_TOKEN: process.env['MILVUS_TOKEN']!};
    if (process.env['OPENAI_API_KEY'])
      env['OPENAI_API_KEY'] = process.env['OPENAI_API_KEY'];
    servers['claude-context'] = {
      command: 'claude-context-mcp',
      args: [] as string[],
      env,
    };
  }

  return Object.keys(servers).length > 0 ? servers : undefined;
}

function buildOptions(
  resume?: string, apiKey?: string, mcpServerUrl?: string,
  systemPromptMode?: string, userDir?: string,
  userId?: string,
  model?: ClaudeModel,
) {
  const systemPrompt = buildSystemPrompt(systemPromptMode);
  const mcpServers = buildMcpServers(apiKey, mcpServerUrl, userId);
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
  index_codebase: (i) => `Index codebase ${i.path ?? ''}`.trim(),
  search_code: (i) => `Search code: ${i.query ?? ''}`,
  get_indexing_status: (i) => `Indexing status ${i.path ?? ''}`.trim(),
  clear_index: (i) => `Clear index ${i.path ?? ''}`.trim(),
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
      emit(ws, {type: 'final', sessionId: sid, content: e.result || '', ...(e.structured_output ? {structured_output: e.structured_output} : {})});
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
  if (active?.pendingInputResolve) {
    const resolve = active.pendingInputResolve;
    active.pendingInputResolve = null;
    resolve(data.value);
  }
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
  const active: ActiveQuery = {abortController, queryHandle: null, pendingInputResolve: null};
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

  let gotResult = false;
  const systemPromptMode = data.systemPromptMode ?? 'datagrok';
  try {
    const existingSession = getSession(sid);
    const opts = buildOptions(existingSession, data.apiKey, mcpUrl, data.systemPromptMode, userDir, userId, data.model);
    const canUseTool = async (toolName: string, input: any) => {
      if (toolName === 'AskUserQuestion' || DB_CLIENT_TOOLS.has(toolName)) {
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
    unregisterActiveQuery(sid);
  }
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

const hasApiKey = !!process.env['ANTHROPIC_API_KEY'];
// Subscription auth requires the host's ~/.claude/.credentials.json to be mounted into the container at this path.
const hasSubscription = fs.existsSync('/home/grok/.claude/.credentials.json');
if (hasApiKey)
  console.log('Claude auth: using ANTHROPIC_API_KEY');
else if (hasSubscription)
  console.log('Claude auth: using subscription credentials at ~/.claude/.credentials.json');
else
  console.warn('Claude auth: no ANTHROPIC_API_KEY and no ~/.claude/.credentials.json — API calls will fail');

const server = serve({fetch: app.fetch, port: PORT});
injectWebSocket(server);
startWorkspaceSync();

console.log(`claude-runtime listening on :${PORT}`);
