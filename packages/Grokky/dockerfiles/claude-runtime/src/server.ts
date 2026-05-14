import * as fs from 'node:fs';
import * as path from 'node:path';
import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {createNodeWebSocket} from '@hono/node-ws';
import {query} from '@anthropic-ai/claude-agent-sdk';
import type {SDKMessage, HookCallback} from '@anthropic-ai/claude-agent-sdk';
import type {UserMessage, AbortMessage, InputResponseMessage, ExecResultMessage, ExecResultEntry, OutgoingMessage, ToolInputs, McpInputs, ToolName, McpName} from './types';
import {ClaudeModel} from './types';
import {WORKSPACE} from './constants';
import {syncUserFiles, generatePackageIndex} from './sync/orchestrator';
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
NEVER use \`/workspace/...\` as an absolute path — that path is outside your scope and access to it
is blocked.

Personal user knowledge files (if any) live in the \`agents/\` directory in your current working
directory. Use Glob or \`ls agents/\` to discover them when relevant.

## How to find a registered function

NEVER guess function names, parameter names, or signatures.

\`PACKAGES.md\` in your current working directory lists the packages installed on this instance.

**Inlined catalogs (do NOT call \`list_functions\` / \`get_function\` for these):** the system prompt below already inlines complete catalogs for the most common packages. If the user's intent maps to a function in one of these catalogs, emit the \`datagrok-exec\` block directly — skip MCP discovery. Inlined-catalog packages:
- \`Chem\`, \`Admetica\` — see \`datagrok-chem-toolkit\`
- \`Chembl\`, \`MolTrack\`, \`HitTriage\`, \`Curves\` — see \`datagrok-chem-data\`
- Other inlined skills (df-and-columns, filtering, selection, viewers, grid-customization, calc-column, projects) cover their respective surfaces.

For packages NOT covered by an inlined catalog: use \`list_functions(keyword?, package?)\` MCP tool to discover the function. It returns \`{name, friendlyName, namespace, description, params}\`. Use \`get_function(name)\` only when you need full parameter details beyond what \`list_functions\` returned.

For packages with a curated knowledge index, \`get_package_knowledge(packageName)\` returns the
authoritative apiRef and docs paths.

Only fall back to filesystem search when the request is about platform JS API surface
(types, classes, methods on DG.DataFrame, DG.Viewer, etc.) — those live in \`workspace/js-api/src/\`.

## Code Execution

When the user asks you to **do** something (add a viewer, modify data, run a script, etc.):
1. If the function is in an inlined catalog → emit the \`datagrok-exec\` block directly.
2. Otherwise resolve the function via \`list_functions\` (or \`get_package_knowledge\` if applicable), THEN emit the block.

NEVER emit a \`\`\`datagrok-exec block with a function name you didn't see in either the inlined catalog or the MCP discovery results. Inventing function names is the #1 cause of silent failures.
Regular \`\`\`javascript blocks are for explanations only and will NOT run.

## Projects: MCP tools + ONE datagrok-exec for current view/layout

Project assembly and sharing go through MCP tools — they are governed and audited. The
"current DataFrame" / "current view layout" only have ids inside the browser, so one
\`datagrok-exec\` block resolves those ids; everything else is MCP.

Recipe for "create a project named X with the current table view (and optionally share)":

1. Call \`create_project(name)\` ONCE. Save the returned id as \`projectId\`.
2. If (and only if) the user said "with the current table / view / layout / dataframe",
   emit ONE \`datagrok-exec\` block that saves the layout and returns the ids:

   \`\`\`datagrok-exec
   const layout = view.saveLayout();
   await grok.dapi.layouts.save(layout);
   return {tableInfoId: t.getTableInfo().id, layoutId: layout.id};
   \`\`\`

   The block's return value is delivered back to you on the next turn as an
   \`<exec-result>\` JSON block. Read \`tableInfoId\` and \`layoutId\` from there.
3. Call \`add_entity_to_project(projectId, tableInfoId)\` AND
   \`add_entity_to_project(projectId, layoutId)\` — both with the ids you just received.
4. If the user asked to share, call \`share_entity\` with the project UUID
   (the \`projectId\` from step 1) and the group name(s). Group names are plain strings;
   do NOT call \`get_group\` first.
5. Call \`get_project_link(projectId)\` ONCE and include the returned url in your reply.

Strict rules:

- Call \`create_project\` EXACTLY ONCE per turn — never twice, never in parallel.
- NEVER parallelize step 1 with step 2, or step 2 with step 3. Each step must complete
  (and you must have read its result) before the next begins. The exec block in step 2
  emits to a different channel than tool calls; if you fire it together with
  \`create_project\` in step 1, you do not yet have \`projectId\`, and you have not yet
  received the ids you need for step 3.
- NEVER use \`grok.dapi.projects.save()\`, \`project.addChild()\`, \`project.addLink()\`,
  or \`grok.dapi.permissions.grant\` inside a datagrok-exec block — those bypass the
  governed MCP path. The only legitimate uses of datagrok-exec in this workflow are:
  saving a NEW layout (\`grok.dapi.layouts.save\`) and reading ids
  (\`t.getTableInfo().id\`).
- If an MCP tool throws, surface the error verbatim — do NOT fall back to in-browser
  workarounds.

## Output formats

The Datagrok UI parses \`\`\`datagrok-*\`\`\` fenced blocks specially and renders them
as interactive elements. Each block has a corresponding skill with the exact JSON shape,
globals, and edge cases — open the skill before emitting the block.

## Clarifying ambiguous requests

You have the AskUserQuestion tool available. When the user's request could reasonably be interpreted
in multiple ways (which plot type, which column, which kind of cleanup), you MUST use AskUserQuestion
to clarify before acting. NEVER guess when there are multiple valid options.`;

// Core skills whose bodies are inlined into the system prompt so Claude does not need to
// load them via the Skill tool. Each saves a 2-3s tool-call round-trip per session.
const INLINED_SKILL_NAMES = [
  'datagrok-exec',
  'datagrok-df-and-columns',
  'datagrok-filtering',
  'datagrok-selection',
  'datagrok-viewers',
  'datagrok-grid-customization',
  'datagrok-calc-column',
  'datagrok-chem-toolkit',
  'datagrok-chem-data',
  'datagrok-projects',
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

const pendingExecResults = new Map<string, ExecResultEntry[]>();

function storeExecResults(clientId: string, results: ExecResultEntry[]): void {
  if (!results || results.length === 0)
    return;
  const existing = pendingExecResults.get(clientId);
  if (existing)
    existing.push(...results);
  else
    pendingExecResults.set(clientId, [...results]);
}

function consumePendingExecResults(clientId: string): ExecResultEntry[] | undefined {
  const r = pendingExecResults.get(clientId);
  if (r && r.length > 0) {
    pendingExecResults.delete(clientId);
    return r;
  }
  return undefined;
}

function formatExecResults(results: ExecResultEntry[]): string {
  const parts = results.map((r) => {
    if (r.ok) {
      const body = r.value === undefined ? '(no return value)' : JSON.stringify(r.value);
      return `block #${r.blockIndex}: ${body}`;
    }
    return `block #${r.blockIndex} FAILED: ${r.error ?? 'unknown error'}`;
  });
  return `<exec-result>\nThe previous datagrok-exec block(s) you emitted ran in the browser. Their return values:\n${parts.join('\n')}\n</exec-result>\n\n`;
}

function handleExecResult(data: ExecResultMessage): void {
  const sid = data.sessionId ?? '';
  if (!sid || !Array.isArray(data.results))
    return;
  storeExecResults(sid, data.results);
}

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

function buildMcpServers(apiKey?: string, mcpServerUrl?: string, userId?: string, apiUrlOverride?: string): Record<string, any> | undefined {
  const servers: Record<string, any> = {};

  servers['datagrok-knowledge'] = createPackageKnowledgeServer(userId);

  const mcpUrl = mcpServerUrl || '';
  if (mcpUrl) {
    const apiUrl = apiUrlOverride ?? apiUrlFromMcpUrl(mcpUrl);
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
  apiUrlOverride?: string,
) {
  const systemPrompt = buildSystemPrompt(systemPromptMode);
  const mcpServers = buildMcpServers(apiKey, mcpServerUrl, userId, apiUrlOverride);
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

  const pendingResults = consumePendingExecResults(sid);
  if (pendingResults)
    message = formatExecResults(pendingResults) + message;

  // Don't block this turn on workspace git pull — it runs every 30 min in the background and
  // a stale read for one turn is fine.
  void awaitWorkspaceSync();

  const mcpUrl = rewriteForDocker(data.mcpServerUrl || '');
  const apiUrlOverride = data.apiUrl ? rewriteForDocker(data.apiUrl) : undefined;
  const abortController = new AbortController();
  const active: ActiveQuery = {abortController, queryHandle: null, pendingInputResolve: null};
  registerActiveQuery(sid, active);

  const userDir = data.apiKey ? await ensureUserDir(data.apiKey) : undefined;
  const userId = userDir ? path.basename(userDir) : undefined;

  const apiUrl = apiUrlOverride ?? apiUrlFromMcpUrl(mcpUrl);
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
    if (userDir && systemPromptMode !== 'bash' && systemPromptMode !== 'none') {
      const packageIndex = await generatePackageIndex(userId);
      if (packageIndex)
        await fs.promises.writeFile(path.join(userDir, 'PACKAGES.md'), packageIndex);
    }
    const existingSession = getSession(sid);
    const opts = buildOptions(existingSession, data.apiKey, mcpUrl, data.systemPromptMode, userDir, userId, data.model, apiUrlOverride);
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

      if (data.type === 'exec_result') {
        handleExecResult(data);
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

async function preWarmSDK(): Promise<void> {
  const abortController = new AbortController();
  try {
    // Warm Anthropic's prompt cache for the *real* (Sonnet + full Datagrok prompt + plugin)
    // shape so the user's first turn doesn't pay the 3-4s cache-creation cost.
    const q = query({
      prompt: (async function*() {
        yield {
          type: 'user' as const, session_id: '',
          message: {role: 'user' as const, content: 'ping'},
          parent_tool_use_id: null,
        };
      })(),
      options: {
        systemPrompt: DATAGROK_SYSTEM_PROMPT,
        allowedTools: ['Read'],
        plugins: [{type: 'local' as const, path: '/app/plugin'}],
        model: ClaudeModel.Sonnet,
        effort: 'low' as const,
        thinking: {type: 'disabled' as const},
        permissionMode: 'acceptEdits' as const,
        abortController,
      } as any,
    });
    for await (const event of q) {
      if (event.type === 'system') {
        abortController.abort();
        try { q.close(); } catch { /* ignore */ }
        return;
      }
    }
  } catch {
  }
}
preWarmSDK();

console.log(`claude-runtime listening on :${PORT}`);
