import {createSdkMcpServer, tool as sdkTool} from '@anthropic-ai/claude-agent-sdk';
import type {HookCallback} from '@anthropic-ai/claude-agent-sdk';
import {z} from 'zod/v4';
import type {ToolInputs, McpInputs, ToolName, McpName, ClientToolDef} from './types';
import {ClaudeModel} from './types';
import {WORKSPACE} from './constants';
import type {Verifier} from './verify';
import type {GroundingGate} from './grounding';
import {buildSystemPrompt} from './prompts';

// ---------------------------------------------------------------------------
// Workspace access guard — PreToolUse hook blocking absolute /workspace paths
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// URL plumbing — container-to-host rewrites and MCP request headers
// ---------------------------------------------------------------------------

export function rewriteForDocker(url: string): string {
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

export function apiUrlFromMcpUrl(mcpUrl: string): string | undefined {
  const idx = mcpUrl.indexOf('/docker/containers/proxy/');
  return idx > 0 ? mcpUrl.substring(0, idx) : undefined;
}

// ---------------------------------------------------------------------------
// Browser MCP tools — in-process server whose tools execute in the Datagrok tab
// ---------------------------------------------------------------------------

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

const EXEC_TIMEOUT_MS = 300_000;
const VERIFY_TIMEOUT_MS = 60_000;
const SHOW_ENTITIES_TIMEOUT_MS = 30_000;
const VIEW_TOOL_TIMEOUT_MS = 120_000;

type AwaitInput = (toolName: string, input: any, timeoutMs: number) => Promise<any>;

const asResult = (o: unknown) => ({content: [{type: 'text' as const, text: JSON.stringify(o)}]});

// In-process MCP server whose datagrok_exec tool round-trips JS to the browser tab (via
// awaitInput) and returns the result. Synchronous, so Claude reports only what ran — fixes B2.
// Created per handleMessage so the ws/sid/active closure behind awaitInput stays current.
export function createBrowserExecServer(awaitInput: AwaitInput) {
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
            return asResult(await awaitInput('datagrok_exec', {code}, EXEC_TIMEOUT_MS));
          } catch (e: any) {
            return asResult({success: false, error: e.message});
          }
        },
      ),
      sdkTool(
        'datagrok_verify',
        'Confirm a previous datagrok_exec / MCP action actually took effect, by re-reading LIVE state. ' +
        'Runs in a fresh scope (globals: grok, ui, DG, view, t) — it CANNOT see variables from your action ' +
        'code, so it must re-derive from t/view/grok. Returns {passed, observed?, error?}. After any action ' +
        'you MUST call this before reporting success; passed:false means it did NOT work — fix it, do not report done.',
        {
          assertion: z.string().describe(
            'Async JS (await works) that MUST `return` the observed result; truthy = verified. Re-read ' +
            'whatever you changed from t/view/grok. Never return a bare literal.',
          ),
          description: z.string().describe('Short human-readable statement of what is being verified.'),
        },
        async ({assertion, description}) => {
          try {
            return asResult(await awaitInput('datagrok_verify', {assertion, description}, VERIFY_TIMEOUT_MS));
          } catch (e: any) {
            return asResult({passed: false, error: e.message});
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
            return asResult(await awaitInput('datagrok_show_entities', {entities}, SHOW_ENTITIES_TIMEOUT_MS));
          } catch (e: any) {
            return asResult({success: false, error: e.message});
          }
        },
      ),
    ],
  });
}

// ---------------------------------------------------------------------------
// View tools — in-process server built from browser-declared tool definitions.
// Each call round-trips to the browser (input_request), where the current view's
// tool implementation runs. Declared fresh per turn, so the tool set follows the view.
// ---------------------------------------------------------------------------

const RESERVED_TOOL_NAMES = new Set(['datagrok_exec', 'datagrok_verify', 'datagrok_show_entities']);
const TOOL_NAME_RE = /^[A-Za-z0-9_-]+$/;

function zodShapeFromJsonSchema(schema: any): Record<string, z.ZodType> {
  const shape: Record<string, z.ZodType> = {};
  const required = new Set<string>(Array.isArray(schema?.required) ? schema.required : []);
  for (const [key, p] of Object.entries<any>(schema?.properties ?? {})) {
    let t: z.ZodType =
      p?.type === 'string' ? (Array.isArray(p.enum) && p.enum.length ? z.enum(p.enum) : z.string()) :
      p?.type === 'number' || p?.type === 'integer' ? z.number() :
      p?.type === 'boolean' ? z.boolean() :
      p?.type === 'array' ? z.array(z.any()) :
      z.any();
    if (p?.description)
      t = t.describe(p.description);
    if (!required.has(key))
      t = t.optional();
    shape[key] = t;
  }
  return shape;
}

export function createViewToolsServer(awaitInput: AwaitInput, tools: ClientToolDef[]) {
  const safe = tools.filter((t) =>
    t?.name && TOOL_NAME_RE.test(t.name) && !RESERVED_TOOL_NAMES.has(t.name) && t.description);
  if (safe.length === 0)
    return undefined;
  return createSdkMcpServer({
    name: 'datagrok-view',
    version: '1.0.0',
    tools: safe.map((t) => sdkTool(
      t.name,
      t.description,
      zodShapeFromJsonSchema(t.inputSchema),
      async (input: any) => {
        try {
          return asResult(await awaitInput(t.name, input, VIEW_TOOL_TIMEOUT_MS));
        } catch (e: any) {
          return asResult({success: false, error: e.message});
        }
      },
    )),
  });
}

// ---------------------------------------------------------------------------
// SDK query options — assembles mcpServers, hooks, and per-mode settings for query()
// ---------------------------------------------------------------------------

function buildMcpServers(
  browserExecServer: ReturnType<typeof createBrowserExecServer>, apiKey?: string, mcpServerUrl?: string,
  viewToolsServer?: ReturnType<typeof createViewToolsServer>,
): Record<string, any> | undefined {
  const servers: Record<string, any> = {};

  servers['datagrok-browser'] = browserExecServer;
  if (viewToolsServer)
    servers['datagrok-view'] = viewToolsServer;

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

export function buildOptions(
  browserExecServer: ReturnType<typeof createBrowserExecServer>,
  resume?: string, apiKey?: string, mcpServerUrl?: string,
  systemPromptMode?: string, userDir?: string,
  model?: ClaudeModel,
  forkSession?: boolean, resumeAt?: string,
  verifier?: Verifier,
  groundingGate?: GroundingGate,
  viewToolsServer?: ReturnType<typeof createViewToolsServer>,
) {
  const systemPrompt = buildSystemPrompt(systemPromptMode);
  const mcpServers = buildMcpServers(browserExecServer, apiKey, mcpServerUrl, viewToolsServer);
  // Bash and 'none' modes are minimal — no output-format skills, no reasoning. Full-prompt turns
  // get thinking + higher effort so the "ground answers in sources, don't answer from memory" rule
  // has a deliberation step to fire in (see DATAGROK_PROMPT) before the model commits an answer.
  const fullMode = !systemPromptMode || (systemPromptMode !== 'bash' && systemPromptMode !== 'none');
  const postToolHooks = [
    ...(verifier ? [verifier.postToolUse] : []),
    ...(groundingGate ? [groundingGate.postToolUse] : []),
  ];
  const stopHooks = [
    ...(verifier ? [verifier.stop] : []),
    ...(groundingGate ? [groundingGate.stop] : []),
  ];
  return {
    systemPrompt,
    allowedTools: ['Read', 'Glob', 'Grep', 'Edit', 'Write', 'Bash', 'WebSearch', 'WebFetch', 'AskUserQuestion'],
    ...(fullMode ? {plugins: [{type: 'local' as const, path: '/app/plugin'}]} : {}),
    ...(mcpServers ? {mcpServers} : {}),
    strictMcpConfig: true,
    permissionMode: 'acceptEdits' as const,
    model: model ?? ClaudeModel.Sonnet,
    effort: fullMode ? 'high' as const : 'low' as const,
    thinking: fullMode ? {type: 'enabled' as const, budgetTokens: 1500} : {type: 'disabled' as const},
    includePartialMessages: true,
    cwd: userDir || WORKSPACE,
    hooks: {
      PreToolUse: [{hooks: [blockWorkspaceAccess]}],
      ...(postToolHooks.length ? {PostToolUse: [{hooks: postToolHooks}]} : {}),
      ...(stopHooks.length ? {Stop: [{hooks: stopHooks}]} : {}),
    },
    // Only the turn after an abort forks (off resumeAt); normal turns resume in place.
    ...(resume ? {resume, ...(forkSession ? {forkSession: true} : {}), ...(resumeAt ? {resumeSessionAt: resumeAt} : {})} : {}),
  };
}

// ---------------------------------------------------------------------------
// Tool-activity summaries — one-line descriptions shown in chat while tools run
// ---------------------------------------------------------------------------

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
  datagrok_exec: (_i) => 'Execute in browser',
  datagrok_verify: (i) => `Verify: ${i.description ?? 'action outcome'}`,
  datagrok_show_entities: (i) => `Show ${(i.entities ?? []).length} entit${(i.entities ?? []).length === 1 ? 'y' : 'ies'}`,
};

export function toolSummary(name: string, input: Record<string, unknown>): string {
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
