import {query} from '@anthropic-ai/claude-agent-sdk';
import type {SDKMessage} from '@anthropic-ai/claude-agent-sdk';

export interface InvokeResult {
  text: string;
  toolUses: Array<{name: string; input: unknown}>;
  usage?: {
    input_tokens?: number;
    output_tokens?: number;
    cache_read_input_tokens?: number;
    cache_creation_input_tokens?: number;
  };
  durationMs: number;
  error?: string;
  finishReason?: string;
}

async function* singleTurnPrompt(message: string) {
  yield {
    type: 'user' as const,
    session_id: '',
    message: {role: 'user' as const, content: message},
    parent_tool_use_id: null,
  };
}

/**
 * Send one prompt through the SDK with the Datagrok system prompt.
 *
 * Mirrors the production `buildOptions()` shape minus what's irrelevant for eval:
 * - No `plugins: [{type: 'local', path: '/app/plugin'}]` — we inline the skills
 *   into the system prompt instead, so the model sees them without round-tripping
 *   through the Skill tool. This is also what production effectively does today.
 * - No MCP servers — eval doesn't need any external state.
 * - `allowedTools: []` — pure prompt-only eval. We want to see what Claude WRITES,
 *   not what it does with Bash/Grep. Tool use would balloon cost and not test
 *   skill content. The skills are tested for their textual recommendations, not
 *   their agent-loop behavior.
 */
export async function invokeClaude(opts: {
  prompt: string;
  systemPrompt: string;
  model?: 'sonnet' | 'opus' | 'haiku';
  timeoutMs?: number;
  // --- Knowledge-layer research extensions. Defaults preserve the original
  // tools-off, no-MCP, no-plugin behavior so the existing skill eval `run.ts`
  // is unaffected. The _research/ arm runner sets these per arm. ---
  allowedTools?: string[];                              // [] => no tools (original behavior)
  mcpServers?: Record<string, unknown>;                 // per-arm MCP (semantic / graph)
  plugins?: Array<{type: 'local'; path: string}>;       // per-arm skill plugin(s)
  cwd?: string;                                          // working dir for Grep/Glob/Read
  permissionMode?: 'acceptEdits' | 'bypassPermissions' | 'default';
}): Promise<InvokeResult> {
  const t0 = Date.now();
  const abortController = new AbortController();
  const timer = setTimeout(() => abortController.abort(), opts.timeoutMs ?? 120_000);

  // The SDK shells out to the bundled Claude Code CLI, which refuses to start when
  // it detects CLAUDECODE=1 (nested-session guard). We're not really nesting — we're
  // a standalone Node script that happens to be launched by a Claude Code agent in
  // some scenarios — so strip these markers from the SDK's environment view.
  const prevClaudeCode = process.env['CLAUDECODE'];
  const prevEntry = process.env['CLAUDE_CODE_ENTRYPOINT'];
  delete process.env['CLAUDECODE'];
  delete process.env['CLAUDE_CODE_ENTRYPOINT'];

  const textChunks: string[] = [];
  const toolUses: Array<{name: string; input: unknown}> = [];
  let usage: InvokeResult['usage'];
  let finishReason: string | undefined;
  let errMsg: string | undefined;

  try {
    const q = query({
      prompt: singleTurnPrompt(opts.prompt),
      options: {
        systemPrompt: opts.systemPrompt,
        // Default [] => no tools (original skill-text eval). Arm runner passes a
        // whitelist (e.g. Grep/Glob/Read [+ Skill, + mcp__* tools]).
        allowedTools: opts.allowedTools ?? [],
        ...(opts.mcpServers ? {mcpServers: opts.mcpServers} : {}),
        ...(opts.plugins ? {plugins: opts.plugins} : {}),
        ...(opts.cwd ? {cwd: opts.cwd} : {}),
        model: opts.model ?? 'sonnet',
        // Match production: low effort, no thinking, accept edits (no-op without tools).
        // @ts-ignore — `effort` is in the runtime types but not always typed
        effort: 'low' as const,
        // @ts-ignore — thinking is in the runtime types
        thinking: {type: 'disabled' as const},
        // bypassPermissions when tools are on, so headless runs don't hang on prompts.
        permissionMode: opts.permissionMode ?? 'acceptEdits',
        includePartialMessages: false,
        abortController,
        strictMcpConfig: true,
      } as any,
    });

    for await (const event of q as AsyncIterable<SDKMessage>) {
      const e = event as any;
      if (event.type === 'assistant') {
        for (const block of e.message?.content ?? []) {
          if (block.type === 'text' && typeof block.text === 'string') {
            textChunks.push(block.text);
          } else if (block.type === 'tool_use') {
            toolUses.push({name: block.name, input: block.input});
          }
        }
      } else if (event.type === 'result') {
        if (e.subtype === 'success') {
          usage = e.usage;
          finishReason = 'success';
        } else {
          finishReason = e.subtype || 'unknown';
          errMsg = (e.errors ?? []).join(', ') || e.subtype || 'unknown';
        }
      }
    }
  } catch (e: any) {
    errMsg = String(e?.message ?? e);
  } finally {
    clearTimeout(timer);
    if (prevClaudeCode !== undefined) process.env['CLAUDECODE'] = prevClaudeCode;
    if (prevEntry !== undefined) process.env['CLAUDE_CODE_ENTRYPOINT'] = prevEntry;
  }

  return {
    text: textChunks.join(''),
    toolUses,
    usage,
    durationMs: Date.now() - t0,
    error: errMsg,
    finishReason,
  };
}
