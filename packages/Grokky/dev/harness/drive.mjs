// T2 driver — a fake browser for the claude-runtime WebSocket protocol.
//
// Speaks exactly what src/claude/runtime-client.ts speaks, so a turn exercises the real
// runtime: system prompt assembly, MCP servers, the verify/grounding gates, the revision
// protocol, session resume. What it does NOT need is a browser, a Datagrok instance, or
// pub serve — browser-side tool calls are answered from a per-case stub table instead.
//
// Everything a turn observably did comes back in the result: which tools ran (in order),
// how many browser round-trips, the code each exec ran, timings, and SDK token/cost metrics.

import {WebSocket} from 'ws';

export const DEFAULT_URL = process.env.GROKKY_DEV_WS ?? 'ws://localhost:5355/ws';

/** Browser-side tool answers. A stub is a literal value or (input) => value. */
const DEFAULT_STUBS = {
  // Claiming success unconditionally would let a turn "pass" the verify gate having done
  // nothing, so the default reports that no browser executed the code. Cases that care
  // override this with a scripted result.
  datagrok_exec: () => ({success: true, returnValue: null}),
  datagrok_verify: () => ({passed: true, observed: true}),
  datagrok_show_entities: () => ({success: true}),
  // Unanswerable headlessly — reply empty so the turn is not left hanging on the gate.
  AskUserQuestion: () => ({answers: {}}),
};

export class RuntimeDriver {
  constructor({url = DEFAULT_URL, apiKey, mcpServerUrl} = {}) {
    this.url = url;
    this.apiKey = apiKey;
    this.mcpServerUrl = mcpServerUrl;
    this.ws = null;
  }

  connect() {
    return new Promise((resolve, reject) => {
      this.ws = new WebSocket(this.url);
      const fail = (e) => reject(new Error(`cannot reach ${this.url} — is \`node dev/runtime.mjs up\` running? (${e.message ?? e})`));
      this.ws.once('open', () => {
        this.ws.removeListener('error', fail);
        resolve(this);
      });
      this.ws.once('error', fail);
    });
  }

  close() {
    this.ws?.close();
    this.ws = null;
  }

  /**
   * Runs one turn and resolves with everything it observably did.
   * @param prompt        user message
   * @param opts.sessionId    reuse to continue a conversation (resume path)
   * @param opts.model        'haiku' | 'sonnet' | 'opus' — omit for the runtime default
   * @param opts.systemPromptMode  'bash' | 'none' — omit for the full Datagrok prompt
   * @param opts.stubs        per-tool overrides, merged over DEFAULT_STUBS
   * @param opts.clientTools  view-tool defs to declare for this turn
   * @param opts.outputSchema forces structured output
   * @param opts.timeoutMs    abort and resolve with an error after this long
   */
  turn(prompt, opts = {}) {
    const sessionId = opts.sessionId ?? `drive-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
    const stubs = {...DEFAULT_STUBS, ...(opts.stubs ?? {})};
    const timeoutMs = opts.timeoutMs ?? 120_000;

    const started = Date.now();
    const out = {
      sessionId, ok: false, content: '', structured: undefined,
      ttftMs: null, totalMs: null,
      tools: [],          // bare tool names, in call order
      execCode: [],       // the JS each datagrok_exec asked the browser to run
      roundTrips: 0,      // browser input_request round-trips
      queued: false, unverified: false, revision: undefined,
      metrics: null, error: null,
    };

    return new Promise((resolve) => {
      let done = false;
      const finish = (patch = {}) => {
        if (done) return;
        done = true;
        clearTimeout(timer);
        this.ws.off('message', onMessage);
        out.totalMs = Date.now() - started;
        resolve(Object.assign(out, patch));
      };

      const timer = setTimeout(() => {
        this.ws.send(JSON.stringify({type: 'abort', sessionId}));
        finish({error: `turn timed out after ${timeoutMs / 1000}s`});
      }, timeoutMs);

      const onMessage = (raw) => {
        let m;
        try {
          m = JSON.parse(raw.toString());
        } catch {
          return;
        }
        // sync_status / auth_* carry no sessionId; everything else must match ours.
        if (m.sessionId && m.sessionId !== sessionId) return;

        switch (m.type) {
        case 'chunk':
          if (out.ttftMs === null) out.ttftMs = Date.now() - started;
          out.content += m.content;
          break;
        case 'tool_activity':
          if (m.name) out.tools.push(m.name);
          break;
        case 'queued':
          out.queued = true;
          break;
        case 'revision_start':
          out.revision = 'started';
          break;
        case 'input_request': {
          out.roundTrips++;
          if (m.toolName === 'datagrok_exec') out.execCode.push(m.input?.code ?? '');
          const stub = stubs[m.toolName];
          const value = typeof stub === 'function' ? stub(m.input) :
            stub !== undefined ? stub :
              {success: false, error: `no stub for tool '${m.toolName}'`};
          Promise.resolve(value).then((v) =>
            this.ws.send(JSON.stringify({type: 'input_response', sessionId, requestId: m.requestId, value: v})));
          break;
        }
        case 'final':
          // A 'replaced' revision supersedes the streamed text; 'kept' means it still stands.
          finish({
            ok: true,
            content: m.revision === 'replaced' ? m.content : (out.content || m.content || ''),
            structured: m.structured_output,
            unverified: !!m.unverified,
            revision: m.revision,
            metrics: m.metrics ?? null,
          });
          break;
        case 'error':
          finish({error: m.message});
          break;
        case 'aborted':
          finish({error: 'aborted'});
          break;
        }
      };

      this.ws.on('message', onMessage);
      this.ws.send(JSON.stringify({
        type: 'user_message', sessionId, message: prompt,
        ...(this.apiKey ? {apiKey: this.apiKey} : {}),
        ...(this.mcpServerUrl ? {mcpServerUrl: this.mcpServerUrl} : {}),
        ...(opts.model ? {model: opts.model} : {}),
        ...(opts.systemPromptMode ? {systemPromptMode: opts.systemPromptMode} : {}),
        ...(opts.clientTools?.length ? {clientTools: opts.clientTools} : {}),
        ...(opts.outputSchema ? {outputSchema: opts.outputSchema} : {}),
      }));
    });
  }
}

/** Convenience: connect, run one turn, disconnect. */
export async function driveOnce(prompt, opts = {}) {
  const d = await new RuntimeDriver(opts).connect();
  try {
    return await d.turn(prompt, opts);
  } finally {
    d.close();
  }
}
