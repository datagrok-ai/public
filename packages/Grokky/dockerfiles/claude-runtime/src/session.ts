import {randomUUID} from 'node:crypto';
import {query} from '@anthropic-ai/claude-agent-sdk';
import type {SDKMessage} from '@anthropic-ai/claude-agent-sdk';
import type {UserMessage, AbortMessage, InputResponseMessage, OutgoingMessage, ImageAttachment} from './types';
import {syncUserFiles} from './sync/orchestrator';
import {ensureUserDir} from './user/user-dir';
import {awaitWorkspaceSync, markQueryStart, markQueryEnd} from './sync/workspace';
import {Verifier} from './verify';
import {GroundingGate} from './grounding';
import {createBrowserExecServer, createViewToolsServer, toolSummary, buildOptions, rewriteForDocker, apiUrlFromMcpUrl} from './query-options';

// ---------------------------------------------------------------------------
// Session registry — SDK session ids and resume/fork points, LRU-capped
// ---------------------------------------------------------------------------

const MAX_SESSIONS = 200;

// forkNext is set on abort: the next turn forks off lastCleanUuid, dropping the aborted turn.
interface SessionRecord {sdkId: string; lastCleanUuid?: string; forkNext?: boolean}

const sessions = new Map<string, SessionRecord>();

// Last assistant uuid of the in-flight turn; committed to the record on a clean result.
const pendingUuid = new Map<string, string>();

function storeSession(clientId: string, sdkId: string, lastCleanUuid?: string): void {
  sessions.delete(clientId);
  sessions.set(clientId, {sdkId, lastCleanUuid});
  if (sessions.size > MAX_SESSIONS)
    sessions.delete(sessions.keys().next().value!);
}

function getSession(clientId: string): SessionRecord | undefined {
  const rec = sessions.get(clientId);
  if (rec !== undefined) {
    sessions.delete(clientId);
    sessions.set(clientId, rec);
  }
  return rec;
}

// ---------------------------------------------------------------------------
// Client streaming — emit plus the fence-aware text-delta filter
// ---------------------------------------------------------------------------

export interface WsSender {
  send(data: string): void;
}

export function emit(ws: WsSender, msg: OutgoingMessage): void {
  ws.send(JSON.stringify(msg));
}

function emitChunk(ws: WsSender, sid: string, content: string): void {
  emit(ws, {type: 'chunk', sessionId: sid, content});
}

type FenceMode = 'prose' | 'other';
interface FenceState { mode: FenceMode; carry: string; lineInProgress: boolean; }
const fenceStates = new Map<string, FenceState>();
const FENCE_RE = /^```([\w-]*)\s*$/;

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
    for (let i = 0; i < lines.length; i++) {
      const couldBeFence = i > 0 || !st.lineInProgress;
      const fence = couldBeFence ? FENCE_RE.exec(lines[i]) : null;
      if (!fence) continue;

      if (i > groupStart)
        emitChunk(ws, sid, lines.slice(groupStart, i).join('\n') + '\n');

      emitChunk(ws, sid, lines[i] + '\n');
      st.mode = st.mode === 'prose' ? 'other' : 'prose';
      groupStart = i + 1;
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

// ---------------------------------------------------------------------------
// SDK event forwarding — translates SDK stream events into the client protocol
// ---------------------------------------------------------------------------

function forwardEvent(ws: WsSender, sid: string, event: SDKMessage, verifier?: Verifier): void {
  const e = event as any;
  switch (event.type) {
  case 'assistant':
    if (e.uuid)
      pendingUuid.set(sid, e.uuid);
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
    if (e.subtype === 'success') {
      // Commit the resume point only on clean completion — aborted turns never reach here.
      if (e.session_id)
        storeSession(sid, e.session_id, pendingUuid.get(sid));
      pendingUuid.delete(sid);
      emit(ws, {
        type: 'final', sessionId: sid, content: e.result || '',
        ...(e.structured_output ? {structured_output: e.structured_output} : {}),
        ...(verifier?.exhausted ? {unverified: true} : {}),
        // Turn metrics the SDK already computed on the result message — surfaced for the
        // latency benchmark harness (see docs/BENCHMARK.md). All optional; older runtimes omit them.
        metrics: {
          inputTokens: e.usage?.input_tokens ?? null,
          outputTokens: e.usage?.output_tokens ?? null,
          cacheReadTokens: e.usage?.cache_read_input_tokens ?? null,
          cacheCreationTokens: e.usage?.cache_creation_input_tokens ?? null,
          costUsd: e.total_cost_usd ?? null,
          numTurns: e.num_turns ?? null,
          durationMs: e.duration_ms ?? null,
          durationApiMs: e.duration_api_ms ?? null,
        },
      });
    } else
      emit(ws, {type: 'error', sessionId: sid, message: (e.errors ?? []).join(', ') || e.subtype || 'unknown'});
    break;
  }
}

// ---------------------------------------------------------------------------
// Browser round-trips — tool calls sent to the Datagrok tab and correlated replies
// ---------------------------------------------------------------------------

interface ActiveQuery {
  abortController: AbortController;
  queryHandle: ReturnType<typeof query> | null;
  // requestId → resolver, so parallel tool calls each await their own browser reply.
  pendingInputs: Map<string, (value: any) => void>;
}

// Round-trips a tool call to the browser: emits input_request, resolves on the matching
// input_response. The requestId keeps parallel calls from crossing wires; rejects on abort.
function awaitBrowserInput(ws: WsSender, sid: string, active: ActiveQuery, toolName: string, input: any, timeoutMs?: number): Promise<any> {
  if (active.abortController.signal.aborted)
    return Promise.reject(new Error('aborted'));
  const requestId = randomUUID();
  emit(ws, {type: 'input_request', sessionId: sid, requestId, toolName, input});
  return new Promise<any>((resolve, reject) => {
    const timer = timeoutMs ? setTimeout(() => {
      if (active.pendingInputs.delete(requestId))
        reject(new Error(`no browser response for ${toolName} after ${timeoutMs / 1000}s`));
    }, timeoutMs) : undefined;
    active.pendingInputs.set(requestId, (value) => {
      if (timer)
        clearTimeout(timer);
      resolve(value);
    });
    active.abortController.signal.addEventListener('abort', () => {
      if (timer)
        clearTimeout(timer);
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

function unregisterActiveQuery(sid: string, q?: ActiveQuery): void {
  if (q && activeQueries.get(sid) !== q)
    return;
  activeQueries.delete(sid);
  markQueryEnd();
}

export function handleInputResponse(ws: WsSender, data: InputResponseMessage): void {
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

// ---------------------------------------------------------------------------
// Turn execution — per-session FIFO queue, the SDK query loop, and abort
// ---------------------------------------------------------------------------

const QUEUED_HEARTBEAT_MS = 60_000;

const sessionChains = new Map<string, Promise<void>>();

async function* promptStream(message: string, images?: ImageAttachment[]) {
  const content = images?.length ?
    [
      ...images.map((img) => ({
        type: 'image' as const,
        source: {type: 'base64' as const, media_type: img.mediaType, data: img.data},
      })),
      ...(message ? [{type: 'text' as const, text: message}] : []),
    ] :
    message;
  yield {
    type: 'user' as const,
    session_id: '',
    message: {role: 'user' as const, content},
    parent_tool_use_id: null,
  };
}

export async function handleMessage(ws: WsSender, data: UserMessage): Promise<void> {
  const sid = data.sessionId ?? '';
  const message = data.message ?? '';
  const images = data.images;
  if (!message && !images?.length)
    return emit(ws, {type: 'error', sessionId: sid, message: 'Empty message'});
  const prev = sessionChains.get(sid);
  let release!: () => void;
  const turn = new Promise<void>((resolve) => { release = resolve; });
  const tail = prev ? prev.then(() => turn) : turn;
  sessionChains.set(sid, tail);
  try {
    if (prev) {
      emit(ws, {type: 'queued', sessionId: sid});
      const heartbeat = setInterval(() => emit(ws, {type: 'queued', sessionId: sid}), QUEUED_HEARTBEAT_MS);
      try {
        await prev;
      } finally {
        clearInterval(heartbeat);
      }
    }
    await runTurn(ws, data, sid, message);
  } finally {
    release();
    if (sessionChains.get(sid) === tail)
      sessionChains.delete(sid);
  }
}

async function runTurn(ws: WsSender, data: UserMessage, sid: string, message: string): Promise<void> {
  const images = data.images;

  // Don't block this turn on workspace git pull — it runs every 30 min in the background and
  // a stale read for one turn is fine.
  void awaitWorkspaceSync();

  const mcpUrl = rewriteForDocker(data.mcpServerUrl || '');
  const abortController = new AbortController();
  const active: ActiveQuery = {abortController, queryHandle: null, pendingInputs: new Map()};
  registerActiveQuery(sid, active);

  let gotResult = false;
  let verifier: Verifier | undefined;
  let groundingGate: GroundingGate | undefined;
  try {
    const userDir = data.apiKey ? await ensureUserDir(data.apiKey) : undefined;

    const apiUrl = apiUrlFromMcpUrl(mcpUrl);
    if (apiUrl && data.apiKey) {
      // Fire-and-forget: file sync writes to disk in userDir; the model reads from disk on demand.
      syncUserFiles(apiUrl, data.apiKey).catch((e: any) =>
        console.warn('handleMessage: failed to sync user files:', e.message));
    }

    const awaitInput = (toolName: string, input: any, timeoutMs: number) =>
      awaitBrowserInput(ws, sid, active, toolName, input, timeoutMs);
    const browserExecServer = createBrowserExecServer(awaitInput);
    const viewToolsServer = data.clientTools?.length ? createViewToolsServer(awaitInput, data.clientTools) : undefined;
    const fullPromptTurn = !data.systemPromptMode || data.systemPromptMode === 'datagrok';
    verifier = fullPromptTurn ? new Verifier() : undefined;
    groundingGate = fullPromptTurn && !data.outputSchema ? new GroundingGate() : undefined;

    const rec = getSession(sid);
    const opts = buildOptions(browserExecServer, rec?.sdkId, data.apiKey, mcpUrl, data.systemPromptMode, userDir, data.model,
      rec?.forkNext, rec?.forkNext ? rec.lastCleanUuid : undefined, verifier, groundingGate, viewToolsServer);
    const canUseTool = async (toolName: string, input: any) => {
      if (toolName === 'AskUserQuestion') {
        const updatedInput = await awaitBrowserInput(ws, sid, active, toolName, input);
        return {behavior: 'allow' as const, updatedInput};
      }
      return {behavior: 'allow' as const, updatedInput: input};
    };
    const outputFormat = data.outputSchema ? {type: 'json_schema' as const, schema: data.outputSchema as Record<string, unknown>} : undefined;
    const q = query({prompt: promptStream(message, images), options: {...opts, canUseTool, abortController, ...(outputFormat ? {outputFormat} : {})}});
    active.queryHandle = q;
    for await (const event of q) {
      if (abortController.signal.aborted)
        break;
      if (event.type === 'result')
        gotResult = true;
      forwardEvent(ws, sid, event, verifier);
    }
  } catch (e: any) {
    if (!abortController.signal.aborted && (!gotResult || !/exited with code/i.test(String(e.message))))
      emit(ws, {type: 'error', sessionId: sid, message: String(e.message || e)});
  } finally {
    if (verifier?.hadActions)
      console.log(`verify[${sid}]: ${verifier.statsLine()}`);
    if (groundingGate)
      console.log(`grounding[${sid}]: ${groundingGate.summary()}`);
    if (activeQueries.get(sid) === active)
      fenceStates.delete(sid);
    unregisterActiveQuery(sid, active);
  }
}

function abortSessionQuery(sid: string, ws?: WsSender): boolean {
  const active = activeQueries.get(sid);
  if (!active)
    return false;
  // Close before abort so in-flight SDK control responses don't write to an aborted transport.
  try {
    if (active.queryHandle)
      active.queryHandle.close();
  } catch { /* query may have already finished */ }
  active.abortController.abort();
  unregisterActiveQuery(sid, active);
  const rec = sessions.get(sid);
  if (rec)
    rec.forkNext = true;
  pendingUuid.delete(sid);
  if (ws)
    flushFenceState(ws, sid);
  else
    fenceStates.delete(sid);
  return true;
}

export function handleAbort(ws: WsSender, data: AbortMessage): void {
  if (abortSessionQuery(data.sessionId, ws))
    emit(ws, {type: 'aborted', sessionId: data.sessionId});
}

// A dropped connection can never receive results — abort its in-flight queries so the SDK
// stops and the workspace sync isn't blocked by a query nobody is waiting for.
export function handleDisconnect(sessionIds: Iterable<string>): void {
  for (const sid of sessionIds)
    abortSessionQuery(sid);
}
