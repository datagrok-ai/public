import {randomUUID} from 'node:crypto';
import type {AbortMessage, AgentEvent, InputResponseMessage, UserMessage} from './types';
import {raceIdle, IDLE} from './watchdog';
import {emit, emitFiltered, flushFenceState, clearFenceState,
  isRevising, beginRevision, endRevision} from './stream-filter';
import type {WsSender} from './stream-filter';
import {syncUserFiles} from './sync/orchestrator';
import {userDirs} from './user/user-dir';
import {markQueryStart, markQueryEnd} from './sync/workspace';
import {Verifier} from './verify';
import {GroundingGate, isSmallTalk} from './grounding';
import {buildSystemPrompt} from './prompts';
import {WORKSPACE, apiUrlFromMcpUrl, rewriteForDocker} from './constants';
import {waitForClaim, endTask} from './tasks';
import {runTurn as runAgentTurn} from './agent/turn';
import {hasSession, parkSession} from './agent/pool';

// Longest a turn may emit nothing before it is treated as wedged. Bounds silence between events,
// not total turn length.
const TURN_IDLE_TIMEOUT_MS = 90000;

// Visible streamed text of the in-flight turn, per session — the grounding gate reads it at Stop
// to decide whether the answer makes platform how-to claims (grounding.ts makesPlatformClaims).
const turnText = new Map<string, string>();

const TURN_LIMITS = {
  full: {maxTurns: 100, maxBudgetUsd: 2.5},
  bash: {maxTurns: 10, maxBudgetUsd: 0.25},
};

const NO_REVISION_RE = /^`?NO_REVISION`?\.?$/;

function forwardEvent(ws: WsSender, sid: string, event: AgentEvent, verifier?: Verifier): void {
  if (event.type === 'text') {
    if (!isRevising(sid))
      turnText.set(sid, (turnText.get(sid) ?? '') + event.text);
    return emitFiltered(ws, sid, event.text);
  }
  if (event.type === 'tool_activity')
    return emit(ws, {type: 'tool_activity', sessionId: sid, name: event.name, summary: event.summary});
  if (event.type === 'activity')
    return;
  const revising = endRevision(sid);
  flushFenceState(ws, sid);
  switch (event.type) {
    case 'result': {
      const noRevision = revising && NO_REVISION_RE.test(event.text.trim());
      return emit(ws, {
        type: 'final', sessionId: sid, content: noRevision ? '' : event.text,
        ...(revising ? {revision: noRevision ? 'kept' as const : 'replaced' as const} : {}),
        ...(event.structuredOutput !== undefined ? {structured_output: event.structuredOutput} : {}),
        ...(verifier?.exhausted ? {unverified: true} : {}),
        metrics: event.metrics,
      });
    }
    case 'limit':
      return emit(ws, {type: 'error', sessionId: sid,
        message: `This request hit its ${event.kind === 'turns' ? 'turn' : 'cost'} limit and was stopped. ` +
          'Try a narrower request.'});
    case 'session_reset':
      return emit(ws, {type: 'session_reset', sessionId: sid});
    case 'error':
      return emit(ws, {type: 'error', sessionId: sid, message: event.message});
  }
}

interface ActiveQuery {
  abortController: AbortController;
  /** The connection that sent the turn; only it can stop the turn (one chat can be open in several tabs). */
  ws: WsSender;
  // requestId → resolver, so parallel tool calls each await their own browser reply.
  pendingInputs: Map<string, (value: unknown) => void>;
}

// Round-trips a tool call to the browser: emits input_request, resolves on the matching
// input_response. The requestId keeps parallel calls from crossing wires; rejects on abort.
function awaitBrowserInput(
  sid: string, active: ActiveQuery, toolName: string, input: unknown, timeoutMs: number,
): Promise<unknown> {
  if (active.abortController.signal.aborted)
    return Promise.reject(new Error('aborted'));
  const requestId = randomUUID();
  emit(active.ws, {type: 'input_request', sessionId: sid, requestId, toolName, input});
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => {
      if (active.pendingInputs.delete(requestId))
        reject(new Error(`no browser response for ${toolName} after ${timeoutMs / 1000}s`));
    }, timeoutMs);
    active.pendingInputs.set(requestId, (value) => {
      clearTimeout(timer);
      resolve(value);
    });
    active.abortController.signal.addEventListener('abort', () => {
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

const QUEUED_HEARTBEAT_MS = 60_000;

const sessionChains = new Map<string, Promise<void>>();

export async function handleMessage(ws: WsSender, data: UserMessage): Promise<void> {
  const sid = data.sessionId ?? '';
  const message = data.message ?? '';
  const images = data.images;
  if (!message && !images?.length) {
    if (data.taskId)
      endTask(data.taskId);
    return emit(ws, {type: 'error', sessionId: sid, message: 'Empty message'});
  }
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
  // Registered before anything is awaited, so an abort that arrives at any point reaches this turn.
  const abortController = new AbortController();
  const active: ActiveQuery = {abortController, ws, pendingInputs: new Map()};
  registerActiveQuery(sid, active);

  let verifier: Verifier | undefined;
  let groundingGate: GroundingGate | undefined;
  try {
    // Reset at once: the browser resends the turn with its transcript, and that one waits for admission.
    if (data.resumeExpected && !hasSession(sid))
      return emit(ws, {type: 'session_reset', sessionId: sid});
    const mcpUrl = rewriteForDocker(data.mcpServerUrl || '');
    // Admission: wait for the celery claim (when the queued route is deployed) — the
    // queue container's concurrency limit is the only cap on concurrent turns.
    if (data.taskId) {
      emit(ws, {type: 'queued', sessionId: sid});
      const heartbeat = setInterval(() => emit(ws, {type: 'queued', sessionId: sid}), QUEUED_HEARTBEAT_MS);
      try {
        if (!await waitForClaim(data.taskId, sid, abortController.signal) && !abortController.signal.aborted)
          console.warn(`task[${sid}]: no claim for ${data.taskId} — failing open`);
      } finally {
        clearInterval(heartbeat);
      }
    }
    if (abortController.signal.aborted)
      return;

    const apiUrl = apiUrlFromMcpUrl(mcpUrl);
    const userDir = data.apiKey ? await userDirs.ensureDir(data.apiKey, apiUrl) : undefined;

    if (apiUrl && data.apiKey) {
      // First sync per user blocks so workspace/ exists for this turn; warm turns return
      // straight from disk (see orchestrator.ts).
      await syncUserFiles(apiUrl, data.apiKey).catch((e: any) =>
        console.warn('handleMessage: failed to sync user files:', e.message));
    }

    const mode = data.systemPromptMode === 'bash' || data.systemPromptMode === 'none' ?
      data.systemPromptMode : 'datagrok';
    const fullPromptTurn = mode === 'datagrok';
    // The visible answer stays until `final` says whether the hidden revision replaces it.
    const onGateBlock = () => {
      if (!isRevising(sid)) {
        flushFenceState(ws, sid);
        beginRevision(sid);
        emit(ws, {type: 'revision_start', sessionId: sid});
      }
    };
    const gates = data.gates ?? {};
    verifier = fullPromptTurn && gates.verify !== false ? new Verifier(onGateBlock) : undefined;
    turnText.delete(sid); // fresh answer buffer — the gate must judge THIS turn's text only
    groundingGate = fullPromptTurn && gates.grounding !== false && !data.outputSchema &&
      !isSmallTalk(data.message ?? '') ?
      new GroundingGate(onGateBlock, () => turnText.get(sid) ?? '') : undefined;

    const it = runAgentTurn({
      sessionKey: sid, text: message, images: data.images, systemPrompt: buildSystemPrompt(mode),
      effort: fullPromptTurn ? 'high' : undefined, model: data.model,
      cwd: userDir || WORKSPACE, userDir, resumeExpected: data.resumeExpected,
      outputSchema: data.outputSchema, clientTools: data.clientTools,
      datagrokMcp: mcpUrl ? {url: mcpUrl, headers: {
        ...(data.apiKey ? {'x-user-api-key': data.apiKey, 'Authorization': data.apiKey} : {}),
        ...(apiUrl ? {'x-datagrok-api-url': apiUrl} : {}),
      }} : undefined,
      limits: TURN_LIMITS[fullPromptTurn ? 'full' : 'bash'],
      gates: [verifier, groundingGate].filter((g) => g !== undefined),
      awaitBrowser: (toolName, input, timeoutMs) => awaitBrowserInput(sid, active, toolName, input, timeoutMs),
      signal: abortController.signal,
    });
    // Each event is raced against the idle timeout, so a turn that goes silent is bounded (see
    // watchdog.ts): a wedged agent emits nothing at all, which a plain await would wait on forever.
    let pendingStep: Promise<IteratorResult<AgentEvent, undefined>> | null = null;
    for (;;) {
      pendingStep = pendingStep ?? it.next();
      const step = await raceIdle(pendingStep, TURN_IDLE_TIMEOUT_MS);
      if (step === IDLE) {
        // Work inside the agent sends activity; a browser round-trip (a tool running in the tab, a question the user
        // has not answered) is silent but carries its own timeout.
        if (active.pendingInputs.size > 0)
          continue;
        console.warn(`watchdog[${sid}]: no agent event for ${TURN_IDLE_TIMEOUT_MS / 1000}s — killing turn`);
        emit(ws, {type: 'error', sessionId: sid,
          message: 'The assistant stopped responding and the turn was cancelled. Please try again.'});
        abortController.abort();
        parkSession(sid, `turn ${sid} went silent`);
        break;
      }
      pendingStep = null;
      if (step.done || abortController.signal.aborted)
        break;
      // A command's answer is its output, the final message: models remark before the tool call despite the prompt.
      if (mode === 'bash' && step.value.type === 'text')
        continue;
      forwardEvent(ws, sid, step.value, verifier);
    }
  } catch (e: any) {
    if (!abortController.signal.aborted)
      emit(ws, {type: 'error', sessionId: sid, message: String(e.message || e)});
  } finally {
    // Every terminal path — final, error, abort, disconnect, watchdog kill — ends the
    // admission task here, so the celery slot can never outlive the turn it admitted.
    if (data.taskId)
      endTask(data.taskId);
    if (verifier?.hadActions)
      console.log(`verify[${sid}]: ${verifier.statsLine()}`);
    if (groundingGate)
      console.log(`grounding[${sid}]: ${groundingGate.summary()}`);
    if (activeQueries.get(sid) === active) {
      clearFenceState(sid);
      endRevision(sid);
      turnText.delete(sid);
    }
    unregisterActiveQuery(sid, active);
  }
}

function abortSessionQuery(sid: string, ws?: WsSender): boolean {
  const active = activeQueries.get(sid);
  if (!active)
    return false;
  active.abortController.abort();
  unregisterActiveQuery(sid, active);
  endRevision(sid);
  turnText.delete(sid);
  if (ws)
    flushFenceState(ws, sid);
  else
    clearFenceState(sid);
  return true;
}

export function handleAbort(ws: WsSender, data: AbortMessage): void {
  if (activeQueries.get(data.sessionId)?.ws === ws && abortSessionQuery(data.sessionId, ws))
    emit(ws, {type: 'aborted', sessionId: data.sessionId});
}

/** Stops the turns of these sessions that `ws` sent (all of them without `ws`: a revoked admission task). A dropped
 * connection can never receive results, and a turn nobody waits for would hold the workspace sync. */
export function handleDisconnect(sessionIds: Iterable<string>, ws?: WsSender): void {
  for (const sid of sessionIds) {
    if (!ws || activeQueries.get(sid)?.ws === ws)
      abortSessionQuery(sid);
  }
}
