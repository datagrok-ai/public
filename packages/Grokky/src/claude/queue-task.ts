/** Queued-task admission for AI chat turns: the queued aiChatTurnTask call holds one celery
 *  slot for exactly the turn's duration while the turn streams over the browser WebSocket;
 *  no conversation content ever enters the queue or the platform logs. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ClaudeRuntimeClient, ClaudeModel} from './runtime-client';

// Bounds a hold whose turn end was never observed (runtime restarted mid-poll).
const HOLD_CAP_MS = 30 * 60_000;

// ---------------------------------------------------------------- worker side

/** Holds one admission slot: long-polls the runtime until it reports the turn ended. */
export async function holdChatTurnTask(sessionId: string, taskId: string): Promise<string> {
  const started = Date.now();
  const client = ClaudeRuntimeClient.getInstance();
  if (!await client.discover())
    throw new Error('claude-runtime container not found');
  const containerId = client.runtimeContainerId!;
  const post = async (path: string) => {
    const r = await grok.dapi.docker.dockerContainers.fetchProxy(containerId, path, {
      // fetchProxy authenticates via browser cookies; headless workers must send the token.
      method: 'POST', headers: {'Content-Type': 'application/json', 'Authorization': grok.dapi.token},
      body: JSON.stringify({taskId, sessionId}),
    });
    if (!r.ok)
      throw new Error(`${path}: HTTP ${r.status}`);
    return r.json();
  };
  for (;;) {
    try {
      (globalThis as any).DG_TASK_PROGRESS?.(null, `holding AI turn slot (${Math.round((Date.now() - started) / 1000)}s)`);
    } catch (_) {
      // Task revoked platform-side (admin cancel): give the slot back and abort the turn.
      await post('/task/release').catch(() => {});
      return JSON.stringify({status: 'released', heldMs: Date.now() - started});
    }
    const out = await post('/task/claim');
    if (out.status !== 'running')
      return JSON.stringify({...out, heldMs: Date.now() - started});
    if (Date.now() - started > HOLD_CAP_MS)
      return JSON.stringify({status: 'timeout', heldMs: Date.now() - started});
  }
}

// --------------------------------------------------------------- browser side

let taskFunc: DG.Func | null | undefined = undefined;

async function taskRoute(): Promise<DG.Func | null> {
  if (taskFunc === undefined) {
    try {
      taskFunc = await grok.dapi.functions.filter('name = "aiChatTurnTask" and options.queue = "true"').first() ?? null;
    } catch (e: any) {
      // Lookup failed (not "function absent") — retry on the next turn instead of caching null.
      console.warn(`Grokky: task route lookup failed: ${e?.message ?? e}`);
      return null;
    }
  }
  return taskFunc;
}

export interface ChatTurnOptions {
  systemPromptMode?: string;
  model?: ClaudeModel;
  clientTools?: {name: string; description: string; inputSchema?: object}[];
}

/** Sends one chat turn, admitted through a queued task when `aiChatTurnTask` is deployed,
 *  falling back to a plain send otherwise; a failed admission task is logged, not fatal. */
export async function sendChatTurn(
  client: ClaudeRuntimeClient, sessionId: string, message: string, options: ChatTurnOptions,
): Promise<void> {
  await client.ensureConnected();
  const func = await taskRoute();
  if (!func)
    return client.send(sessionId, message, options);
  // Send first: a task admitting a message that never reached the runtime would hold its slot until HOLD_CAP_MS.
  const taskId = crypto.randomUUID();
  client.send(sessionId, message, {...options, taskId});
  func.prepare({sessionId: sessionId, taskId: taskId}).call()
    .catch((e: any) => console.warn(`Grokky: turn admission task failed (turn proceeds unmetered): ${e?.message ?? e}`));
}
