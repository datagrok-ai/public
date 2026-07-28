import type {OutgoingMessage} from './types';

// Client streaming: message emission plus the fence-aware text-delta filter.
//
// Deliberately free of any Claude Agent SDK import so it stays pure text processing —
// that keeps it unit-testable under any Node version (the SDK is ESM-only, which a
// CommonJS require() cannot load before Node 22).

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

// Sessions whose turn is in a gate-demanded revision: the already-streamed answer stays visible in
// the browser, the revision text is NOT streamed, and the `result` decides kept vs replaced.
const revisingSids = new Set<string>();

export const isRevising = (sid: string): boolean => revisingSids.has(sid);
export const beginRevision = (sid: string): void => void revisingSids.add(sid);
/** Ends a revision, reporting whether one was in progress. */
export const endRevision = (sid: string): boolean => revisingSids.delete(sid);
export const clearFenceState = (sid: string): void => void fenceStates.delete(sid);

// Streams text from a Claude text_delta. Holds only the partial trailing line when it could
// still become a fence marker (starts with `` ` `` at line start); everything else is emitted
// immediately under the current mode. Complete lines are batched per same-mode group so a
// 50-line delta yields ~3 emits, not 50.
export function emitFiltered(ws: WsSender, sid: string, text: string): void {
  if (revisingSids.has(sid))
    return;
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

export function flushFenceState(ws: WsSender, sid: string): void {
  const st = fenceStates.get(sid);
  if (st?.carry)
    emitChunk(ws, sid, st.carry);
  fenceStates.delete(sid);
}
