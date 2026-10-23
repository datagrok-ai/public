import type {AgentEvent, AgentProcessEvent, AgentTurnRequest} from '../types';
import {toolSummary} from '../tool-summary';
import {bareToolName} from '../verify';
import {killDescendants} from '../watchdog';
import {acquire, dropSession, hasExited, hasSession, processForSession, SessionLost} from './pool';
import type {AgentProcess} from './pool';

class Channel<T> {
  private queue: T[] = [];
  private waiters: ((r: IteratorResult<T, undefined>) => void)[] = [];
  private closed = false;

  push(v: T): void {
    if (this.closed)
      return;
    const w = this.waiters.shift();
    if (w)
      w({value: v, done: false});
    else
      this.queue.push(v);
  }

  end(): void {
    this.closed = true;
    for (const w of this.waiters.splice(0))
      w({value: undefined, done: true});
  }

  next(): Promise<IteratorResult<T, undefined>> {
    if (this.queue.length)
      return Promise.resolve({value: this.queue.shift()!, done: false});
    if (this.closed)
      return Promise.resolve({value: undefined, done: true});
    return new Promise((resolve) => this.waiters.push(resolve));
  }
}

/** Runs one turn on the session's agent process: gate feedback continues it, browser calls go to the tab.
 * The iterator ends after the terminal event, or silently once `req.signal` aborts. */
export function runTurn(req: AgentTurnRequest): AsyncIterator<AgentEvent, undefined> {
  const chan = new Channel<AgentEvent>();
  let p: AgentProcess | undefined;
  let finish!: () => void;
  const settled = new Promise<void>((resolve) => { finish = resolve; });
  const onExit = (): void => {
    const last = p?.stderrTail.at(-1);
    chan.push({type: 'error', message: `The assistant stopped unexpectedly${last ? `: ${last}` : '.'}`});
    finish();
  };
  const sessionLost = (): void => {
    dropSession(req.sessionKey, 'stored session lost');
    chan.push({type: 'session_reset'});
  };

  const onEvent = (ev: AgentProcessEvent): void => {
    switch (ev.type) {
      case 'text':
      case 'activity':
        return chan.push(ev);
      case 'tool_start':
        return chan.push({type: 'tool_activity', name: bareToolName(ev.name), summary: toolSummary(ev.name, ev.args)});
      case 'tool_end':
        for (const g of req.gates)
          g.afterTool(ev.name, ev.args, ev.result);
        return chan.push({type: 'activity'});
      case 'browser_call':
        req.awaitBrowser(ev.name, ev.input, ev.timeoutMs).then(
          (value) => p?.send({type: 'browser_result', id: ev.id, value}),
          (e) => p?.send({type: 'browser_result', id: ev.id, error: String(e?.message ?? e)}));
        return chan.push({type: 'activity'});
      case 'stopped': {
        const feedback = req.gates.map((g) => g.atStop()).filter((f): f is string => !!f);
        if (feedback.length && !req.signal.aborted)
          return p?.send({type: 'continue', text: feedback.join('\n\n')});
        chan.push({type: 'result', text: ev.text, structuredOutput: ev.structuredOutput, metrics: ev.metrics});
        return finish();
      }
      case 'limit':
      case 'error':
        chan.push(ev);
        return finish();
      case 'session_lost':
        sessionLost();
        return finish();
      case 'cancelled':
        return finish();
    }
  };

  const drive = async (): Promise<void> => {
    if (req.signal.aborted)
      return;
    for (;;) {
      // Checked with nothing awaited before `acquire` takes the session, so it cannot vanish in between.
      if (req.resumeExpected && !hasSession(req.sessionKey))
        return chan.push({type: 'session_reset'});
      p = await acquire(req.sessionKey);
      // A cancelled turn ends a moment after the abort that ended it on this side.
      for (let prev = p.turn; prev; prev = p.turn)
        await prev.settled;
      if (req.signal.aborted)
        return;
      // The process may have exited as the turn before this one ended; the next acquire starts another.
      if (!hasExited(p))
        break;
    }
    p.turn = {events: onEvent, settled};
    p.proc.once('exit', onExit);
    p.turns++;
    p.send({type: 'turn', text: req.text, images: req.images, systemPrompt: req.systemPrompt, effort: req.effort,
      cwd: req.cwd, userDir: req.userDir, outputSchema: req.outputSchema,
      clientTools: req.clientTools, datagrokMcp: req.datagrokMcp, limits: req.limits, restore: p.restore});
    await settled;
  };

  req.signal.addEventListener('abort', () => {
    const proc = processForSession(req.sessionKey);
    if (proc?.turn?.settled === settled) {
      proc.send({type: 'cancel'});
      // Hermes stops the tool it is running, but not what that tool started in a session of its own.
      void settled.then(() => {
        if (!proc.disposed && proc.proc.pid)
          killDescendants(proc.proc.pid);
      });
    }
    chan.end();
  }, {once: true});

  drive().catch((e: any) => {
    if (e instanceof SessionLost)
      return sessionLost();
    if (!req.signal.aborted)
      chan.push({type: 'error', message: String(e?.message ?? e)});
  }).finally(() => {
    finish();
    p?.proc.off('exit', onExit);
    if (p?.turn?.settled === settled) {
      p.turn = null;
      p.lastUsed = Date.now();
    }
    chan.end();
  });
  return {next: () => chan.next()};
}
