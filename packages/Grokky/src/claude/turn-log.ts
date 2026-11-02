import {_package} from '../package';
import {featureOf} from './ai-feature';
import type {FinalEvent} from './runtime-client';

const MAX_TEXT_LENGTH = 60000;

type Turn = {prompt?: string, chunks: string[], visibleText: string, toolCalls: string[]};

/** What administrators chose to log besides usage, from the package settings. */
export function contentLogging(): {prompts: boolean, responses: boolean, toolCalls: boolean} {
  const settings = _package.settings ?? {};
  return {
    prompts: String(settings.logPrompts) === 'true',
    responses: String(settings.logResponses) === 'true',
    toolCalls: String(settings.logToolCalls) === 'true',
  };
}

/** Logs each completed turn as an `ai-turn` usage event, and the text the settings allow as `ai-turn-content`.
 * Separate events, because retention deletes whole events: the text can expire before the usage does.
 * Parameter names avoid `token`: the platform's log redactor blanks any parameter named so.
 * A session's turns are a queue: the runtime runs them in order, and a late `aborted` ends the older one. */
export class TurnLog {
  private turns = new Map<string, Turn[]>();

  start(sessionId: string, prompt?: string): void {
    const turn = {prompt, chunks: [], visibleText: '', toolCalls: []};
    const queue = this.turns.get(sessionId);
    if (queue)
      queue.push(turn);
    else
      this.turns.set(sessionId, [turn]);
  }

  chunk(sessionId: string, content: string): void {
    this.turns.get(sessionId)?.[0].chunks.push(content);
  }

  toolActivity(sessionId: string, name: string, summary: string): void {
    this.turns.get(sessionId)?.[0].toolCalls.push(`${name}: ${summary}`);
  }

  /** The chunks after this are the hidden revision; a final that keeps the original answer has no content. */
  revisionStart(sessionId: string): void {
    const turn = this.turns.get(sessionId)?.[0];
    if (!turn)
      return;
    turn.visibleText = turn.chunks.join('');
    turn.chunks = [];
  }

  drop(sessionId: string): Turn | undefined {
    const queue = this.turns.get(sessionId);
    const turn = queue?.shift();
    if (queue && !queue.length)
      this.turns.delete(sessionId);
    return turn;
  }

  clear(): void {
    this.turns.clear();
  }

  finish(final: FinalEvent, status: 'ok' | 'error' | 'aborted'): void {
    if (!final.sessionId)
      return;
    const turn = this.drop(final.sessionId);
    const metrics = final.metrics;
    const ids = {turnId: crypto.randomUUID(), chatId: final.sessionId, feature: featureOf(final.sessionId)};
    const usage: {[key: string]: string | number} = {...ids, status};
    const values: {[key: string]: string | number | null | undefined} = {
      model: metrics?.model,
      inputCount: metrics?.inputTokens, outputCount: metrics?.outputTokens,
      cacheReadCount: metrics?.cacheReadTokens, cacheCreationCount: metrics?.cacheCreationTokens,
      costUsd: metrics?.costUsd, numTurns: metrics?.numTurns,
      durationMs: metrics?.durationMs, apiDurationMs: metrics?.durationApiMs,
    };
    for (const [name, value] of Object.entries(values)) {
      if (value != null)
        usage[name] = value;
    }
    _package.logger.usage('ai-turn', usage, 'ai-turn');
    if (turn)
      this.logContent(ids, turn, final);
  }

  private logContent(ids: {[key: string]: string}, turn: Turn, final: FinalEvent): void {
    const logging = contentLogging();
    const content: {[key: string]: string} = {};
    if (logging.prompts && turn.prompt)
      content.prompt = turn.prompt.slice(0, MAX_TEXT_LENGTH);
    if (logging.responses) {
      const structured = final.structured_output ? JSON.stringify(final.structured_output) : '';
      const response = final.content || turn.visibleText || structured || turn.chunks.join('');
      if (response)
        content.response = response.slice(0, MAX_TEXT_LENGTH);
    }
    if (logging.toolCalls && turn.toolCalls.length)
      content.toolCalls = turn.toolCalls.join('\n').slice(0, MAX_TEXT_LENGTH);
    if (Object.keys(content).length)
      _package.logger.usage('ai-turn-content', {...ids, ...content}, 'ai-turn-content');
  }
}
