import * as ui from 'datagrok-api/ui';
import {_package} from '../package';
import {featureOf} from './ai-feature';
import type {FinalEvent} from './runtime-client';

const MAX_TEXT_LENGTH = 60000;

type Turn = {prompt?: string, logResponses: boolean, logToolCalls: boolean, chunks: string[], visibleText: string, toolCalls: string[]};

export function contentLogging(): {prompts: boolean, responses: boolean, toolCalls: boolean} {
  const settings = _package.settings ?? {};
  return {
    prompts: String(settings.logPrompts) === 'true',
    responses: String(settings.logResponses) === 'true',
    toolCalls: String(settings.logToolCalls) === 'true',
  };
}

export function contentLoggingNotice(): HTMLElement | null {
  const logging = contentLogging();
  return logging.prompts || logging.responses || logging.toolCalls ?
    ui.divText('Your administrators log the text of your AI requests.', 'grokky-content-log-notice') : null;
}

export class TurnLog {
  private turns = new Map<string, Turn[]>();

  start(sessionId: string, prompt?: string): void {
    const logging = contentLogging();
    const turn: Turn = {prompt: logging.prompts ? prompt : undefined, logResponses: logging.responses,
      logToolCalls: logging.toolCalls, chunks: [], visibleText: '', toolCalls: []};
    const queue = this.turns.get(sessionId);
    if (queue)
      queue.push(turn);
    else
      this.turns.set(sessionId, [turn]);
  }

  chunk(sessionId: string, content: string): void {
    const turn = this.turns.get(sessionId)?.[0];
    if (turn?.logResponses)
      turn.chunks.push(content);
  }

  toolActivity(sessionId: string, name: string, summary: string): void {
    const turn = this.turns.get(sessionId)?.[0];
    if (turn?.logToolCalls)
      turn.toolCalls.push(`${name}: ${summary}`);
  }

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
    // Names avoid `token`: the log redactor blanks such parameters.
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
    const content: {[key: string]: string} = {};
    if (turn.prompt)
      content.prompt = turn.prompt.slice(0, MAX_TEXT_LENGTH);
    if (turn.logResponses) {
      const structured = final.structured_output ? JSON.stringify(final.structured_output) : '';
      const response = final.content || turn.visibleText || structured || turn.chunks.join('');
      if (response)
        content.response = response.slice(0, MAX_TEXT_LENGTH);
    }
    if (turn.toolCalls.length)
      content.toolCalls = turn.toolCalls.join('\n').slice(0, MAX_TEXT_LENGTH);
    if (Object.keys(content).length)
      _package.logger.usage('ai-turn-content', {...ids, ...content}, 'ai-turn-content');
  }
}
