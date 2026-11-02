import * as ui from 'datagrok-api/ui';
import {_package} from '../package';
import {featureOf} from './ai-feature';
import type {FinalEvent} from './runtime-client';

// A longer value fails the btree index on event_parameter_values, and the log batch with it.
const MAX_TEXT_BYTES = 2500;

type Turn = {turnId: string, logResponses: boolean, logToolCalls: boolean, chunks: string[], visibleText: string,
  toolCalls: string[]};

function contentLogging(): {prompts: boolean, responses: boolean, toolCalls: boolean} {
  const settings = _package.settings ?? {};
  return {
    prompts: String(settings.logPrompts) === 'true',
    responses: String(settings.logResponses) === 'true',
    toolCalls: String(settings.logToolCalls) === 'true',
  };
}

function capped(text: string): string {
  return new TextDecoder().decode(new TextEncoder().encode(text).slice(0, MAX_TEXT_BYTES));
}

export function contentLoggingNotice(): HTMLElement | null {
  const logging = contentLogging();
  const logged = [logging.prompts && 'prompts', logging.responses && 'responses', logging.toolCalls && 'tool calls']
    .filter((name) => name);
  return logged.length ? ui.divText(`Your administrators log your AI ${logged.join(', ')}.`,
    'grokky-content-log-notice') : null;
}

export class TurnLog {
  private turns = new Map<string, Turn[]>();

  start(sessionId: string, prompt?: string): void {
    const logging = contentLogging();
    const logsText = prompt != null;
    const turn: Turn = {turnId: crypto.randomUUID(), logResponses: logsText && logging.responses,
      logToolCalls: logsText && logging.toolCalls, chunks: [], visibleText: '', toolCalls: []};
    const queue = this.turns.get(sessionId);
    if (queue)
      queue.push(turn);
    else
      this.turns.set(sessionId, [turn]);
    if (prompt && logging.prompts)
      this.logContent(sessionId, turn.turnId, {prompt: capped(prompt)});
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
    const turnId = turn?.turnId ?? crypto.randomUUID();
    const usage: {[key: string]: string | number} = {turnId, chatId: final.sessionId,
      feature: featureOf(final.sessionId), status};
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
    if (!turn)
      return;
    const content: {[key: string]: string} = {};
    if (turn.logResponses) {
      const structured = final.structured_output ? JSON.stringify(final.structured_output) : '';
      const response = final.revision === 'replaced' ? final.content :
        turn.visibleText || turn.chunks.join('') || final.content || structured;
      if (response)
        content.response = capped(response);
    }
    if (turn.toolCalls.length)
      content.toolCalls = capped(turn.toolCalls.join('\n'));
    if (Object.keys(content).length)
      this.logContent(final.sessionId, turnId, content);
  }

  private logContent(sessionId: string, turnId: string, content: {[key: string]: string}): void {
    _package.logger.usage('ai-turn-content', {turnId, chatId: sessionId, feature: featureOf(sessionId), ...content},
      'ai-turn-content');
  }
}
