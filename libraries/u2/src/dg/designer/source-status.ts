/* What a data source is doing right now, and the explicit run that changes it. A source keeps its
   own clock — nothing on the canvas belongs to it — so without this the designer shows a form full
   of empty fields whether the source is loading, was never run, or failed. Read through the DD5
   BindSource protocol, so every source answers the same three questions. */
import * as grok from 'datagrok-api/grok';
import {Component} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {isBindSource} from '../../spec/bind-source.js';
import type {BindSource} from '../../spec/bind-source.js';

export interface SourceStatus {
  /** The source's own vocabulary: idle, loading, ready or done, error. */
  state: string;
  /** How much it holds, or null where it holds nothing countable. */
  count: number | null;
  /** What {@link count} counts, for a message that reads like English. */
  counts: 'rows' | 'entities';
  error: string | null;
}

const LOADING = 'loading';

/** What `component` reports now, or null when it is not a data source (or has nothing to report).
 * Reads signals: called inside an effect, it follows the source. */
export function sourceStatus(component: Component): SourceStatus | null {
  if (!isBindSource(component))
    return null;
  const state = step(component, 'state');
  const rows = step(component, 'rowCount');
  const total = step(component, 'total');
  const error = step(component, 'error');
  if (state === undefined && rows === undefined && total === undefined)
    return null;
  const count = typeof rows === 'number' ? rows : typeof total === 'number' ? total : null;
  return {
    state: state === undefined || state === null ? '' : String(state),
    count,
    counts: typeof rows === 'number' ? 'rows' : 'entities',
    error: error === undefined || error === null ? null : String(error),
  };
}

/** The status as one line — the balloon a Refresh answers with, and the chip's hover text. */
export function statusText(at: SourceStatus): string {
  if (at.error !== null)
    return at.error;
  if (at.state === 'idle')
    return 'not run yet — use Refresh';
  const held = at.count === null ? '' : `${at.count} ${at.counts}`;
  if (at.state === '' || held === '')
    return at.state === '' ? held : at.state;
  return `${at.state} · ${held}`;
}

/** Runs the source once and says what came of it. This is the one run design time never makes by
 * itself (DD9), so it must be visible: a refresh that changes nothing on screen is the same as no
 * refresh at all — which is exactly how it read. */
export async function refreshSource(component: Component, name: string): Promise<void> {
  const refresh = component.getFunctions().find((f) => f.name === 'refresh');
  if (refresh === undefined)
    return;
  await refresh.apply({});
  const at = await settled(component);
  if (at === null)
    grok.shell.info(`${name}: refreshed`);
  else if (at.error !== null)
    grok.shell.error(`${name}: ${at.error}`);
  else
    grok.shell.info(`${name}: ${statusText(at)}`);
}

/** A refresh whose function answers nothing (a pager's) still lands eventually: one effect watches
 * until the source stops loading and disposes itself on the reading that resolves this. */
function settled(component: Component): Promise<SourceStatus | null> {
  const now = sourceStatus(component);
  if (now === null || now.state !== LOADING)
    return Promise.resolve(now);
  return new Promise<SourceStatus | null>((resolve) => {
    const scope = new Scope();
    scope.effect(() => {
      const at = sourceStatus(component);
      if (at === null || at.state !== LOADING) {
        resolve(at);
        scope.dispose();
      }
    });
  });
}

function step(source: BindSource, name: string): unknown {
  const at = source.bindStep(name);
  return at === null || isBindSource(at) ? undefined : at.value;
}
