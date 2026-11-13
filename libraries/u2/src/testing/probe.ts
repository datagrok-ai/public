import {operability} from 'datagrok-api/u2core';
import type {OperabilityRule} from 'datagrok-api/u2core';
import {Names, Query} from './query.js';

/** Options of {@link expectOperable}. */
export interface ExpectOperableOptions {
  /** More selectors of what a user operates there, beyond native controls and interactive roles. */
  interactive?: string[];
}

/** Fails, listing every finding, when anything under `within` a user is offered cannot be operated as offered
 * (`DG.U2.operability`): covered or taking no pointer, inert or disabled with no reason, an icon button whose
 * accessible name is not its tooltip, an option taller than its rail. Hit-tests where the page has layout.
 * @example
 * expectOperable(dialog); */
export function expectOperable(within: Element | {root: Element}, options?: ExpectOperableOptions): void {
  const root = Query.element(within);
  const off = Query.hidden(root);
  if (off !== undefined)
    throw new Error(`expectOperable: ${off}: nothing in it can be operated`);
  const findings = operability(root, {hitTest: Reach.layout(), interactive: options?.interactive});
  if (findings.length > 0) {
    throw new Error(`expectOperable: ${findings.length} control(s) in ${Names.describe(root)} cannot be operated ` +
      `as offered:\n${findings.map((f) => `  - ${f.rule}: ${f.detail}`).join('\n')}`);
  }
}

const BLOCKING: OperabilityRule[] = ['covered', 'no-pointer', 'inert', 'disabled-without-reason'];
const POINTER: OperabilityRule[] = ['covered', 'no-pointer'];

export class Reach {
  static layout(): boolean {
    return typeof document.elementFromPoint === 'function';
  }

  static blocked(el: HTMLElement, pointer: boolean): string | undefined {
    const finding = operability(el, {hitTest: pointer && Reach.layout()}).find((f) => BLOCKING.includes(f.rule) &&
      (pointer || !POINTER.includes(f.rule)) && (f.element === el || f.element.contains(el)));
    if (finding !== undefined)
      return finding.detail;
    const hidden = Query.hidden(el);
    if (hidden !== undefined)
      return hidden;
    const who = Names.called(el);
    const inert = el.closest('[inert]');
    if (inert !== null)
      return `${who} is in ${Names.describe(inert)}, which is inert${Reach.because(inert)}`;
    const off = el.closest('[aria-disabled="true"], [disabled]');
    if (off !== null)
      return `${who} is disabled${Reach.because(off) || Reach.because(el)}`;
    const source = pointer ? Reach.noPointer(el) : undefined;
    return source === undefined ? undefined : `${who} takes no pointer: pointer-events: none` +
      (source === el ? '' : ` on ${Names.describe(source)}`);
  }

  static noPointer(el: Element): Element | undefined {
    let source: Element | undefined;
    for (let at: Element | null = el; at !== null; at = at.parentElement) {
      const value = getComputedStyle(at).pointerEvents;
      if (value === 'none')
        source = at;
      else if (value)
        break;
    }
    return source;
  }

  static because(el: Element): string {
    const reason = Names.squash(el.getAttribute('aria-description') ?? '') ||
      Names.squash(el.getAttribute('title') ?? '');
    return reason === '' ? '' : `: ${reason}`;
  }
}
