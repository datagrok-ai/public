/* What a Section and an Accordion pane share about folding: the user's own open/close, remembered per
   browser under a key, and the header line that says what a folded body holds — or what is wrong in it. */
import {Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import type {Scope} from '../../core/scope.js';
import type {Text} from '../../core/elements.js';

export class Fold {
  /** The open/close stored under `key`, or null when the user never chose one. */
  static stored(kind: string, key: string | undefined): boolean | null {
    if (key === undefined)
      return null;
    try {
      const value = globalThis.localStorage?.getItem(`u2.${kind}.${key}`);
      return value === 'true' ? true : value === 'false' ? false : null;
    }
    catch {
      return null;
    }
  }

  static store(kind: string, key: string | undefined, expanded: boolean): void {
    if (key === undefined)
      return;
    try {
      globalThis.localStorage?.setItem(`u2.${kind}.${key}`, String(expanded));
    }
    catch {
    }
  }

  /** Writes `text` into `hint` and follows `validity`: `<cls>-invalid` on `root`, aria-invalid and the
   * message as the header's tooltip (else `tooltip`), and — while folded — the message in place of the
   * hint, so a fold never hides an error. An element hint (a summary with swatches) is kept, hidden
   * meanwhile. */
  static bindHint(scope: Scope, fold: FoldHint): void {
    const {root, header, hint, cls, text, validity, folded} = fold;
    let message = hint;
    if (text instanceof HTMLElement) {
      message = document.createElement('span');
      message.className = `${cls}-message`;
      hint.append(text, message);
    }
    scope.effect(() => {
      const problem = validity?.value ?? null;
      const shown = problem !== null && folded();
      root.classList.toggle(`${cls}-invalid`, problem !== null);
      if (problem === null)
        header.removeAttribute('aria-invalid');
      else
        header.setAttribute('aria-invalid', 'true');
      if (validity !== undefined)
        header.title = problem ?? fold.tooltip?.() ?? '';
      if (text instanceof HTMLElement) {
        text.hidden = shown;
        message.textContent = shown ? problem : '';
      }
      else
        hint.textContent = shown ? problem :
          text instanceof Signal ? String(text.value ?? '') : (text as string | undefined) ?? '';
    });
  }
}

export interface FoldHint {
  root: HTMLElement;
  header: HTMLElement;
  hint: HTMLElement;
  /** The fold's class prefix: `u2-section` marks the root `u2-section-invalid`. */
  cls: string;
  text: Text | HTMLElement | undefined;
  validity: ReadonlySignal<string | null> | undefined;
  /** Read in the effect, so the hint follows the fold. */
  folded: () => boolean;
  /** The header's tooltip while nothing is wrong — asked each time, so a renamed fold keeps its own. */
  tooltip?: () => string;
}
