import {Signal, ReadonlySignal} from './signals.js';
import {Ownership} from './owner.js';
import type {Owner} from './owner.js';

/** Why something is off — the text the user reads ("Pick a source file first") — or `null` (or `''`) while it is
 * on; a signal of either is followed. Never a boolean: a control that is off says why. */
export type Reason = string | null | ReadonlySignal<string | null>;

const isSignal = (x: unknown): x is ReadonlySignal<string | null> => x instanceof Signal;
const text = (x: string | null | undefined): string | null => x === undefined || x === '' ? null : x;

export const Reason = {
  /** What `reason` says now, not followed: `null` while on. */
  now(reason: Reason | undefined): string | null {
    return text(isSignal(reason) ? reason.peek() : reason);
  },

  /** What `reason` says, followed inside an effect or a `computed`: `null` while on. */
  read(reason: Reason | undefined): string | null {
    return text(isSignal(reason) ? reason.value : reason);
  },

  /** Hands `apply` what `reason` says now and, for a signal, on every change — owned by `owner`, else by whoever
   * is building. Nothing given, nothing runs. */
  follow(reason: Reason | undefined, apply: (reason: string | null) => void, owner?: Owner): void {
    if (reason === undefined)
      return;
    if (isSignal(reason))
      Ownership.scopeOf(owner, 'a live `disabled` reason').effect(() => apply(text(reason.value)));
    else
      apply(text(reason));
  },

  /** `el` off with `reason`, or on with `null`: `aria-disabled`, and the reason as its `aria-description` — what a
   * screen reader says on focus and the control state reports as `disabledReason`. */
  apply(el: HTMLElement, reason: string | null): void {
    if (!reason) {
      el.removeAttribute('aria-disabled');
      el.removeAttribute('aria-description');
    }
    else {
      el.setAttribute('aria-disabled', 'true');
      el.setAttribute('aria-description', reason);
    }
  },

  /** A region off with `reason`: inert — nothing in it takes a press or the focus — and dimmed, the reason its
   * `aria-description`, where the operability probe and the testing kit read it; on again with `null`. Say the
   * reason where the user sees it too: a region takes no pointer, so no tooltip of its own shows. */
  region(el: HTMLElement, reason: string | null): void {
    const off = text(reason);
    el.classList.toggle('u2-region-off', off !== null);
    if (off === null) {
      el.removeAttribute('inert');
      el.removeAttribute('aria-description');
    }
    else {
      el.setAttribute('inert', '');
      el.setAttribute('aria-description', off);
    }
  },

  /** The reason `el` is off for: its `aria-description` while `aria-disabled`, else `null`. */
  of(el: Element): string | null {
    return el.getAttribute('aria-disabled') === 'true' ? text(el.getAttribute('aria-description')) : null;
  },

  /** A u2 button is off by `aria-disabled`, never by the native attribute: it keeps its tab stop, so the keyboard
   * reaches it and its tooltip says why, while a click, Enter or Space on it does nothing. `el.disabled` reads and
   * writes that state; a native `disabled` (an input sets one on the controls inside it) still reads as off. */
  focusable(el: HTMLButtonElement): HTMLButtonElement {
    Object.defineProperty(el, 'disabled', {
      configurable: true,
      get: () => el.getAttribute('aria-disabled') === 'true' || el.hasAttribute('disabled'),
      set: (x: boolean) => {
        if (x)
          el.setAttribute('aria-disabled', 'true');
        else
          el.removeAttribute('aria-disabled');
      },
    });
    // registered first: the button's own listeners and the ones it bubbles to never hear the press (a
    // capture-phase listener on an ancestor still does)
    el.addEventListener('click', (e) => {
      if (el.getAttribute('aria-disabled') !== 'true')
        return;
      e.preventDefault();
      e.stopImmediatePropagation();
    }, true);
    return el;
  },
};
