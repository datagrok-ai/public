/* Notifications — the u2 counterpart of the platform's `Balloon` (d4/widgets/balloon): a
   top-right stack of per-type balloons; on a platform page the platform layer routes them into the
   platform's own stack instead (`presenter`, src/dg). An info or a warning auto-hides, an error
   stays until the user closes it; hovering pauses the countdown; `oneTimeKey` persists across sessions,
   `singleKey` dedupes concurrent balloons; copy and close reveal on hover; a click anywhere
   on the balloon closes it. Deliberately not ported: `appendHtml` for string content — a
   string renders as text, never as markup. */
import {Overlay} from '../../core/overlay.js';
import {Scope} from '../../core/scope.js';
import {batch} from '../../core/signals.js';
import {div} from '../../core/elements.js';
import {icon} from './icon.js';
import {iconButton} from '../actions/buttons.js';

export type NotifyType = 'info' | 'warning' | 'error';

export interface NotifyOptions {
  /** Close after {@link timeout}; on for an info or a warning, off for an error, which stays until closed. */
  autoHide?: boolean;
  /** Auto-hide delay in milliseconds; 10000 by default. */
  timeout?: number;
  /** Shown once ever — remembered in localStorage under `u2-notify-<key>`. */
  oneTimeKey?: string;
  /** At most one balloon with this key at a time; released when it closes. */
  singleKey?: string;
  /** Appended to the balloon text when the user copies it. */
  copyText?: string;
  /** A click anywhere on the balloon closes it; on by default. */
  closeOnClick?: boolean;
  /** Runs once the balloon is gone — its timer, a click, ✕, `close()` — owned by whoever showed it. */
  onClose?: () => void;
}

/** Text; an element; or a function building the balloon's content in a scope of its own, ended with the balloon:
 * an Undo link raised from a row's handler keeps working after the row goes and leaves nothing behind. */
export type NotifyContent = string | HTMLElement | (() => HTMLElement);

export interface NotifyHandle {
  root: HTMLElement;
  close(): void;
}

/** Shows a balloon somewhere else than u2's own stack — the platform's, inside the platform
 * (`src/dg`) — and calls `closed` once it is gone; undefined leaves it to u2's stack. Runs after the
 * `oneTimeKey` and `singleKey` checks, which u2 keeps; `options` come with u2's `autoHide` and `timeout` applied. */
export type NotifyPresenter = (content: string | HTMLElement, type: NotifyType,
  options: NotifyOptions & {autoHide: boolean, timeout: number}, closed: () => void) => NotifyHandle | undefined;

const TYPE_ICONS: Record<NotifyType, string> = {
  info: 'info-circle', warning: 'exclamation-triangle', error: 'times-circle',
};

class Notify {
  /** Set by the platform layer: one balloon stack, the platform's, on a platform page. */
  static presenter: NotifyPresenter | undefined;
  private static _container: HTMLElement | undefined;
  private static readonly _singleKeys = new Set<string>();
  private static readonly _open = new Set<() => void>();

  static get host(): HTMLElement {
    if (!Notify._container || !Notify._container.isConnected) {
      const el = div([], 'u2-notify-container');
      el.dataset.u2 = 'notify';
      Overlay.host.append(el);
      Notify._container = el;
    }
    return Notify._container;
  }

  static info(content: NotifyContent, options: NotifyOptions = {}): NotifyHandle | null {
    return Notify.show(content, 'info', options);
  }

  static warning(content: NotifyContent, options: NotifyOptions = {}): NotifyHandle | null {
    return Notify.show(content, 'warning', options);
  }

  static error(content: NotifyContent, options: NotifyOptions = {}): NotifyHandle | null {
    return Notify.show(content, 'error', options);
  }

  /** Returns null when the balloon is suppressed by `oneTimeKey` or `singleKey`: then nothing is built and
   * `onClose` never runs. */
  static show(content: NotifyContent, type: NotifyType = 'info',
    options: NotifyOptions = {}): NotifyHandle | null {
    const {oneTimeKey, singleKey} = options;
    if (singleKey !== undefined && Notify._singleKeys.has(singleKey))
      return null;
    if (oneTimeKey !== undefined && !Notify._claimOneTime(oneTimeKey))
      return null;
    if (singleKey !== undefined)
      Notify._singleKeys.add(singleKey);

    const shower = Scope.ambient;
    const own = new Scope(undefined, {label: 'notify', lazy: true});
    const built = typeof content === 'function' ? own.build(content) : content;
    let ended = false;
    const end = () => {
      if (ended)
        return;
      ended = true;
      own.dispose();
      const onClose = options.onClose;
      if (onClose !== undefined)
        batch(() => Scope.runWith(shower?.isDisposed ? undefined : shower, onClose));
    };

    let routedClose: (() => void) | undefined;
    let gone = false;
    const released = () => {
      gone = true;
      if (routedClose !== undefined)
        Notify._open.delete(routedClose);
      if (singleKey !== undefined)
        Notify._singleKeys.delete(singleKey);
      end();
    };
    const autoHide = options.autoHide ?? type !== 'error';
    const timeout = options.timeout ?? 10000;
    const routed = Notify.presenter?.(built, type, {...options, autoHide, timeout}, released);
    if (routed !== undefined) {
      routedClose = routed.close;
      // a presenter that could not show it released it already
      if (!gone)
        Notify._open.add(routed.close);
      return routed;
    }

    const body = div([], 'u2-notify-content');
    if (typeof built === 'string')
      body.textContent = built;
    else
      body.append(built);

    // the balloon outlives whoever raised it (a handler, an effect run): its controls have no owner
    const balloon = Scope.runWith(undefined, () => div([icon(TYPE_ICONS[type], {cls: 'u2-notify-type-icon'}), body],
      `u2-notify u2-notify-${type}`));
    balloon.setAttribute('role', type === 'info' ? 'status' : 'alert');

    let closed = false;
    const close = () => {
      if (closed)
        return;
      closed = true;
      window.clearTimeout(timer);
      balloon.remove();
      Notify._open.delete(close);
      if (singleKey !== undefined)
        Notify._singleKeys.delete(singleKey);
      end();
    };
    Notify._open.add(close);

    const copy = Scope.runWith(undefined, () => iconButton('copy', () => {
      const text = options.copyText === undefined ?
        balloon.textContent ?? '' : `${balloon.textContent}\n\n${options.copyText}`;
      try {
        void navigator.clipboard?.writeText(text);
      } catch {
        // no clipboard in this host — the copy affordance is best-effort
      }
    }, {tooltip: 'Copy to clipboard'}));
    copy.addEventListener('click', (e) => e.stopPropagation());
    balloon.append(div([copy, Scope.runWith(undefined, () => iconButton('times', close, {tooltip: 'Close'}))],
      'u2-notify-icons'));

    if (options.closeOnClick ?? true)
      balloon.addEventListener('click', close);

    let timer = 0;
    if (autoHide) {
      const bar = div([], 'u2-notify-timer');
      bar.style.animationDuration = `${timeout}ms`;
      balloon.append(bar);
      let startedAt = Date.now();
      let remaining = timeout;
      timer = window.setTimeout(close, remaining);
      balloon.addEventListener('pointerenter', () => {
        window.clearTimeout(timer);
        remaining -= Date.now() - startedAt;
      });
      balloon.addEventListener('pointerleave', () => {
        startedAt = Date.now();
        timer = window.setTimeout(close, Math.max(remaining, 0));
      });
    }

    Notify.host.append(balloon);
    return {root: balloon, close};
  }

  static closeAll(): void {
    for (const close of [...Notify._open])
      close();
  }

  /** localStorage is unavailable in some hosts (headless tests, sandboxed frames) — then every
   * show is a first time. */
  private static _claimOneTime(key: string): boolean {
    try {
      const storageKey = `u2-notify-${key}`;
      if (window.localStorage.getItem(storageKey) !== null)
        return false;
      window.localStorage.setItem(storageKey, new Date().toISOString());
      return true;
    } catch {
      return true;
    }
  }
}

/** `notify.info('Saved')`, `notify.error(...)` — the u2 toast service. */
export const notify = Notify;
