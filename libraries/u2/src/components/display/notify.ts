/* Notifications — the u2 counterpart of the platform's `Balloon` (d4/widgets/balloon): a
   top-right stack of per-type balloons. Ported semantics: info auto-hides, warning and error
   stay until closed; hovering pauses the countdown; `oneTimeKey` persists across sessions,
   `singleKey` dedupes concurrent balloons; copy and close reveal on hover; a click anywhere
   on the balloon closes it. Deliberately not ported: `appendHtml` for string content — a
   string renders as text, never as markup. */
import {Overlay} from '../../core/overlay.js';
import {div} from '../../core/elements.js';
import {icon} from './icon.js';
import {iconButton} from '../actions/buttons.js';

export type NotifyType = 'info' | 'warning' | 'error';

export interface NotifyOptions {
  /** Close after {@link timeout}; info balloons by default, warnings and errors stay until closed. */
  autoHide?: boolean;
  /** Auto-hide delay in milliseconds; 5000 by default. */
  timeout?: number;
  /** Shown once ever — remembered in localStorage under `u2-notify-<key>`. */
  oneTimeKey?: string;
  /** At most one balloon with this key at a time; released when it closes. */
  singleKey?: string;
  /** Appended to the balloon text when the user copies it. */
  copyText?: string;
  /** A click anywhere on the balloon closes it; on by default. */
  closeOnClick?: boolean;
}

export interface NotifyHandle {
  root: HTMLElement;
  close(): void;
}

const TYPE_ICONS: Record<NotifyType, string> = {
  info: 'info-circle', warning: 'exclamation-triangle', error: 'times-circle',
};

class Notify {
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

  static info(content: string | HTMLElement, options: NotifyOptions = {}): NotifyHandle | null {
    return Notify.show(content, 'info', options);
  }

  static warning(content: string | HTMLElement, options: NotifyOptions = {}): NotifyHandle | null {
    return Notify.show(content, 'warning', options);
  }

  static error(content: string | HTMLElement, options: NotifyOptions = {}): NotifyHandle | null {
    return Notify.show(content, 'error', options);
  }

  /** Returns null when the balloon is suppressed by `oneTimeKey` or `singleKey`. */
  static show(content: string | HTMLElement, type: NotifyType = 'info',
    options: NotifyOptions = {}): NotifyHandle | null {
    const {oneTimeKey, singleKey} = options;
    if (singleKey !== undefined && Notify._singleKeys.has(singleKey))
      return null;
    if (oneTimeKey !== undefined && !Notify._claimOneTime(oneTimeKey))
      return null;
    if (singleKey !== undefined)
      Notify._singleKeys.add(singleKey);

    const body = div([], 'u2-notify-content');
    if (typeof content === 'string')
      body.textContent = content;
    else
      body.append(content);

    const balloon = div([icon(TYPE_ICONS[type], {cls: 'u2-notify-type-icon'}), body],
      `u2-notify u2-notify-${type}`);
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
    };
    Notify._open.add(close);

    const copy = iconButton('copy', () => {
      const text = options.copyText === undefined ?
        balloon.textContent ?? '' : `${balloon.textContent}\n\n${options.copyText}`;
      try {
        void navigator.clipboard?.writeText(text);
      } catch {
        // no clipboard in this host — the copy affordance is best-effort
      }
    }, {tooltip: 'Copy to clipboard'});
    copy.addEventListener('click', (e) => e.stopPropagation());
    balloon.append(div([copy, iconButton('times', close, {tooltip: 'Close'})], 'u2-notify-icons'));

    if (options.closeOnClick ?? true)
      balloon.addEventListener('click', close);

    let timer = 0;
    if (options.autoHide ?? type === 'info') {
      const timeout = options.timeout ?? 5000;
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
