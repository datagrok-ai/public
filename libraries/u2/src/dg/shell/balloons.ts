/* On a platform page u2's notifications join the platform's balloon stack — one stack, one look —
   through `grok.shell.info / warning / error`; on any other page (a gallery, a test, a standalone
   host) the platform is not there and u2 keeps its own toasts. Installed at dg-index import. */
import * as grok from 'datagrok-api/grok';
import {notify} from '../../components/display/notify.js';
import type {NotifyHandle, NotifyOptions, NotifyType} from '../../components/display/notify.js';

/** Where the platform stacks its balloons (`balloon.dart`). */
const CONTAINER = '.d4-balloon-container';

export class PlatformBalloons {
  /** Whether the page runs the platform: its Dart side registers `grok_Balloon`. */
  static get available(): boolean {
    return typeof (globalThis as {grok_Balloon?: unknown}).grok_Balloon === 'function';
  }

  /** u2's three types are the platform's; `autoHide` goes over as u2 resolved it (an error stays),
   * `timeout` in seconds and `copyText` as it is. `oneTimeKey` and `singleKey` stay u2's.
   * The platform closes a balloon on any click, so `closeOnClick: false` is lost. */
  static show(content: string | HTMLElement, type: NotifyType,
    options: NotifyOptions & {autoHide: boolean, timeout: number}, closed: () => void): NotifyHandle | undefined {
    if (!PlatformBalloons.available)
      return undefined;
    const before = document.querySelector(CONTAINER)?.lastElementChild ?? null;
    // the platform reads a string as markup (an error's as text); u2's string is always text
    const shown = typeof content === 'string' && type !== 'error' ? PlatformBalloons._text(content) : content;
    const balloon = {autoHide: options.autoHide, copyText: options.copyText,
      timeout: Math.max(1, Math.round(options.timeout / 1000))};
    if (type === 'error')
      grok.shell.error(shown, balloon);
    else if (type === 'warning')
      grok.shell.warning(shown, balloon);
    else
      grok.shell.info(shown, balloon);
    const container = document.querySelector(CONTAINER);
    const root = container?.lastElementChild;
    if (container === null || !(root instanceof HTMLElement) || root === before) {
      closed();
      return {root: document.createElement('div'), close: () => {}};
    }
    // the platform closes it on its own — its timer, a click, its ✕
    const observer = new MutationObserver(() => {
      if (root.isConnected)
        return;
      observer.disconnect();
      closed();
    });
    observer.observe(container, {childList: true});
    return {root, close: () => root.remove()};
  }

  private static _text(text: string): HTMLElement {
    const el = document.createElement('span');
    el.textContent = text;
    return el;
  }
}

notify.presenter = (content, type, options, closed) => PlatformBalloons.show(content, type, options, closed);
