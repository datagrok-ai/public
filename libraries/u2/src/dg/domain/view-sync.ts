/* The platform half of a path-bearing u2 app (Astra A7's second boundary): one history entry per
   user MOVE (typing is not navigating), the address bar as the way back in, and the shell glue that
   ties an app to the view it was docked in. Nothing here knows about domain tables — the second u2
   app with a path reuses it instead of copying `DomainApp`. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import type {Control} from '../../core/component.js';
import type {ReadonlySignal} from '../../core/signals.js';
import type {DomainSession} from '../../sources/session.js';
import {confirmDiscard} from '../../sources/session.js';
import {DomainAddress} from './address.js';

/** What a move is judged by: two states differing only in `query`/`search` are typing. */
export interface ViewState {
  page: string;
  entity: string | null;
  query: string;
  search: string;
  mode: string;
}

export interface ViewSyncHost {
  /** The address the app is at now. */
  readonly path: ReadonlySignal<string>;
  /** The state the current path stands for. */
  state(): ViewState;
  /** The entity value that means "a draft" — a move off it is not a move. */
  readonly draftEntity: string;
  /** Restores an address (the app's `open`). */
  open(address: string): Promise<unknown>;
  /** Whether an address is one of this app's. */
  owns(path: string): boolean;
  /** Whether the shell is showing this app's view. */
  isCurrentView(): boolean;
}

export class ViewSync {
  /** How long the query and the search text must stand still before the entry they opened settles. */
  static readonly SETTLE_MS = 600;

  private _restoring = 0;
  /** Set while a burst of typing is rewriting the entry it opened ({@link _arm}). */
  private _typing: ReturnType<typeof setTimeout> | undefined;

  constructor(private readonly _host: ViewSyncHost, owner: Control) {
    // Back restores the app itself: the platform router re-parses the URL on `popstate` but never
    // asks a JS-hosted view (`handlePath` is not called on the stand — probed 2026-09-15), so the
    // address the entry carries would only move the URL. A router that does call `handlePath`
    // lands on the same `open`, which is idempotent.
    const onPop = () => this._restoreFromUrl();
    window.addEventListener('popstate', onPop);
    owner.own(() => window.removeEventListener('popstate', onPop));
    // One history entry per user move. The shell mirrors `view.path` onto the address bar with
    // `replaceState` (`routing.dart` onViewUrlChanged, a synchronous stream), so a move costs no
    // entry of its own — the app pushes it, and the shell's replace then rewrites that same entry.
    // The effect reads `path` and nothing else: a dependency on the state signals themselves would
    // put it behind `appView`'s mirror in the flush order, too late to push.
    // Typing is not navigating: five keystrokes in the search box were five entries and five
    // Backs. A move that only changed the text rewrites the entry in place and arms a timer; the
    // entry is pushed once the text has settled (or when the next real move needs one).
    let state = _host.state();
    owner.effect(() => {
      const path = _host.path.value;
      const previous = state;
      state = _host.state();
      if (this._restoring !== 0 || typeof history === 'undefined' || !this._moved(previous, state))
        return;
      if (ViewSync._typed(previous, state)) {
        // the first keystroke opens ONE entry for the search — Back still reaches the list as it
        // was before it — and every keystroke after it rewrites that entry in place
        if (this._typing === undefined)
          history.pushState(null, '', path);
        else
          history.replaceState(null, '', path);
        this._arm();
        return;
      }
      this._disarm();
      history.pushState(null, '', path);
    });
    owner.own(() => this._disarm());
  }

  /** Nonzero while a move is being restored from the address bar: nothing is pushed. */
  get restoring(): boolean {
    return this._restoring !== 0;
  }

  /** Runs `fn` as a RESTORE — the entry exists already, so no history entry is made. */
  async restore<T>(fn: () => Promise<T>): Promise<T> {
    this._restoring++;
    try {
      return await fn();
    }
    finally {
      this._restoring--;
    }
  }

  /** A move that only rewrote the query or the search text: the page the user is on did not
   * change, so it is typing, not navigating. */
  private static _typed(from: ViewState, to: ViewState): boolean {
    return from.page === to.page && from.entity === to.entity && from.mode === to.mode;
  }

  /** A move rather than a rewrite of the page in place: the business key replacing the id the page
   * was addressed by, and the draft that became a row, are the same page. */
  private _moved(from: ViewState, to: ViewState): boolean {
    return from.page !== to.page || from.query !== to.query || from.search !== to.search ||
      from.mode !== to.mode || (from.entity !== to.entity && from.entity !== this._host.draftEntity);
  }

  /** The address bar back into the app — only while the shell shows this app's view, and only for
   * an address under one of its routes: another view's URL is not the app's to read. */
  private _restoreFromUrl(): void {
    if (this._host.isCurrentView() && this._host.owns(location.pathname))
      void this._host.open(`${location.pathname}${location.search}`);
  }

  /** Holds the typing burst open; the entry settles once the text has stood still. */
  private _arm(): void {
    this._disarm();
    this._typing = setTimeout(() => this._typing = undefined, ViewSync.SETTLE_MS);
  }

  private _disarm(): void {
    if (this._typing !== undefined)
      clearTimeout(this._typing);
    this._typing = undefined;
  }
}

export interface ViewMountOptions {
  view: DG.ViewBase;
  /** The base the app was built with, and whether the caller pinned it (`path` was given). */
  base: string;
  pinned: boolean;
  /** The app's base now, read after every rebase. */
  baseOf(): string;
  /** Re-points the app at the route the shell mounted its view at. */
  rebase(prefix: string): void;
  open(address: string): Promise<unknown>;
  /** The query keys whose presence makes a cold URL worth replaying. */
  params: readonly string[];
  /** The session whose pending changes gate the ✕ and the unload. */
  session: DomainSession;
  /** The app's own scope — every subscription is unsubscribed through it. */
  own(dispose: () => void): void;
}

/** Ties an app to the view the shell docked it in: find-or-activate, the route the shell really
 * mounted it at, a cold deep link replayed once, the pane's ✕ gate and the browser's unload gate.
 * Everything `DomainTable.app()` did around `appView`, in one place. */
export function mountView(options: ViewMountOptions): void {
  const {view, base, pinned, session} = options;
  // The shell mounts a package app at its own route and prepends it to every path the view
  // reports (`View.path` = the app call's prefix + the view's own): the app rebases onto that
  // route once the view is docked, so a zero-code `table.app()` lives at `/apps/<App>` and its
  // deep links are the shell's. `/domains/<schema>/<table>` is what a view outside an app keeps.
  const mounted = (): string => {
    const full = view.path ?? '';
    const at = full.indexOf('?');
    const here = at < 0 ? full : full.slice(0, at);
    const own = options.baseOf();
    if (!pinned && here.endsWith(own) && here.length > own.length)
      options.rebase(here.slice(0, here.length - own.length));
    return options.baseOf();
  };
  // both routes: the one the shell mounted the view at, and the address a view outside an app
  // keeps — a `/domains/…` link must still reach an app that has rebased onto `/apps/…`
  view.acceptsPath = (p) => {
    const here = p.toLowerCase();
    return DomainAddress.under(here, mounted().toLowerCase()) || DomainAddress.under(here, base.toLowerCase());
  };
  // the router has updated the address bar before it calls the handler (`ViewBase.path`), and
  // hands over the path alone — the row `/domains/…` carries as a segment is in it
  view.handlePath = (p) => {
    mounted();
    void options.open(`${p}${location.search}`);
  };
  // A cold deep link (`/apps/Stockroom?entity=…`) reaches the app func, never a path handler —
  // and the func is handed the path under the app root, never the query. The app opens the
  // address bar itself, from the snapshot taken here: docking the view rewrites the URL to the
  // view's own path (`routing.dart` setViewPath) before any event of ours runs. Only the app's
  // own address counts — opened from the tree, it is another view's URL.
  const from = {pathname: location.pathname.toLowerCase(), search: location.search};
  let replayed = false;
  const replay = (): void => {
    const at = mounted().toLowerCase();
    // only the app's own parameters: a URL carrying nothing but the platform's (`browse=`) has
    // no page to restore, and opening it would drop the query the app was built with
    const deep = new URLSearchParams(from.search);
    if (replayed || !options.params.some((key) => deep.has(key)) || !DomainAddress.under(from.pathname, at))
      return;
    replayed = true;
    void options.open(from.search);
  };
  const added = grok.events.onViewAdded.subscribe((v) => {
    if (v.dart !== view.dart)
      return;
    added.unsubscribe();
    replay();
  });
  // The view's path is the shell's only once it has docked, which is not guaranteed to be before
  // onViewAdded: the route is derived again every time the view becomes current, so a rebase that
  // lost the race still lands and the snapshot is still replayed (once).
  const current = grok.events.onCurrentViewChanged.subscribe((e) => {
    if (e.args?.current?.dart === view.dart)
      replay();
  });
  // the pane's ✕: cancelled here, then closed for real once the user has decided
  const removing = grok.events.onViewRemoving.subscribe((e) => {
    // saving too: the rows read clean for the whole write-back, and closing through it would
    // drop the batch's own re-read
    if (e.args.view.dart !== view.dart || !(session.isDirty.peek() || session.isSaving.peek()))
      return;
    e.preventDefault();
    void confirmDiscard(session, {action: 'close the view'}).then((ok) => ok && view.close());
  });
  options.own(() => {
    added.unsubscribe();
    current.unsubscribe();
    removing.unsubscribe();
  });
}
