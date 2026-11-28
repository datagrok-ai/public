// What runs only while a page of the app is open (the history hook, the panels' watch, the session's listeners):
// started with the first page, stopped once the last one closes, so nothing of the app stays on the platform's
// globals while the user is elsewhere. No imports of the app's modules: stages.ts, handlers.ts and colors.ts use it.
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** A page of the app: a stage's view (pmxStage) or Compare (pmxAlone); see ui.ts stageOf and stages.ts STANDS_ALONE. */
export const isAppView = (v: DG.ViewBase | null | undefined): boolean =>
  !!v?.root?.dataset.pmxStage || !!v?.root?.dataset.pmxAlone;

type Stop = () => void;
const starters: (() => Stop)[] = [];
/** The stops of what runs now; null while the app is closed. */
let running: Stop[] | null = null;

/** `start` runs with every opening of the app (now, when it is open) and the stop it answers as the last page
 *  closes. */
export function whileAppOpen(start: () => Stop): void {
  starters.push(start);
  running?.push(start());
}

/** Whether the app is open (what whileAppOpen registered runs). */
export const appIsOpen = () => running !== null;

/** The app opened (its function, a file handler, a stage shown): starts what runs while it is open, once. */
export function appOpened(): void {
  if (running)
    return;
  running = [];
  // (a page of the app removed and none left, a tick later: a move adds the new page before the old one closes)
  const removed = grok.events.onViewRemoved.subscribe((gone) => {
    if (!isAppView(gone))
      return;
    setTimeout(() => {
      if (running && ![...grok.shell.views].some((v) => isAppView(v) && v.dart !== gone.dart))
        appClosed();
    });
  });
  running.push(() => removed.unsubscribe());
  for (const s of starters)
    running.push(s());
}

/** Stops what runs while the app is open (the last page closed). */
export function appClosed(): void {
  const stops = running;
  running = null;
  // (the last started stops first)
  for (const s of (stops ?? []).reverse()) {
    try {
      s();
    } catch (e) {
      console.warn('PMX Studio: a listener could not be stopped', e);
    }
  }
}
