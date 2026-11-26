/* What makes a stage save as the platform's domain apps do, with domain storage on (store.ts): Ctrl+S saves the
   view's plan or study, and closing the view or the browser tab with changes not saved asks first (Save / Discard /
   Cancel), as EMS's own apps ask. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {ViewHost} from '@datagrok-libraries/u2/dg';

import {okNamed} from './ui';

/** How a study table is saved, by the Monitor view it is followed in: Analyze, which changes the same table (the
 *  patients' exposure, the effect it keeps, exclusions), saves it the same way. */
type Saver = () => Promise<void>;
export const tableSavers = new WeakMap<DG.DataFrame, Saver>();

/** A save refused because the row is shared with the user to read only: says so, and offers to save their version as
 *  a copy of their own (`copy`, which says where it went). Nothing else is shown: no error, no stack. */
export function offerCopy(copy: () => Promise<unknown>): void {
  const d = ui.dialog('Shared with you to read')
    .add(ui.divText('You can open it but not change it. Save your changes as your own copy?'))
    .onOK(() => void copy().catch((e) => grok.shell.error(`The copy could not be saved: ${e?.message ?? e}`)));
  okNamed(d, 'Save as my copy');
  d.show();
}

/** The gates not yet released (each lets go when its view closes, however it closes): for the tests. */
export const armedGates = new Set<string>();

/** Ctrl+S (Cmd+S) and the unsaved-changes gate for the host's view, until it closes (the host's end releases it, and
 *  so do the view's removal and a tab closed after it). `save` resolves once saved (or rejects, having said why);
 *  `dirty` says whether there is anything to lose; `unsaved`, whether it was never saved at all (said so: it has no
 *  "changes" to speak of). */
export function saveGate(host: ViewHost, what: string, save: () => Promise<unknown>, dirty: () => boolean,
  unsaved: () => boolean = () => false): void {
  const view = host.view;
  const id = `${what} ${Math.random().toString(36).slice(2)}`;
  let leaving = false;
  const onKey = (e: KeyboardEvent) => {
    // (the shell hands out a new wrapper of the current view each time: compared by the view behind it)
    if ((e.ctrlKey || e.metaKey) && !e.altKey && e.key.toLowerCase() === 's' && grok.shell.v?.dart === view.dart) {
      e.preventDefault();
      e.stopPropagation();
      void save().catch(() => {});
    }
  };
  // (the browser asks "Leave site?" only while the view is open and something in it is not saved: a view closed
  // without its removal being seen lets go here)
  const open = () => [...grok.shell.views].some((v) => v.dart === view.dart);
  const onUnload = (e: BeforeUnloadEvent) => {
    if (!open()) {
      release();
      return;
    }
    if (!dirty())
      return;
    // (said in the console, so a "Leave site?" can be told from the platform's own)
    console.warn(`PMX Studio: unsaved ${what} in "${view.name}"`);
    e.preventDefault();
    e.returnValue = '';
  };
  const removing = grok.events.onViewRemoving.subscribe((data: DG.EventData) => {
    if (data.args?.view?.dart !== view.dart || leaving || !dirty())
      return;
    data.preventDefault();
    // (let go at once: a removal the shell reports for another wrapper of the view would leave the tab's "Leave
    // site?" armed after a Discard)
    const close = () => {
      leaving = true;
      release();
      view.close();
    };
    const d = ui.dialog('Unsaved changes').add(ui.divText(unsaved() ? `This ${what} is not saved yet.` :
      `This ${what} has changes that are not saved.`))
      .onOK(() => void save().then(close, () => {}))
      .addButton('Discard', () => {
        d.close();
        close();
      });
    okNamed(d, 'Save');
    d.show();
  });
  const release = () => {
    armedGates.delete(id);
    document.removeEventListener('keydown', onKey, true);
    window.removeEventListener('beforeunload', onUnload);
    removing.unsubscribe();
    removed.unsubscribe();
  };
  const removed = grok.events.onViewRemoved.subscribe((v: DG.ViewBase) => {
    if (v?.dart === view.dart)
      release();
  });
  document.addEventListener('keydown', onKey, true);
  window.addEventListener('beforeunload', onUnload);
  armedGates.add(id);
  host.own(release);
}
