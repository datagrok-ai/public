/* What makes a stage save as the platform's domain apps do, with domain storage on (store.ts): Ctrl+S saves the
   view's plan or study, and closing the view or the browser tab with changes not saved asks first (Save / Discard /
   Cancel), as EMS's own apps ask. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {ViewHost} from '@datagrok-libraries/u2/dg';
import {Emitter} from '@datagrok-libraries/u2';

import {dialogLine, okNamed, modal, showCentred} from './ui';
import {PLAN_TAG, parsePlan} from './plan-file';
import {ROW_TAG, inTrash, markTableFromFile, rowOfTable} from './store';
import {studyNameOf} from './study-file';
import {StudyMeta, emptyMeta, metaOf} from './study-meta';
import {noBidi} from './arm-name';
import {forgetForm} from './form-memory';

/** How a study table is saved, by the Monitor view it is followed in: Analyze, which changes the same table (the
 *  patients' exposure, the effect it keeps, exclusions), saves it the same way. */
type Saver = () => Promise<void>;
export const tableSavers = new WeakMap<DG.DataFrame, Saver>();

/** A save refused because the row is shared with the user to read only: says so, and offers to save their version as
 *  a copy of their own (`copy`, which says where it went). Nothing else is shown: no error, no stack. */
export function offerCopy(copy: () => Promise<unknown>, refused?: any): void {
  // (one deleted since it was opened, from another tab or by someone else: said so, and saved as a new one, QA 0.1.3
  // round 3 chaos 11)
  // every save of a study or plan in the trash comes here, whichever page or key asked (QA 0.1.3 round 4 chaos N2)
  const gone = refused?.code === 'deleted';
  const what = refused?.kind === 'plan' ? 'plan' : 'study';
  const d = modal('offer-copy', gone ? 'In the trash' : 'Shared with you to read');
  if (!d)
    return;
  d.add(ui.divText(gone ? `This ${what} is in the trash. Save your changes as a new ${what}?` :
    'You can open it but not change it. Save your changes as your own copy?'))
    .onOK(() => void copy().catch((e) => grok.shell.error(`The copy could not be saved: ${e?.message ?? e}`)));
  okNamed(d, gone ? 'Save as new' : 'Save as my copy');
  d.show();
}

/** The gates not yet released (each lets go when its view closes, however it closes): for the tests. */
export const armedGates = new Set<string>();

/** Whether a view (by the platform's view behind it) has changes not saved: one moving to another stage stays open
 *  rather than be closed (stages.ts). */
type Dirty = () => boolean;
export const dirtyViews = new Map<unknown, Dirty>();

/** Asks, for a view (by the platform's view behind it) with changes not saved, Save / Discard / Cancel: true once
 *  saved or discarded (the view may then close without asking again), false on Cancel. */
export const leaveAsks = new Map<unknown, (stay?: boolean, after?: string) => Promise<boolean>>();

/** The tables each view tracks as saved (savedTables), by the platform's view behind it, and whether each changed. */
const tracked = new Map<unknown, {tables: Set<DG.DataFrame>, changed: (df: DG.DataFrame) => boolean}>();
/** The study tables each view shows (holdTable), by the platform's view behind it. */
const held = new Map<unknown, Set<DG.DataFrame>>();
/** The tables the app itself put in the workspace (a study opened from a row, a file, the demo): only these does it
 *  ever take out again. One the user opened outside PMX Studio (a CSV from Files, then picked) stays theirs. */
const openedHere = new WeakSet<object>(); // (by the platform's table behind it: the shell hands out other wrappers)
/** The views whose changes were discarded: what they show is not kept to come back to (stages.ts leaveOthers). */
export const discardedViews = new Set<unknown>();

/** `df` shown by the host's view: put in the workspace when it is not there (the app opened it), and held while the
 *  view is open. The view gone, its tables go too, unless still used (dropUnused): each closed page left its study
 *  in the workspace, and the Study field listed "X", "X (2)", "X (3)" (QA 0.1.3 round 3 technician N10, analyst 14,
 *  chaos 16). `named`: what shows the study again by its name, once a name the platform suffixed is given back. */
export function holdTable(host: ViewHost, df: DG.DataFrame, named?: () => void): void {
  if (!inWorkspace(df)) {
    // (its name as the Study field lists it: no right-to-left marks to show it reversed, QA 0.1.3 round 4 chaos N5)
    if (noBidi(df.name) !== df.name)
      df.name = noBidi(df.name).trim() || 'Study';
    const name = df.name;
    grok.shell.addTable(df);
    openedHere.add(df.dart);
    if (df.name !== name)
      meant.set(df.dart, {df, name, named});
  }
  // (a tick later: a page opening with its study holds it while its content is built, before the host has its view)
  void Promise.resolve().then(() => {
    const view = host.view.dart;
    if (!held.has(view)) {
      held.set(view, new Set());
      host.own(() => {
        const mine = [...held.get(view) ?? []];
        held.delete(view);
        // (its changes read as the view goes: a changed table is never taken out without asking)
        const t = tracked.get(view);
        const unchanged = mine.filter((x) => !t?.changed(x));
        setTimeout(() => dropUnused(unchanged), 2000);
      });
    }
    held.get(view)!.add(df);
  });
}

/** Whether `df` is still used: shown by an open view of the app (held or tracked), or by any table view of the
 *  platform's (Open as table, the user's own grid). */
function used(df: DG.DataFrame, except?: unknown): boolean {
  const open = new Set([...grok.shell.views].map((v) => v.dart));
  open.delete(except);
  // (by the platform's table behind the wrapper: the shell may hand out another wrapper of the same table)
  const same = (ts: Iterable<DG.DataFrame>) => [...ts].some((x) => x.dart === df.dart);
  return [...held].some(([v, ts]) => open.has(v) && same(ts)) ||
    [...tracked].some(([v, t]) => open.has(v) && same(t.tables)) ||
    [...grok.shell.tableViews].some((v) => v.dataFrame?.id === df.id);
}
/** Whether an open page (or table view) still shows `df`: a study row opened again takes up only such a table, else
 *  reads the row as saved (store.ts tableOfRow). */
export const inUse = (df: DG.DataFrame): boolean => used(df);
const inWorkspace = (df: DG.DataFrame) => grok.shell.tables.some((t) => t.dart === df.dart);

/** The tables the platform put in the workspace as "X (2)" because the page they replace still held its own "X" (two
 *  saved studies of one name, Ctrl+O from one to the other): the name each was meant to have, given back as that
 *  table goes, and what shows it again (QA 0.1.3 r5: the only study open read "X (2)"). */
const meant = new Map<unknown, {df: DG.DataFrame, name: string, named?: () => void}>();

/** Takes `df` out of the workspace when the app put it there and nothing uses it now (see used). */
function takeOut(df: DG.DataFrame, except?: unknown): void {
  if (!openedHere.has(df.dart) || !inWorkspace(df) || used(df, except))
    return;
  // (renamed before the close: the platform's field lists the tables again as one goes, by their names then, and
  // is told its pick once the list has it)
  const freed = [...meant].filter(([, m]) => m.df.dart !== df.dart && inWorkspace(m.df) &&
    !grok.shell.tables.some((t) => t.dart !== df.dart && t.name === m.name));
  freed.forEach(([k, m]) => {
    meant.delete(k);
    m.df.name = m.name;
  });
  meant.delete(df.dart);
  grok.shell.closeTable(df);
  freed.forEach(([, m]) => m.named?.());
}

/** A page let go before it was shown (Cancel on the unsaved prompt it asked): the tables it put in the workspace go
 *  at once, none left behind as "X (2)" in the Study field (QA 0.1.3 round 4 technician R4-5). */
export function letGoTables(view: DG.ViewBase): void {
  for (const df of held.get(view.dart) ?? [])
    takeOut(df, view.dart);
}

/** The Study field's list (Monitor's and Analyze's): no study in the trash, and each saved study once, the one shown
 *  first (QA 0.1.3 round 4 chaos N2, N5, 16: a deleted study still listed and opened; "X" beside its cleaned copy). */
export function studyPickList(names: string[], current: string | null | undefined): string[] {
  const tables = grok.shell.tables;
  const byName = (n: string | null | undefined) => tables.find((t) => t.name === n);
  const seen = new Set<string>();
  const shown = rowOfTable(byName(current))?.id;
  if (shown)
    seen.add(shown);
  return names.filter((n) => {
    const t = byName(n);
    if (n === current || !t)
      return true;
    if (inTrash(t))
      return false;
    const id = rowOfTable(t)?.id;
    if (!id)
      return true;
    if (seen.has(id))
      return false;
    seen.add(id);
    return true;
  });
}

/** Each of `tables` taken out of the workspace unless used (a page opened since may have taken it up again). */
const dropUnused = (tables: DG.DataFrame[]) => tables.forEach((df) => takeOut(df));

/** Discard: the tables `view` changed are taken out (the app's own, unused elsewhere), so the study opens again as
 *  saved (QA 0.1.3 round 3 chaos 1: a discarded study stayed in the workspace and its row reopened it, changes and
 *  all). A table of the user's own, opened outside the app, keeps its values: the app never closes it. */
function dropChanged(view: DG.ViewBase): void {
  const mine = tracked.get(view.dart);
  for (const df of mine?.tables ?? []) {
    if (mine!.changed(df)) {
      // (nor the form values kept for it: the study opens again as saved, QA 0.1.3 round 4 chaos N1)
      forgetForm(df);
      takeOut(df, view.dart);
    }
  }
}

/** The class of a view with changes not saved (its header says so). */
export const DIRTY_CLASS = 'pmx-dirty';

/** What reports a change: a table's events (savedTables), a form's inputs. */
export type Changes = {subscribe: (next: () => void) => {unsubscribe(): void}};

/** The tables of a view as last saved, and whether each changed since, without serialising one on every check. A
 *  counter per table follows the platform's events (its values, rows, columns, names, metadata) while the host's
 *  view is open; `content` (the table's text, say) is taken again only once the counter moved since the last
 *  comparison, so a change undone back to the saved state reads as saved again and a check with nothing changed
 *  costs nothing. `tags`: the small parts a table is saved with, compared on every check. `changes` fires on each
 *  change of a table followed and on each `set` (a save, wherever it was asked). */
export function savedTables(host: ViewHost, content: (df: DG.DataFrame) => string,
  tags: (df: DG.DataFrame) => string = () => ''): {set(df: DG.DataFrame): void, has(df: DG.DataFrame): boolean,
    changed(df: DG.DataFrame): boolean, changes: Changes} {
  const changes = new Emitter();
  const versions = new WeakMap<DG.DataFrame, number>();
  const saved = new WeakMap<DG.DataFrame, {content: string, tags: string, same: number, differs: number}>();
  const version = (df: DG.DataFrame): number => {
    let v = versions.get(df);
    if (v === undefined) {
      versions.set(df, v = 0);
      const bump = () => {
        versions.set(df, versions.get(df)! + 1);
        changes.fire();
      };
      for (const e of [df.onDataChanged, df.onColumnNameChanged, df.onMetadataChanged])
        host.subscribe(e, bump);
    }
    return v;
  };
  const self = {
    changes,
    has: (df: DG.DataFrame) => saved.has(df),
    set: (df: DG.DataFrame) => {
      saved.set(df, {content: content(df), tags: tags(df), same: version(df), differs: -1});
      // (a tick later: a page opening with its study sets it while its content is built, before the host has its view)
      void Promise.resolve().then(() => {
        const view = host.view.dart;
        if (!tracked.has(view)) {
          tracked.set(view, {tables: new Set(), changed: (x) => self.changed(x)});
          // (let go a moment after the view: holdTable reads its changes as the view goes)
          host.own(() => setTimeout(() => tracked.delete(view)));
        }
        tracked.get(view)!.tables.add(df);
      });
      changes.fire();
    },
    changed: (df: DG.DataFrame): boolean => {
      const s = saved.get(df);
      if (!s)
        return false;
      if (tags(df) !== s.tags)
        return true;
      const v = version(df);
      if (v === s.same || v === s.differs)
        return v === s.differs;
      const same = content(df) === s.content;
      if (same)
        s.same = v;
      else
        s.differs = v;
      return !same;
    },
  };
  return self;
}

/** Whether a table is one the app opened and saves itself: a study (study-file.ts STUDY_TAG), one with its plan
 *  (plan-file.ts PLAN_TAG), a saved row (store.ts ROW_TAG) or one it opened from its own files (APP_TABLE_TAG: the
 *  demo's mouse PK, which held the platform's "Leave site?" on an untouched demo, QA 0.1.3 r4). Anything else, the
 *  platform's prompt stays for. */
export const appTable = (df: DG.DataFrame): boolean =>
  ['pmx.study', 'pmx.plan', 'pmx.row', 'pmx.app'].some((t) => !!df.getTag(t));
/** The tag of a table the app opened from its own files (nothing of the user's to lose in it): appTable. */
export const APP_TABLE_TAG = 'pmx.app';

/** Ctrl+S (Cmd+S) and the unsaved-changes gate for the host's view, until it closes (the host's end releases it, and
 *  so do the view's removal and a tab closed after it). `save` resolves once saved (or rejects, having said why);
 *  `dirty` says whether there is anything to lose; `unsaved`, whether it was never saved at all (said so: it has no
 *  "changes" to speak of); `changes`, what reports a change (the marker is checked a moment after each). */
export function saveGate(host: ViewHost, what: string, save: () => Promise<unknown>, dirty: () => boolean,
  unsaved: () => boolean = () => false, changes?: Changes): void {
  const view = host.view;
  const id = `${what} ${Math.random().toString(36).slice(2)}`;
  let leaving = false;
  const onKey = (e: KeyboardEvent) => {
    // (the shell hands out a new wrapper of the current view each time: compared by the view behind it)
    if ((e.ctrlKey || e.metaKey) && !e.altKey && e.key.toLowerCase() === 's' && grok.shell.v?.dart === view.dart) {
      e.preventDefault();
      e.stopPropagation();
      // (from a field being typed in: its value committed first, as leaving the field does, then saved; QA 0.1.3
      // designer D16, Ctrl+S in a field did nothing)
      const field = document.activeElement;
      if (field instanceof HTMLElement && view.root?.contains(field) && field.matches('input, textarea, select')) {
        field.dispatchEvent(new Event('change', {bubbles: true}));
        field.blur();
        setTimeout(() => {
          void save().then(mark, mark);
          field.focus();
        }, 50);
        return;
      }
      void save().then(mark, mark);
    }
  };
  // the marker of changes not saved ("Unsaved changes" after the header's study ID; the view's class pmx-dirty):
  // checked a moment after each change the view reports (`changes`: its tables' events, its form's inputs, a save
  // from a button) and after each save through here; nothing is checked while nothing changes
  let timer: ReturnType<typeof setTimeout> | undefined;
  const mark = () => {
    const on = open() && dirty();
    view.root?.classList.toggle(DIRTY_CLASS, on);
    let tag = view.root?.querySelector<HTMLElement>(`.${DIRTY_CLASS}-mark`);
    const id = view.root?.querySelector('.pmx-head-id');
    if (on && !tag && id) {
      tag = document.createElement('span');
      tag.className = `${DIRTY_CLASS}-mark`;
      tag.textContent = 'Unsaved changes';
      Object.assign(tag.style, {flex: 'none', fontSize: '12px', fontStyle: 'italic', color: 'var(--dg-orange, #c77c02)'});
      id.after(tag);
    }
    // (a study or plan never saved says so already, "Not saved yet": not "Unsaved changes" beside it, QA 0.1.3 round 4
    // designer U-d)
    if (tag)
      tag.hidden = !on || !!id?.textContent?.startsWith('Not saved yet');
    // the platform's own "Leave site?" follows ours: the study's tables, opened or written by the view, leave the
    // workspace "changed" though nothing of the study is unsaved (its documented clear, while this view is current).
    // The platform keeps one flag for the whole workspace (no table's own): cleared only while every open table is
    // the app's, else a table of the user's own keeps its prompt
    if (!on && open() && grok.shell.v?.dart === view.dart && grok.shell.tables.every(appTable))
      grok.shell.clearDirtyFlag();
  };
  // (a burst of changes, a pasted column: checked once, after the last)
  const soon = () => {
    clearTimeout(timer);
    timer = setTimeout(mark, 300);
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
  // Save / Discard / Cancel: true once saved or discarded (discarded, the gate lets go: the view closes without
  // asking again), false on Cancel or a save that failed (it said why)
  // `stay`: asked for another study picked in the view's Study field, the view staying open: Discard then only
  // drops the changes, and the gate stays
  // `after`: what saving comes before, said in the question ("before following it": designer wording 9)
  const ask = (stay = false, after = ''): Promise<boolean> => new Promise((resolve) => {
    const d = modal(`unsaved:${id}`, 'Unsaved changes');
    if (!d)
      return resolve(false);
    // the study or plan named as its header names it (QA 0.1.3 round 3: technician N16, chaos 25, analyst W10)
    const name = view.root?.querySelector('.pmx-head-name')?.textContent?.trim();
    const it = name ? `“${name}”` : `this ${what}`;
    // (OK resolves once its save is done; Discard at once; Cancel, ✕ and Esc as Cancel)
    let answer = false;
    let saving = false;
    d.add(ui.divText(unsaved() ? `${it[0].toUpperCase()}${it.slice(1)} is not saved yet. Save it?` :
      `Save your changes to ${it}${after ? ` ${after}` : ''}?`))
      .onOK(() => {
        saving = true;
        void save().then(() => resolve(true), () => {
          mark();
          resolve(false);
        });
      })
      .addButton('Discard', () => {
        answer = true;
        if (stay) {
          dropChanged(view);
          d.close();
          return;
        }
        leaving = true;
        // (what was discarded is not kept anywhere: the saved version is what opens next)
        discardedViews.add(view.dart);
        dropChanged(view);
        // (let go at once: a removal the shell reports for another wrapper of the view would leave the tab's "Leave
        // site?" armed after a Discard)
        release();
        d.close();
      });
    okNamed(d, 'Save');
    d.onClose.subscribe(() => answer ? resolve(true) : !saving && resolve(false));
    // in the middle of the window, whatever asked (Back's opened by the toolbox), the buttons in the platform's one
    // order, Cancel · Discard · Save, as every dialog of the app (the reorder held only some of the time: two orders,
    // QA 0.1.3 round 4 designer U2)
    showCentred(d);
  });
  const removing = grok.events.onViewRemoving.subscribe((data: DG.EventData) => {
    if (data.args?.view?.dart !== view.dart || leaving || !dirty())
      return;
    data.preventDefault();
    void ask().then((ok) => {
      if (!ok)
        return;
      leaving = true;
      release();
      view.close();
    });
  });
  const release = () => {
    armedGates.delete(id);
    setTimeout(() => discardedViews.delete(view.dart), 5000);
    dirtyViews.delete(view.dart);
    leaveAsks.delete(view.dart);
    document.removeEventListener('keydown', onKey, true);
    changed?.unsubscribe();
    clearTimeout(timer);
    window.removeEventListener('beforeunload', onUnload);
    removing.unsubscribe();
    removed.unsubscribe();
  };
  const removed = grok.events.onViewRemoved.subscribe((v: DG.ViewBase) => {
    if (v?.dart === view.dart)
      release();
  });
  document.addEventListener('keydown', onKey, true);
  const changed = changes?.subscribe(soon);
  // (the view made current again: the platform's workspace flag is cleared from the current view only, see mark)
  host.subscribe(grok.events.onCurrentViewChanged, () => grok.shell.v?.dart === view.dart && soon());
  window.addEventListener('beforeunload', onUnload);
  // (and checked a few times while the view's tables and first computation arrive: a reload asks nothing of a view
  // opened and left as it was)
  for (const ms of [1500, 4000, 10000])
    setTimeout(() => armedGates.has(id) && mark(), ms);
  armedGates.add(id);
  dirtyViews.set(view.dart, dirty);
  leaveAsks.set(view.dart, ask);
  host.own(release);
}

/** Save as… for a study table: a name asked (said under the field when it cannot be saved, Save waiting), then the
 *  table renamed as a new study (no row, no ID, not a copy; its plan renamed with it) and `save` run: the one opened
 *  stays as it was saved. Monitor's and Analyze's File › Save as… (QA 0.1.3 round 3 designer U3). */
export function saveStudyAs(df: DG.DataFrame, setMeta: (m: StudyMeta) => void, save: () => Promise<unknown>): void {
  const d = modal('save-study-as', 'Save study as');
  if (!d)
    return;
  const name = ui.input.string('Name', {value: studyNameOf(df)});
  const nameError = dialogLine('');
  nameError.style.color = 'var(--dg-red, #c0504d)';
  d.add(name).add(nameError).onOK(() => {
    const chosen = (name.value ?? '').trim();
    if (!chosen || chosen.length > 100)
      return;
    // a study followed from a plan goes by its plan's name (studyNameOf): the plan it keeps is renamed with it
    try {
      const p = df.getTag(PLAN_TAG) ? parsePlan(df.getTag(PLAN_TAG)!) : null;
      if (p)
        df.setTag(PLAN_TAG, JSON.stringify({...p, name: chosen}));
    } catch {} // a damaged plan: saving says what is wrong
    df.name = chosen;
    df.setTag(ROW_TAG, '');
    markTableFromFile(df, '');
    setMeta({...metaOf(df) ?? emptyMeta(), id: '', copyOf: undefined});
    void save().catch(() => {});
  });
  const ok = okNamed(d, 'Save');
  const sync = () => {
    const n = (name.value ?? '').trim().length;
    nameError.textContent = n === 0 ? 'Give the study a name' : n > 100 ? `A name of 100 characters at most (this ` +
      `one has ${n})` : '';
    nameError.hidden = !nameError.textContent;
    ok.disabled = !nameError.hidden;
  };
  name.input.addEventListener('input', sync);
  sync();
  d.show();
}
