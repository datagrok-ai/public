/* The Monitor view, built as Analyze is: the form generated from pmxMonitorStudy's signature, the results from its
   outputs. The groups are shown, unless the study is blinded (its details): then they stay hidden until someone
   reveals them, and that goes in the study's log. A study's measurements
   are an ordinary table (paste into its grid): Monitor only reads it, and saving the study is a gesture of its own.
   Tabs by question: who needs attention today, how the groups will be formed, how the study is growing, and whether to
   go on, only on the reviews fixed at the start (MONITOR_DESIGN.md). */
import {maxOf} from './format';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {badge, BasicTable, button, buttonWithIcon, computed, div, divH, divV, signal,
  span} from '@datagrok-libraries/u2';
import {funcForm, ViewHost} from '@datagrok-libraries/u2/dg';
import type {TableInput, ViewContent} from '@datagrok-libraries/u2/dg';

import {NOT_A_STUDY, analyzeView, readsAsStudy} from './analyze-view';
import {BandChart} from './band-chart';
import {weekdaysOf} from './frames';
import {isFiniteNumber} from './engine';
import {StudyLink, groupColors, hasOwnColors, linkStudy, paint, resetColors} from './colors';
import {allocateGroups, balanceNote, day0Weekday, MonitorFrames, monitorFrames, MonitorInputs, latestReadings,
  mm3Text, outsideBand, naturalCompare, parseLooks, promptsByChance, removedAnimals, shortOfAnimals,
  weightMissing, nextMeasurement, whyAndStatus, animalIds} from './monitor';
import {FROZEN_ATLAS_NOTE, PLAN_TAG, SavedPlan, atlasNote, parsePlan, planKey, planKeyOf, planOfStudy,
  randomizedGroups, shownPath} from './plan-file';
import {planInput} from './plan-vs-observed';
import {canReturn, keepOnLeaving, linkRow, openRow, openStage, placeOf, showStage} from './stages';
import {LIMITS_TAG, LINE_TAG, LOG_TAG, LOOKS_TAG, REVEALED, VOLUME_TAG, WEIGHT_LABEL, WEIGHT_TAG, addLog,
  groupsShown, logOf,
  LogEntry, newStudyTable, parseStudy, EXCLUDED, EXPOSURE_TAG, exposureOfTag, stamped, lineOf, setHistoryOn,
  restoreFields, roleColumns, saveStudy, studyNameOf, studyOf, tableOf, ticked} from './study-file';
import {exportName} from './study-file';
import {TWIN_NOTE, savedHistory} from './history-models';
import {markLine, modal} from './ui';
import {WeightRow, groupWeights, readSchedule} from './schedule';
import {customizePage, panelTabs, Tab} from './panels';
import {stagePage, studyHeader} from './study-header';
import {openStudyDialog, stageLinks, takesFiles} from './hub-view';
import {animalMarks, animalReadings, blindColumn, emptyStudy, openStudyTable, studyDataLine,
  studyOf as studyNamesOf,
  keepStudyField, readingColumns, tagStudy, withoutExcluded} from './study-table';
import {analysisDeclarationOf, evaluationDay as evaluationDayOf} from './analyze';
import {EmptyState, NoControl} from './crux';
import {META_TAG, StudyMeta, armDoses, dateOf, dateText, decisionEntry, emptyMeta, metaOf, renamedGroups, setMetaOn,
  stopsOf, withGroups, withId, weightChangeText} from './study-meta';
import {renamedExposure} from './exposure';
import {ROW_TAG, Saved, domainStorageOn, isReadOnly, markTableFromFile, newRowNote, ownerOf, rowOfTable, save, saveCopy,
  saveError, setRowOn, tableFileOwner, tableOfRow, withCode, rowsOf as studyRowsOf} from './store';
import {holdTable, inUse, leaveAsks, offerCopy, saveGate, saveStudyAs, savedTables, studyPickList,
  tableSavers} from './native';
import {applyKept, keptOfCall, recallForm, rememberForm} from './form-memory';
import {CAL, DayLook, dayActions, dayToday, icsOf, mark, marked, startMismatch, weeksGrid} from './calendar';
import {LineInput} from './line-input';
import {esc, printSheet, sheetHtml} from './sheet';
import {Exports, download, reasonOf, stageMenu, reportHtml, stageToolbox} from './report';
import {PmxAnimal, PmxGroup, markReading, popupFor} from './handlers';
import {AnimalPanels, AnimalView, growthFrame, placeWords} from './animal-panel';
import {EntryAnimal, LOG_VALUES, enterDayDialog, loadDraft} from './day-entry';
import {studyInputsOf} from './explore';
import {rowsOf, tolerabilityTab} from './tolerability';
import {findCombinations, parseArmName} from './arm-name';
import {linePrior, modelName} from './priors';
import {IndependenceBand, bandWhy, bandWords, independenceBand, notCombined} from './synergy';
import {askedOf, autoRun, block, callInputs, cardGrid, dayText, exampleLine, infoNote, onTableEdits, TABLE_EDITS,
  dialogLine, isBodyWeight, okNamed, showCentred,
  keepNamedColumns, noteLine, notes, notesBox, pct, plainError, settleForm, sized, markStage, stageOf,
  statCard, swatch, volumeAxis, wrapCells, stagePanes, whenSettled, percentRange} from './ui';

/** Monitor's tabs and their panels, in order and width (twelfths): today's work and how the groups grow; the
 *  measurements and the groups formed; how the animals bear it; whether each group is worth continuing. */
const MONITOR_TABS: Tab[] = [
  // (what needs attention first: the day's work; the calendar and the data checks folded to a line)
  {name: 'Today', layout: [{id: 'today', span: 12}, {id: 'growth', span: 12}, {id: 'calendar', span: 12, folded: true},
    {id: 'data', span: 12, folded: true}]},
  // (enrollment first: forming the groups is the study's first step; the sheet is every day's after it)
  {name: 'Measurement sheet', layout: [{id: 'enrollment', span: 12}, {id: 'sheet', span: 12}]},
  {name: 'Welfare', layout: [{id: 'welfare', span: 12}]},
  {name: 'Decisions', layout: [{id: 'decisions', span: 12}]}];

/** Below this chance of success, on a planned review, Monitor suggests stopping the arm. A person decides. */
export const FUTILITY = 0.1;
/** An arm whose TGI interval starts at this (the effect well clear of none) has a clear result already in sight. */
export const CLEAR_TGI = 50;

/** The "chance of a clear result" in words. Only the low end is a figure: on the atlas's real arms a low chance held,
 *  a high one was too hopeful (MONITOR_DESIGN.md §8), so no high figure is shown. Above the low end it reads "not
 *  low", except for an arm whose TGI interval already starts far from no effect (`tgiLower` at CLEAR_TGI or more),
 *  where "not low" would read as doubt about a clearly effective arm: "High", a word from its interval, not the
 *  model's figure. */
export function chanceWords(success: number, tgiLower: number): string {
  if (success < FUTILITY)
    return `Low (under ${pct(FUTILITY)})`;
  return tgiLower >= CLEAR_TGI ? `High (TGI at least ${Math.round(tgiLower)}%)` : `Not low (${pct(FUTILITY)} or more)`;
}

interface AnimalRow {[c: string]: any}

/** The details as a blinded study shows them before its groups are revealed: each group's doses under "Group 1",
 *  "Group 2"… (the calendar and its .ics say when to dose, not what) */
const unnamed = (m: StudyMeta, revealed: boolean): StudyMeta => {
  if (revealed)
    return m;
  const name = blindName(m);
  const known = (x: {group: string}) => m.groups.some((g) => g.group === x.group);
  // the schedules a name gives (no schedule typed) read before the name goes, one per distinct schedule; the
  // exceptions and doses given under the same "Group k" (those of a group no longer in the study left out)
  return {...m, groups: m.groups.flatMap((g) => [...new Set<string>(g.schedule.trim() ? [g.schedule] :
    parseArmName(g.group).map((p) => p.schedule ?? ''))].map((schedule) => ({...g, group: name(g.group),
    schedule, compound: '', dose: null}))),
  exceptions: m.exceptions.filter(known).map((e) => ({...e, group: name(e.group)})),
  given: m.given.filter(known).map((g) => ({...g, group: name(g.group)})),
  ...m.stopped ? {stopped: Object.fromEntries(Object.entries(m.stopped).map(([g, d]) => [name(g), d]))} : {}};
};
/** A group as a blinded study names it before it is revealed: "Group k", numbered by group, not by row (a
 *  combination has a row per drug). */
const blindName = (m: StudyMeta) => (group: string) =>
  `Group ${[...new Set(m.groups.map((g) => g.group))].indexOf(group) + 1}`;
/** What "Enter today's values" asks before its dialog opens: the day to enter, when there is one (null: no dated
 *  day 0, or day 0's date not this table's: the dialog asks for the day number) and the plan's measurement days
 *  (null: no plan). A day off the plan's calendar asks first ('confirm'); any other opens at once. A day entered off
 *  the calendar is then measured like any other. */
export const entryStep = (day: number | null, planned: ((d: number) => boolean) | null): 'open' | 'confirm' =>
  day !== null && day >= 0 && planned !== null && !planned(day) ? 'confirm' : 'open';

/** The question a day off the plan's calendar asks before its dialog opens: the day, and the plan's next measuring
 *  day (QA 0.1.3 round 3: technician N20). */
export function offPlanAsk(day: number, planned: (d: number) => boolean): string {
  const next = Array.from({length: 14}, (_, i) => day + i + 1).find(planned);
  const then = next === undefined ? '' : ` (the next one is day ${next})`;
  return `Day ${day} is not a measurement day in the plan${then}: enter anyway?`;
}

/** A log's time (kept in UTC) as the reader's clock reads it: 2026-11-27 14:05. */
export function localTime(iso: string): string {
  const t = new Date(iso);
  const two = (n: number) => String(n).padStart(2, '0');
  return Number.isNaN(t.getTime()) ? iso : `${t.getFullYear()}-${two(t.getMonth() + 1)}-${two(t.getDate())} ` +
    `${two(t.getHours())}:${two(t.getMinutes())}`;
}
const grams = (x: unknown) => isFiniteNumber(x) ? x.toFixed(1) : '—';
const mm3 = (x: unknown) => isFiniteNumber(x) ? mm3Text(x) : '—';
function fill<R>(host: HTMLElement, items: R[], key: (r: R) => string,
  columns: [string, (r: R) => string | HTMLElement][], onRowClick?: (r: R) => void): void {
  const t = new BasicTable<R>({items, keyOf: key, columns: columns.map(([header, render], j) =>
    ({header, render, align: j === 0 || header === 'Group' ? 'left' : 'center'})), onRowClick});
  if (onRowClick)
    t.root.style.cursor = 'pointer';
  t.root.style.width = '100%';
  wrapCells(t.root);
  host.replaceChildren(t.root);
}

/** The tags a study is saved from (studyOf), besides its measurements. */
const SAVED_TAGS = [PLAN_TAG, META_TAG, LIMITS_TAG, LINE_TAG, VOLUME_TAG, EXPOSURE_TAG, LOOKS_TAG, LOG_TAG];

/** A file opened in place of a study just followed from Plan (still empty) is that study: it takes the plan, the
 *  details and the reviews the followed table carried, unless it carries a plan of its own. True when it took them. */
export function takeOverFollowed(followed: DG.DataFrame, opened: DG.DataFrame): boolean {
  if (opened === followed || followed.rowCount > 0 || !followed.getTag(PLAN_TAG) || opened.getTag(PLAN_TAG))
    return false;
  for (const t of SAVED_TAGS) {
    const v = followed.getTag(t);
    if (v && !opened.getTag(t))
      opened.setTag(t, v);
  }
  return true;
}

/** A log entry as the log shows it: one with its values listed after LOG_VALUES (a day removed) as its line, and
 *  "show values ▸" that unfolds them (QA 0.1.3 round 3: technician N19, chaos 26). */
function logWhat(text: string): string | HTMLElement {
  const at = text.indexOf(LOG_VALUES);
  if (at < 0)
    return text;
  const list = span(text.slice(at + LOG_VALUES.length));
  list.style.display = 'none';
  const more = ui.link('show values ▸', () => {
    const open = list.style.display === 'none';
    list.style.display = open ? 'block' : 'none';
    more.textContent = open ? 'hide values ▾' : 'show values ▸';
  });
  return ui.divV([ui.span([`${text.slice(0, at)} · `, more]), list]);
}

/** The step of a followed study with day 0 entered and no groups yet. */
const DAY0_DONE = 'Day 0 entered. Next: Propose groups (random groups with similar tumor sizes).';

/** The Monitor view; `initial`: the table of a study to follow. */
export async function monitorView(initial?: DG.DataFrame): Promise<DG.ViewBase> {
  const exports: Exports = [];
  const host = ViewHost.create({name: 'PMX Studio', content: (self) => monitorContent(self, initial,
    exports)});
  const view = markStage(host.view, 'monitor');
  // a saved study's URL names it from the start: a reload before its first forecast opens it again
  linkRow(view, 'study', rowOfTable(initial)?.id);
  return view;
}

/** Each study table's groups and their animals as they stood when Monitor first read it or last saved it: a group
 *  renamed since is found by its animals (study-meta.ts renamedGroups). */
const groupsBefore = new WeakMap<DG.DataFrame, Map<string, string[]>>();

/** What the view shows, built under its host, which owns it all (the form, the charts, the panes). */
function monitorContent(self: ViewHost, initial: DG.DataFrame | undefined, exports: Exports): ViewContent {
  const call = DG.Func.byName('PmxStudio:pmxMonitorStudy').prepare();
  const result = signal<MonitorFrames | undefined>(undefined);
  const error = signal<string | null>(null);
  // the error is an empty state (nothing measured yet), not a failure
  const empty = signal(false);
  const stale = signal(false);
  const busy = signal(false);
  // the log's version: what it records (a reveal, the blinding switched) changes what is shown
  const logVersion = signal(0);
  // the study table's values: a flag or an exclusion changes what the lists show without a new forecast
  const dataVersion = signal(0);
  let latest = 0;
  // no Run: each change forecasts again once the changes stop
  let preparing = false;
  const auto = autoRun(() => form, () => preparing ? undefined : run());
  // the study table the view last took in (setStudy): its own reviews and model, not the last study's
  let studySet: DG.DataFrame | null = null;
  // (the form settled: a study field changing from then on is the user's pick, not the form filling itself in)
  let ready = false;
  const form = funcForm(call, {layout: 'tall', overrides: {line: {input: new LineInput({label: 'Tumor model'})}},
    onInputChanged: (name: string) => {
      // another study picked in the field: taken in as an opened one is (setStudy), never with the last one's fields
      const df = name === 'study' && ready ? studyTable() : null;
      if (df && df !== studySet && !preparing)
        return void switchTo(df);
      // (a value the user set for this study: kept for it, given back when the page shows it again)
      if (name !== 'study' && !preparing && studySet && studySet === studyTable()) {
        rememberForm(studySet, 'monitor', keptOfCall(call));
        // the ones saved with the study (its body weight column, limits, model) on its table at once: a change of
        // theirs is a change of the study, marked unsaved; the others are the page's only
        if (['bodyWeight', 'humane', 'weightLimit', 'line'].includes(name))
          remember();
      }
      stale.value = result.peek() !== undefined;
      auto.changed();
    }});
  settleForm(form, {control: 'Auto (the vehicle group)', evaluationDay: 'Auto (endpoint reached)'});
  keepNamedColumns(form, {bodyWeight: isBodyWeight});

  const studyTable = (): DG.DataFrame | null => {
    const v = call.getParamValue('study');
    return v instanceof DG.DataFrame ? v : null;
  };
  // the plan's text: the one picked in the form, else the study table's own
  const planText = () => call.getParamValue('plan') || studyTable()?.getTag(PLAN_TAG);
  // the study as it is saved: its table with the form's limits
  const savedStudy = (df: DG.DataFrame) =>
    studyOf(df, call.getParamValue('humane') ?? null, call.getParamValue('weightLimit') ?? null);
  // the groups, unless the study is blinded and no one has revealed them
  const revealed = computed(() => {
    void result.value;
    void logVersion.value;
    const df = studyTable();
    return !df || groupsShown(df);
  });
  // the Group input is not the user's to change while the groups are hidden
  form.effect(() => {
    const field = form.getInput('group');
    if (field)
      field.disabled = revealed.value ? null : 'Hidden until the groups are revealed';
  });
  // the planned reviews can be changed (a typo fixed) until the study is saved with them or its first review has
  // been measured; then they are the study's, and why says so
  const looksLocked = (df: DG.DataFrame | null): string | null => {
    const fixed = df?.getTag(LOOKS_TAG);
    if (!df || !fixed)
      return null;
    if (rowOfTable(df))
      return 'The reviews were fixed when this study was saved with them';
    let first: number | undefined;
    try {
      first = parseLooks(fixed)[0];
    } catch {
      return null;
    }
    const [day, volume] = [roleCol(df, 'day'), roleCol(df, 'volume')];
    for (let i = 0; first !== undefined && day && volume && i < df.rowCount; i++) {
      if (!day.isNone(i) && !volume.isNone(i) && day.get(i) >= first)
        return `The reviews were fixed once the first one (day ${first}) was measured`;
    }
    return null;
  };
  const lockLooks = () => {
    const field = form.getInput('looks');
    if (field)
      field.disabled = looksLocked(studyTable());
  };
  // what goes with the table: its body weight column and the protocol's limits, for the study file and for Analyze
  const remember = () => {
    const df = studyTable();
    if (!df)
      return;
    const weight = call.getParamValue('bodyWeight');
    df.setTag(WEIGHT_TAG, weight instanceof DG.Column ? weight.name : typeof weight === 'string' ? weight : '');
    df.setTag(LIMITS_TAG, JSON.stringify({humane: call.getParamValue('humane') ?? null,
      weightLimit: call.getParamValue('weightLimit') ?? null}));
    // an empty field leaves the model the study was saved with
    const line = call.getParamValue('line');
    if (typeof line === 'string' && line) {
      df.setTag(LINE_TAG, line);
      setHistoryOn(df, line);
    }
  };
  const planState = {stale, error, result: computed(() => result.value && {analysis: result.value.status}),
    reading: Promise.resolve() as Promise<unknown>, changed: () => auto.changed()};
  // the study table the result on screen is for; a switch to another clears what was shown of the last (true)
  let resultTable: DG.DataFrame | null = null;
  const shownFor = (df: DG.DataFrame | null): boolean => {
    if (df === resultTable)
      return false;
    resultTable = df;
    allocated.value = null;
    return true;
  };
  const run = async (prepare?: () => Promise<void>) => {
    busy.value = true;
    const mine = ++latest;
    try {
      // a study the view fills in changes the form as it goes: none of those changes runs by itself
      preparing = prepare !== undefined;
      await prepare?.();
      preparing = false;
      auto.cancel();
      await planState.reading;
      // a table picked that is no study (a PK table): said so, not "No control group found"
      const picked = studyTable();
      if (picked && !readsAsStudy(picked))
        throw new Error(NOT_A_STUDY);
      // a measurements file opened in place of the study just followed: it is that study, with its plan
      const opened = studyTable();
      if (initial && opened && opened !== initial && takeOverFollowed(initial, opened))
        lockLooks();
      // the reviews once fixed hold, whatever the field says
      if (looksLocked(studyTable()))
        call.setParamValue('looks', studyTable()!.getTag(LOOKS_TAG));
      // a blinded study's group column is hidden from the form's list, which would drop it: given directly
      const df0 = studyTable();
      if (df0 && !groupsShown(df0))
        call.setParamValue('group', roleColumns(df0).group);
      const sent = askedOf(call);
      tagStudy(studyInputsOf(call));
      const frames = await monitorFrames(callInputs<MonitorInputs>(call));
      if (mine !== latest)
        return;
      shownFor(studyTable());
      error.value = null;
      empty.value = false;
      stale.value = askedOf(call) !== sent;
      // the reviews go with the study from its first run, as the field has them until they are fixed
      const df = studyTable();
      if (df && !looksLocked(df))
        df.setTag(LOOKS_TAG, String(call.getParamValue('looks') ?? ''));
      lockLooks();
      remember();
      result.value = frames;
    } catch (e: any) {
      // only the last run speaks: an older one failing late would write over a newer result's place
      if (mine === latest) {
        // another study (a new empty one): nothing of the last one stays on screen, its header and figures included
        if (shownFor(studyTable()))
          result.value = undefined;
        // a study followed from its plan, before day 0 or before its groups: the step to take, not a red error nor
        // "No control group found" (QA 0.1.3 round 3: technician N9, designer B8)
        const planned = !!planText() && !hasGroups(studyTable());
        const firstStep = planned && (e instanceof EmptyState || e instanceof NoControl);
        empty.value = e instanceof EmptyState || firstStep;
        const step = hasVolumes(studyTable()) ? DAY0_DONE : 'Enter day-0 volumes, then Propose groups.';
        error.value = firstStep ? step : e instanceof EmptyState ? e.message : e instanceof NoControl ?
          (revealed.peek() ? e.message : `There is no group "${e.control}" in the study table: check its spelling.`) :
          plainError(e, 'The forecast could not be computed. Try again; if it fails again, reload the page.');
        // what is shown is no longer what the table or the form says
        stale.value = result.peek() !== undefined;
      }
    } finally {
      if (mine === latest) {
        preparing = false;
        busy.value = false;
      }
    }
  };

  // a group given to some animal; a volume in some row (a study before its groups, before its day 0)
  const anyIn = (df: DG.DataFrame | null, role: string) => {
    const c = df ? call.getParamValue(role) instanceof DG.Column ? call.getParamValue(role) as DG.Column :
      df.col(role) ?? df.col('~' + role) : null;
    return !!c && Array.from({length: c.length}, (_, i) => i)
      .some((i) => !c.isNone(i) && String(c.get(i)).trim() !== '');
  };
  const hasGroups = (df: DG.DataFrame | null) => anyIn(df, 'group');
  const hasVolumes = (df: DG.DataFrame | null) => anyIn(df, 'volume');
  // the form reads the group column by name; after a rename it follows
  const followGroup = (df: DG.DataFrame) => {
    const input = form.getInput('group');
    let name: string | null = null;
    try {
      name = roleColumns(df).group.name;
    } catch {
      // no group column: nothing to follow
    }
    if (input && name && input.peek() !== name)
      input.value = name;
  };

  // another study picked in the Study field while this one has changes not saved: the field goes back to it and asks
  // first (Save / Discard / Cancel), before the form takes the new one; Cancel keeps it, its fields as they were (QA
  // 0.1.3 round 3 analyst 8). ponytail: the form gives no hook before a change, so the field is set back at once
  const switchTo = async (df: DG.DataFrame) => {
    const old = studySet;
    // (the view read only with a study already shown: the first one comes while the page is built, before its view)
    const ask = old ? leaveAsks.get(self.view.dart) : undefined;
    if (old && ask && grok.shell.tables.includes(old) && asSaved.changed(old)) {
      (form.getInput('study') as unknown as TableInput | null)!.table = old;
      if (!await ask(true))
        return void run(() => setStudy(old));
    }
    return void run(() => setStudy(df));
  };
  const setStudy = async (df: DG.DataFrame) => {
    studySet = df;
    // a blinded study opens with its group column hidden (`~group`), before any view of the table exists
    blindColumn(df);
    // in the workspace, for the Study field; its table view only when asked for (Data › Measurements, a new study):
    // a view of its own per study opened piled up beside the app's (QA 0.1.3 designer #18)
    // (named "X (2)" while the page it replaces held its own "X": its name back as that one goes, picked again)
    holdTable(self, df, () => {
      const field = form.getInput('study');
      if (field)
        field.value = df.name;
    });
    // back to Monitor once the table is open; a view not yet shown is left to showStage: making it current would add
    // it, and showStage would add it a second time
    if ([...grok.shell.views].some((v) => v.dart === self.view.dart)) {
      grok.shell.v = self.view;
      // (the URL names this study at once: a file opened in place of a saved one no longer reloads the one before,
      // QA 0.1.3 round 4 chaos N6)
      linkRow(self.view, 'study', rowOfTable(df)?.id);
    }
    const study = form.getInput('study');
    if (study)
      study.value = df.name;
    // (a study with no reviews of its own has none: not the last study's, shown in this view before)
    const looks = df.getTag(LOOKS_TAG) || null;
    const field = form.getInput('looks');
    if (field && (looks || field.value))
      field.value = looks;
    await form.settled;
    // the form's table change resets the call's other fields as it settles: the reviews go in after it
    if (looks || call.getParamValue('looks')) {
      if (field)
        field.value = looks;
      call.setParamValue('looks', looks);
    }
    restoreFields(form, df);
    followGroup(df);
    // (what the user set for it on this page before: Back, Forward, the field picking it again)
    applyKept(form, recallForm(df, 'monitor'));
    await form.settled;
    // the study's details read again: they are cached against the result, which a plan followed with nothing
    // measured never changes (its day 0 date stayed null: no date beside the day, no question on Save day)
    metaVersion.value++;
    lockLooks();
  };
  const newStudy = () => run(async () => {
    const df = emptyStudy();
    await setStudy(df);
    tagStudy({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: 'body_weight'});
    openStudyTable(df);
  });
  whenSettled(form, () => {
    keepStudyField(form.getInput('study') as unknown as TableInput, () => call.getParamValue('study') ?? null);
    const study = form.getInput('study') as unknown as {actions?: {id: string, icon: string, tooltip: string,
      run: () => void}[]} | null;
    // one way in for new data, as in Analyze: File > Open (Ctrl+O) and a drop, through the importer's checks (the
    // separator, decimal commas, dates as days, units); the platform's own open-file, Files and database icons read a
    // table raw (a European CSV's dates taken for its volumes). The field picks among the studies open, or a new one
    if (study?.actions) {
      study.actions = [{id: 'new', icon: 'plus',
        tooltip: 'Start an empty study table, to type or paste the measurements into', run: () => void newStudy()}];
    }
    const picker = form.getInput('study') as unknown as TableInput | null;
    if (picker) {
      const offer = picker.setItems.bind(picker);
      picker.setItems = (items: string[]) => offer(studyPickList(items.filter((name) => {
        const t = grok.shell.tables.find((x) => x.name === name);
        return !t || name === picker.value || readsAsStudy(t);
      }), picker.value));
      picker.setItems(grok.shell.tables.map((t) => t.name));
    }
    ready = true;
    // (no table given: the study the field picked by itself is taken in as an opened one, its own model and fields
    // with it, not the form's defaults: QA 0.1.3 round 3, a CT26 plan via the Monitor step read MC38)
    const first = initial ?? studyTable();
    if (first)
      void run(() => setStudy(first));
  });
  planInput(form, call, planState);

  const start = divV([
    span('Open the study\'s measurements (one row per animal and day, day 0 the randomization), check the columns ' +
      'and the forecasts follow. Paste or type each session\'s values into the table: each is checked as you enter ' +
      'it, and the forecasts update.'),
    exampleLine(),
  ], 'pmx-monitor-start');
  Object.assign(start.style, {gap: 'var(--dg-space-m)', maxWidth: '760px'});

  const status = computed(() => result.value?.status);
  const cell = <T = number>(name: string) => computed(() => {
    const s = status.value;
    return (s && !s.col(name)?.isNone(0) ? s.get(name, 0) ?? undefined : undefined) as T | undefined;
  });
  const card = statCard.bind(null, form);
  const metaVersion = signal(0);
  const meta = computed(() => {
    void result.value;
    void metaVersion.value;
    void logVersion.value; // the arms stopped are read from the log
    return metaOf(studyTable()) ?? emptyMeta();
  });
  // saved (an ID given) or measured on: the reviews' lock follows
  form.effect(() => {
    void meta.value;
    lockLooks();
  });
  // the expected end, how it stands against the plan's
  const [needing, endMedian, plannedEnd] = [cell('attention'), cell('end_median'), cell('planned_end')];
  const versusPlan = () => {
    const [e, p] = [endMedian.value, plannedEnd.value];
    if (e === undefined || p === undefined)
      return '';
    const d = Math.round(e - p);
    return d === 0 ? 'as planned' : `${Math.abs(d)} ${Math.abs(d) === 1 ? 'day' : 'days'} ` +
      `${d > 0 ? 'later' : 'earlier'} than planned`;
  };
  // (a change is shown against the previous result of the same study only)
  const attentionCard = card('Animals needing attention', needing,
    () => 'Animals listed under Today, each with its reason', undefined, true, studyTable);
  form.effect(() => attentionCard.classList.toggle('pmx-stat-attn', (needing.value ?? 0) > 0));
  // the next measurement from today on, never a day already past: a planned day passed with no values is said
  // missed beside it; an ended study has none (QA 0.1.3 round 3: technician N12, designer W1, chaos 7, analyst 14)
  const nextDue = computed(() => {
    const s = status.value;
    if (!s || s.get('ended', 0))
      return {next: undefined, missed: null, dated: false};
    const latest = s.get('day', 0) as number;
    const today = dayToday(meta.value.start);
    const future = String(s.get('future_days', 0) ?? '').split(',').filter((x) => x).map(Number);
    const at = today !== null && today >= latest ? today : null;
    return {...nextMeasurement(future.length ? future : [s.get('next_day', 0)], at), dated: at !== null};
  });
  const nextValue = computed(() => nextDue.value.next);
  const nextLabel = computed(() => {
    const {next, missed, dated} = nextDue.value;
    const on = (d: number, weekday = false) => !dated ? '' : dateOf(meta.value.start, d)?.toLocaleDateString('en-GB',
      {...weekday ? {weekday: 'short'} : {}, day: 'numeric', month: 'short', timeZone: 'UTC'}) ?? '';
    // (the day "Next: day 14", its date beside: "Day 14" alone read as today, QA 0.1.3 round 4 chaos wording)
    return [next === undefined ? '' : on(next, true) || 'Next measurement',
      missed === null ? '' : `day ${missed}${on(missed) ? ` (${on(missed)})` : ''} missed`]
      .filter((x) => x).join(' · ') || 'Next measurement';
  });
  const nextCard = card(nextLabel, nextValue, () => 'The next measuring day from today on (the plan\'s weekdays, ' +
    'else the study\'s usual interval); a planned day passed with no values is said missed',
  (x) => `Next: ${dayText(x).toLowerCase()}`, true, studyTable);
  form.effect(() => nextCard.hidden = !!status.value?.get('ended', 0));
  // two figures, what to do now; the expected end (against the plan's) is in the study's status line
  const cards = cardGrid([attentionCard, nextCard], 'pmx-monitor-cards', 150);
  const endStatus = computed(() => {
    const e = endMedian.value;
    // a study every animal has left (a removal reason each), or one measured past its expected end, says so: not an
    // end to expect
    const latest = status.value?.get('day', 0) as number | undefined;
    if (status.value?.get('ended', 0))
      return 'ended';
    if (e === undefined)
      return '';
    if (latest !== undefined && latest > Math.round(e))
      return `past the expected end (~${dayText(e).toLowerCase()})`;
    // ("ends ~day 14"; against the plan only when it differs: designer wording 13)
    const vs = versusPlan();
    return `ends ~${dayText(e).toLowerCase()}` + (vs && vs !== 'as planned' ? ` (${vs})` : '');
  });
  const endTip = () => {
    const [a, b] = [cell('end_lower').peek(), cell('end_upper').peek()];
    // not "8 in 10": halfway through, the atlas's MC38 studies fell in this range about half the time (§8)
    return `The day the median control reaches the endpoint, likely between ${dayText(a)} and ` +
      `${b === undefined ? 'later than the days forecast' : dayText(b)}; the plan's last day: ` +
      `${plannedEnd.peek() === undefined ? 'no plan' : dayText(plannedEnd.peek())}`;
  };

  // Today
  const attention = div([], 'pmx-monitor-attention');
  const sheet = div([], 'pmx-monitor-sheet');
  // an animal's menu on its row (its name, the first cell), in the attention list and the measurement sheet alike
  for (const list of [attention, sheet]) {
    list.addEventListener('contextmenu', (e) => {
      const row = (e.target as HTMLElement).closest('tbody tr');
      const name = (row?.querySelector('[data-animal]') as HTMLElement | null)?.dataset.animal ??
        row?.querySelector('td')?.textContent?.trim();
      if (!name)
        return;
      e.preventDefault();
      popupFor(animalPanels.objectOf(name), {causedBy: e});
    });
  }
  const noteBox = notesBox();
  // each value as it goes in, against what the last run expected: no run needed to hear "measure again"
  const entered = div([], 'pmx-monitor-entered');
  const checkEntered = () => {
    const [r, df] = [result.peek(), studyTable()];
    // the columns as the last run read them (the table's own names, whatever they are)
    let cols: ReturnType<typeof roleColumns> | null = null;
    try {
      cols = df ? roleColumns(df) : null;
    } catch {
      // not a study yet
    }
    if (!r || !df || !cols) {
      entered.replaceChildren();
      return;
    }
    const {animal, day, volume} = cols;
    const nextDay = r.status.get('next_day', 0);
    const values = Array.from({length: df.rowCount}, (_, i) => ({animal: String(animal.get(i)),
      day: day.get(i), volume: volume.isNone(i) ? NaN : volume.get(i)}));
    const off = outsideBand(values, nextDay, new Map(rowsOf(r.animals).map((a) => [a.animal, a as any])));
    // a negative volume is impossible on any day: corrected, not measured again
    const bad = values.filter((v) => v.volume < 0);
    entered.replaceChildren(...bad.length === 0 ? [] : [noteLine('Impossible, correct it: ' + bad.map((o) =>
      `${o.animal} at ${mm3(o.volume)} mm³ on day ${o.day}`).join('; ') + ' (a volume cannot be negative)', 'error')],
    ...off.length === 0 ? [] : [noteLine(`Measure again (${dayText(nextDay)}): ` +
      off.map((o) => `${o.animal} at ${mm3(o.volume)} mm³, expected ${mm3(o.lower)}–${mm3(o.upper)}`).join('; '),
    'warning')]);
  };
  let watching: {unsubscribe(): void}[] = [];
  // edits come in bursts (a paste, values typed one after another): the checks and the lists redone once after them
  const edited = () => {
    checkEntered();
    dataVersion.value++;
  };
  // (unsubscribed: a pending check goes with it)
  const stopWatching = () => {
    watching.forEach((w) => w.unsubscribe());
    watching = [];
  };
  form.effect(() => {
    void result.value;
    stopWatching();
    const df = studyTable();
    if (df)
      watching = [DG.debounce(df.onDataChanged, 200).subscribe(edited)]; // (rows added are among its events)
    checkEntered();
  });
  const groupCols = (): [string, (r: AnimalRow) => string][] => revealed.value ? [['Group', (a) => a.group]] : [];
  // an animal's expected next value: the range its flag reads ("outside the expected …: measure again"), the one
  // range the lists, the sheet, the day's entry and the flags say (QA 0.1.3 round 3: technician N7)
  // (the forecast's own narrower range, the old "Next expected" column, and its median on hover: nothing lost)
  const expectedRange = (a: AnimalRow): [number, number] | null =>
    a.next_median === null || a.check_lower === null || a.check_upper === null ? null : [a.check_lower, a.check_upper];
  const expectedText = (a: AnimalRow, unit = ''): string | HTMLElement => {
    const range = expectedRange(a);
    if (!range)
      return '—';
    // the likely value and its narrow range shown, the wide range that asks to measure again on hover (designer
    // wording 16: "47–405 mm³" read as a guess, not a forecast)
    const likely = a.next_lower !== null && a.next_upper !== null ?
      ` (${mm3(a.next_lower)}–${mm3(a.next_upper)})` : '';
    const el = span(`~${mm3(a.next_median)}${unit}${likely}`);
    ui.tooltip.bind(el, `Most likely ${mm3(a.next_median)} mm³${likely}. Outside ${mm3(a.check_lower)}–` +
      `${mm3(a.check_upper)} mm³: measure again.`);
    return el;
  };
  // an animal clicked: Analyze's panel with what Monitor has (its growth against the forecast, its place in its group
  // and its colour only once the groups are revealed, its weight, its details), cleared with the view
  let shownAnimal: object | null = null;
  self.own(() => {
    if (shownAnimal && grok.shell.o === shownAnimal)
      grok.shell.o = null;
  });
  const animalRow = (animal: string) => {
    const r = result.peek();
    return r ? rowsOf(r.animals).find((a) => String(a.animal) === animal) : undefined;
  };
  const animalView = (animal: string): AnimalView => {
    const n = studyNamesOf(studyTable());
    const r = result.peek();
    const a = animalRow(animal);
    if (!n || !r || !a)
      return {};
    const color = revealed.peek() ? colorsOf(r)[String(a.group)] : undefined;
    const readings = animalReadings(n, animal);
    const next = r.status.get('next_day', 0);
    const growth = growthFrame(readings.days, readings.volumes, color, typeof a.next_median === 'number' &&
      typeof next === 'number' ? {day: next, lower: a.next_lower, median: a.next_median, upper: a.next_upper} : null);
    // its place in its group by its last volume, the smallest last (a blinded study: none, it would show the groups)
    const bars = revealed.peek() ? rowsOf(r.animals).filter((b) => b.group === a.group &&
      typeof b.last_volume === 'number').map((b) => ({animal: String(b.animal), group: String(b.group),
      value: b.last_volume as number})) : [];
    const place = placeWords(bars, animal);
    // a treated animal: Analyze's answer (against its twin) is not Monitor's, said where Analyze puts it (A11)
    const treated = revealed.peek() && String(a.group) !== String(r.status.get('control', 0));
    return {verdict: treated ? {text: 'Twin comparison: open in Analyze', tone: 'neutral'} : null,
      curve: growth ? {title: 'Growth', frame: growth} : null,
      group: place ? {bars, order: 'desc', summary: `${place}, ${mm3(a.last_volume)} mm³`} : null,
      weightLimit: call.getParamValue('weightLimit') ?? 20};
  };
  // the animal shown, marked in the lists as Analyze's list marks it (an accent bar at the row's start)
  const markShown = (animal: string) => [attention, sheet].forEach((list) => list.querySelectorAll('tbody tr')
    .forEach((row) => {
      const on = !!row.querySelector(`[data-animal="${CSS.escape(animal)}"]`);
      const first = row.querySelector<HTMLElement>('td');
      if (first)
        first.style.boxShadow = on ? 'inset 3px 0 0 var(--dg-accent, #2083d5)' : '';
    }));
  const growthAxis = volumeAxis();
  const animalPanels = new AnimalPanels({names: () => studyNamesOf(studyTable()),
    color: (animal) => {
      const r = result.peek();
      const a = animalRow(animal);
      return r && a && revealed.peek() ? colorsOf(r)[String(a.group)] : undefined;
    },
    view: animalView, log: growthAxis.log, onShow: markShown,
    // another study opened, a value corrected: the panel shown follows (A1, T11)
    changed: () => result.value});
  self.own(() => animalPanels.dispose());
  const showAnimal = (a: AnimalRow) => {
    if (studyNamesOf(studyTable()) && result.peek())
      shownAnimal = animalPanels.show(String(a.animal));
  };
  // only for a study with body weights
  const weightCols = (r: MonitorFrames): [string, (a: AnimalRow) => string | HTMLElement][] =>
    r.animals.col('weight_change') ?
      [['Body weight vs day 0', (a) => {
        // a gain says so (+4%), as a loss does (−4%); near the limit as the reason it needs attention says it
        const el = span(typeof a.weight_change === 'number' ?
          weightChangeText(a.weight_change, call.getParamValue('weightLimit') ?? 20) : '—');
        if (typeof a.last_weight === 'number')
          ui.tooltip.bind(el, `Last weighed: ${grams(a.last_weight)} g`);
        return el;
      }]] : [];
  form.effect(() => {
    const r = result.value;
    void dataVersion.value;
    if (!r)
      return;
    // an animal flagged or excluded says so by its name, as the study table colours its rows
    const n = studyNamesOf(studyTable());
    const marks = n ? animalMarks(n) : new Map<string, string>();
    const animalCell = (a: AnimalRow) => {
      const m = marks.get(String(a.animal));
      const el = span(String(a.animal));
      if (m) {
        el.append(' ', badge(m === 'excluded' ? 'Excluded' : m === 'partly' ? 'Some values excluded' : 'Flagged',
          {variant: m === 'flagged' ? 'warning' : 'default'}));
      }
      el.dataset.animal = String(a.animal);
      // its details on a click on its name only: a click elsewhere in a row leaves the context panel as it is
      el.style.cursor = 'pointer';
      el.title = 'Show this animal\'s details';
      el.addEventListener('click', () => showAnimal(a));
      return el;
    };
    // by animal: the engine's order, group by group, would show the groups
    const all = rowsOf(r.animals).sort((a, b) => naturalCompare(String(a.animal), String(b.animal)));
    const listed = all.filter((a) => a.attention);
    // a day under way: the groups not reached yet in one line under the real alerts, not an animal each
    const waiting = r.status.get('waiting', 0) as number;
    const groupsLeft = String(r.status.get('waiting_groups', 0) ?? '').split(', ').filter((g) => g);
    const waitLine = waiting ? ui.divText(`${waiting} ${waiting === 1 ? 'animal' : 'animals'} to measure today ` +
      `(day ${r.status.get('day', 0)}): ${revealed.value ? groupsLeft.join(', ') : `${groupsLeft.length} ` +
      `${groupsLeft.length === 1 ? 'group' : 'groups'}`} not measured yet.`) : null;
    if (listed.length === 0) {
      attention.replaceChildren(noteLine(waiting ? 'No animal measured so far needs attention.' :
        'No animal needs attention today.', 'success'));
    } else {
      // why apart from what to do now (QA 0.1.3 round 3, technician: "Needs attention split Why / Status")
      fill(attention, listed, (a) => a.animal, [['Animal', animalCell], ...groupCols(),
        ['Why', (a) => whyAndStatus(a.attention).why],
        ['Status', (a) => {
          const said = whyAndStatus(a.attention).status;
          if (!/Likely a typo|Impossible volume/.test(a.attention) || !n)
            return said;
          const fix = ui.link('Correct…', () => markReading(new PmxAnimal(n, String(a.animal)), 'correct'),
            'Correct this value: the old one stays in the log');
          fix.addEventListener('click', (e) => e.stopPropagation());
          return ui.span([...said ? [`${said} `] : [], fix]);
        }], ['Last volume', (a) => `${mm3(a.last_volume)} mm³, day ${a.last_day}`],
        ...weightCols(r), ['Next expected', (a) => expectedText(a, ' mm³')]]);
    }
    if (waitLine)
      attention.append(waitLine);
    fill(sheet, all, (a) => a.animal, [['Animal', animalCell], ...groupCols(),
      ['Last volume, mm³', (a) => mm3(a.last_volume)], ...weightCols(r),
      ['Next expected, mm³', (a) => expectedText(a)]]);
    // a prompt to measure again, not a verdict: said with what it costs by chance (no correction over animals and
    // sessions, which would let a mistyped value through)
    const chance: number | null = r.status.col('check_chance')?.isNone(0) === false ?
      r.status.get('check_chance', 0) : null;
    const above: number = chance && r.status.col('check_above')?.isNone(0) === false ?
      r.status.get('check_above', 0) : 0.5;
    const side = above > 2 / 3 ? ', mostly above it' : above < 1 / 3 ? ', mostly below it' : '';
    // one line; what it costs by chance on hover
    const again = infoNote('A value outside its range: measure again; a second reading that agrees stands.');
    ui.tooltip.bind(again, chance ? `In the atlas's published controls a correct value fell outside its range about ` +
      `1 time in ${Math.round(1 / chance)} (the average over studies)${side}: with ${all.length} animals, ` +
      `${promptsByChance(all.length, chance)} by chance.` :
      `How often a correct value falls outside its range is not measured for ${r.status.get('line', 0)} yet.`);
    sheet.append(again);
    noteBox.replaceChildren(...notes(r.status.get('notes', 0)));
  });

  // Enrollment
  const allocationTable = div([], 'pmx-monitor-allocation');
  // `onlyNew` (a study under way): only the animals without a group are placed, in the room the groups have left
  const propose = async (onlyNew = false) => {
    const df = studyTable();
    if (!df)
      throw new Error('Open the study\'s measurements first');
    // the plan the form gives, else the table's own
    const text = planText();
    if (!text)
      throw new Error('Propose groups needs the study\'s plan, for its groups and animals: add it in the Plan field.');
    const plan: SavedPlan = parsePlan(text);
    const [animal, day, volume] = ['animal', 'day', 'volume'].map((role) => roleCol(df, role));
    if (!animal || !day || !volume)
      throw new Error('Propose groups needs the animals, their days and volumes: pick those columns in the form');
    // each animal's latest volume
    const wName = df.getTag(WEIGHT_TAG);
    const weightCol = roleCol(df, 'bodyWeight') ?? (wName ? df.col(wName) : null);
    const had = roleCol(df, 'group');
    const groupOf = (i: number) => had && !had.isNone(i) ? String(had.get(i)).trim() : '';
    const placed = new Map<string, string>();
    for (let i = 0; i < df.rowCount; i++) {
      if (groupOf(i))
        placed.set(String(animal.get(i)), groupOf(i));
    }
    // (ponytail: a group is matched to the plan's by name; a study that renamed its control places none in it)
    const left = (g: string, n: number) => Math.max(0, n - [...placed.values()].filter((x) => x === g).length);
    // each animal's latest volume (under way: of those without a group)
    const latest = latestReadings(Array.from({length: df.rowCount}, (_, i) => i)
      .filter((i) => !onlyNew || !placed.has(String(animal.get(i)))).map((i) => ({animal: String(animal.get(i)),
        day: day.get(i), volume: volume.isNone(i) ? NaN : volume.get(i),
        weight: weightCol && !weightCol.isNone(i) ? Number(weightCol.get(i)) : null})));
    const planned = randomizedGroups(plan);
    const groups = onlyNew ? planned.map((g) => ({...g, n: left(g.name, g.n)})).filter((g) => g.n > 0) : planned;
    if (onlyNew && !groups.length)
      throw new Error('Every group already has its planned animals: none left to place the new ones in.');
    // what the engine would refuse, said as what to do (a retry or a reload changes neither)
    const unsized = groups.filter((g) => !(g.n > 0)).map((g) => g.name);
    if (unsized.length) {
      throw new Error(`Every group needs at least one animal, and the plan gives none to ${unsized.join(', ')}: size ` +
        `${unsized.length > 1 ? 'them' : 'it'} in Plan, save the plan and follow it again.`);
    }
    if (latest.animal.length === 0)
      throw new Error('No animal has a tumor volume yet: enter each animal\'s volume, then propose the groups.');
    const needed = groups.reduce((t, g) => t + g.n, 0);
    if (!onlyNew && latest.animal.length < needed) {
      throw new Error(shortOfAnimals(needed, latest.animal, Array.from({length: df.rowCount}, (_, i) => i)
        .filter((i) => !animal.isNone(i)).map((i) => String(animal.get(i)))));
    }
    // a new seed per study: it decides who gets what
    const seed = crypto.getRandomValues(new Uint32Array(1))[0];
    const {allocation, balance, weightUsed} = await allocateGroups(latest.animal, latest.volume, groups, seed,
      weightCol ? latest.weight : undefined);
    const note = balanceNote(!!weightCol, weightMissing(latest.volume, latest.weight), weightUsed);
    const group = new Map(rowsOf(allocation).map((r) => [r.animal, r.group]));
    const col = roleCol(df, 'group') ?? df.columns.addNewString(groupsShown(df) ? 'group' : '~group');
    for (let i = 0; i < df.rowCount; i++) {
      if (!onlyNew || !groupOf(i))
        col.set(i, group.get(String(animal.get(i))) ?? '');
    }
    // more animals than the plan's places: those the balance left out are named, not dropped in silence
    const out = latest.animal.filter((a) => !group.get(a)).sort(naturalCompare);
    const outNote = out.length ? `${out.length} ${out.length === 1 ? 'animal' : 'animals'} left without a group ` +
      `(the plan has ${needed} places for ${latest.animal.length}): ${out.join(', ')}` : '';
    record(`${onlyNew ? `New animals placed (${latest.animal.length - out.length})` : 'Groups allocated'} ` +
      `(seed ${seed})${weightUsed ? '; balanced on body weight too' : ''}${outNote ? `; ${outNote}` : ''}`);
    allocated.value = {bal: rowsOf(balance), weightUsed, note, outNote};
    proposeProblem.value = '';
  };
  // the balance of the groups just proposed; while the groups are hidden (blinded then or since), named as the
  // calendar names them ("Group 1"), not by their treatment, and the animals left out without their names
  const allocated = signal<{bal: AnimalRow[], weightUsed: boolean, note: string | null, outNote: string} | null>(null);
  form.effect(() => {
    const a = allocated.value;
    const shown = revealed.value;
    if (!a)
      return void allocationTable.replaceChildren();
    const named = (b: AnimalRow) => shown ? b.group : `Group ${a.bal.indexOf(b) + 1}`;
    fill(allocationTable, a.bal, (b) => b.group, [['Group', named], ['Animals', (b) => String(b.n)],
      ['Geometric mean volume', (b) => `${mm3(b.geo_mean)} mm³`],
      ['Volume spread (CV)', (b) => `${Math.round(100 * Math.sqrt(Math.exp(b.sd_log * b.sd_log) - 1))}%`],
      ...a.weightUsed ? [['Body weight, g (mean ± sd)',
        (b: AnimalRow) => `${grams(b.weight_mean)} ± ${grams(b.weight_sd)}`] as [string, (b: AnimalRow) => string]] :
        []]);
    if (a.note)
      allocationTable.append(ui.divText(a.note));
    if (a.outNote)
      allocationTable.append(noteLine(`${a.outNote}.`));
  });
  // an action's failure is said where it happened and goes away (the run's error line stays the run's)
  const failed = (what: string) => (e: any) =>
    grok.shell.error(plainError(e, `${what}. Try again; if it fails again, reload the page.`));
  // Propose groups' refusal said on the page, beside the buttons that asked, as well as in a balloon (QA 0.1.3 round 4
  // designer D4: only a fading balloon said why nothing happened)
  const proposeProblem = signal('');
  const proposeFailed = (what: string) => (e: any) => {
    proposeProblem.value = plainError(e, `${what}. Try again; if it fails again, reload the page.`);
    failed(what)(e);
  };
  const problemLine = () => {
    const line = span('');
    line.style.color = CAL.red;
    form.effect(() => {
      line.textContent = proposeProblem.value;
      line.hidden = !proposeProblem.value;
    });
    return line;
  };
  // groups already in the table are someone's allocation, perhaps the one the animals are being dosed by, so they
  // are replaced only when asked twice
  const proposeAsked = () => {
    const df = studyTable();
    const groups = df ? roleCol(df, 'group') : null;
    const allocated = !!df && !!groups &&
      Array.from({length: df.rowCount}, (_, i) => i).some((i) => !groups.isNone(i) && String(groups.get(i)).trim());
    if (!allocated)
      return void propose().catch(proposeFailed('The groups could not be proposed'));
    // under way (measured after day 0): the animals keep their groups; only those without one are placed
    const day = roleCol(df, 'day');
    if (day && Array.from({length: df.rowCount}, (_, i) => day.get(i)).some((x) => typeof x === 'number' && x > 0)) {
      const animal = roleCol(df, 'animal');
      const without = new Set(Array.from({length: df.rowCount}, (_, i) => i)
        .filter((i) => groups!.isNone(i) || !String(groups!.get(i)).trim()).map((i) => String(animal?.get(i))));
      if (!animal || !without.size) {
        grok.shell.info('The study is under way and every animal has its group: groups are proposed before ' +
          'dosing starts. A new animal without a group would be placed in its room.');
        return;
      }
      return void propose(true).catch(proposeFailed('The new animals could not be placed'));
    }
    const d = modal('replace-groups', 'Replace the groups?');
    if (!d)
      return;
    d.add(dialogLine('Every animal gets a new random group', 'From its latest volume, in place of the one it has. ' +
        'The study\'s log keeps the new seed, not the old groups.'))
      .onOK(() => void propose().catch(proposeFailed('The groups could not be proposed')));
    okNamed(d, 'Replace');
    showCentred(d);
  };
  const proposeButton = buttonWithIcon('Propose groups', 'random', proposeAsked);
  proposeButton.title = 'Fill the table\'s group column with random groups alike in their tumors';
  // the animal room's sheet, from the groups in the table and the seed the log kept
  const printRandomization = () => {
    const df = studyTable();
    if (!df || !roleCol(df, 'group')) {
      grok.shell.warning('The sheet needs the study\'s groups: propose them, or open a table with a group column.');
      return;
    }
    const seed = [...logOf(df)].reverse().map((e) => /seed (\d+)/.exec(e.what)?.[1]).find((s) => s) ?? null;
    try {
      printSheet(sheetHtml(df.name, metaOf(df), df, seed));
      addLog(df, 'Randomization sheet printed (it shows the groups)');
      logVersion.value++;
    } catch (e: any) {
      failed('The sheet could not be printed')(e);
    }
  };
  const sheetButton = buttonWithIcon('Print randomization sheet', 'print', printRandomization);
  sheetButton.title = 'A page for the animal room: each animal, its group and its volume at randomization, ' +
    'what each group receives, and lines to sign. It shows the groups.';

  // Growth: by group, once revealed
  const colorsOf = (r: MonitorFrames) => {
    const names = [...new Set(r.days.col('group')!.toList() as string[])];
    const m = metaOf(studyTable());
    const colors = groupColors(names, r.status.get('control', 0), m?.colors, m?.slots);
    return Object.fromEntries(names.map((n) => [n, colors[n]]));
  };
  // each combination found in the groups' names (unless the plan unticked it), its drugs started together: where A+B
  // would run if its two drugs acted independently, dashed beside it, and whether it slows growth more or less on the
  // last day measured (synergy.ts independenceBand); computed after each reading, shown when it arrives
  const independent = signal(new Map<string, {band: IndependenceBand, a: string, b: string}>());
  let bandedFor: MonitorFrames | undefined;
  form.effect(() => {
    const r = result.value;
    const s = studyInputsOf(call);
    if (!revealed.value || r === bandedFor)
      return;
    const mine = bandedFor = r;
    if (independent.peek().size)
      independent.value = new Map();
    if (!r || !s)
      return;
    const control = r.status.get('control', 0) as string;
    const treated = [...new Set(r.days.col('group')!.toList() as string[])].filter((g) => g !== control);
    const off = notCombined(planText());
    // the band reads A and B alone as they act inside A+B: not for drugs started on other days (a staggered design)
    const starts = (g: string) => parseArmName(g).map((p) => readSchedule(p.schedule)?.first ?? 0);
    const combos = findCombinations(treated).filter((c) => !off.includes(treated[c.ab]) &&
      new Set([c.a, c.b, c.ab].flatMap((i) => starts(treated[i]))).size === 1);
    if (!combos.length)
      return;
    void linePrior(r.status.get('line', 0)).then(async (lp) => {
      const logBeta = lp.study_prior?.mean_log_beta;
      if (logBeta === undefined)
        return;
      const found = new Map<string, {band: IndependenceBand, a: string, b: string}>();
      for (const c of combos) {
        const [a, b, ab] = [treated[c.a], treated[c.b], treated[c.ab]];
        const band = await independenceBand(s, {control, a, b, ab}, Math.exp(logBeta), r.status.get('day', 0))
          .catch(() => null);
        if (band)
          found.set(ab, {band, a, b});
      }
      if (bandedFor === mine && found.size)
        independent.value = found;
    }).catch(() => {});
  });
  // (in words, designer wording 14; with two combinations, each named)
  const independentName = (ab: string) => independent.peek().size > 1 ? `${ab}: expected if its two drugs simply ` +
    'add up' : 'Expected if the two drugs simply add up';
  // the dashed categories of the growth chart: the bands of independence (BandChart reads them each time it draws)
  const dashed: string[] = [];
  const growthData = computed(() => {
    const r = result.value;
    if (!r || !revealed.value)
      return undefined;
    // the bands of independence as more groups, in their combination's colour (a copy: the frame stays Monitor's)
    const bands = independent.value;
    const t = bands.size ? r.days.clone() : r.days;
    dashed.splice(0, dashed.length, ...[...bands.keys()].map(independentName));
    for (const [ab, {band: x}] of bands) {
      x.days.forEach((d, j) => {
        const i = t.rowCount;
        t.rows.addNew();
        t.set('group', i, independentName(ab));
        t.set('day', i, d);
        t.set('lower', i, Math.exp(x.base + x.lower[j]));
        t.set('median', i, Math.exp(x.base + x.mid[j]));
        t.set('upper', i, Math.exp(x.base + x.upper[j]));
      });
    }
    const colors = colorsOf(r);
    for (const ab of bands.keys())
      colors[independentName(ab)] = colors[ab];
    paint(t, 'group', colors);
    // what the study is read against: the endpoint, today and the expected end (the chart's scale follows them)
    const lines = t.meta.formulaLines;
    lines.clear();
    const endpoint = Number(call.getParamValue('endpoint'));
    if (endpoint > 0)
      lines.addLine({title: `Endpoint ${endpoint} mm³`, formula: `\${median} = ${endpoint}`});
    // today as the calendar has it: from day 0's date, else the last day measured; none before day 0
    const today = dayToday(meta.value.start) ?? r.status.get('day', 0) as number;
    if (today >= 0)
      lines.addLine({title: `Today, day ${today}`, formula: `\${day} = ${today}`, color: '#2083d5'});
    const end = r.status.get('end_median', 0);
    if (Number.isFinite(end) && end !== null)
      lines.addLine({title: `Expected end, day ${end}`, formula: `\${day} = ${end}`, color: '#9aa5b1'});
    // the notebook on the days it speaks of: what happened stays next to how the groups grew
    const notes = new Map<number, string[]>();
    for (const n of meta.value.notes)
      notes.set(n.day, [...notes.get(n.day) ?? [], n.animal ? `${n.animal}: ${n.text}` : n.text]);
    for (const [day, texts] of notes) {
      const text = texts.join('; ');
      lines.addLine({title: text.length > 32 ? `${text.slice(0, 31)}…` : text, formula: `\${day} = ${day}`,
        color: '#8b5cc4'});
    }
    return t;
  });
  const growth = new BandChart(growthData, {x: 'day', lower: 'lower', upper: 'upper', line: 'median', split: 'group',
    xTitle: 'Day after randomization', yTitle: 'Median volume (mm³)', log: growthAxis.log, markers: true, fit: true,
    redrawOn: meta, dashed});
  growth.root.dataset.pmx = 'monitor-growth';
  // the reveal where it is asked for, as in Decisions and Welfare (`reveal` is defined with the log, below)
  const growthHidden = divV([infoNote('<b>The groups are hidden.</b> Reveal them to see each group\'s growth.'),
    button('Reveal the groups', () => reveal())]);
  growthHidden.style.gap = 'var(--dg-space-m)';
  // under the chart, a line per combination: where it stands against its band on the last day measured
  const bandLines = div([], 'pmx-band-lines');
  bandLines.dataset.pmx = 'monitor-independence';
  form.effect(() => {
    bandLines.replaceChildren(...[...independent.value].map(([ab, x]) => {
      const line = ui.divText(bandWords(ab, x.band));
      line.title = bandWhy(x.a, x.b);
      return line;
    }));
  });
  form.effect(() => {
    growth.root.hidden = !revealed.value;
    growthAxis.row.hidden = !revealed.value;
    growthHidden.hidden = revealed.value;
    bandLines.hidden = !revealed.value;
  });

  // Calendar: the study by its days, today among them once day 0 has a date
  const setMeta = (m: StudyMeta) => {
    const df = studyTable();
    if (df) {
      setMetaOn(df, m);
      metaVersion.value++;
    }
  };
  // a column the form maps a role to, else the study table's own of that name
  const roleCol = (df: DG.DataFrame, role: string): DG.Column | null => {
    const c = call.getParamValue(role);
    return c instanceof DG.Column ? c : df.col(role) ?? df.col('~' + role);
  };
  // each group's animals, by ID
  const animalsByGroup = (df: DG.DataFrame): Map<string, string[]> => {
    const [animal, group] = [roleCol(df, 'animal'), roleCol(df, 'group')];
    const out = new Map<string, string[]>();
    for (let i = 0; animal && group && i < df.rowCount; i++) {
      if (!group.isNone(i) && !animal.isNone(i))
        out.set(String(group.get(i)), [...out.get(String(group.get(i))) ?? [], String(animal.get(i))]);
    }
    return out;
  };
  form.effect(() => {
    void result.value;
    const df = studyTable();
    if (df && !groupsBefore.has(df))
      groupsBefore.set(df, animalsByGroup(df));
  });
  const studyPlan = (): SavedPlan | null => {
    const text = planText();
    try {
      return text ? parsePlan(text) : null;
    } catch {
      return null;
    }
  };
  const calendarDays = dayActions(() => meta.peek(), setMeta, () => {
    const df = studyTable();
    const a = df ? roleCol(df, 'animal') : null;
    return a ? [...new Set(a.toList().map(String))].sort(naturalCompare) : [];
  }, {logGiven: true, last: () => monitorEvents.last, measured: () => monitorEvents.measureDays ?? [],
    name: (g) => revealed.peek() ? g : blindName(meta.peek())(g)});
  const calendarSummary = divV([], 'pmx-calendar-summary');
  Object.assign(calendarSummary.style, {gap: '4px', lineHeight: '1.45', marginBottom: 'var(--dg-space-m)'});
  const calendarBox = div([]);
  // the folded calendar's one line
  const calendarFolded = signal('');
  let monitorEvents: {last: number, events: [number, string][], end?: number,
    /** the days the calendar marks as measured or to measure (its day panel's Measure) */
    measureDays?: number[]} = {last: 0, events: []};
  form.effect(() => {
    const r = result.value;
    const m = meta.value;
    const df = studyTable();
    if (!r || !df) {
      calendarSummary.replaceChildren();
      calendarFolded.value = '';
      calendarBox.replaceChildren();
      return;
    }
    const dayCol = roleCol(df, 'day');
    const volCol = roleCol(df, 'volume');
    const measured = new Set<number>();
    for (let i = 0; i < df.rowCount; i++) {
      if (dayCol && volCol && !dayCol.isNone(i) && !volCol.isNone(i))
        measured.add(dayCol.get(i));
    }
    const latest = r.status.get('day', 0) as number;
    const plan = studyPlan();
    const startWeekday = day0Weekday(m.start, plan);
    const weekdays = new Set(weekdaysOf(plan?.inputs.measurementDays as string | undefined));
    // the plan's measurement days: day 0 and its weekdays (no plan: only the days measured count)
    const planned = (d: number) => plan !== null && (d === 0 || weekdays.has((startWeekday + d) % 7));
    const today = dayToday(m.start);
    dayState = {planned: plan ? planned : () => true, measured, plan: plan !== null};
    // always there: off the plan's calendar it asks first (entryStep)
    enterButton.hidden = false;
    const round = (x: number | undefined) => x === undefined ? undefined : Math.round(x);
    const [end, lo, hi] = [round(cell('end_median').value), round(cell('end_lower').value),
      round(cell('end_upper').value)];
    const planEnd = plan?.predicted.studyEnd ?? null;
    // the calendar's last day: the study's own span, and today while it is no more than two weeks past it (a study
    // opened again a year on is not a year of empty weeks, each a measurement missed); the likely end's range to a
    // week past the expected end at most (halfway through it can run weeks on, each an empty week to scroll)
    const span = Math.max(latest, end === undefined ? hi ?? 0 : Math.min(hi ?? end, end + 7), planEnd ?? 0);
    const last = Math.max(span, Math.min(today ?? 0, span + 14));
    // the day the study is on: today's, or without a date the last day measured
    const now = Math.min(today ?? latest, last);
    const looks = new Set(String(df.getTag(LOOKS_TAG) ?? '').split(/[,;\s]+/).filter((x) => x).map(Number));
    const missing = Array.from({length: Math.max(0, now)}, (_, d) => d).filter((d) => planned(d) && !measured.has(d));

    const lines: HTMLElement[] = [];
    const todayMark = mark('Today', CAL.blue, true);
    if (today === null) {
      lines.push(marked(todayMark, `Day 0 has no date yet: today is taken as the last day measured, day ` +
        `${latest}. Set the date from any day's menu (right click) and it follows today.`));
    } else if (measured.size && startMismatch(m.start, latest)) {
      // measurements past today's study day: the date given to day 0 is not this table's
      lines.push(marked(todayMark, `${startMismatch(m.start, latest)} Set day 0's date from any day's menu (right ` +
        'click).'));
    } else if (today < 0)
      lines.push(marked(todayMark, `Day 0 is in ${-today} days, on ${dateText(dateOf(m.start, 0)!, true)}.`));
    else {
      lines.push(marked(todayMark, `Day ${today}, ${dateText(dateOf(m.start, today)!, true)}: ${today} days in` +
        (end === undefined ? '.' : today <= end ? `, about ${end - today} to the expected end (day ${end}` +
          (lo !== undefined && hi !== undefined ? `, likely ${Math.max(lo, today)} to ${hi}).` : ').') :
          `, past the expected end (day ${end}).`)));
    }
    if (today !== null && planned(today) && !measured.has(today))
      lines.push(marked(mark('Measure today', CAL.blue, true), 'Today is a measurement day and has no values yet.'));
    lines.push(marked(mark('Measured', CAL.slate, true), `${measured.size} days with measurements, the last ` +
      `day ${latest}.`));
    if (missing.length > 0) {
      lines.push(marked(mark('No data', CAL.red, false), `The plan measures on day${missing.length > 1 ? 's' : ''} ` +
        `${missing.join(', ')}, with no values in the table.`));
    }
    // each dosed group: the doses its schedule gives (with the study's exceptions) up to today against the study's,
    // for the animals it has and their last body weights; once doses are logged, how many dosing days so far are
    if (revealed.value) {
      const groupCol = roleCol(df, 'group');
      const animalCol = roleCol(df, 'animal');
      const weightCol = df.col(df.getTag(WEIGHT_TAG) ?? '') ?? null;
      const excluded = df.col(EXCLUDED);
      const rows: WeightRow[] = [];
      for (let i = 0; groupCol && animalCol && i < df.rowCount; i++) {
        if (groupCol.isNone(i) || animalCol.isNone(i))
          continue;
        rows.push({group: String(groupCol.get(i)), animal: String(animalCol.get(i)),
          day: dayCol && !dayCol.isNone(i) ? Number(dayCol.get(i)) : -Infinity,
          weight: weightCol && !weightCol.isNone(i) ? Number(weightCol.get(i)) : null, excluded: ticked(excluded, i)});
      }
      for (const [group, {animals, kg}] of groupWeights(rows, weightCol?.meta.units ?? null)) {
        const d = armDoses(group, m, planEnd ?? span, animals, kg, Math.max(-1, now), false);
        if (d) {
          const line = marked(mark('Dose', CAL.green, false), `${group}: ${d.text}.`);
          line.title = `${d.tip} The doses so far are the schedule's up to today, less the skipped and with the ` +
            'added; log the doses given, skip or add one from a day\'s menu in the calendar.';
          lines.push(line);
        }
      }
    }
    // the days after the expected end stay on the calendar (the end is a forecast: the study may run on to its
    // likely range), but what falls on them is said to be conditional, here, on each day and in the .ics
    const runsOn = (d: number) => end !== undefined && d > end && d >= now && !measured.has(d);
    if (end !== undefined && last > end && Math.max(now, end + 1) <= last) {
      lines.push(marked(mark('If it runs on', CAL.amber, false), `Days ${Math.max(now, end + 1)}–${last} are past ` +
        `the expected end (day ${end}): their doses and measurements happen only if the study runs on.`));
    }
    calendarSummary.replaceChildren(...lines);
    calendarFolded.value = `${measured.size} days measured, the last day ${latest}` +
      (missing.length ? ` · ${missing.length} planned without data` : '') + (end === undefined ? '' :
      ` · expected end day ${end}`);

    const look = (d: number, own: HTMLElement[]): DayLook => {
      const marks: HTMLElement[] = [];
      if (d === today)
        marks.push(mark('Today', CAL.blue, true));
      if (d === 0)
        marks.push(mark('Randomize', CAL.blue, false));
      if (measured.has(d))
        marks.push(mark('Measured', CAL.slate, true));
      else if (planned(d)) {
        marks.push(d < now ? mark('No data', CAL.red, false) : d === today ? mark('Measure today', CAL.blue, true) :
          mark(runsOn(d) ? 'Measure if running' : 'Measure', CAL.slate, false));
      }
      if (looks.has(d))
        marks.push(mark('Review', CAL.dark, false));
      marks.push(...own);
      if (d === end)
        marks.push(mark('Expected end', CAL.dark, true));
      if (planEnd !== null && d === planEnd && d !== end)
        marks.push(mark('Planned end', CAL.dark, false));
      const detail = runsOn(d) && d !== today ? `past the expected end (day ${end}): its doses and measurements ` +
        'only if the study runs on' : d === today ? 'today' : measured.has(d) ? 'measured' : planned(d) ? (d < now ?
        'a measurement day with no values' : 'a measurement day') : undefined;
      // the days gone by shaded, so where the study stands reads at a glance; the likely end's days in amber
      return {marks, detail, border: d === today ? `2px solid ${CAL.blue}` : undefined,
        tint: d < now ? `color-mix(in srgb, ${CAL.slate} 16%, transparent)` :
          lo !== undefined && hi !== undefined && d >= lo && d <= hi ?
            `color-mix(in srgb, ${CAL.amber} 10%, transparent)` : undefined};
    };
    // from today (or the last day measured) on: what is still to do, for the .ics export
    monitorEvents = {last, events: Array.from({length: last - Math.max(0, now) + 1}, (_, i) => Math.max(0, now) + i)
      .flatMap((d): [number, string][] => [...planned(d) && !measured.has(d) ? [[d, 'Measure tumors'] as [number,
        string]] : [], ...looks.has(d) ? [[d, 'Planned review'] as [number, string]] : [],
      ...d === end ? [[d, 'Expected end of study'] as [number, string]] : []]), end,
    measureDays: Array.from({length: last + 1}, (_, d) => d).filter((d) => measured.has(d) || planned(d))};
    calendarBox.replaceChildren(weeksGrid({meta: unnamed(m, revealed.value), startWeekday, last, day: look,
      control: String(r.status.get('control', 0)), onMenu: calendarDays.menu,
      onNote: calendarDays.addNote, onSelect: calendarDays.select}));
  });

  // Enter today's values: a day typed in a dialog, written to the table only on Save day (day-entry.ts)
  const weightColOf = (df: DG.DataFrame) => roleCol(df, 'bodyWeight') ?? df.col(df.getTag(WEIGHT_TAG) ?? '') ?? null;
  // every animal of the study on `day`, in the measuring order (blinded: by animal, as Monitor lists them): one that
  // missed the last session says when it was last measured; one noted dead or removed before the day, or given a
  // removal reason in the table, is listed greyed and not counted
  const entryAnimals = (day: number): EntryAnimal[] => {
    const [r, df] = [result.peek(), studyTable()];
    if (!r || !df)
      return [];
    const {animal, day: dayCol, volume} = roleColumns(df);
    const weight = weightColOf(df);
    const last = new Map<string, {day: number, volume: number | null, weight: number | null}>();
    const [onDay, days] = [new Set<string>(), new Set<number>()];
    for (let i = 0; i < df.rowCount; i++) {
      if (dayCol.isNone(i) || animal.isNone(i))
        continue;
      const [a, dd] = [String(animal.get(i)), dayCol.get(i) as number];
      days.add(dd);
      if (dd === day)
        onDay.add(a);
      const was = last.get(a);
      if (dd >= day || (was && was.day > dd))
        continue;
      last.set(a, {day: dd, volume: volume.isNone(i) ? was?.volume ?? null : volume.get(i),
        weight: weight && !weight.isNone(i) ? Number(weight.get(i)) : was?.weight ?? null});
    }
    // the session before the day: an animal last measured before it missed it
    const previous = maxOf([...days].filter((x) => x < day));
    const gone = new Map<string, string>();
    for (const x of meta.peek().notes) {
      const m = /^Not measured: (dead|removed)/.exec(x.text);
      if (m && x.day < day)
        gone.set(x.animal, `Off the study: ${m[1]} on day ${x.day}`);
    }
    for (const [a, x] of removedAnimals(df, animal, dayCol)) {
      if (x.day < day)
        gone.set(a, `Off the study on day ${x.day}: ${x.why}`);
    }
    const nextDay = r.status.get('next_day', 0) as number;
    const order = [...new Set(r.days.col('group')!.toList() as string[])];
    const shown = revealed.peek();
    return rowsOf(r.animals).sort((a, b) => (shown ? order.indexOf(a.group) - order.indexOf(b.group) : 0) ||
      naturalCompare(String(a.animal), String(b.animal)))
      .map((a) => {
        const name = String(a.animal);
        const l = last.get(name);
        // the range: the next measuring day's; a day before it, after the animal's last value, between that value and
        // the next day's range (a Saturday's 720 mm³ for a 1,913 mm³ tumor is outside both)
        // any day after its last value has the range Needs attention shows (expectedRange; an off-calendar day 12
        // read "—" for every animal, QA 0.1.3 round 4 technician R4-6, chaos N10)
        const own = day > a.last_day ? expectedRange(a) : null;
        const around = !!own && day !== nextDay && l?.volume != null && a.last_day === l.day && day < nextDay;
        const [lower, upper] = !own ? [null, null] : around ?
          [Math.min(own[0], l!.volume!), Math.max(own[1], l!.volume!)] : own;
        const missedLast = l && Number.isFinite(previous) && l.day < previous;
        // counted removed early by Welfare (not measured on its group's last day) and not measured since: off the
        // study for the day's list too, greyed and not counted (QA 0.1.3 round 3: technician F22)
        if (!gone.has(name) && a.removed_early && day > a.last_day)
          gone.set(name, `Off the study: counted as removed early (not measured since day ${a.last_day})`);
        return {animal: name, group: String(a.group), lastDay: l?.day ?? null,
          lastVolume: l?.volume ?? null, lastWeight: l?.weight ?? null, lower, upper,
          expected: lower !== null ? `${mm3(lower)}–${mm3(upper)}` : '—',
          ...!onDay.has(name) && gone.has(name) ? {gone: gone.get(name)} : {},
          ...missedLast ? {note: `Last measured day ${l!.day}`} : {}};
      });
  };
  // the calendar's planned and measured days, as its last reading had them
  let dayState: {planned: (d: number) => boolean, measured: Set<number>, plan: boolean} =
    {planned: () => true, measured: new Set(), plan: false};
  // before any measurement: the animals the table lists with no volume yet (their day 0 is still to enter)
  const unmeasured = (df: DG.DataFrame): Map<string, string> => {
    const [animal, volume, group] = [roleCol(df, 'animal'), roleCol(df, 'volume'), roleCol(df, 'group')];
    const out = new Map<string, string>();
    if (!animal || !volume)
      return out;
    const measuredOnes = new Set<string>();
    for (let i = 0; i < df.rowCount; i++) {
      if (animal.isNone(i))
        continue;
      const a = String(animal.get(i)).trim();
      if (!volume.isNone(i))
        measuredOnes.add(a);
      else if (a && !out.has(a))
        out.set(a, group && !group.isNone(i) ? String(group.get(i)) : '');
    }
    measuredOnes.forEach((a) => out.delete(a));
    return out;
  };
  // a planned study with no animals in its table yet: their IDs, as many as the plan has, typed or pasted here (the
  // rows are written with their day-0 values, on Save day)
  const askAnimals = (df: DG.DataFrame, plan: SavedPlan) => {
    const n = randomizedGroups(plan).reduce((t, g) => t + g.n, 0);
    const d = modal('day0-animals', 'The study\'s animals');
    if (!d)
      return;
    const ids = ui.input.textArea('Animal IDs', {value: Array.from({length: n}, (_, i) =>
      `M${String(i + 1).padStart(String(n).length, '0')}`).join('\n')});
    const count = dialogLine('');
    d.add(dialogLine(`The plan randomizes ${n} animals. Paste one ID per line, as on the cages.`))
      .add(ids).add(count).onOK(() => {
        const {ids: list, twice} = animalIds(ids.value ?? '', n);
        if (list.length && !twice.length)
          enterDayZero(df, new Map(list.map((a) => [a, ''])));
      });
    const next = okNamed(d, 'Next');
    // the IDs as typed, counted against the plan's; one typed twice is named, and Next waits for it to be one
    // (QA 0.1.3 round 3: technician N5)
    const sync = () => {
      const {ids: list, twice, line} = animalIds(ids.value ?? '', n);
      count.textContent = line;
      count.style.color = twice.length ? CAL.red : list.length !== n ? CAL.amber : '';
      next.disabled = !list.length || twice.length > 0;
      next.title = !list.length ? 'Type the animals\' IDs first' : twice.length ? `Each ID once: ${twice.join(', ')} ` +
        `${twice.length === 1 ? 'appears' : 'appear'} more than once` : '';
    };
    ids.input.addEventListener('input', sync);
    sync();
    showCentred(d);
  };
  // day 0 (randomization) of a study with no measurement yet: its animals as the table lists them, else as given,
  // else as its draft kept them (the IDs typed in the animals' step, never the defaults again: QA 0.1.3 r5)
  const enterDayZero = (df: DG.DataFrame, given?: Map<string, string>) => {
    const study = rowOfTable(df)?.id ?? (metaOf(df)?.id || df.name);
    const listed = unmeasured(df);
    const drafted = new Map((loadDraft(study, 0)?.animals ?? []).map((a): [string, string] => [a, '']));
    const ids = listed.size ? listed : given ?? drafted;
    const [animal, day, volume] = [roleCol(df, 'animal'), roleCol(df, 'day'), roleCol(df, 'volume')];
    const plan = studyPlan();
    if (!ids.size && plan && animal && day && volume)
      return askAnimals(df, plan);
    if (!ids.size || !animal || !day || !volume) {
      return void grok.shell.warning('Day 0 needs the animals: add their IDs to the study table, a row each (or ' +
        'paste their day-0 volumes there), then enter the day here.');
    }
    const group = roleCol(df, 'group') ?? df.columns.addNewString(groupsShown(df) ? 'group' : '~group');
    const animals = [...ids].sort(([a], [b]) => naturalCompare(a, b)).map(([a, g]): EntryAnimal => ({animal: a,
      group: g, lastDay: null, lastVolume: null, lastWeight: null, lower: null, upper: null, expected: '—'}));
    // the day-0 weights are each animal's baseline: a study with no weight column yet gets one to type them into
    // (QA 0.1.3 round 3: technician N5)
    let weight = weightColOf(df);
    if (!weight) {
      weight = df.columns.addNewFloat('body_weight');
      weight.meta.units = 'g';
      weight.meta.friendlyName = WEIGHT_LABEL;
      df.setTag(WEIGHT_TAG, weight.name);
      // (the form's field too, as restoreFields gives it: the study is saved with the column the field names)
      const field = form.getInput('bodyWeight');
      if (field)
        field.value = weight.name;
    }
    enterDayDialog({df, cols: {animal, group, day, readings: readingColumns(df, volume, weight)},
      study, day: 0, animals: (d) => d === 0 ? animals : [],
      warning: (d) => d === 0 ? null : 'Before the first measurements, only day 0 (randomization) can be entered',
      date: (d) => longDate(d, null),
      showGroup: revealed.peek(), saved: () => {
        logVersion.value++;
        metaVersion.value++;
        auto.changed();
      }});
  };
  // a study day's date, as the day's entry says it beside the day (null: day 0 has no date, or one the data past
  // `latest` contradicts)
  const longDate = (d: number, latest: number | null): string | null => {
    const at = latest !== null && startMismatch(meta.peek().start, latest) ? null : dateOf(meta.peek().start, d);
    return at ? at.toLocaleDateString('en-GB', {weekday: 'short', day: 'numeric', month: 'short', year: 'numeric',
      timeZone: 'UTC'}) : null;
  };
  const enterDay = () => {
    const [r, df] = [result.peek(), studyTable()];
    if (!df)
      return void grok.shell.warning('Open a study first: the day\'s list is its animals');
    if (!r)
      return enterDayZero(df);
    const today = dayToday(meta.peek().start);
    const latest = r.status.get('day', 0) as number;
    const planEnd = studyPlan()?.predicted.studyEnd;
    // today's study day, unless day 0 has no date or its date is not this table's (measured past it)
    const day = today !== null && today >= 0 && today >= latest ? today : null;
    if (entryStep(day, dayState.plan ? dayState.planned : null) === 'confirm') {
      const d = modal('off-plan day', 'Enter today\'s values');
      d?.add(ui.divText(offPlanAsk(day!, dayState.planned))).onOK(() => openEntry(df, r, day, latest, planEnd)).show();
      return;
    }
    openEntry(df, r, day, latest, planEnd);
  };
  const openEntry = (df: DG.DataFrame, r: MonitorFrames, day: number | null, latest: number,
    planEnd: number | undefined) => {
    const c = roleColumns(df);
    enterDayDialog({df, cols: {animal: c.animal, group: c.group, day: c.day,
      readings: readingColumns(df, c.volume, weightColOf(df))},
    study: rowOfTable(df)?.id ?? (metaOf(df)?.id || df.name),
    day: day ?? r.status.get('next_day', 0), animals: entryAnimals,
    warning: (d) => dayState.measured.has(d) ? `Day ${d} already has values: you are editing them` :
      dayState.planned(d) ? null : offPlanAsk(d, dayState.planned).replace(': enter anyway?', ''),
    limits: {last: latest, end: Number.isFinite(planEnd) ? Math.round(planEnd!) : null, today: day},
    date: (d) => longDate(d, latest),
    showGroup: revealed.peek(), saved: () => {
      logVersion.value++;
      metaVersion.value++;
    }});
  };
  const enterButton = buttonWithIcon('Enter today\'s values', 'edit', enterDay);
  enterButton.classList.add('pmx-go');
  enterButton.title = 'Type the day\'s volumes animal by animal, each checked as you go; written to the study on ' +
    'Save day';
  enterButton.hidden = true;
  // before the first measurements: shown while the table lists animals with no volume yet (their day 0)
  form.effect(() => {
    void dataVersion.value;
    const df = studyTable();
    if (!result.value)
      enterButton.hidden = !(df && (unmeasured(df).size || (!df.rowCount && studyPlan())));
  });

  // Decisions
  const decisionBody = div([], 'pmx-monitor-decisions');
  const logBox = div([], 'pmx-monitor-log');
  const record = (what: string) => {
    const df = studyTable();
    if (df) {
      addLog(df, what);
      blindColumn(df);
      followGroup(df);
      logVersion.value++;
    }
  };
  // a stop against the suggestion (Continue): asked once more, with a short reason the log keeps
  const stopAgainst = (arm: string, review: number, chance: string) => {
    const d = modal('stop-against', `Stop ${arm}?`);
    if (!d)
      return;
    const why = ui.input.string('Reason', {nullable: false});
    why.setTooltip('Kept in the study\'s log with the decision');
    d.add(dialogLine('The suggestion is to continue', `Chance of a clear result: ${chance}. Stopping ends this ` +
      'group\'s doses after this review.')).add(why)
      .onOK(() => {
        if (!why.value?.trim())
          return;
        // one decision, one line
        record(decisionEntry(arm, review, 'stop', `against the suggestion to continue: ${why.value.trim()}`));
      });
    const ok = okNamed(d, 'Stop');
    const sync = () => {
      ok.disabled = !why.value?.trim();
      ok.title = ok.disabled ? 'Write the reason first: the log keeps it with the decision' : '';
    };
    why.input.addEventListener('input', sync);
    sync();
    showCentred(d);
  };
  // revealing is for good (the log keeps who and when): asked once more, then `then` (Analyze, which shows them)
  const reveal = (then?: () => void) => {
    const d = modal('reveal-groups', 'Reveal the groups?');
    if (!d)
      return;
    d.add(dialogLine('Everyone who opens this study sees its groups from now on',
      'The log records who revealed them and when; blinding cannot be restored for this study.'))
      .onOK(() => {
        record(REVEALED);
        then?.();
      });
    okNamed(d, 'Reveal');
    showCentred(d);
  };
  form.effect(() => {
    const r = result.value;
    void logVersion.value;
    // (no study analyzed: nothing to decide on, not the last study's table)
    if (!r)
      return void decisionBody.replaceChildren();
    // no reviews fixed (a study without a plan): nothing to decide on, a line saying where they come from
    if (!String(call.getParamValue('looks') ?? '').trim()) {
      decisionBody.replaceChildren(infoNote('No planned reviews: add them in <b>Planned reviews</b>, or follow a ' +
        'plan.'));
      return;
    }
    if (!revealed.value) {
      decisionBody.replaceChildren(infoNote('<b>The groups are hidden.</b> The chance of success tells how the ' +
        'groups ' +
        'differ: reveal them only when someone must decide.'), button('Reveal the groups', () => reveal()));
      return;
    }
    const d = rowsOf(r.decisions).filter((a) => a.planned);
    // a review on a day under way waits for its last group: a chance read on half the groups would be read again
    if (r.status.get('look_waits', 0)) {
      decisionBody.replaceChildren(infoNote(`A planned review falls on ${dayText(r.status.get('day', 0))}, ` +
        `which is not complete: ${esc(String(r.status.get('waiting_groups', 0)))} not measured yet. ` +
        'The chance of success shows once every group has its values.'));
      return;
    }
    if (!r.status.get('at_look', 0) || d.length === 0) {
      const next = r.status.get('next_look', 0);
      decisionBody.replaceChildren(infoNote(r.status.col('next_look')?.isNone(0) || !next ?
        'No planned review ahead. The chance of success is shown only on the reviews fixed at the start.' :
        `The next planned review is on ${dayText(next)}. The chance of success is shown only on the reviews fixed at ` +
        'the start: looking between them would make it easier to stop on chance.'));
      return;
    }
    const review = r.status.get('day', 0) as number;
    const log = studyTable() ? logOf(studyTable()!) : [];
    // a stop holds from the review it was taken on, whichever review this is, until undone; a continue is this
    // review's
    const stopped = stopsOf(log);
    // (the last word on this review stands: a continue undone is no longer one)
    const continued = (arm: string) => {
      const [yes, undone] = [decisionEntry(arm, review, 'continue'), decisionEntry(arm, review, 'continue undone')];
      return log.map((l) => l.what).filter((w) => w === yes || w === undone).pop() === yes;
    };
    fill(decisionBody, d, (a) => a.arm, [['Group', (a) => swatch(a.arm, colorsOf(r)[a.arm])],
      // only the low end as a figure (MONITOR_DESIGN.md §8); a clearly effective arm reads "High" (chanceWords)
      ['Chance of a clear result', (a) => chanceWords(a.success, a.tgi_lower)],
      // (a range of one value, as at the study's end, is that value: not "68% (68–68%)")
      ['TGI if it continues', (a) => Math.round(a.tgi_lower) === Math.round(a.tgi_upper) ?
        `${Math.round(a.tgi_median)}%` : `${Math.round(a.tgi_median)}% (${percentRange(a.tgi_lower, a.tgi_upper)})`],
      ['Suggests', (a) => badge(a.success < FUTILITY ? 'Consider stopping' : 'Continue',
        {variant: a.success < FUTILITY ? 'warning' : 'success'})],
      // once decided, the decision in place of the buttons (the log keeps it, with the review's day); a stop ends
      // the arm's doses after its review (the calendar's too), so it can be undone
      ['You decide', (a) => {
        const stop = stopped[a.arm];
        if (stop !== undefined) {
          return divH([span(`Stopped on day ${stop}, no more doses`),
            button('Undo', () => record(decisionEntry(a.arm, stop, 'stop undone')))]);
        }
        return continued(a.arm) ? divH([span('Continued'),
          button('Undo', () => record(decisionEntry(a.arm, review, 'continue undone')))]) : divH([
          button('Continue', () => record(decisionEntry(a.arm, review, 'continue'))),
          button('Stop', () => a.success < FUTILITY ? record(decisionEntry(a.arm, review, 'stop')) :
            stopAgainst(a.arm, review, chanceWords(a.success, a.tgi_lower)))]);
      }]]);
  });
  form.effect(() => {
    void logVersion.value;
    void result.value;
    const df = studyTable();
    const log = df ? logOf(df) : [];
    // while the groups are hidden, an entry names them as the calendar does ("Group k": a decision names its arm);
    // the longest names first, so a combination is not taken for one of its drugs
    const m = meta.value;
    const order = [...new Set([...m.groups.map((g) => g.group), ...result.value?.days.col('group')?.toList() ?? []])];
    const named = revealed.value ? [] : [...order].sort((a, b) => b.length - a.length);
    const what = (text: string) => named.reduce((t, g) => t.split(g).join(`Group ${order.indexOf(g) + 1}`), text);
    if (log.length === 0)
      logBox.replaceChildren(span('Nothing recorded yet.', 'pmx-muted'));
    else {
      // the time and who in narrow columns, the text left and wrapping in the room left
      const t = new BasicTable<LogEntry>({items: log, keyOf: (l) => l.when + l.what, columns: [
        {header: 'When (your time)', render: (l) => localTime(l.when), width: '130px', align: 'left'},
        {header: 'Who', render: (l) => l.who, width: '110px', align: 'left'},
        {header: 'What', render: (l) => logWhat(what(l.what)), align: 'left'}]});
      t.root.style.width = '100%';
      wrapCells(t.root);
      logBox.replaceChildren(t.root);
    }
  });
  // Welfare: Analyze's Tolerability on the weights so far; by group, so once revealed
  const welfareFrames = computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    const control = r.status.get('control', 0) as string;
    const today = r.status.get('day', 0) as number;
    const arms = r.tolerability.rowCount === 0 ? [] :
      [...new Set(r.tolerability.col('group')!.toList() as string[])].filter((g) => g !== control);
    return {...r, analysis: DG.DataFrame.fromColumns([DG.Column.fromStrings('control', [control])]),
      comparisons: DG.DataFrame.fromColumns([DG.Column.fromStrings('arm', arms)]),
      volumes: r.days.clone(DG.BitSet.create(r.days.rowCount, (i) => r.days.get('day', i) <= today))};
  });
  const welfare = tolerabilityTab(form, welfareFrames, colorsOf);
  const welfareBody = divV(welfare.content);
  const welfareHidden = divV([infoNote('<b>The groups are hidden.</b> Body weight by group tells how the groups ' +
    'differ. The animals to look at are under Today, without their group.'),
  button('Reveal the groups', () => reveal())]);
  welfareHidden.style.gap = 'var(--dg-space-m)';
  form.effect(() => {
    const has = (result.value?.tolerability.rowCount ?? 0) > 0;
    welfareHidden.hidden = !has || revealed.value;
    welfareBody.hidden = has && !revealed.value;
  });

  // what a study row had when opened or last saved: its measurements and the tags it is saved from (domain storage:
  // closing the view with anything else asks first)
  // (the table's text taken only once it changed: native.ts savedTables)
  const asSaved = savedTables(self, (df) => df.toCsv(), (df) => SAVED_TAGS.map((t) => df.getTag(t) ?? '').join('|'));
  const saveStudyNow = (): Promise<void> => {
    const df = studyTable();
    if (!df)
      return Promise.reject(new Error('No study table'));
    remember();
    // through setMeta: the header shows the ID the study has just been given
    const groups = roleCol(df, 'group');
    const names = groups ? groups.categories.filter((g) => g) : [];
    // a group renamed in the table since it was read or saved (the same animals under a new name) takes its details,
    // kept exposure and analysis's TGI with it
    const now = animalsByGroup(df);
    const renames = renamedGroups(groupsBefore.get(df) ?? now, now);
    groupsBefore.set(df, now);
    if (renames.size && df.getTag(EXPOSURE_TAG))
      df.setTag(EXPOSURE_TAG, JSON.stringify(renamedExposure(exposureOfTag(df), renames) ?? {}));
    const details = withGroups(metaOf(df) ?? emptyMeta(), names, renames);
    if (!domainStorageOn()) {
      setMeta(withId(details));
      const said = (p: string) => grok.shell.info(`Study saved: ${shownPath(p)}`);
      void Promise.resolve().then(() => saveStudy(savedStudy(df), df.name))
        .then(said, failed('The study could not be saved'));
      return Promise.resolve();
    }
    // a row: the one the table was opened from or saved to, at its version (another's save since is refused)
    return withCode(details).then(async (m) => {
      setMeta(m);
      const s = await stamped(savedStudy(df));
      const at = rowOfTable(df);
      try {
        const saved = await save('study', s, at, {fileOwner: tableFileOwner(df)});
        setRowOn(df, saved);
        linkRow(self.view, 'study', saved.id);
        head.refresh(); // its ID is now its own: no longer "Not saved yet"
        markTableFromFile(df, '');
        // (its ID was another study's: saved as a new one under a new ID)
        if (saved.meta) {
          setMeta(saved.meta);
          grok.shell.info(newRowNote('study', saved, m.id));
        }
      } catch (e) {
        // a study shared with the user to read: their version goes to a copy of their own, if they say so
        if (!at || !isReadOnly(e))
          throw e;
        return offerCopy(async () => {
          const {saved, doc} = await saveCopy('study', s, at);
          setRowOn(df, saved);
          linkRow(self.view, 'study', saved.id);
          df.name = doc.name;
          setMeta(doc.meta!);
          asSaved.set(df);
          grok.shell.info(`Saved as your copy: ${esc(doc.name)}`);
        }, e);
      }
      asSaved.set(df);
      grok.shell.info('Study saved');
    }).catch((e) => {
      // (kept until dismissed: a save refused, another's version saved since, is read and acted on, QA 0.1.3 F12)
      grok.shell.error(saveError(e, 'study'), {autoHide: false});
      throw e;
    });
  };
  const saveButton = buttonWithIcon('Save study', 'save', () => void saveStudyNow().catch(() => {}));
  // Save as…: the study under a name of its own, a new study with its own ID (the one opened stays as it was)
  const saveAs = () => {
    const df = studyTable();
    if (df)
      saveStudyAs(df, setMeta, () => saveStudyNow().then(() => head.refresh()));
  };
  if (domainStorageOn())
    exports.push(['Save as…', saveAs, 'File', () => studyTable() ? null : 'Open a study first']);
  if (domainStorageOn()) {
    // the table as first shown, once the view has read it (what it writes on its first reading is not a change)
    form.effect(() => {
      // (read once a run has ended, with a forecast or without: a followed study with nothing measured yet has none,
      // and its URL still names it, QA 0.1.3 designer N3)
      const read = !!result.value || !!error.value;
      if (!read)
        return;
      const df = studyTable();
      // (the URL names the study's row: a reload opens it again; a table not saved yet, none)
      const id = rowOfTable(df)?.id;
      // (the view is made once the content is: the first run waits for it)
      queueMicrotask(() => linkRow(self.view, 'study', id));
      // (any table, saved or not: a demo, a file or a plan just followed, its values entered since, asks too: QA 0.1.3
      // round 3 technician N1, N6, designer B3)
      if (df && !asSaved.has(df))
        asSaved.set(df);
    });
    queueMicrotask(() => {
      saveGate(self, 'study', saveStudyNow, () => {
        const now = studyTable();
        return !!now && asSaved.changed(now);
      }, undefined, asSaved.changes);
    });
  }
  saveButton.title = domainStorageOn() ? 'Save the plan, the measurements and the log in PMX Studio (yours until ' +
    'you share it)' : 'Save the plan, the measurements and the log to your files';
  let toAnalyze = false;
  // a day under way (groups not measured yet): Analyze on the last day every group has, when its own day (the
  // control at the endpoint, else the last day measured) would be the day under way; said, not done in silence
  const completeDay = (): number | undefined => {
    const r = result.peek();
    if (!r || !r.status.get('waiting', 0))
      return undefined;
    const complete = r.status.get('complete_day', 0) as number;
    try {
      if (evaluationDayOf(analysisDeclarationOf(withoutExcluded(callInputs<MonitorInputs>(call)))).day <= complete)
        return undefined;
    } catch {
      // not read: the last complete day all the same
    }
    grok.shell.info(`Day ${r.status.get('day', 0)} is not complete (${r.status.get('waiting_groups', 0)} not ` +
      `measured yet): Analyze reads day ${complete}, the last day every group was measured.`);
    return complete;
  };
  const analyzeThis = () => {
    const df = studyTable();
    // (a second click while Analyze is being built would open it twice)
    if (!df || toAnalyze)
      return;
    // Analyze shows the groups: that is revealing them, asked first
    if (!revealed.peek())
      return reveal(analyzeThis);
    toAnalyze = true;
    // in place of an Analyze view already open
    [...grok.shell.views].find((v) => stageOf(v) === 'analyze')?.close();
    // Analyze's Save (and Ctrl+S) saves this table as this view does, while it is this view's table
    tableSavers.set(df, () => studyTable() === df ? saveStudyNow() :
      Promise.reject(new Error('Open the study in Monitor to save it')));
    // (the same table on: its changes go with it, and this page stays open with them)
    void analyzeView({study: df, evaluationDay: completeDay()}).then((v) => showStage('analyze', v, false, true))
      .catch(failed('Analyze could not be opened'))
      .finally(() => toAnalyze = false);
  };
  const analyzeButton = buttonWithIcon('Analyze this study', 'chart-line', analyzeThis);
  analyzeButton.classList.add('pmx-go');
  analyzeButton.title = 'Analyze this table as planned; in a blinded study this reveals the groups, and the log ' +
    'says so';
  // (in the Study menu too: on a narrow window the button gives way first, QA 0.1.3 designer D15)
  exports.push(['Analyze this study', () => analyzeButton.click(), 'Study', () => reasonOf(analyzeButton)]);

  // the three limits an animal is checked against, as the last run read them: in the attention panel's (i)
  const TODAY_HINT = 'Who needs attention today; a person decides what to do';
  const todayHint = span(TODAY_HINT);
  form.effect(() => {
    if (!result.value)
      return;
    const v = (k: string) => call.getParamValue(k);
    todayHint.textContent = `${TODAY_HINT}. Limits: endpoint ${v('endpoint')} mm³ (median control), humane limit ` +
      `${v('humane')} mm³, weight loss ${v('weightLimit') ?? 20}%`;
  });
  // the study's data in Datagrok's own table view; here, one line on it and the way there (a folded panel: its line
  // is the panel's summary)
  const dataLine = studyDataLine(() => studyInputsOf(call));
  self.own(dataLine.dispose);
  const dataSummary = signal('');
  form.effect(() => {
    void result.value;
    void error.value;
    dataLine.refresh();
    dataSummary.value = dataLine.root.querySelector('span')?.textContent?.replace(/ Right-click.*$/, '') ?? '';
  });
  // a study with animals and volumes but no groups yet has no forecast, so no panels: forming the groups is here
  // (the step itself is the error line's, said once: DAY0_DONE)
  const firstGroups = divH([buttonWithIcon('Propose groups', 'random', proposeAsked), problemLine()]);
  Object.assign(firstGroups.style, {gap: 'var(--dg-space-m)', alignItems: 'center'});
  // before day 0: Enter today's values in the page itself, not only in the toolbox (QA 0.1.3 round 3: technician N9,
  // designer B8)
  const firstDay = divH([buttonWithIcon('Enter today\'s values', 'edit', enterDay),
    span('Day 0: each animal\'s tumor volume and body weight, then Propose groups.')]);
  firstDay.querySelector('button')?.classList.add('pmx-go');
  Object.assign(firstDay.style, {gap: 'var(--dg-space-m)', alignItems: 'center'});
  form.effect(() => {
    void dataVersion.value;
    void error.value;
    const df = studyTable();
    firstDay.hidden = !!result.value || !df || !(unmeasured(df).size || (!df.rowCount && studyPlan()));
    firstGroups.hidden = !!result.value || !error.value || !planText() || !firstDay.hidden;
  });
  const panels = panelTabs([
    {id: 'today', title: 'Needs attention', hint: todayHint, content: [enterButton, entered, attention, noteBox]},
    {id: 'data', title: 'Data checks', hint: 'The measurements in Datagrok\'s table view: values the outlier rules ' +
      'find, flagged and excluded. Right-click a value in the table to flag or exclude it', content: [dataLine.root],
    summary: dataSummary},
    {id: 'growth', title: 'Growth', hint: 'Each group so far and what is expected up to the end',
      content: [growthHidden, growthAxis.row, sized(growth, 340), bandLines]},
    {id: 'calendar', title: 'Study calendar', hint: 'The study by its days: today, the days measured and to measure, ' +
      'doses, reviews and the expected end. Double-click a day to note what happened, right-click it for more',
    content: [calendarSummary, calendarBox], summary: calendarFolded},
    {id: 'sheet', title: 'Measurement sheet', hint: 'The range the next value of each animal should fall in: check ' +
      'each value as you measure it', content: [sheet]},
    {id: 'decisions', title: 'Decisions', hint: 'Whether each treated group is still worth continuing, and what was ' +
      'decided', content: [decisionBody, block('Log', 'Who did what, and when', logBox)]},
    {id: 'welfare', title: 'Welfare', hint: 'Whether the treatments are being tolerated, from the body weights so far',
      content: [welfareHidden, welfareBody]},
    {id: 'enrollment', title: 'Enrollment', hint: 'Random groups alike in their starting tumors',
      content: [divH([proposeButton, sheetButton]), problemLine(), allocationTable]},
  ], MONITOR_TABS, 'pmx-monitor-panels', 'monitor', form);
  exports.push(['Customize this page', () => customizePage(panels.root), 'View'],
    // the study's first steps, at the top of the toolbox and the menu as well as on the Enrollment panel
    ['Enter today\'s values…', enterDay, 'Study', () => {
      const df = studyTable();
      return result.value || (df && (unmeasured(df).size || (!df.rowCount && studyPlan()))) ? null :
        'Open a study with its animals first';
    }],
    ['Propose groups…', proposeAsked, 'Study'], ['Print randomization sheet', printRandomization, 'Study']);
  // the study's colours: a colour changed in a chart is the study's, in every stage open on it
  const link: StudyLink = {meta: () => metaOf(studyTable()), setMeta, colors: () => {
    const r = result.peek();
    return r ? colorsOf(r) : {};
  }, table: studyTable, plan: () => planKeyOf(studyTable()),
  control: () => result.peek()?.status.get('control', 0) ?? null};
  self.own(linkStudy(link));
  exports.push(['Reset colors', () => resetColors(link), 'View',
    () => hasOwnColors(link) ? null : 'The colors are the automatic ones']);

  // a study whose plan was made with another atlas says so above its numbers (the plan the form gives, else the
  // table's); the latest answer stands
  const atlasLine = div([], 'ui-label-details');
  atlasLine.hidden = true;
  let atlasAsked = 0;
  form.effect(() => {
    void result.value;
    const asked = ++atlasAsked;
    void atlasNote(studyPlan()?.provenance, FROZEN_ATLAS_NOTE).then((t) => {
      if (asked === atlasAsked)
        atlasLine.hidden = !(atlasLine.textContent = t ?? '');
    });
  });
  // a model planned from its history: its twin is not validated for it (D694)
  const twinLine = markLine(TWIN_NOTE);
  form.effect(() => {
    void result.value;
    const df = studyTable();
    twinLine.hidden = !savedHistory(String(call.getParamValue('line') || (df ? lineOf(df) : '') || ''));
  });
  const content = divV([atlasLine, twinLine, cards, panels.root], 'pmx-monitor-content');
  const needStudy = () => studyTable() ?? (() => {
    throw new Error('Open a study first');
  })();
  exports.push(
    ['Report (PDF)…', async () => {
      const df = needStudy();
      // the report is what is on screen: not while that is behind the table or the form
      if (busy.peek() || stale.peek() || !result.peek())
        throw new Error('The forecasts are not up to date yet: export the report once they are');
      printSheet(await panels.all(() => reportHtml({title: df.name, stage: 'Monitor', meta: metaOf(df), content,
        groups: roleCol(df, 'group')?.categories.filter((g) => g),
        planAtlas: studyPlan() ? studyPlan()!.provenance?.atlas ?? null : undefined, log: logOf(df),
        lastDay: result.peek()!.status.get('day', 0) as number})));
    }],
    ['Study file (.pmxstudy)', () => {
      const df = needStudy();
      remember();
      void stamped(savedStudy(df))
        .then((s) => download(`${exportName(df)}.pmxstudy`, JSON.stringify(s), 'application/json'),
          failed('The study file could not be made'));
    }],
    ['Measurements (CSV)', () => {
      const df = needStudy();
      // (a table with no measurement yet would download its header alone)
      if (!df.rowCount)
        throw new Error('No measurements yet: enter a day first, then export them');
      download(`${exportName(df)}.csv`, csvOf(df), 'text/csv');
    }],
    ['Measurements', () => openStudyTable(needStudy()), 'Data > Open as table'],
    ['Calendar (.ics)', () => {
      const df = needStudy();
      download(`${exportName(df)}.ics`, icsOf(exportName(df), unnamed(meta.peek(), revealed.peek()), monitorEvents.last,
        monitorEvents.events, monitorEvents.end),
      'text/calendar');
    }]);
  Object.assign(content.style, {gap: 'var(--dg-space-l)'});
  form.effect(() => {
    content.hidden = !result.value;
    // a study open with nothing measured yet: its empty state says what to do, not this too
    start.hidden = !!result.value || empty.value;
  });
  const panes = stagePanes(form, content, {busy, result, error, stale, empty,
    before: [start, firstDay, firstGroups],
    busyText: 'Forecasting the study…', staleText: 'the forecasts shown are for the last complete inputs'});
  // a session's values pasted or typed into the table: the forecasts follow, after a pause (one value, then the next)
  self.own(onTableEdits(form, computed(() => [result.value, error.value]),
    studyTable, () => auto.changed(TABLE_EDITS)));
  // the table's own subscriptions are not built under the host: dropped with it
  self.own(stopWatching);
  const head = studyHeader(form, 'monitor', computed(() => {
    const r = result.value;
    // the details changed (an ID given on saving, a note, day 0's date): the header shows them as they are
    void metaVersion.value;
    void error.value;
    const df = studyTable();
    // a study followed from its plan and not measured yet has its header too: its name, model and where it stands
    // (QA 0.1.3 round 3: technician N9, designer B8)
    if (!r && !(df && planText()))
      return undefined;
    const colors = r ? colorsOf(r) : {};
    const day = cell('day').value;
    const line = r ? String(r.status.get('line', 0) ?? '') : df ? lineOf(df) ?? '' : '';
    const start = metaOf(df)?.start ?? null;
    // today's study day beside the last measured, when they differ ("today day 12 · last measured day 9")
    const today = dayToday(start);
    const where = !r ? hasVolumes(df) ? '' : 'day 0 not entered yet' : day === undefined ? '' :
      today !== null && today > day && !startMismatch(start, day) ? `today day ${today} · last measured day ${day}` :
        dayText(day);
    // day 0 dated after the data, said with its date (the calendar says how to set it) (QA 0.1.3 round 3: W13)
    const mismatch = r && startMismatch(start, day) ?
      `day 0 dated ${dateText(dateOf(start, 0)!, true)}, after the data allow (day ${day})` : '';
    return {name: df ? studyNameOf(df) : '', groups: revealed.value ?
      Object.entries(colors).map(([name, color]) => ({name, color})) : [],
    status: [line ? modelName(line) : '', where, r ? endStatus.value : '', revealed.value ? '' : 'groups hidden',
      mismatch].filter((x) => x).join(' · '),
    statusTip: () => r ? [endTip(), startMismatch(start, day) ?? ''].filter((x) => x).join(' ') : '',
    groupsHidden: !revealed.value, meta: () => metaOf(df),
    // with the store, an ID is a study's once saved (a plan's or a file's may be taken: the store gives another)
    unsaved: () => domainStorageOn() && !!df && !rowOfTable(df), setMeta: df ? (m) => {
      const was = !!metaOf(df)?.blinded;
      setMeta(m);
      if (m.blinded !== was)
        record(m.blinded ? 'Blinded entry on' : 'Blinded entry off: the groups are shown');
    } : undefined, groupObject: (g, color) => new PmxGroup(studyNamesOf(df), g, color)};
    // the header's Analyze step is this study's analysis (unless that would reveal a blinded study's groups: the
    // button does that, and says so); Plan opens with the study's plan, but not on a blinded study's hidden groups
  }), (s) => s === 'analyze' && studyTable() && revealed.peek() ? analyzeThis() : void openStage(s,
    revealed.peek() ? planOfStudy(planText(), metaOf(studyTable())?.id) : null), computed(() => ({analyze:
    !result.value && !canReturn('analyze') ? 'Needs a study with measurements: open one, or follow a plan and ' +
      'measure' : null})));
  // saving is always in reach, as in Plan (File › Save): a study's measurements, log and notes live in the table until
  // saved
  form.effect(() => analyzeButton.hidden = !result.value);
  const links = stageLinks({stage: 'monitor'}, () => head.details());
  // another stage shown in this view's place: this study, to come back to
  keepOnLeaving(self, async () => {
    const df = studyTable();
    return df ? () => monitorView(df) : null;
  });
  queueMicrotask(() => stageMenu(self.view, exports, {open: () => void openStudyDialog(undefined, {stage: 'monitor'}),
    save: {text: () => saveButton.textContent?.trim() || 'Save study', run: () => saveButton.click(),
      why: () => !studyTable() ? 'Open a study first' : reasonOf(saveButton)}}, links, (d) => self.own(d)));
  const toolbox = stageToolbox(exports, links.study, links.help, panels.root);
  // its "On this page" only while the page shows its sections (QA 0.1.3 round 3: technician N9, designer B8)
  form.effect(() => {
    const pane = toolbox.querySelector('.pmx-index')?.closest<HTMLElement>('.d4-accordion-pane');
    if (pane)
      pane.style.display = result.value ? '' : 'none';
  });
  return {
    root: takesFiles(stagePage(head.row, panes), {stage: 'monitor'}),
    ribbon: [[head.steps], [analyzeButton]],
    toolbox,
    path: '/monitor',
    status: computed(() => busy.value ? 'Forecasting the study…' :
      result.value ? modelName(String(status.value?.get('line', 0))) : ''),
  };
}

/** Whether `df` is the study followed from `plan`: by the plan's study ID (a plan and the study followed from it
 *  share one), never by its content, so a plan edited since still finds its study (QA 0.1.3 round 4 designer D1: an
 *  endpoint changed gave a new "(2)" study); a plan with no ID yet, by its name and inputs. */
export function followedFrom(df: DG.DataFrame, plan: SavedPlan): boolean {
  if (!df.getTag(PLAN_TAG))
    return false;
  const id = plan.study?.id;
  return id ? metaOf(df)?.id === id : planKeyOf(df) === planKey(plan);
}

/** Plan's "Follow this study": the study already followed from the plan (open, else its saved row), else a new
 *  measurements table carrying the plan (`fresh`), opened in Monitor (the plan's changes are asked about before:
 *  plan-view). The one way a plan goes to Monitor: Follow, and a plan file opened there (hub-view openBrought: a
 *  plan whose study is followed made a second, unsaved study of the same ID, QA 0.1.3 r4). */
export async function followPlan(plan: SavedPlan, fresh = () => newStudyTable(plan)): Promise<void> {
  const open = grok.shell.tables.find((t) => followedFrom(t, plan));
  if (open)
    return void await showStage('monitor', await monitorView(open));
  const id = plan.study?.id;
  const [row] = id && domainStorageOn() ? await studyRowsOf('study', id).catch(() => []) : [];
  if (row) {
    const v = await openRow('study', row.id);
    return void await showStage(placeOf(v) === 'analyze' ? 'analyze' : 'monitor', v);
  }
  await showStage('monitor', await monitorView(fresh()));
}

/** A .pmxstudy file opened from Files or a drop; `stage`: the one picked (else by whether it is finished). */
export async function openStudyFile(content: string, row: Saved | null = null,
  stage?: 'monitor' | 'analyze'): Promise<DG.DataFrame[]> {
  const s = parseStudy(content);
  // (a study row's document: its next save updates that row)
  const df = row ? tableOfRow(s, row, undefined, inUse) : tableOf(s);
  if (!row) // (a file's: its first save is the user's copy when its owner is someone else)
    markTableFromFile(df, ownerOf('study', s));
  // a finished study (saved from Analyze, the Overview's stage 3) opens in Analyze: in Monitor its last animals would
  // read as needing attention and its next measurement day as due
  if ((stage ?? (s.finished === true ? 'analyze' : 'monitor')) === 'analyze')
    await showStage('analyze', await analyzeView({study: df}));
  else
    await showStage('monitor', await monitorView(df));
  return [];
}

/** The measurements as CSV, without the group column while the groups are hidden. */
const csvOf = (df: DG.DataFrame): string => df.toCsv(groupsShown(df) ? undefined :
  {columns: df.columns.names().filter((c) => !c.startsWith('~'))});
