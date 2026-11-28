import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Scope, signal, untracked} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {inMm3} from './analyze';
import {BandChart} from './band-chart';
import {Bar, WaterfallChart} from './spider-chart';
import {PmxAnimal} from './pmx-objects';
import type {StudyInputs as StudyNames} from './explore';
import {animalFacts, animalReadings} from './study-table';
import {PanelFact, PanelSection, PanelVerdict, factList, objectPanel, setUserSetting, showDepth, sized,
  userSetting} from './ui';

const MEASURED = 'Measured';
const EXPECTED = 'Expected';
const WEIGHED = 'Body weight';
const OTHERS = 'Others';
const GREY = '#9aa5b1';
const FOLDS = 'animalPanelFolds';

/** What a stage knows of an animal beyond its readings (each: none where the stage has none). */
export interface AnimalView {
  /** its answer, and under it its effect bar (Analyze, a treated animal) */
  verdict?: PanelVerdict | null;
  bar?: HTMLElement | null;
  /** its curve: against its twin, or its own growth (against the forecast in Monitor); open by default */
  curve?: {title: string, frame: DG.DataFrame} | null;
  /** its group's animals by one value (`order` as the page's waterfall), and its place there in words */
  group?: {bars: Bar[], order: 'asc' | 'desc', summary: string} | null;
  /** the day it was removed early */
  removedDay?: number | null;
  /** the study's weight-loss limit (%) */
  weightLimit?: number | null;
  /** the stage's own facts, after the readings (best response, doubling, endpoint) */
  extra?: PanelFact[];
}

let live = 0;
/** The animal panels' charts alive now, across views. */
export const liveAnimalCharts = () => live;

/** 1st, 2nd, 3rd, 4th, 11th, 21st. */
export function ordinal(k: number): string {
  const s = ['th', 'st', 'nd', 'rd'];
  const v = k % 100;
  return `${k}${s[(v - 20) % 10] ?? s[v] ?? s[0]}`;
}

/** An animal's place among `bars`, the smallest value first (1: the smallest); 0: not among them. */
export const rankOf = (bars: Bar[], animal: string) => {
  const me = bars.find((b) => b.animal === animal);
  return me ? bars.filter((b) => b.value < me.value).length + 1 : 0;
};

/** An animal's place among `bars` in words, counted from the nearer end: "the smallest of 10", "3rd smallest of 10",
 *  "2nd largest of 10", "the largest of 10" (QA 0.1.3 A11: the largest read "10th smallest of 10"); '': not among
 *  them. */
export function placeWords(bars: Bar[], animal: string): string {
  const me = bars.find((b) => b.animal === animal);
  if (!me)
    return '';
  const n = bars.length;
  if (n === 1)
    return 'the only one';
  const low = bars.filter((b) => b.value < me.value).length + 1;
  const high = bars.filter((b) => b.value > me.value).length + 1;
  return low <= high ? `${low === 1 ? 'the' : ordinal(low)} smallest of ${n}` :
    `${high === 1 ? 'the' : ordinal(high)} largest of ${n}`;
}

/** An animal's measured volumes as the band chart draws them, and with `expected` the forecast's range on its next
 *  day from its last reading (Monitor's); undefined: nothing measured. */
export function growthFrame(days: number[], volumes: number[], color: string | undefined,
  expected?: {day: number, lower: number, median: number, upper: number} | null): DG.DataFrame | undefined {
  const rows = days.map((day, i) => ({day, series: MEASURED, lower: volumes[i], upper: volumes[i],
    median: volumes[i]}));
  if (!rows.length)
    return undefined;
  const last = rows[rows.length - 1];
  if (expected && expected.day > last.day)
    rows.push({...last, series: EXPECTED}, {series: EXPECTED, ...expected});
  const df = DG.DataFrame.fromObjects(rows)!;
  df.col('series')!.meta.colors.setCategorical({[MEASURED]: color ?? '#2083d5', [EXPECTED]: GREY});
  df.col('day')!.meta.units = 'd';
  inMm3(df, ['lower', 'upper', 'median']);
  return df;
}

/** An animal's body weight as a % change from its first weighing, its limit and the day it was removed as reference
 *  lines; its worst change in words, to the tenth as the tables give it, with the limit ("Lost 19.6% by day 9 (limit
 *  20%)": "−20% at most" beside a table's 19.6%, round 3 designer 4). Null: never weighed. */
export function weightFrame(weighed: [number, number][], color: string | undefined, limit?: number | null,
  removedDay?: number | null): {frame: DG.DataFrame, summary: string} | null {
  if (!weighed.length || !(weighed[0][1] > 0))
    return null;
  const rows = weighed.map(([day, w]) => {
    const change = Math.round(1000 * (w / weighed[0][1] - 1)) / 10;
    return {day, series: WEIGHED, lower: change, upper: change, median: change};
  });
  const df = DG.DataFrame.fromObjects(rows)!;
  df.col('series')!.meta.colors.setCategorical({[WEIGHED]: color ?? '#2083d5'});
  df.col('day')!.meta.units = 'd';
  for (const c of ['lower', 'upper', 'median']) {
    df.col(c)!.meta.units = '%';
    df.col(c)!.meta.format = '0.0';
  }
  const line = (title: string, formula: string) => df.meta.formulaLines.addLine({title, formula, color: GREY});
  if (limit)
    line(`−${limit}%, humane endpoint`, `\${median} = ${-limit}`);
  if (removedDay !== null && removedDay !== undefined)
    line(`Removed, day ${removedDay}`, `\${day} = ${removedDay}`);
  const worst = rows.reduce((a, b) => b.median < a.median ? b : a);
  // (a loss that rounds to 0% is none: "0% at most" read as a contradiction, QA 0.1.3)
  const loss = Math.abs(worst.median).toFixed(1);
  // the worst, then where it stands on the last day weighed when that is another ("Lost at most 2.1% (day 2); day 9:
  // +1%", designer wording 18)
  const last = rows[rows.length - 1];
  const now = Math.round(last.median);
  const since = last === worst ? '' : `; day ${last.day}: ${now > 0 ? '+' : now < 0 ? '−' : ''}${Math.abs(now)}%`;
  return {frame: df, summary: worst.median < 0 && Math.round(-worst.median) > 0 ? `Lost at most ${loss}% (day ` +
    `${worst.day}${limit ? `, limit ${limit}%` : ''})${since}` : 'No loss'};
}

type Which = 'curve' | 'group' | 'weight';

/** One table, whichever JS wrapper holds it. */
const sameTable = (a?: DG.DataFrame | null, b?: DG.DataFrame | null) => !!a && !!b && (a === b || a.dart === b.dart);

/** A view's animal panels: `show` (every entry point's) makes the animal the current object; its handler renders
 *  `render`'s panel and the Actions. */
export class AnimalPanels {
  private readonly curveData = signal<DG.DataFrame | undefined>(undefined);
  private readonly groupData = signal<{bars: Bar[], colors: {[g: string]: string}, title: string,
    order: 'asc' | 'desc', groups: string[]} | undefined>(undefined);
  private readonly weightData = signal<DG.DataFrame | undefined>(undefined);
  private scope: Scope | null = null;
  private charts: Partial<Record<Which, {root: HTMLElement}>> = {};
  // the last panel's sections, by title: their folds read when the next is built (no event tells a fold)
  private panes: [string, DG.AccordionPane][] = [];
  private readonly folds: {[title: string]: boolean};
  // the details of the animals this view made: another object current lets the charts go
  private readonly mine = new WeakSet<object>();
  // owns the subscription and, as its child, the charts
  private readonly owner = Scope.root('animal panels');
  // the study's data as last seen (`changed`), and as the panel shown was built: a change while something else was
  // current (the Correct dialog) builds it again once it is back (QA 0.1.3 round 3, analyst 6)
  private version = 0;
  private built = -1;
  // each panel's charts go to the latest panel built: an earlier one opening later does not take them (designer B4)
  private generation = 0;
  private disposed = false;

  constructor(private readonly o: {
    names: () => StudyNames | null,
    /** the animal's group colour; undefined: the groups are hidden (a blinded study) */
    color: (animal: string) => string | undefined,
    view: (animal: string) => AnimalView,
    /** the page's Linear / Log switch, which the curve follows */
    log?: ReadonlySignal<boolean>,
    /** what else showing an animal does: the page marks it */
    onShow?: (animal: string) => void,
    /** read in an effect: when it changes (another study, a correction, a renamed group), the shown panel follows */
    changed?: () => unknown}) {
    let kept: {[title: string]: boolean} = {};
    try {
      kept = JSON.parse(userSetting(FOLDS) ?? '{}') ?? {};
    } catch {
      // none kept
    }
    this.folds = kept;
    this.owner.subscribe(grok.events.onCurrentObjectChanged, () => {
      if (!this.shown())
        this.release();
      else if (this.built !== this.version)
        this.refresh();
    });
    const changed = o.changed;
    if (changed) {
      this.owner.effect(() => {
        changed();
        untracked(() => {
          this.version++;
          this.refresh();
        });
      });
    }
  }

  /** The animal shown, as its study is now: of another study, or gone from it (its group renamed away, the study
   *  replaced), the panel is cleared with its actions (Correct, Exclude on a study no longer open: QA 0.1.3 A1); else
   *  built again (a correction, a renamed group: T11, D18). */
  refresh(): void {
    const x = this.shown();
    if (!x)
      return;
    // (before the panel is built: the platform's change of object comes back here first)
    this.built = this.version;
    const n = this.o.names();
    const there = sameTable(n?.study, x.names?.study) && n!.study.col(n!.animal)?.toList()
      .some((a) => String(a) === x.animal);
    // (not frozen, forced: a follow, not a pick; `o =` freezes the panel for a second and keeps an equal object)
    grok.shell.setCurrentObject(there ? this.objectOf(x.animal) : null, false, true);
  }

  /** The animal as its object, its panel this view's. */
  objectOf(animal: string): PmxAnimal {
    const details = () => this.render(animal);
    this.mine.add(details);
    const names = this.o.names();
    // its actions only while its study is the one on screen: a panel left from another study changes nothing
    // (QA 0.1.3 round 3, analyst 1: a correction landed in a study no longer open)
    const live = () => !this.disposed && sameTable(this.o.names()?.study, names?.study);
    return new PmxAnimal(names, animal, {color: this.o.color(animal), details, live});
  }

  // this view's animal in the context panel; null: something else is there
  private shown(): PmxAnimal | null {
    const x = grok.shell.o;
    return x instanceof PmxAnimal && x.o.details && this.mine.has(x.o.details) ? x : null;
  }

  /** Every entry point's: the animal in the context panel, marked on the page. */
  show(animal: string): PmxAnimal {
    // (built from the data as it is now: not built again as it arrives)
    this.built = this.version;
    this.o.onShow?.(animal);
    return showDepth(this.objectOf(animal));
  }

  /** The animal's panel: the charts already made given its data. */
  render(animal: string): HTMLElement {
    this.keepFolds();
    this.built = this.version;
    const gen = ++this.generation;
    const chart = (which: Which) => gen === this.generation ? this.chart(which) : ui.div();
    const n = this.o.names();
    if (!n)
      return ui.divText(animal);
    const v = this.o.view(animal);
    const color = this.o.color(animal);
    const r = animalReadings(n, animal);
    const sections: (PanelSection & {open: boolean})[] = [];
    this.curveData.value = v.curve?.frame;
    if (v.curve)
      sections.push({title: v.curve.title, open: true, content: () => chart('curve')});
    this.groupData.value = v.group ? {title: '', order: v.group.order, groups: [OTHERS, r.group],
      bars: v.group.bars.map((b) => ({...b, group: b.animal === animal ? r.group : OTHERS})),
      colors: {[OTHERS]: '#c4cbd3', [r.group]: color ?? GREY}} : undefined;
    if (v.group) {
      sections.push({title: 'Among its group', summary: v.group.summary, open: false,
        content: () => chart('group')});
    }
    const w = n.bodyWeight ? weightFrame(r.weighed, color, v.weightLimit, v.removedDay) : null;
    this.weightData.value = w?.frame;
    if (w)
      sections.push({title: 'Body weight', summary: w.summary, open: false, content: () => chart('weight')});
    sections.push({title: 'Details', open: false, content: () => factList(animalFacts(n, animal,
      {color, extra: v.extra, weightLimit: v.weightLimit}))});
    for (const s of sections)
      s.open = this.folds[s.title] ?? s.open;
    return objectPanel({title: animal, color, subtitle: color && r.group ? r.group : undefined,
      verdict: v.verdict, chart: v.bar ?? undefined, facts: [], sections,
      accordion: (acc) => this.panes = acc.panes.map((p, i) => [sections[i].title, p])});
  }

  /** The charts let go (the next animal makes them again); the folds kept. */
  release(): void {
    this.keepFolds();
    this.panes = [];
    live -= Object.keys(this.charts).length;
    this.charts = {};
    this.scope?.dispose();
    this.scope = null;
  }

  dispose(): void {
    this.disposed = true;
    // the view gone (another study opened in its place): its animal leaves the context panel with it
    if (this.shown())
      grok.shell.o = null;
    this.release();
    this.owner.dispose();
  }

  private chart(which: Which): HTMLElement {
    let c = this.charts[which];
    if (!c) {
      this.scope ??= this.owner.child('animal panel charts');
      const band = (data: ReadonlySignal<DG.DataFrame | undefined>, yTitle: string, log?: ReadonlySignal<boolean>) =>
        new BandChart(data, {x: 'day', lower: 'lower', upper: 'upper', line: 'median', split: 'series',
          xTitle: 'Day after randomization', yTitle, log, markers: true, legend: which === 'curve'});
      c = this.scope.build(() => which === 'curve' ? band(this.curveData, 'Tumor volume', this.o.log) :
        which === 'weight' ? band(this.weightData, 'Weight change') :
          new WaterfallChart(this.groupData, undefined, undefined, {compact: true}));
      c.root.classList.add(`pmx-animal-${which}`);
      this.charts[which] = c;
      live++;
    }
    return sized(c, which === 'group' ? 80 : which === 'curve' ? 180 : 140);
  }

  private keepFolds(): void {
    let changed = false;
    for (const [title, pane] of this.panes) {
      const open = pane.expanded;
      if (this.folds[title] !== open) {
        this.folds[title] = open;
        changed = true;
      }
    }
    if (changed)
      setUserSetting(FOLDS, JSON.stringify(this.folds));
  }
}
