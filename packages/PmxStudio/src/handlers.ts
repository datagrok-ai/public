import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as echarts from 'echarts/core';
import {divV} from '@datagrok-libraries/u2';

import type {StudyInputs as StudyNames} from './explore';
import {PmxAnimal, PmxDay, PmxFigure, PmxGroup, PmxSavedStudy} from './pmx-objects';
import {animalMarks, animalPanel, markRows} from './study-table';
import {recolorOf} from './colors';
import {plural} from './wide-table';
import type {MenuAt} from './calendar';
import {download, failed, pngOf} from './report';
import {helpPage, objectPanel} from './ui';

// the platform's way (ObjectHandler + registerParamFunc): one list of actions for the right click and the Actions pane

export {PmxAnimal, PmxGroup, PmxFigure, PmxSavedStudy};

/** `apart`: starts a group of its own in the Actions pane (a line above it); `menuOnly`: in the right-click menu only
 *  (copying a name is not what the panel is for). */
type Action<T> = {name: string, run: (x: T) => void, when?: (x: T) => boolean, apart?: boolean, menuOnly?: boolean,
  /** the link's words in the panel, when they depend on the object (the name is the registered function's) */
  label?: (x: T) => string};

/** A study's rows of one animal or one group, and the words a mark is logged with. */
function rowsOf(n: StudyNames, role: 'animal' | 'group', value: string): {rows: number[], what: string} {
  const col = n.study.col(n[role])!;
  const rows = Array.from({length: n.study.rowCount}, (_, i) => i).filter((i) => String(col.get(i)) === value);
  return {rows, what: `${role} ${value} (${plural(rows.length, 'row')})`};
}

const inStudy = (n: StudyNames | null, role: 'animal' | 'group'): n is StudyNames => !!n?.study.col(n[role]);
const markOf = (x: PmxAnimal) => inStudy(x.names, 'animal') ? animalMarks(x.names).get(x.animal) : undefined;
/** `x` shown again in the context panel once marked, when it was there as the action started (as a new object: the
 *  platform renders a new one), so its actions follow its mark: an animal just excluded offers Include. Read before:
 *  while the mark's dialog is open, the dialog is the current object. */
const again = <T extends PmxAnimal | PmxGroup>(x: T, copy: () => T) => {
  const o = grok.shell.o;
  const shown = o === x || (o?.constructor === x.constructor && String(o) === String(x));
  return () => {
    if (shown)
      grok.shell.o = copy();
  };
};
const markRole = (role: 'animal' | 'group', value: string, n: StudyNames, action: 'flag' | 'exclude' | 'clear',
  done?: () => void) => {
  const {rows, what} = rowsOf(n, role, value);
  markRows(n.study, rows, action, what, done);
};
const markAnimal = (x: PmxAnimal, action: 'flag' | 'exclude' | 'clear') => markRole('animal', x.animal, x.names!,
  action, again(x, () => new PmxAnimal(x.names, x.animal, x.o)));
const markGroup = (x: PmxGroup, action: 'flag' | 'exclude' | 'clear') => markRole('group', x.group, x.names!,
  action, again(x, () => new PmxGroup(x.names, x.group, x.color, x.o)));
const selectRows = (role: 'animal' | 'group', value: string, n: StudyNames) => {
  const these = new Set(rowsOf(n, role, value).rows);
  n.study.selection.init((i) => these.has(i));
};
const copyName = (name: string) => void navigator.clipboard.writeText(name);

/** An animal's actions, in the menu's order. */
export const ANIMAL_ACTIONS: Action<PmxAnimal>[] = [
  {name: 'Show against its twin', run: (x) => x.o.twin!(), when: (x) => !!x.o.twin},
  {name: 'Flag this animal…', run: (x) => markAnimal(x, 'flag'), apart: true,
    when: (x) => inStudy(x.names, 'animal')},
  {name: 'Exclude from the analysis…', run: (x) => markAnimal(x, 'exclude'),
    when: (x) => inStudy(x.names, 'animal') && markOf(x) !== 'excluded'},
  {name: 'Include in the analysis', run: (x) => markAnimal(x, 'clear'),
    when: (x) => markOf(x) === 'excluded' || markOf(x) === 'partly'},
  {name: 'Clear the flag', run: (x) => markAnimal(x, 'clear'),
    when: (x) => markOf(x) === 'flagged'},
  {name: 'Select the animal’s rows', run: (x) => selectRows('animal', x.animal, x.names!), apart: true,
    when: (x) => inStudy(x.names, 'animal')},
  {name: 'Copy the animal’s name', run: (x) => copyName(x.animal), menuOnly: true},
];

/** A group's actions. */
export const GROUP_ACTIONS: Action<PmxGroup>[] = [
  {name: 'Flag every animal of the group…', run: (x) => markGroup(x, 'flag'),
    when: (x) => inStudy(x.names, 'group')},
  {name: 'Exclude the group from the analysis…', run: (x) => markGroup(x, 'exclude'),
    when: (x) => inStudy(x.names, 'group')},
  {name: 'Clear its marks', run: (x) => markGroup(x, 'clear'),
    when: (x) => inStudy(x.names, 'group')},
  // (named apart from the animal's: the platform keys its functions by name)
  {name: 'Select the group’s rows', run: (x) => selectRows('group', x.group, x.names!), apart: true,
    when: (x) => inStudy(x.names, 'group')},
  {name: 'Change color…', run: (x) => {
    const recolor = recolorOf(x.group, x.color, x.names?.study)!;
    let picked = x.color;
    ui.showColorPicker(DG.Color.fromHtml(x.color), (c) => picked = DG.Color.toHtml(c), () => recolor(picked));
  }, when: (x) => !!recolorOf(x.group, x.color, x.names?.study)},
  {name: 'Copy the group’s name', run: (x) => copyName(x.group), menuOnly: true},
];

const chartOf = (x: PmxFigure) => x.o.chart ? echarts.getInstanceByDom(x.o.chart) ?? null : null;

/** A figure's actions. */
export const FIGURE_ACTIONS: Action<PmxFigure>[] = [
  {name: 'Pin these figures to compare', run: (x) => x.o.pin!.set(true), when: (x) => !!x.o.pin && !x.o.pin.on()},
  {name: 'Unpin: compare with the previous result', run: (x) => x.o.pin!.set(false),
    when: (x) => !!x.o.pin?.on()},
  {name: 'Save as PNG', run: (x) => download(`${x.name}.png`, pngOf(chartOf(x)!)), when: (x) => !!chartOf(x),
    apart: true},
  {name: 'Copy as image', run: (x) => void fetch(pngOf(chartOf(x)!)).then((r) => r.blob())
    .then((b) => navigator.clipboard.write([new ClipboardItem({'image/png': b})])).catch(failed),
  when: (x) => !!chartOf(x)},
];

/** The details, then an Actions section with every action that applies, as the right click lists them (the platform
 *  adds no Actions pane for an object of a package's own type). */
export function withActions<T>(details: HTMLElement, x: T, actions: Action<T>[]): HTMLElement {
  const these = actions.filter((a) => !a.menuOnly && (!a.when || a.when(x)));
  if (!these.length)
    return details;
  // a quiet line between the groups of actions (marking apart from the utilities), none above the first
  const items = these.flatMap((a, i) => [...a.apart && i > 0 ? [ui.div([], 'pmx-actions-sep')] : [],
    ui.link(a.label?.(x) ?? a.name, () => a.run(x))]);
  const acc = ui.accordion('pmx-actions');
  acc.addPane('Actions', () => ui.divV(items, {style: {gap: '4px'}}), true);
  return divV([details, acc.root]);
}

class AnimalHandler extends DG.ObjectHandler<PmxAnimal> {
  get type() {
    return 'pmx-animal';
  }
  isApplicable(x: any) {
    return x instanceof PmxAnimal;
  }
  // the app's help in the help panel while an animal is the current object (none: it stayed "Help is loading…")
  get helpUrl() {
    return helpPage();
  }
  getCaption(x: PmxAnimal) {
    return x.animal;
  }
  renderTooltip(_x: PmxAnimal) {
    // the element's own tooltip is the one shown: the platform's beside it doubled it (checked in the browser)
    return null as unknown as HTMLElement;
  }
  renderProperties(x: PmxAnimal) {
    const details = x.o.details?.() ?? (x.names ? animalPanel(x.names, x.animal, {color: x.o.color,
      showGroup: true}) : ui.divText(x.animal));
    return withActions(details, x, ANIMAL_ACTIONS);
  }
}

class GroupHandler extends DG.ObjectHandler<PmxGroup> {
  get type() {
    return 'pmx-group';
  }
  isApplicable(x: any) {
    return x instanceof PmxGroup;
  }
  get helpUrl() {
    return helpPage();
  }
  getCaption(x: PmxGroup) {
    return x.group;
  }
  renderTooltip(_x: PmxGroup) {
    // the element's own tooltip is the one shown: the platform's beside it doubled it (checked in the browser)
    return null as unknown as HTMLElement;
  }
  renderProperties(x: PmxGroup) {
    return withActions(x.o.details?.() ?? ui.divText(x.group), x, GROUP_ACTIONS);
  }
}

class FigureHandler extends DG.ObjectHandler<PmxFigure> {
  get type() {
    return 'pmx-figure';
  }
  isApplicable(x: any) {
    return x instanceof PmxFigure;
  }
  get helpUrl() {
    return helpPage();
  }
  getCaption(x: PmxFigure) {
    return x.name;
  }
  renderTooltip(_x: PmxFigure) {
    // the element's own tooltip is the one shown: the platform's beside it doubled it (checked in the browser)
    return null as unknown as HTMLElement;
  }
  renderProperties(x: PmxFigure) {
    return withActions(figurePanel(x), x, FIGURE_ACTIONS);
  }
}

/** A saved study's actions (its panel only: its row's click opens it, there is no right click on the list). */
export const SAVED_STUDY_ACTIONS: Action<PmxSavedStudy>[] = [
  {name: 'Open', label: (x) => `Open in ${x.o.stage}`, run: (x) => x.o.open()},
  {name: 'Compare with the picked studies', run: (x) => x.o.compare()?.(), when: (x) => !!x.o.compare()},
];

class SavedStudyHandler extends DG.ObjectHandler<PmxSavedStudy> {
  get type() {
    return 'pmx-saved-study';
  }
  isApplicable(x: any) {
    return x instanceof PmxSavedStudy;
  }
  get helpUrl() {
    return helpPage();
  }
  getCaption(x: PmxSavedStudy) {
    return x.name;
  }
  renderTooltip(_x: PmxSavedStudy) {
    return null as unknown as HTMLElement;
  }
  renderProperties(x: PmxSavedStudy) {
    return withActions(x.o.details(), x, SAVED_STUDY_ACTIONS);
  }
}

class DayHandler extends DG.ObjectHandler<PmxDay> {
  get type() {
    return 'pmx-day';
  }
  isApplicable(x: any) {
    return x instanceof PmxDay;
  }
  get helpUrl() {
    return helpPage();
  }
  getCaption(x: PmxDay) {
    return x.name;
  }
  renderTooltip(_x: PmxDay) {
    return null as unknown as HTMLElement;
  }
  renderProperties(x: PmxDay) {
    return x.details();
  }
}

/** A figure's panel: its name, whether it is pinned, and in one sentence what pinning does. */
export function figurePanel(x: PmxFigure): HTMLElement {
  const on = x.o.pin?.on();
  const el = objectPanel({title: x.name, subtitle: x.o.pin ? on ? 'Pinned' : 'Not pinned' : undefined, facts: []});
  if (x.o.pin) {
    el.append(ui.divText(on ? 'New results are drawn beside these.' :
      'Pin them, and each new result is drawn beside these, so you can see what a change does.'));
  }
  return el;
}

let registered = false;

/** The three handlers and their actions, registered once (the first object bound or shown registers them). Each
 *  action is a function of its own, applicable only when it applies to that object (an excluded animal offers
 *  Include, not Exclude). */
export function registerHandlers(): void {
  if (registered)
    return;
  registered = true;
  const add = <T>(h: DG.ObjectHandler<T>, actions: Action<T>[]) => {
    DG.ObjectHandler.register(h);
    for (const a of actions) {
      grok.functions.registerParamFunc(a.name, h.type, (x: T) => a.run(x),
        (x: any) => h.isApplicable(x) && (!a.when || a.when(x)));
    }
  };
  add(new AnimalHandler(), ANIMAL_ACTIONS);
  add(new GroupHandler(), GROUP_ACTIONS);
  add(new FigureHandler(), FIGURE_ACTIONS);
  // (its actions in its panel only: no right click on the Overview's list)
  DG.ObjectHandler.register(new SavedStudyHandler());
  DG.ObjectHandler.register(new DayHandler());
  // an animal, group or figure is its view's: another view made current lets it go (the platform's ribbon named it
  // over the Overview's stages)
  // a stage's object (an animal, a group, a figure) left behind: not the current object of another view; on the
  // Overview and Compare, which have none of their own, nothing is (a study table's cell or row from Analyze either)
  grok.events.onCurrentViewChanged.subscribe(() => {
    // (Compare: stages.ts STANDS_ALONE, read here as its data attribute; stages.ts imports the views)
    const d = grok.shell.v?.root?.dataset;
    if (grok.shell.o != null && (isPmxObject(grok.shell.o) || d?.pmxStage === 'home' || d?.pmxAlone))
      grok.shell.o = null;
  });
}

export const isPmxObject = (x: unknown) => x instanceof PmxAnimal || x instanceof PmxGroup || x instanceof PmxFigure ||
  x instanceof PmxSavedStudy || x instanceof PmxDay;

// registered as the module loads: every view reaches it (report.ts), Plan's figures too
registerHandlers();

/** `el` stands for `x`: its right click offers x's actions, as the platform's objects do. */
export function bindObject<T extends HTMLElement>(x: PmxAnimal | PmxGroup | PmxFigure, el: T): T {
  registerHandlers();
  ui.bind(x, el, {contextMenu: true});
  return el;
}

/** The actions that apply to `x`, in their order. */
export function actionsOf(x: PmxAnimal | PmxGroup | PmxFigure): Action<any>[] {
  const all: Action<any>[] = x instanceof PmxAnimal ? ANIMAL_ACTIONS : x instanceof PmxGroup ? GROUP_ACTIONS :
    FIGURE_ACTIONS;
  return all.filter((a) => !a.when || a.when(x));
}

/** A menu of `x`'s actions where nothing can be bound (a chart's bar, a list built elsewhere); `first`: items of the
 *  place before them (a legend's colour). */
export function popupFor(x: PmxAnimal | PmxGroup | PmxFigure, at: MenuAt, first?: (m: DG.Menu) => void): void {
  const menu = DG.Menu.popup();
  first?.(menu);
  for (const a of actionsOf(x))
    menu.item(a.name, () => a.run(x));
  menu.show(at);
}
