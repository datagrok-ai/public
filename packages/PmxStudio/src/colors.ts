/* A study's colours, the same in every chart, table and chip of Plan, Monitor and Analyze: the control in a slate,
   each treated group in the platform's palette at the place it was given when it first appeared (StudyMeta.slots;
   without one, the free places by its name's place among the names sorted), unless the user gave the group another
   colour. A colour changed in any viewer of a study (its column's categorical colours) becomes the
   study's (StudyMeta.colors, saved with the plan and the study) and every open view of the study follows. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {signal} from '@datagrok-libraries/u2';
import type {StudyMeta} from './study-meta';
import {sameJson, maxOf, minOf} from './format';
import {whileAppOpen} from './app-open';

export type GroupColors = {[group: string]: string};

export const CONTROL_COLOR = '#5b6b7f';
/** The atlas's predicted control, in a lighter grey than the study's own (Analyze's Predicted controls, Compare's
 *  atlas). */
export const PREDICTED_COLOR = '#9aa5b1';

/** A colour with next to no hue (its channels within 24 of each other): it reads as the control's slate. */
const greyish = (argb: number) => {
  const c = [(argb >> 16) & 255, (argb >> 8) & 255, argb & 255];
  return maxOf(c) - minOf(c) < 24;
};

/** The platform's categorical colour `k` (0-based) among its coloured ones: its greys are skipped, so 15 treated arms
 *  take 15 hues (the default palette's 20, less 2 greys, keeps them apart) and none is mistaken for the control. */
export const armColor = (k: number) => {
  const hues = DG.Color.categoricalPalette.filter((c) => !greyish(c));
  return DG.Color.toHtml(hues.length ? hues[k % hues.length] : DG.Color.getCategoricalColor(k));
};

export type Slots = {[group: string]: number};

const byName = (a: string, b: string) =>
  a.localeCompare(b, 'en', {numeric: true, sensitivity: 'base'}) || (a < b ? -1 : a > b ? 1 : 0);

/** `slots` with each of `arms` placed: a group keeps its place; the rest take the free places in their names' order.
 *  Not the order listed: Plan, Monitor and Analyze list a study's groups in different orders, and a group keeps one
 *  colour in all three. */
function placed(arms: string[], slots: Slots): Slots {
  const out: Slots = {};
  for (const a of arms) {
    if (a in slots)
      out[a] = slots[a];
  }
  const used = new Set(Object.values(out));
  let k = 0;
  for (const a of [...new Set(arms)].filter((x) => !(x in out)).sort(byName)) {
    while (used.has(k))
      k++;
    used.add(k);
    out[a] = k;
  }
  return out;
}

/** Each treated arm's automatic colour, from its place (`slots`, StudyMeta.slots). The map keeps `arms`' order. */
export function treatedColors(arms: string[], slots: Slots = {}): GroupColors {
  const at = placed(arms, slots);
  return Object.fromEntries(arms.map((a) => [a, armColor(at[a])]));
}

/** Plan's arms' places after its arms went from `before` to `names`: an arm keeps its place for life (one added,
 *  removed or moved recolours no other), an arm renamed where it stands keeps its own, a new one takes a free place;
 *  the places of arms no longer there are given back. */
export function slotsAfter(slots: Slots, before: string[], names: string[]): Slots {
  const kept: Slots = {};
  for (const n of names) {
    if (n in slots)
      kept[n] = slots[n];
  }
  names.forEach((n, i) => {
    const old = before[i];
    if (!(n in kept) && old !== undefined && old in slots && !names.includes(old) &&
      !Object.values(kept).includes(slots[old]))
      kept[n] = slots[old];
  });
  return placed(names, kept);
}

/** Bumped when a study's colours change: whatever reads armColors or groupColors in a signal's scope draws again. */
export const colorTick = signal(0);

/** The colours of Plan's arms ('Control' and the treated `arms`): the user's (`own`, StudyMeta.colors), else the
 *  automatic ones from their places (`slots`, StudyMeta.slots). The one place a group gets its colour. */
export function armColors(arms: string[], own: GroupColors = {}, slots: Slots = {}): GroupColors {
  void colorTick.value;
  const colors: GroupColors = {'Control': CONTROL_COLOR, ...treatedColors(arms, slots)};
  for (const g of Object.keys(colors))
    colors[g] = own[g] ?? colors[g];
  return colors;
}

/** A study's groups in their colours (armColors), the control named `control` in the control's, whether or not
 *  `groups` lists it. */
export function groupColors(groups: string[], control: string, own: GroupColors = {}, slots: Slots = {}):
GroupColors {
  const colors = armColors(groups.filter((g) => g !== control), own, slots);
  colors[control] = own[control] ?? colors['Control'];
  return colors;
}

/** An open view's study: its details as kept, how to keep them, and its groups in their colours now. `table`: the
 *  measurements table it reads (Monitor and Analyze of one study share it). */
export interface StudyLink {
  meta(): StudyMeta | null;
  setMeta(m: StudyMeta): void;
  colors(): GroupColors;
  table?(): DG.DataFrame | null;
  /** The plan it follows or shows, as one text for all stages (Plan names its control 'Control', the study may not:
   *  a plan opened from a study is the same study even before either has an ID). */
  plan?(): string | null;
  /** The study's control group: its colour is kept as the control's ('Control'), the one Plan's control arm reads. */
  control?(): string | null;
}

const links = new Set<StudyLink>();

/** Follows `link` while its view is open: returns what stops it. A view opened on a study already open elsewhere
 *  takes the colours shown there. */
export function linkStudy(link: StudyLink): () => void {
  const open = [...links].find((l) => sameStudy(link, l) && hasOwnColors(l));
  const m = link.meta();
  const theirs = open?.meta()?.colors;
  if (m && theirs && !sameJson(theirs, m.colors ?? {})) {
    link.setMeta({...m, colors: {...theirs}});
    colorTick.value++;
  }
  links.add(link);
  return () => void links.delete(link);
}

function sameStudy(a: StudyLink, b: StudyLink): boolean {
  if (a === b)
    return true;
  const id = a.meta()?.id;
  if (id && id === b.meta()?.id)
    return true;
  const [pa, pb] = [a.plan?.(), b.plan?.()];
  if (pa && pa === pb)
    return true;
  const [ta, tb] = [a.table?.(), b.table?.()];
  return !!ta && !!tb && ta.dart === tb.dart;
}

/** `link`'s own colours with `group` in `to`: the control's under 'Control', the key every stage reads it by. */
function ownWith(link: StudyLink, own: GroupColors, group: string, to: string): GroupColors {
  const next = {...own};
  const key = group === link.control?.() ? 'Control' : group;
  if (key !== group)
    delete next[group];
  next[key] = to;
  return next;
}

/** The colours of the study read from `df`, as its open stages show them (empty: no stage open on it). */
export function colorsOfTable(df: DG.DataFrame): GroupColors {
  return [...links].find((l) => l.table?.()?.dart === df.dart)?.colors() ?? {};
}

export const hasOwnColors = (link: StudyLink) => Object.keys(link.meta()?.colors ?? {}).length > 0;

/** The study's own colours set to `own` (empty: the automatic ones), in every open view of it. */
export function setColors(link: StudyLink, own: GroupColors): void {
  for (const l of [...links, link].filter((x, i, all) => all.indexOf(x) === i && sameStudy(link, x))) {
    const m = l.meta();
    if (m)
      l.setMeta({...m, colors: {...own}});
  }
  colorTick.value++;
}

export const resetColors = (link: StudyLink) => setColors(link, {});

/** A colour as the platform writes it in a column's categorical colours (an ARGB number or a hex), as #rrggbb. */
function norm(c: unknown): string {
  if (typeof c === 'number')
    return DG.Color.toHtml(c).toLowerCase();
  // (a short hex, #888, reads back as null: kept as written)
  return typeof c === 'string' ? (DG.Color.toHtml(DG.Color.fromHtml(c)) ?? c).toLowerCase() : '';
}

const parse = (s: string | null | undefined): {[k: string]: unknown} => {
  try {
    const x = JSON.parse(s || '{}');
    return x && typeof x === 'object' ? x : {};
  } catch {
    return {};
  }
};

/** A column coloured by group: what was painted on it (to tell the user's change from ours) and, for a column whose
 *  categories only name a group ('Power, Cpd-A'), the group each names. */
interface Painted {df: DG.DataFrame, cols: Map<string, {was: string, groupOf: {[category: string]: string}}>,
  sub: {unsubscribe(): void}}
const painted = new Map<unknown, Painted>();
/** ponytail: the last frames painted are watched (a chart's data is built again on each change); a frame painted
 *  longer ago than that stops following the user's changes */
const WATCHED = 200;
// grid header edits come through this event, not the table's (while a page of the app is open: its grids)
whileAppOpen(() => {
  const sub = grok.events.onEvent('d4-grid-color-coding-changed').subscribe(() => painted.forEach(check));
  return () => sub.unsubscribe();
});

/** Colours `df`'s column `col` by `colors` (from armColors / groupColors), and watches it: a colour the user then
 *  changes in a viewer becomes the study's. `groupOf`: the group each category names, when not itself. */
export function paint(df: DG.DataFrame, col: string, colors: GroupColors, groupOf: {[category: string]: string} = {}):
  void {
  const c = df.col(col);
  if (!c)
    return;
  let p = painted.get(df.dart);
  if (!p) {
    p = {df, cols: new Map(), sub: df.onMetadataChanged.subscribe(() => check(painted.get(df.dart)))};
    painted.set(df.dart, p);
    if (painted.size > WATCHED) {
      const [oldest, q] = painted.entries().next().value!;
      q.sub.unsubscribe();
      painted.delete(oldest);
    }
  }
  p.cols.set(col, {was: JSON.stringify(colors), groupOf});
  c.meta.colors.setCategorical(colors);
}

/** A study's frames (`dfs`) in its colours now: each column whose categorical colours name one of its groups. */
export function follow(link: StudyLink, dfs: DG.DataFrame[]): void {
  const now = link.colors();
  for (const df of dfs) {
    for (const c of df.columns.toList()) {
      if (c.meta.colors.getType() !== DG.COLOR_CODING_TYPE.CATEGORICAL)
        continue;
      const cur = parse(c.getTag(DG.TAGS.COLOR_CODING_CATEGORICAL));
      const mine = Object.keys(cur).filter((k) => k in now);
      if (mine.length)
        paint(df, c.name, {...cur as GroupColors, ...Object.fromEntries(mine.map((k) => [k, now[k]]))});
    }
  }
}

/** A watched frame's columns read again: a group's colour the user changed goes to the study whose colour it had. */
function check(p: Painted | undefined): void {
  if (!p)
    return;
  for (const [name, x] of p.cols) {
    const tag = p.df.col(name)?.getTag(DG.TAGS.COLOR_CODING_CATEGORICAL) ?? '';
    if (tag === x.was)
      continue;
    const [before, after] = [parse(x.was), parse(tag)];
    x.was = tag;
    const changed = new Map<StudyLink, GroupColors>();
    for (const [k, v] of Object.entries(after)) {
      const [from, to] = [norm(before[k]), norm(v)];
      if (!from || from === to)
        continue;
      const g = x.groupOf[k] ?? k;
      const link = linkOfColor(g, from);
      if (link)
        changed.set(link, ownWith(link, changed.get(link) ?? {...link.meta()?.colors}, g, to));
    }
    changed.forEach((own, link) => setColors(link, own));
  }
}

/** The open study whose group `group` is drawn in `color` (ponytail: two open studies with a group of one name and
 *  colour: the first). */
const linkOfColor = (group: string, color: string, df?: DG.DataFrame | null) => {
  const drawing = [...links].filter((l) => norm(l.colors()[group]) === norm(color));
  // two open studies with a group of one name and colour: the one whose table it is
  return (df ? drawing.find((l) => l.table?.()?.dart === df.dart) : undefined) ?? drawing[0];
};

/** What sets the group `group`, drawn in `color`, to another colour in its study (`df`'s, when given); null: no
 *  open study draws it so. */
export function recolorOf(group: string, color: string, df?: DG.DataFrame | null): ((to: string) => void) | null {
  const link = linkOfColor(group, color, df);
  return link ? (to) => setColors(link, ownWith(link, {...link.meta()?.colors}, group, norm(to))) : null;
}

/** Adds Change color… to `menu` when an open study draws `group` in `color`: the platform's colour picker, the
 *  colour picked then the study's. Whether it was added. */
export function addColorItem(menu: DG.Menu, group: string, color: string): boolean {
  const recolor = recolorOf(group, color);
  if (!recolor)
    return false;
  menu.item('Change color…', () => {
    let picked = color;
    ui.showColorPicker(DG.Color.fromHtml(color), (c) => picked = DG.Color.toHtml(c), () => recolor(picked));
  });
  return true;
}
