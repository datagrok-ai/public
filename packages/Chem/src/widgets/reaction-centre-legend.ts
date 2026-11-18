/**
 * The reaction centre's legend (the reaction-centre colours, green made, red broken, blue order changed, mean nothing
 * to a reader without one): one short line under a reaction picture that colours its centre, a swatch and a few words
 * for each colour the picture has, as Crux counts the bonds it colours (`ChemReactionScheme.get_centre_counts`).
 * Nothing at all under a picture that colours none (an unmapped reaction). The Reaction panel and its full-screen view
 * put it right under their large picture, centred on it ({@link pictureWithLegend}); a route's steps, coloured as the
 * route is, share it. The Reaction Sketcher shows the panel. Grid cells have none.
 */
import * as ui from 'datagrok-api/ui';
import {cruxCentreCounts} from '../rendering/crux-reaction-drawing';

/** The colours as Crux draws them (crux-core scheme-draw layout.rs CENTRE_MADE, CENTRE_BROKEN, CENTRE_CHANGED), in
 * the order of its counts: made, broken, order changed. */
export const CENTRE_COLOURS = ['#009E4C', '#DB2828', '#1E6BEA'];
/** The words: the line's name, then each colour's. */
export const LEGEND_TITLE = 'Reaction centre:';
export const LEGEND_WORDS = ['bond made', 'bond broken', 'bond order changed'];
/** What the line's tooltip says. */
export const LEGEND_TOOLTIP = 'Coloured from the atom map numbers of each step: the bonds the step makes, breaks, or ' +
  'changes the order of. A reaction without map numbers is not coloured.';
/** The legend's class, and the data attribute that lists the colours shown (the tests read them). */
export const LEGEND_CLASS = 'chem-reaction-centre-legend';
/** The swatch's side and the gaps (CSS px), the words' size. */
const SWATCH_PX = 10;
const GAP_PX = 4;
const ITEM_GAP_PX = 12;
const FONT_SIZE = '11px';

/** What the pictures of `reactions` colour together, as Crux counts it: per colour, whether any has it. Nothing while
 * Crux is not loaded, or for a text it reads no reaction from. */
function coloured(reactions: string[]): boolean[] {
  const any = [false, false, false];
  for (const r of reactions)
    (cruxCentreCounts(r) ?? []).forEach((n, i) => any[i] = any[i] || n > 0);
  return any;
}

/** A legend for the pictures of `reactions` (one picture, or a route's steps): a swatch and words for each colour they
 * show, or an empty, hidden element when they show none. {@link setCentreLegend} fills it again. */
export function centreLegend(reactions: string[]): HTMLElement {
  const legend = ui.div([], LEGEND_CLASS);
  ui.tooltip.bind(legend, LEGEND_TOOLTIP);
  setCentreLegend(legend, reactions);
  return legend;
}

/** Fills `legend` for the pictures of `reactions` (a picture swapped for its mapped one, the Reaction panel's Map). */
export function setCentreLegend(legend: HTMLElement, reactions: string[]): void {
  const shown = coloured(reactions);
  ui.empty(legend);
  legend.dataset.colours = shown.map((s, i) => s ? CENTRE_COLOURS[i] : '').filter((c) => c).join(' ');
  legend.style.display = shown.some((s) => s) ? 'flex' : 'none';
  if (!shown.some((s) => s))
    return;
  Object.assign(legend.style, {flexWrap: 'wrap', alignItems: 'center', justifyContent: 'center',
    columnGap: `${ITEM_GAP_PX}px`, rowGap: `${GAP_PX}px`, fontSize: FONT_SIZE, color: 'var(--grey-5)'});
  legend.append(ui.span([LEGEND_TITLE]));
  shown.forEach((s, i) => {
    if (!s)
      return;
    const swatch = ui.div([], {style: {width: `${SWATCH_PX}px`, height: `${SWATCH_PX}px`,
      background: CENTRE_COLOURS[i], borderRadius: '2px', flex: 'none'}});
    legend.append(ui.divH([swatch, ui.span([LEGEND_WORDS[i]])], {style: {alignItems: 'center', gap: `${GAP_PX}px`}}));
  });
}

/** The room between a picture and its legend right under it (CSS px). */
export const PICTURE_LEGEND_GAP = 8;

/** A reaction picture and its legend as one group, the legend right under the picture, centred on it: the Reaction
 * panel's large picture and its full-screen view's (a route's legend sat at the panel's end, under its steps and
 * left-aligned, and at the end of the full-screen view's list of steps, while a single reaction's sat under its
 * picture). `picture` may be replaced in place (swapped for its mapped one). */
export function pictureWithLegend(picture: HTMLElement, legend: HTMLElement): HTMLElement {
  legend.style.flex = 'none';
  const group = ui.div([picture, legend]);
  Object.assign(group.style, {display: 'flex', flexDirection: 'column', alignItems: 'center',
    gap: `${PICTURE_LEGEND_GAP}px`});
  return group;
}
