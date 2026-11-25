/**
 * The Export panel for a reaction: a dialog that previews the reaction at print size and saves what it previews.
 * - Two presets:
 *   - "ACS Document 1996": the publication look, black, bonds 0.508 cm, labels Arial 10 pt; its units are points;
 *   - "Screen": the cells' look, Chem's molecule draw options; its units are CSS px, 96 an inch.
 * - Widths:
 *   - a journal's single column (8.3 cm) or double column (17.8 cm). The scheme keeps its print size and takes as many
 *     rows as the width needs (the drawing's `column` option, crux-core spike scheme-column). A Chem drawing engine
 *     without that option fits the scheme into the width on up to three rows, smaller when it must;
 *   - free: the scheme's natural width, on one row.
 *   "Fit each step on one line" (a column breaks a route at its arrows, so no step of the demo route was on one line):
 *   each step drawn as its own reaction on a row of its own, its number left of it, every step at one scale, the
 *   largest at which the widest fits the width (the column's; free: print size), so the bonds are one length across the
 *   figure; the intermediate ends one row and starts the next.
 *   Under half of print size (a column too), a warning under the note says the text's size, too small to read in print,
 *   and the wider width to choose.
 *   "As drawn", offered with a drawing (the Reaction Sketcher's canvas, a cell's ChemDraw page
 *   or RXN block with coordinates): the reaction as it is drawn, not arranged, the sketcher's own picture at the
 *   style's print size.
 * - Formats:
 *   - SVG, sized in pt or px so that a document places it at print size;
 *   - PNG at 300 or 600 dpi, its pHYs chunk saying so: a single column at 300 dpi is 980 px wide;
 *   - CDXML, the drawing's page as Crux Sketch's reaction writer writes it (rendering/cdxml-reaction-cells.ts
 *     reactionCdxml), saved as reaction.cdxml; said "not yet" while the vendored sketcher writes no reaction as CDXML;
 *   - Copy as picture: the PNG on the clipboard.
 * - The background white or transparent.
 * - The reaction centre's colours on or off, off for ACS 1996. A drawing exported with them carries their legend,
 *   drawn into the picture under the scheme.
 * - The map numbers on or off: off for ACS 1996; for Screen as the caller shows them, else off;
 *   a preset chosen sets them again. The box governs CDXML too: a mapped reaction's page keeps its atom mapping
 *   (ChemDraw's ReactionStepAtomMap) with the box on and leaves it out with the box off, and the note says which.
 * - No PDF.
 *
 * Opened from the Reaction panel (its Export icon), the Reaction Sketcher's ribbon (Export...) and a reaction cell's
 * menu ("Export reaction...", in both modes: the reaction features are Crux's whatever the Chemical Engine; RDKit mode
 * loads Chem's Crux engine for it).
 */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChemWasm, cruxLoadFailed, cruxModuleSync, noteCruxError} from '../engine/crux-module';
import {CHEM_CENTRE, StepTexts, chemDrawOptions, cruxCentreCounts, cruxDrawsStepTexts, cruxDrawsSteps,
  cruxForReactions, stepTextsOfStep, tooLargeToDraw} from '../rendering/crux-reaction-drawing';
import {CRUX_OFF, TOO_LARGE, inkRows, withResolution} from '../utils/reactions/reaction-document';
import {CENTRE_COLOURS, LEGEND_TITLE, LEGEND_WORDS} from './reaction-centre-legend';
import {DrawnPicture, DrawnPictureOptions, cdxmlReactionSmiles, cellPicture, isCdxmlText, moleculesOnNoStep,
  noStepWords, reactionCdxml} from '../rendering/cdxml-reaction-cells';
import {reactionsHelpUrl} from '../utils/reactions/reactions-help';

type Scheme = InstanceType<ChemWasm['ChemReactionScheme']>;

/** A preset: its name, the drawing's look, its units (an inch's worth, and the SVG's unit), and whether the reaction
 * centre's colours are on when it is chosen. */
export interface ExportPreset {
  label: string;
  look: 'publication' | 'host';
  unitsPerInch: number;
  unit: 'pt' | 'px';
  centre: boolean;
}

export const ACS_1996 = 'ACS Document 1996';
export const SCREEN = 'Screen';
export const EXPORT_PRESETS: ExportPreset[] = [
  {label: ACS_1996, look: 'publication', unitsPerInch: 72, unit: 'pt', centre: false},
  {label: SCREEN, look: 'host', unitsPerInch: 96, unit: 'px', centre: true},
];

/** The widths: a label and its centimetres; `null`, the scheme's own width. */
export const SINGLE_COLUMN = 'single column (8.3 cm)';
export const DOUBLE_COLUMN = 'double column (17.8 cm)';
export const FREE_WIDTH = 'free';
export const EXPORT_WIDTHS: {label: string, cm: number | null}[] = [
  {label: SINGLE_COLUMN, cm: 8.3}, {label: DOUBLE_COLUMN, cm: 17.8}, {label: FREE_WIDTH, cm: null},
];

/** The map numbers' box, and what the note says of a mapped reaction's CDXML page, which the box governs too (one
 * rule for every format; the page keeps the mapping as ChemDraw's ReactionStepAtomMap). */
export const MAP_NUMBERS = 'Map numbers';
/** What the Map numbers and Reaction centre colours boxes say, off, of a reaction without map numbers. */
export const NOT_MAPPED_TOOLTIP = 'The reaction has no atom map numbers: no centre to colour, no numbers to show';
export const CDXML_WITH_MAPPING = 'CDXML with the atom mapping (ChemDraw\'s ReactionStepAtomMap)';
export const CDXML_WITHOUT_MAPPING = 'CDXML without the atom mapping';

/** A CDXML page's reaction steps' atom maps (their ReactionStepAtomMap attribute). */
const STEP_ATOM_MAP = /\s+ReactionStepAtomMap\s*=\s*"[^"]*"/g;

/** Whether CDXML `page` maps any atom of its reaction steps. */
export function cdxmlMapped(page: string): boolean {
  return /ReactionStepAtomMap\s*=\s*"\s*\d/.test(page);
}

/** CDXML `page` as Save writes it: as written with `keep`, else its steps' atom maps left out. */
export function cdxmlMapping(page: string, keep: boolean): string {
  return keep ? page : page.replace(STEP_ATOM_MAP, '');
}

/** Whether the map numbers are on when preset `preset` is chosen: off for ACS 1996 (a journal's scheme); for Screen
 * (the cells' look) as the caller shows them (`caller`), off when it does not say. */
export function presetMapNumbers(preset: string, caller?: 'show' | 'hide'): boolean {
  return presetOf(preset).look === 'host' && caller === 'show';
}

/** The formats Save writes; Copy as picture is a button of its own (the PNG at the chosen dpi, 300 for SVG). */
export const SVG = 'SVG';
export const PNG_300 = 'PNG 300 dpi';
export const PNG_600 = 'PNG 600 dpi';
export const CDXML = 'CDXML';
export const EXPORT_FORMATS = [SVG, PNG_300, PNG_600, CDXML];
export const COPY_AS_PICTURE = 'Copy as picture';
/** What the panel says of CDXML: being written, not written by this version, written; and the file it saves. */
export const CDXML_NOT_YET = 'CDXML: Chem cannot write a reaction as CDXML yet';
export const CDXML_WRITING = 'Writing the CDXML…';
export const CDXML_WRITTEN = 'CDXML: the drawing as a page ChemDraw opens (its style and width are set in ChemDraw)';
export const CDXML_FILE = 'reaction.cdxml';
/** What the panel says when the drawing has texts (conditions, captions) that the SVG and PNG, drawn from the reaction,
 * leave out and its CDXML keeps (Export once dropped the text over the arrow without a word; it names the formats
 * that leave them out). */
export const TEXTS_CDXML_ONLY = 'Texts over the arrow are kept in CDXML only; SVG and PNG leave them out';

/** What one export is: a preset, a width, the background, the colours, the map numbers (left out: the drawing's own
 * choice, as the cells draw them), and whether each step is fitted on a line of its own ({@link FIT_STEPS}). */
export interface ExportChoice {
  preset: string;
  width: string;
  transparent: boolean;
  centre: boolean;
  mapNumbers?: 'show' | 'hide';
  fitSteps?: boolean;
  /** The reaction as it is drawn ({@link AS_DRAWN}), in place of the arrangements. */
  asDrawn?: boolean;
}

/** The choice that exports the reaction as it is drawn on the canvas, not arranged, offered
 * when the caller has a drawing ({@link ExportOptions.asDrawn}: the Reaction Sketcher's canvas, a cell's ChemDraw page
 * or RXN block with coordinates). Its picture is the sketcher's (crux-sketch's picture API): at the style's
 * print size (ACS 1996's 14.4 pt bond), in its colours (ACS 1996 black, Screen in element colours), the map numbers
 * as the box says; the width and the arrangements do not apply, nor the reaction centre's colours, which the canvas
 * does not draw. */
export const AS_DRAWN = 'As drawn';
const AS_DRAWN_TOOLTIP = 'The reaction as it is drawn on the canvas, not arranged: at the style\'s print size, as ' +
  'wide as it is drawn';
/** What the note says while the drawing's picture is made, and when it is not made. */
export const DRAWN_DRAWING = 'Drawing the picture…';
export const DRAWN_NOT_MADE = 'The drawing could not be pictured';
/** What the note says of a picture as drawn by a sketcher that takes no options (a crux-sketch before its picture API):
 * the canvas's own size, colours and map numbers. */
export const DRAWN_OWN_LOOK = 'at the canvas\'s own size and colours';

/** The choice that puts each step on a line of its own, and what it says of itself. */
export const FIT_STEPS = 'Fit each step on one line';
const FIT_STEPS_TOOLTIP = 'Each step drawn on a row of its own, its number left of it, all at one size: the largest ' +
  'at which the widest step fits the width';
/** A figure drawn under this percent of print size: a warning under the note says its text's size, too small to read
 * in print ("drawn at 37% of print size; its text 3.7 pt", grey in the note, was read as no warning). */
const SMALL_STEP_PERCENT = 50;
/** The warning's class (the tests find it by it). */
export const NOTE_WARNING_CLASS = 'chem-reaction-export-warning';
/** The text of a step at print size in each preset (units): ACS 1996's labels and agent names are 10 pt Arial; the
 * Screen look's labels are 8.7 CSS px (an H 6.2 px tall, crux-core 3330cc1e6). */
const TEXT_SIZE: {[preset: string]: number} = {[ACS_1996]: 10, [SCREEN]: 8.7};
/** Each look's bond length at print size (units: ACS 1996's 0.508 cm is 14.4 pt; the cells' look, Crux's own 22 px).
 * The steps of {@link FIT_STEPS} are drawn at a fraction of it through Crux's `bondLength`: a box smaller than a
 * step's natural size would shrink its bonds more than the box (its texts and margins shrink less), unevenly from
 * step to step. */
const PRINT_BOND: {[look: string]: number} = {publication: 14.4, host: 22};
/** The least bond, in a look's units, at which Crux draws map numbers: a scheme drawn smaller leaves them out, as they
 * would cover the molecules (crux-core scheme-draw MAP_NUMBERS_MIN_BOND_PX, 10). In the Export panel's units that is
 * 70% of ACS 1996's print size, 46% of Screen's (the default route in a single column, at 37%, showed none with Map
 * numbers on, and nothing said why). */
const MAP_NUMBERS_MIN_BOND = 10;
/** The least scale a step is drawn at; the passes that find the largest scale at which the widest step fits, and how
 * close they get (a half percent). */
const MIN_STEP_SCALE = 0.05;
const STEP_SCALE_PASSES = 8;
const STEP_SCALE_CLOSE = 0.005;
/** The blank between two steps' rows, and between a step's number and its row (in the number's font size). */
const STEP_ROW_GAP = 0.6;
const STEP_NUMBER_GAP = 0.6;
/** A step's number: black, as the route's arrow numbers; its SVG class (the tests read the numbers by it). */
const STEP_NUMBER_COLOUR = '#000000';
export const STEP_NUMBER_CLASS = 'chem-reaction-export-step';

/** A step on its row ({@link FIT_STEPS}): its index in the route, its number, its box (units: its natural size at
 * `bond`), the bond length it is drawn with (units), and its size against print size (percent, 100 at most: the
 * steps' one scale). */
export interface StepRow {
  step: number;
  label: string;
  x: number;
  y: number;
  width: number;
  height: number;
  bond: number;
  percent: number;
}

/** The dialog's title, its classes (the tests find the preview and its note by them) and the file names. */
export const EXPORT_TITLE = 'Export reaction';
export const PREVIEW_CLASS = 'chem-reaction-export-preview';
export const NOTE_CLASS = 'chem-reaction-export-note';
export const SVG_FILE = 'reaction.svg';
export const PNG_FILE = 'reaction.png';

/** The preview's box as the dialog opens (CSS px, its padding inside), the least it is made, and the most a picture is
 * enlarged in it: the box fills the dialog as the dialog is resized, the picture redrawn to fit it (the preview
 * stayed 720 x 420 in a dialog made larger). What Save writes keeps its print size. */
const PREVIEW_OPEN_WIDTH = 720;
const PREVIEW_OPEN_HEIGHT = 420;
/** What the dialog takes besides the preview's box (its title, the six choices, the note, the buttons: CSS px), so that
 * it opens within a short window. */
const PREVIEW_DIALOG_CHROME = 400;
const PREVIEW_MIN_HEIGHT = 120;
const PREVIEW_PAD = 12;
const PREVIEW_MAX_ZOOM = 4;
/** The class of the box of the choices (Style ... Map numbers), and its columns' least width and gap (CSS px): one
 * column in the dialog as it opens, two once it is made about 830 px wide. */
export const CHOICES_CLASS = 'chem-reaction-export-choices';
const CHOICE_COLUMN_MIN = 380;
const CHOICE_COLUMN_GAP = 16;

/** The zoom that fits a picture `width` x `height` (CSS px at print size) into a box `boxWidth` x `boxHeight`: at most
 * {@link PREVIEW_MAX_ZOOM}; 1 for a box or a picture of no size. */
export function previewZoom(width: number, height: number, boxWidth: number, boxHeight: number): number {
  if (width <= 0 || height <= 0 || boxWidth <= 0 || boxHeight <= 0)
    return 1;
  return Math.min(PREVIEW_MAX_ZOOM, boxWidth / width, boxHeight / height);
}
/** The largest side of a canvas, device px; a PNG larger at its dpi is written at less, and says so. */
const CANVAS_MAX_SIDE = 16000;
/** Without the drawing's `column` option: the room Crux is given for each of up to three rows (in the scheme's
 * natural height), and the blank kept around the ink (units), as the Reactions app's export does. */
const FALLBACK_MAX_ROWS = 3;
const FALLBACK_ROW_ROOM = 1.5;
const FALLBACK_MARGIN = 6;
/** The steps' numbers' font in each preset ({@link FIT_STEPS}; units), the legend's family too; the swatch's side and
 * the gaps in the legend's font size. */
const LEGEND_FONT: {[preset: string]: {size: number, family: string}} = {
  [ACS_1996]: {size: 8, family: 'Arial, Helvetica, sans-serif'},
  [SCREEN]: {size: 11, family: 'Roboto, Arial, sans-serif'},
};
/** The legend's size against the figure's (it was oversized: Screen's 11 px under its 8.7 px labels, and full size
 * under a scheme drawn smaller): this share of the look's label size ({@link TEXT_SIZE}: ACS 1996's 8 pt under its
 * 10 pt labels), times the scheme's scale, never under {@link LEGEND_MIN_SCALE} of it. */
const LEGEND_TEXT_SHARE = 0.8;
const LEGEND_MIN_SCALE = 0.6;

/** The legend's font in preset `preset` under a scheme drawn at `percent` of print size (null: print size). */
function legendFontOf(preset: string, percent: number | null): {size: number, family: string} {
  const scale = Math.max(LEGEND_MIN_SCALE, Math.min(1, (percent ?? 100) / 100));
  return {size: TEXT_SIZE[preset] * LEGEND_TEXT_SHARE * scale, family: LEGEND_FONT[preset].family};
}
const LEGEND_COLOUR = '#555555';
const LEGEND_SWATCH = 0.8;
const LEGEND_GAP = 0.4;
const LEGEND_ITEM_GAP = 1.2;
const LEGEND_LINE = 1.5;
/** The margin left and right of the legend, and above it (units). */
const LEGEND_PAD = 4;

export function presetOf(label: string): ExportPreset {
  return EXPORT_PRESETS.find((p) => p.label === label) ?? EXPORT_PRESETS[0];
}

export function widthOf(label: string): number | null {
  return EXPORT_WIDTHS.find((w) => w.label === label)?.cm ?? null;
}

/** `cm` in a preset's units, whole (8.3 cm: 235 pt, 314 CSS px). */
export function unitsOf(cm: number, unitsPerInch: number): number {
  return Math.round(cm / 2.54 * unitsPerInch);
}

/** A swatch or a text of the legend, at its place (units, from the legend's top left). */
interface LegendPart {x: number, y: number, text: string, colour?: string}

/** Where an export's parts are, in its preset's units: the picture's size, Crux's box in it, the legend under it. */
export interface ExportLayout {
  width: number;
  height: number;
  /** Crux's box, and its top relative to the picture's (the blank above the ink cut off without `column`). */
  boxWidth: number;
  boxHeight: number;
  boxTop: number;
  /** The scheme's part of the height (the legend's is the rest). */
  schemeHeight: number;
  /** The bonds against print size, percent: 100 at print size, less when the scheme was made smaller to fit; null when
   * the drawing does not say. */
  percent: number | null;
  /** Whether the drawing kept print size by taking rows (its `column` option). */
  column: boolean;
  /** {@link FIT_STEPS}: each step's row, in the route's order (the scheme's part is theirs; none otherwise). The texts
   * over a step's arrow a host adds later map onto these by `step`. */
  rows: StepRow[];
  legend: LegendPart[];
  /** The legend's font, sized to the figure; and the steps' numbers' ({@link FIT_STEPS}). */
  legendFont: {size: number, family: string};
  numberFont: {size: number, family: string};
  preset: ExportPreset;
  /** Crux's options besides the box. */
  options: {[key: string]: any};
}

/** A reaction read once for the Export panel, laid out and drawn for each choice. Every Crux call that throws is
 * reported to the module (`noteCruxError`: a trap it could not recover from turns Crux off) and leaves the export with
 * {@link error}. */
export class ReactionExport {
  readonly reaction: string;
  /** The texts over and under each step's arrow, drawn with the reaction (Crux's `stepTexts`); null, none. */
  readonly stepTexts: StepTexts | null;
  private scheme: Scheme | null = null;
  // Each step read as a reaction of its own ({@link FIT_STEPS}), once first asked for.
  private steps: Scheme[] | null = null;
  // Whether the steps are drawn from the route's own scheme ({@link stepSchemes}); `steps` then holds the route's
  // scheme once a step, nothing of their own to free.
  private ownSteps = false;
  error: string | null = null;
  private readonly layouts = new Map<string, ExportLayout>();

  constructor(reaction: string, stepTexts: StepTexts | null = null) {
    this.reaction = reaction;
    this.stepTexts = stepTexts;
    const crux = cruxLoadFailed() ? null : cruxModuleSync() as unknown as ChemWasm | null;
    if (!reaction.trim())
      this.error = 'There is no reaction to export';
    else if (tooLargeToDraw(reaction))
      this.error = TOO_LARGE;
    else if (!crux)
      this.error = CRUX_OFF;
    else {
      try {
        this.scheme = new crux.ChemReactionScheme(reaction, '');
      } catch (e) {
        noteCruxError(e);
        this.error = cruxLoadFailed() ? CRUX_OFF : 'This text does not read as a reaction';
      }
    }
  }

  private live(): Scheme | null {
    if (this.scheme !== null && cruxLoadFailed()) {
      this.scheme = null; // what Crux held dropped, never freed on a broken instance
      this.steps = null;
      this.error = CRUX_OFF;
    }
    return this.scheme;
  }

  // The steps' schemes, read once: the route's own for each step where Crux draws a step alone from the route's
  // layout (chem-wasm's `step`, so that an intermediate is drawn the same way in both its steps), else each
  // step read as a reaction of its own (the scheme's `get_steps`).
  private stepSchemes(scheme: Scheme): Scheme[] {
    if (this.steps === null && cruxDrawsSteps()) {
      this.ownSteps = true;
      this.steps = scheme.get_steps().map(() => scheme);
    }
    if (this.steps === null) {
      const crux = cruxModuleSync() as unknown as ChemWasm;
      const read: Scheme[] = [];
      try {
        for (const text of scheme.get_steps())
          read.push(new crux.ChemReactionScheme(text, ''));
      } catch (e) {
        read.forEach((s) => s.free());
        throw e;
      }
      this.steps = read;
    }
    return this.steps;
  }

  private lost(e: unknown): void {
    noteCruxError(e);
    const scheme = this.scheme;
    const steps = this.ownSteps ? [] : this.steps ?? [];
    this.scheme = null;
    this.steps = null;
    this.layouts.clear();
    this.error = cruxLoadFailed() ? CRUX_OFF : 'The reaction could not be drawn';
    if (!cruxLoadFailed()) {
      scheme?.free();
      steps.forEach((s) => s.free());
    }
  }

  /** Crux's options for `choice`, besides the box. The background is the export's own (white or none). */
  static options(choice: ExportChoice): {[key: string]: any} {
    const preset = presetOf(choice.preset);
    // a page, not a cell: no levels of detail (a cell's: the agents left out under 9 px bonds; crux-core's `levels`
    // option, which an engine before it ignores: there a single column drew the default route without its agents)
    const o: {[key: string]: any} = {look: preset.look, centre: choice.centre, clearBackground: false, levels: false};
    // the host's look: Chem's molecule options and its reaction centre's widths, as its cells (ACS 1996: Crux's own)
    if (preset.look === 'host')
      Object.assign(o, {drawOptions: chemDrawOptions(), ...CHEM_CENTRE});
    if (choice.mapNumbers)
      o.mapNumbers = choice.mapNumbers;
    // the box on: map numbers at any size, as asked (chem-wasm's `mapNumbersMinBond`; a file before it
    // ignores the key and leaves them out under 10 units a bond, which the warning says)
    if (choice.mapNumbers === 'show')
      o.mapNumbersMinBond = 0;
    return o;
  }

  /** The layout of `choice`, or null when the reaction cannot be drawn ({@link error} says why). */
  layout(choice: ExportChoice): ExportLayout | null {
    const key = JSON.stringify(choice);
    const kept = this.layouts.get(key);
    if (kept)
      return kept;
    const scheme = this.live();
    if (!scheme)
      return null;
    const preset = presetOf(choice.preset);
    const options = {...ReactionExport.options(choice), ...(this.stepTexts ? {stepTexts: this.stepTexts} : {})};
    const cm = widthOf(choice.width);
    let l: ExportLayout;
    try {
      const natural = Array.from(scheme.get_size(JSON.stringify(options)));
      const base = {preset, options, legend: [], legendFont: LEGEND_FONT[preset.label],
        numberFont: LEGEND_FONT[preset.label], boxTop: 0, rows: []};
      if (choice.fitSteps) {
        const fitted = stepRows(this.stepSchemes(scheme), options, (k) => this.stepOptions(k), cm === null ? null :
          unitsOf(cm, preset.unitsPerInch), LEGEND_FONT[preset.label]);
        l = {...base, ...fitted, boxWidth: fitted.width, boxHeight: fitted.height, schemeHeight: fitted.height,
          percent: Math.min(100, ...fitted.rows.map((r) => r.percent)), column: false};
      } else if (cm === null) {
        l = {...base, width: natural[0], height: natural[1], boxWidth: natural[0], boxHeight: natural[1],
          schemeHeight: natural[1], percent: 100, column: false};
      } else {
        const width = unitsOf(cm, preset.unitsPerInch);
        const asked = {...options, width, column: true};
        const column = Array.from(scheme.get_size(JSON.stringify(asked)));
        if (column.length > 2 && column[0] === width) {
          l = {...base, width, height: column[1], boxWidth: width, boxHeight: column[1], schemeHeight: column[1],
            percent: column[2], column: true};
          l.options = asked;
        } else
          l = {...base, ...fallbackFrame(scheme, options, natural, width), width, column: false};
      }
    } catch (e) {
      this.lost(e);
      return null;
    }
    const counts = choice.centre ? cruxCentreCounts(this.reaction) : null;
    if (counts && counts.some((n) => n > 0)) {
      l.legendFont = legendFontOf(preset.label, l.percent);
      const {parts, height} = legendLayout(counts, l.width, l.legendFont);
      l.legend = parts.map((p) => ({...p, y: p.y + l.schemeHeight}));
      l.height = l.schemeHeight + height;
    }
    this.layouts.set(key, l);
    return l;
  }

  /** Tests: whether {@link FIT_STEPS}'s steps are drawn from the route's own scheme (chem-wasm's `step`), once laid
   * out. */
  get drawsRouteSteps(): boolean {
    return this.ownSteps;
  }

  // Step `k`'s options ({@link FIT_STEPS}): from the route's scheme, step k fitted to its own box as a one-step
  // scheme is (`oneScale: false`; the route's texts are the options' already); else, a reaction of its own, its texts
  // alone.
  private stepOptions(k: number): {stepTexts?: StepTexts, step?: number, oneScale?: boolean} {
    return this.ownSteps ? {step: k, oneScale: false} : {stepTexts: stepTextsOfStep(this.stepTexts, k)};
  }

  // What `l` draws, box by box: the scheme in its box, or each step on its row at the steps' bond length.
  private boxes(scheme: Scheme, l: ExportLayout):
    {scheme: Scheme, options: object, x: number, y: number, width: number, height: number}[] {
    return l.rows.length ?
      l.rows.map((r) => ({scheme: this.steps![r.step], options: {...l.options, ...this.stepOptions(r.step),
        bondLength: r.bond}, x: r.x, y: r.y, width: r.width, height: r.height})) :
      [{scheme, options: l.options, x: 0, y: l.boxTop, width: l.boxWidth, height: l.boxHeight}];
  }

  /** Draws the export on `g` at `ratio` device px per unit, its top left at the canvas's; false when it cannot. */
  paint(g: CanvasRenderingContext2D, l: ExportLayout, ratio: number): boolean {
    const scheme = this.live();
    if (!scheme)
      return false;
    const boxes = this.boxes(scheme, l);
    for (const box of boxes) {
      g.save();
      try {
        g.setTransform(ratio, 0, 0, ratio, 0, 0);
        g.beginPath();
        g.rect(0, 0, l.width, l.schemeHeight);
        g.clip();
        box.scheme.draw_to_canvas(g.canvas, JSON.stringify({...box.options, width: box.width, height: box.height,
          offsetx: box.x, offsety: box.y}));
      } catch (e) {
        this.lost(e);
        return false;
      } finally {
        g.restore();
      }
    }
    if (l.rows.length) {
      g.save();
      g.setTransform(ratio, 0, 0, ratio, 0, 0);
      g.font = `${l.numberFont.size}px ${l.numberFont.family}`;
      g.textBaseline = 'middle';
      g.fillStyle = STEP_NUMBER_COLOUR;
      for (const r of l.rows)
        g.fillText(r.label, 0, r.y + r.height / 2);
      g.restore();
    }
    if (l.legend.length) {
      g.save();
      g.setTransform(ratio, 0, 0, ratio, 0, 0);
      const size = l.legendFont.size;
      g.font = `${size}px ${l.legendFont.family}`;
      g.textBaseline = 'middle';
      for (const p of l.legend) {
        if (p.colour) {
          g.fillStyle = p.colour;
          g.fillRect(p.x, p.y - size * LEGEND_SWATCH / 2, size * LEGEND_SWATCH, size * LEGEND_SWATCH);
        } else {
          g.fillStyle = LEGEND_COLOUR;
          g.fillText(p.text, p.x, p.y);
        }
      }
      g.restore();
    }
    return true;
  }

  /** The export as an SVG document, sized in its preset's unit (pt for ACS 1996), so that it is placed at print size;
   * null when it cannot be drawn. */
  svg(choice: ExportChoice): string | null {
    const l = this.layout(choice);
    const scheme = this.live();
    if (!l || !scheme)
      return null;
    const boxes = this.boxes(scheme, l);
    let inner = '';
    try {
      for (const box of boxes) {
        inner += box.scheme.get_svg(JSON.stringify({...box.options, width: box.width, height: box.height}))
          .replace(/^<\?xml[^>]*>\s*/, '').replace('<svg', `<svg x='${box.x}' y='${box.y}'`);
      }
    } catch (e) {
      this.lost(e);
      return null;
    }
    // each step's number left of its row ({@link FIT_STEPS}), in the legend's font
    inner += l.rows.map((r) => `<text class='${STEP_NUMBER_CLASS}' x='0' y='${r.y + r.height / 2}' ` +
      `dominant-baseline='middle' font-family='${l.numberFont.family}' font-size='${l.numberFont.size}' ` +
      `fill='${STEP_NUMBER_COLOUR}'>${r.label}</text>\n`).join('');
    const {width, height} = l;
    const unit = l.preset.unit;
    const background = choice.transparent ? '' :
      `<rect x='0' y='0' width='${width}' height='${height}' fill='#FFFFFF'/>\n`;
    const size = l.legendFont.size;
    const legend = l.legend.map((p) => p.colour ?
      `<rect x='${p.x}' y='${p.y - size * LEGEND_SWATCH / 2}' width='${size * LEGEND_SWATCH}' ` +
        `height='${size * LEGEND_SWATCH}' fill='${p.colour}'/>` :
      `<text x='${p.x}' y='${p.y}' dominant-baseline='middle' font-family='${l.legendFont.family}' ` +
        `font-size='${size}' fill='${LEGEND_COLOUR}'>${p.text}</text>`).join('\n');
    return `<?xml version='1.0' encoding='UTF-8'?>\n<svg xmlns='http://www.w3.org/2000/svg' version='1.1' ` +
      `width='${width}${unit}' height='${height}${unit}' viewBox='0 0 ${width} ${height}'>\n${background}` +
      `<svg x='0' y='0' width='${width}' height='${l.schemeHeight}' viewBox='0 0 ${width} ${l.schemeHeight}' ` +
      `overflow='hidden'>\n${inner}</svg>\n${legend}${legend ? '\n' : ''}</svg>\n`;
  }

  /** The export as a PNG at `dpi` (less when that is too large for a canvas: `dpi` in the answer says what it is). A
   * column's PNG is the column's width at that dpi: a single column at 300 dpi is 980 px wide. */
  async png(choice: ExportChoice, dpi: number):
    Promise<{png: Blob, dpi: number, width: number, height: number} | null> {
    const l = this.layout(choice);
    if (!l)
      return null;
    const cm = widthOf(choice.width);
    let px = cm === null ? Math.round(l.width / l.preset.unitsPerInch * dpi) : Math.round(cm / 2.54 * dpi);
    let ratio = px / l.width;
    const side = Math.max(px, l.height * ratio);
    if (side > CANVAS_MAX_SIDE) {
      ratio *= CANVAS_MAX_SIDE / side;
      px = Math.floor(l.width * ratio);
    }
    const canvas = document.createElement('canvas');
    canvas.width = px;
    canvas.height = Math.round(l.height * ratio);
    const g = canvas.getContext('2d')!;
    if (!choice.transparent) {
      g.fillStyle = '#FFFFFF';
      g.fillRect(0, 0, canvas.width, canvas.height);
    }
    if (!this.paint(g, l, ratio))
      return null;
    const png = await new Promise<Blob>((resolve, reject) =>
      canvas.toBlob((b) => b ? resolve(b) : reject(new Error('the PNG was not written')), 'image/png'));
    const written = ratio * l.preset.unitsPerInch;
    return {png: await withResolution(png, written), dpi: Math.round(written), width: canvas.width,
      height: canvas.height};
  }

  /** The export's size on the screen at print size (CSS px, 96 an inch), or null when it cannot be drawn. */
  printSizeOnScreen(choice: ExportChoice): {width: number, height: number} | null {
    const l = this.layout(choice);
    const css = l ? 96 / l.preset.unitsPerInch : 0;
    return l ? {width: l.width * css, height: l.height * css} : null;
  }

  /** The preview: the export on the screen at `zoom` times print size (96 CSS px an inch), as a canvas drawn anew at
   * that size (not a picture stretched). */
  preview(choice: ExportChoice, zoom = 1): HTMLCanvasElement | null {
    const l = this.layout(choice);
    if (!l)
      return null;
    const css = 96 / l.preset.unitsPerInch * zoom;
    const ratio = Math.min(window.devicePixelRatio * css, CANVAS_MAX_SIDE / Math.max(l.width, l.height));
    const canvas = ui.canvas(Math.round(l.width * ratio), Math.round(l.height * ratio));
    canvas.style.width = `${l.width * css}px`;
    canvas.style.height = `${l.height * css}px`;
    canvas.style.display = 'block';
    const g = canvas.getContext('2d')!;
    if (!choice.transparent) {
      g.fillStyle = '#FFFFFF';
      g.fillRect(0, 0, canvas.width, canvas.height);
    }
    return this.paint(g, l, ratio) ? canvas : null;
  }

  free(): void {
    const scheme = this.scheme;
    const steps = this.ownSteps ? [] : this.steps ?? [];
    this.scheme = null;
    this.steps = null;
    this.layouts.clear();
    if (!cruxLoadFailed()) {
      scheme?.free();
      steps.forEach((s) => s.free());
    }
  }
}

/** Without the drawing's `column` option: Crux's box `width` wide with room for up to three rows, the blank around its
 * ink cut off (measured on a canvas at one device px a unit). */
function fallbackFrame(scheme: Scheme, options: object, natural: number[], width: number):
  {height: number, boxWidth: number, boxHeight: number, boxTop: number, schemeHeight: number, percent: number | null} {
  const rows = Math.min(FALLBACK_MAX_ROWS, Math.ceil(natural[0] / width));
  const boxHeight = Math.ceil(rows === 1 ? natural[1] : natural[1] * rows * FALLBACK_ROW_ROOM);
  const canvas = new OffscreenCanvas(width, boxHeight);
  const g = canvas.getContext('2d', {willReadFrequently: true})! as unknown as CanvasRenderingContext2D;
  scheme.draw_to_canvas(canvas, JSON.stringify({...options, width, height: boxHeight}));
  const [top, bottom] = inkRows(g.getImageData(0, 0, width, boxHeight));
  const percent = natural[0] <= width ? 100 : null;
  if (bottom < top)
    return {height: boxHeight, boxWidth: width, boxHeight, boxTop: 0, schemeHeight: boxHeight, percent};
  const height = bottom - top + 1 + 2 * FALLBACK_MARGIN;
  return {height, boxWidth: width, boxHeight, boxTop: FALLBACK_MARGIN - top, schemeHeight: height, percent};
}

/** {@link FIT_STEPS}: each step on a row of its own, its number left of it in `font`, every step at one scale: the
 * largest at which the widest fits the room the `line` (units; null: no limit) leaves, print size at most (one
 * bond length across the figure, as journals expect; rows at scales of their own looked wrong). The
 * scale is a bond length ({@link PRINT_BOND}) each step is drawn with, its box its natural size at it; the rows one
 * under the other. The picture is the line wide (with no line, as its widest row). */
function stepRows(steps: Scheme[], options: {[key: string]: any}, stepOptions: (k: number) => object,
  line: number | null, font: {size: number, family: string}): {width: number, height: number, rows: StepRow[]} {
  const g = new OffscreenCanvas(1, 1).getContext('2d')! as unknown as CanvasRenderingContext2D;
  g.font = `${font.size}px ${font.family}`;
  const numberWidth = Math.ceil(g.measureText(String(steps.length)).width + font.size * STEP_NUMBER_GAP);
  const print = PRINT_BOND[options.look] ?? PRINT_BOND.publication;
  const sizesAt = (scale: number) =>
    steps.map((s, k) =>
      Array.from(s.get_size(JSON.stringify({...options, ...stepOptions(k), bondLength: print * scale}))));
  const widestAt = (sizes: number[][]) => Math.max(1, ...sizes.map(([w]) => w));
  let scale = 1;
  let sizes = sizesAt(scale);
  const room = line === null ? null : Math.max(1, line - numberWidth);
  if (room !== null && widestAt(sizes) > room) {
    // the largest scale at which the widest step fits, halved in on from the proportional guess (a step's width
    // follows the bond length nearly, not exactly, in proportion); a per-step scale would be min(1, room / its width)
    let [lo, hi] = [MIN_STEP_SCALE, 1];
    let fits: number[][] | null = null;
    let guess = Math.max(lo, room / widestAt(sizes));
    for (let pass = 0; pass < STEP_SCALE_PASSES && hi - lo > STEP_SCALE_CLOSE; pass++) {
      const at = sizesAt(guess);
      if (widestAt(at) <= room)
        [lo, fits] = [guess, at];
      else
        hi = guess;
      guess = (lo + hi) / 2;
    }
    scale = lo;
    sizes = fits ?? sizesAt(lo);
  }
  const rows: StepRow[] = [];
  let y = 0;
  sizes.forEach(([w, h], step) => {
    const row = {step, label: `${step + 1}`, x: numberWidth, y, width: Math.max(1, w), height: Math.max(1, h),
      bond: print * scale, percent: Math.round(scale * 100)};
    rows.push(row);
    y += row.height + font.size * STEP_ROW_GAP;
  });
  const height = Math.ceil(rows.length ? y - font.size * STEP_ROW_GAP : 1);
  return {width: line ?? numberWidth + widestAt(sizes), height, rows};
}

/** The legend under an export `width` units wide, for Crux's counts (`[made, broken, changed]`): the title, then a
 * swatch and words for each colour drawn, on as many lines as the width needs; and its height. */
function legendLayout(counts: number[], width: number, font: {size: number, family: string}):
  {parts: LegendPart[], height: number} {
  const g = new OffscreenCanvas(1, 1).getContext('2d')! as unknown as CanvasRenderingContext2D;
  g.font = `${font.size}px ${font.family}`;
  const s = font.size;
  const line = s * LEGEND_LINE;
  const items: {colour?: string, text: string}[] = [{text: LEGEND_TITLE}];
  counts.forEach((n, i) => {
    if (n > 0)
      items.push({colour: CENTRE_COLOURS[i], text: LEGEND_WORDS[i]});
  });
  const parts: LegendPart[] = [];
  let x = LEGEND_PAD;
  let y = LEGEND_PAD + line / 2;
  for (const item of items) {
    const swatch = item.colour ? s * (LEGEND_SWATCH + LEGEND_GAP) : 0;
    const w = swatch + g.measureText(item.text).width;
    if (x > LEGEND_PAD && x + w > width - LEGEND_PAD) {
      x = LEGEND_PAD;
      y += line;
    }
    if (item.colour)
      parts.push({x, y, text: '', colour: item.colour});
    parts.push({x: x + swatch, y, text: item.text});
    x += w + s * LEGEND_ITEM_GAP;
  }
  return {parts, height: y + line / 2 + LEGEND_PAD / 2};
}

/** What the note under the preview says of `l`: its size on paper, and whether it is at print size. */
export function exportNote(l: ExportLayout): string {
  const cm = (units: number) => (units / l.preset.unitsPerInch * 2.54).toFixed(1);
  const size = `${cm(l.width)} × ${cm(l.height)} cm`;
  if (l.rows.length) {
    // the steps' one scale; under half size, the warning says how small their text is
    const percent = Math.min(...l.rows.map((r) => r.percent));
    if (percent >= 100)
      return `${size}, each step on one line, at print size`;
    return `${size}; each step on one line, drawn at ${percent}% of print size`;
  }
  if (l.percent === null)
    return `${size}; made smaller to fit the width`;
  if (l.percent < 100)
    return `${size}; drawn at ${l.percent}% of print size: a part of it is wider than the width`;
  return `${size}, at print size`;
}

/** What the warning under the note says of `l`, drawn under half its print size: its text's size in print, too small
 * to read, and the wider width to choose; null at half size or more. */
export function exportWarning(l: ExportLayout): string | null {
  const percent = l.rows.length ? Math.min(...l.rows.map((r) => r.percent)) : l.percent;
  if (percent === null || percent >= SMALL_STEP_PERCENT)
    return null;
  const pt = (TEXT_SIZE[l.preset.label] * percent / 100 * 72 / l.preset.unitsPerInch).toFixed(1);
  const cm = l.width / l.preset.unitsPerInch * 2.54;
  const wider = cm < (widthOf(DOUBLE_COLUMN) ?? 0) - 0.5 ? 'the double column or the free width' : 'the free width';
  return `Too small to read in print: its text is ${pt} pt, at ${percent}% of print size. Choose ${wider}.`;
}

/** The least percent of print size at which `l`'s look draws map numbers ({@link MAP_NUMBERS_MIN_BOND}). */
export function mapNumbersMinPercent(l: ExportLayout): number {
  return Math.ceil(MAP_NUMBERS_MIN_BOND / (PRINT_BOND[l.preset.look] ?? PRINT_BOND.publication) * 100);
}

/** Whether `l` is drawn too small for Crux to draw its map numbers: its bonds (each fitted step's, else the scheme's)
 * under {@link MAP_NUMBERS_MIN_BOND}. False when its size is not known. */
export function mapNumbersLeftOut(l: ExportLayout): boolean {
  const percent = l.rows.length ? Math.min(...l.rows.map((r) => r.percent)) : l.percent;
  return percent !== null && percent < mapNumbersMinPercent(l);
}

/** What the warning under the note says when Map numbers is on, the reaction has map numbers, and `l` is drawn too
 * small for them ({@link mapNumbersLeftOut}): that they are left out, from what size they show, and the width to
 * choose (`doubleShows`: whether the double column shows them); null otherwise, and null where Crux draws them at any
 * size ({@link cruxDrawsSmallMapNumbers}: the note says they are small, {@link mapNumbersNote}). */
export function mapNumbersWarning(l: ExportLayout, choice: ExportChoice, reaction: string,
  doubleShows: boolean): string | null {
  if (choice.mapNumbers !== 'show' || !hasMapNumbers(reaction) || !mapNumbersLeftOut(l) || cruxDrawsSmallMapNumbers())
    return null;
  const percent = l.rows.length ? Math.min(...l.rows.map((r) => r.percent)) : l.percent;
  const wider = choice.width === SINGLE_COLUMN && doubleShows ? 'the double column or the free width' :
    'the free width';
  return `Map numbers are left out: they are drawn from ${mapNumbersMinPercent(l)}% of print size, and this is ` +
    `${percent}%. Choose ${wider}.`;
}

/** What the note says of map numbers drawn small (an engine with `mapNumbersMinBond` draws them at any size). */
export const MAP_NUMBERS_SMALL = 'Map numbers are small at this width';

/** What the note says where Crux draws map numbers at any size ({@link cruxDrawsSmallMapNumbers}) and `l` is drawn
 * under the size at which they used to be left out: they are small, and the width that draws them larger; null
 * otherwise. */
export function mapNumbersNote(l: ExportLayout, choice: ExportChoice, reaction: string): string | null {
  if (choice.mapNumbers !== 'show' || !hasMapNumbers(reaction) || !mapNumbersLeftOut(l) || !cruxDrawsSmallMapNumbers())
    return null;
  const wider = choice.width === SINGLE_COLUMN ? 'the double column or the free width' : 'the free width';
  return `${MAP_NUMBERS_SMALL}: ${wider} draws them larger`;
}

/** Whether Chem's Crux file draws map numbers at any size when asked (chem-wasm's `mapNumbersMinBond`): a
 * mapped reaction drawn with 4-unit bonds and the key at 0 has map numbers, which a file before it leaves out (it
 * ignores the key). Asked once a module. */
export function cruxDrawsSmallMapNumbers(): boolean {
  const crux = cruxLoadFailed() ? null : cruxModuleSync() as unknown as ChemWasm | null;
  if (!crux)
    return false;
  if (smallMapNumbersOf !== crux) {
    smallMapNumbersOf = crux;
    let scheme: Scheme | null = null;
    try {
      scheme = new crux.ChemReactionScheme('[CH3:1][OH:2]>>[CH2:1]=[O:2]', '');
      smallMapNumbers = scheme.get_svg(JSON.stringify({look: 'publication', width: 200, height: 60, bondLength: 4,
        mapNumbers: 'show', mapNumbersMinBond: 0, levels: false})).includes('class=\'map-number');
    } catch {
      smallMapNumbers = false;
    } finally {
      scheme?.free();
    }
  }
  return smallMapNumbers;
}
/** {@link cruxDrawsSmallMapNumbers}'s answer, and the module it is of. */
let smallMapNumbers = false;
let smallMapNumbersOf: ChemWasm | null = null;

/** Whether a reaction's text holds map numbers: `[C:1]` in reaction SMILES or SMARTS; in an RXN block, an atom line's
 * map number (V2000: the tenth number after the symbol; V3000: the number after the coordinates). */
function hasMapNumbers(reaction: string): boolean {
  if (/:\d+\]/.test(reaction))
    return true;
  if (!reaction.includes('$RXN'))
    return false;
  const coordinate = /^-?\d+\.\d+$/;
  return reaction.split('\n').some((line) => {
    const f = line.trim().split(/\s+/);
    if (f[0] === 'M' && f[1] === 'V30')
      return f.length > 7 && coordinate.test(f[4]) && coordinate.test(f[6]) && +f[7] > 0;
    return f.length > 13 && coordinate.test(f[0]) && coordinate.test(f[2]) && +f[13] > 0;
  });
}

/** {@link AS_DRAWN}'s picture options for `choice`: the style's print size as a bond in points (ACS 1996's 14.4 pt;
 * Screen's 22 CSS px, 16.5 pt) and its colours (ACS 1996 black, Screen in element colours), the map numbers as the
 * box says (the sketcher's own when the choice does not say). */
export function drawnOptions(choice: ExportChoice): DrawnPictureOptions {
  const preset = presetOf(choice.preset);
  // in points, to a thousandth (14.4 * 72 / 72 is 14.399999999999999 in floating point)
  const bondLength = Math.round(PRINT_BOND[preset.look] * 72 / preset.unitsPerInch * 1000) / 1000;
  const options: DrawnPictureOptions = {bondLength, colour: preset.look === 'host'};
  if (choice.mapNumbers)
    options.mapNumbers = choice.mapNumbers === 'show';
  return options;
}

/** The root `<svg>` tag of `svg`. */
const svgRoot = (svg: string) => /<svg\b[^>]*>/.exec(svg)?.[0] ?? '';

/** The width and height `svg`'s root says, in CSS px (a length in px or with no unit), or null. */
export function svgPixels(svg: string): {width: number, height: number} | null {
  const root = svgRoot(svg);
  const n = (name: string) => {
    const m = new RegExp(`\\s${name}\\s*=\\s*['"]([\\d.]+)(px)?['"]`).exec(root);
    return m ? Number(m[1]) : NaN;
  };
  const [width, height] = [n('width'), n('height')];
  return width > 0 && height > 0 ? {width, height} : null;
}

/** `svg` with its background left out (crux-sketch's `<g data-layer='background'>`, a white rectangle). */
export function withoutBackground(svg: string): string {
  return svg.replace(/(<g\s+data-layer\s*=\s*['"]background['"][^>]*>)[\s\S]*?(<\/g>)/, '$1$2');
}

/** A drawing's SVG (CSS px) as Save writes it: sized in `preset`'s unit, so that a document places it at the size it
 * was drawn (ACS 1996: points, 0.75 a CSS px), its view box kept, an XML prolog first. */
export function drawnSvgFile(svg: string, preset: ExportPreset): string {
  const size = svgPixels(svg);
  const root = svgRoot(svg);
  if (!size || !root)
    return svg;
  const k = preset.unitsPerInch / 96;
  const round = (v: number) => Math.round(v * 100) / 100;
  let sized = root.replace(/\swidth\s*=\s*['"][^'"]*['"]/, ` width='${round(size.width * k)}${preset.unit}'`)
    .replace(/\sheight\s*=\s*['"][^'"]*['"]/, ` height='${round(size.height * k)}${preset.unit}'`);
  if (!/\sviewBox\s*=/.test(sized))
    sized = sized.replace(/^<svg\b/, `<svg viewBox='0 0 ${size.width} ${size.height}'`);
  const body = svg.replace(/^\s*<\?xml[^>]*>\s*/, '').replace(root, sized);
  return `<?xml version='1.0' encoding='UTF-8'?>\n${body}`;
}

/** `svg` as an image, once it is decoded. */
function svgImage(svg: string): Promise<HTMLImageElement> {
  return new Promise((resolve, reject) => {
    const image = new Image();
    image.onload = () => resolve(image);
    image.onerror = () => reject(new Error('the picture could not be decoded'));
    image.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
  });
}

/** What the note says of a picture as drawn, `width` x `height` CSS px: its size on paper, and that it is at print
 * size (`takesOptions`), else at the canvas's own size. */
export function drawnNote(width: number, height: number, takesOptions: boolean): string {
  const cm = (px: number) => (px / 96 * 2.54).toFixed(1);
  return `${cm(width)} × ${cm(height)} cm, as drawn, ${takesOptions ? 'at print size' : DRAWN_OWN_LOOK}`;
}

/** What the warning under the note says of a picture as drawn `width` CSS px wide: wider than a double column, it is
 * printed smaller; null when it fits one. */
export function drawnWarning(width: number): string | null {
  const cm = width / 96 * 2.54;
  const double = widthOf(DOUBLE_COLUMN) ?? 0;
  if (cm <= double + 0.05)
    return null;
  return `${cm.toFixed(1)} cm wide, wider than a double column (${double} cm): a journal prints it smaller. ` +
    `Draw it narrower, or turn ${AS_DRAWN} off to arrange it in a column.`;
}

/** The Export panel's inputs (the tests set them as a person would). */
export interface ExportInputs {
  style: DG.InputBase<string | null>;
  width: DG.InputBase<string | null>;
  format: DG.InputBase<string | null>;
  transparent: DG.InputBase<boolean | null>;
  colours: DG.InputBase<boolean | null>;
  fitSteps: DG.InputBase<boolean | null>;
  mapNumbers: DG.InputBase<boolean | null>;
  /** {@link AS_DRAWN}, offered when the caller has a drawing; else null. */
  asDrawn: DG.InputBase<boolean | null> | null;
}

/** What a caller adds to the Export panel of a reaction. */
export interface ExportOptions {
  /** The map numbers as the caller shows them. */
  mapNumbers?: 'show' | 'hide';
  /** The drawing's own CDXML page (the Reactions app: its sketcher's document as drawn, the texts kept), asked once
   * CDXML is chosen, in place of the page written from the reaction text; null when it is not written. */
  cdxml?: () => Promise<string | null>;
  /** The drawing has texts that the SVG and PNG leave out: the panel says so ({@link TEXTS_CDXML_ONLY}), unless they
   * are its steps' texts and Crux draws them ({@link stepTexts}). */
  textsLeftOut?: boolean;
  /** The texts over and under each step's arrow, as the drawing holds them (the Reactions app: crux-sketch's
   * `reactionSteps`), drawn in the SVG, the PNG and the preview as in the drawing ("HATU, DIPEA" typed over the arrow
   * was once in no picture). */
  stepTexts?: StepTexts | null;
  /** How many molecules of the ChemDraw page the reaction was read from are on no step (a CDXML cell's): the SVG
   * and PNG, drawn from the reaction, leave them out, and the panel says so; its CDXML is the
   * page. */
  moleculesLeftOut?: number;
  /** The drawing as it is drawn (the Reaction Sketcher's canvas, a cell's page or RXN block with coordinates): the
   * panel offers {@link AS_DRAWN} with it. */
  asDrawn?: DrawnPicture | null;
}

/** Opens the Export panel for `reaction` (see {@link ExportOptions}). */
export function openReactionExport(reaction: string, options: ExportOptions = {}): DG.Dialog {
  return exportPanel(reaction, options).dialog;
}

/** A picture as drawn ({@link AS_DRAWN}): its SVG (CSS px; its background left out for a transparent one), its size,
 * and the SVG decoded, which the preview draws. */
interface Drawn {
  svg: string;
  width: number;
  height: number;
  image: HTMLImageElement;
}

/** The open panels' inputs, by their dialog's root (the tests reach a panel a cell's menu opened by it). */
const panelInputs = new WeakMap<HTMLElement, ExportInputs>();

/** Tests: the inputs of the Export panel whose dialog's root is `root`. */
export function exportInputsOf(root: HTMLElement): ExportInputs | undefined {
  return panelInputs.get(root);
}

/** {@link openReactionExport}'s dialog, and its inputs. */
export function exportPanel(reaction: string, options: ExportOptions = {}):
  {dialog: DG.Dialog, inputs: ExportInputs} {
  const picture = new ReactionExport(reaction, options.stepTexts ?? null);
  // the drawing's texts drawn with it, where Crux draws a host's step texts; else CDXML alone keeps them
  const textsDrawn = !!options.stepTexts && cruxDrawsStepTexts();
  /** The reaction as CDXML: asked once, when CDXML is first chosen; undefined until written, null when it is not. */
  let cdxml: Promise<string | null> | null = null;
  let written: string | null | undefined = undefined;
  /** The drawing as drawn ({@link AS_DRAWN}), when the caller has one: its pictures, made once per options; the one
   * the preview shows; the turn of the last asked for (a later choice overtakes it). */
  const drawn = options.asDrawn ?? null;
  const drawnPictures = new Map<string, Promise<Drawn | null>>();
  let drawnShown: Drawn | null = null;
  let drawnTurn = 0;
  const shows = (on: boolean): 'show' | 'hide' => on ? 'show' : 'hide';
  const choice: ExportChoice = {preset: ACS_1996, width: SINGLE_COLUMN, transparent: false,
    centre: presetOf(ACS_1996).centre, mapNumbers: shows(presetMapNumbers(ACS_1996, options.mapNumbers))};
  let format = SVG;
  // the box fills the dialog's room (its contents are a flex column), the picture centred in it, fitted to it
  const openHeight = Math.max(PREVIEW_MIN_HEIGHT,
    Math.min(PREVIEW_OPEN_HEIGHT, window.innerHeight - PREVIEW_DIALOG_CHROME));
  const preview = ui.div([], {classes: PREVIEW_CLASS, style: {flex: '1 1 auto', height: `${openHeight}px`,
    minHeight: `${PREVIEW_MIN_HEIGHT}px`, position: 'relative', overflow: 'hidden', padding: `${PREVIEW_PAD}px`,
    background: 'var(--grey-1)'}});
  /** The box's size the picture was last fitted to. */
  let fittedTo = '';
  const note = ui.divText('', NOTE_CLASS);
  note.style.color = 'var(--grey-5)';
  // under the note, a warning's sign in Datagrok's warning orange and its words in the dark orange that reads on the
  // light preview (the light orange alone is faint as text): a figure too small to read in print
  const warningText = ui.divText('');
  const warningSign = ui.iconFA('exclamation-triangle');
  warningSign.style.color = 'var(--orange-2)';
  const warning = ui.divH([warningSign, warningText], NOTE_WARNING_CLASS);
  Object.assign(warning.style, {display: 'none', color: 'var(--orange-3)', gap: '6px', alignItems: 'baseline',
    fontWeight: '500'});
  const colours = ui.input.bool('Reaction centre colours', {value: choice.centre, onValueChanged: (v) => {
    choice.centre = !!v;
    refresh();
  }});
  colours.setTooltip('The bonds each mapped step makes, breaks or changes the order of, with their legend');
  const mapNumbers = ui.input.bool(MAP_NUMBERS, {value: choice.mapNumbers === 'show', onValueChanged: (v) => {
    choice.mapNumbers = shows(!!v);
    refresh();
  }});
  mapNumbers.setTooltip('The atoms\' map numbers, as the reaction is mapped; in CDXML, the atom mapping');
  // a reaction without map numbers has no centre to colour and no numbers to show: both boxes off, and why (they were
  // on and did nothing for a route the Reaction panel had mapped step by step and passed on unmapped)
  const mapped = hasMapNumbers(reaction);
  if (!mapped) {
    for (const input of [colours, mapNumbers]) {
      input.enabled = false;
      input.setTooltip(NOT_MAPPED_TOOLTIP);
    }
  }
  const preset = ui.input.choice('Style', {items: EXPORT_PRESETS.map((p) => p.label), value: choice.preset,
    onValueChanged: (v) => {
      choice.preset = v ?? ACS_1996;
      // each preset's own colours (off for ACS 1996) and map numbers (off for ACS 1996)
      choice.centre = presetOf(choice.preset).centre;
      colours.value = choice.centre;
      choice.mapNumbers = shows(presetMapNumbers(choice.preset, options.mapNumbers));
      mapNumbers.value = choice.mapNumbers === 'show';
      refresh();
    }});
  const width = ui.input.choice('Width', {items: EXPORT_WIDTHS.map((w) => w.label), value: choice.width,
    onValueChanged: (v) => {
      choice.width = v ?? SINGLE_COLUMN;
      refresh();
    }});
  const fitSteps = ui.input.bool(FIT_STEPS, {value: false, onValueChanged: (v) => {
    choice.fitSteps = !!v;
    refresh();
  }});
  fitSteps.setTooltip(FIT_STEPS_TOOLTIP);
  // as drawn: the width, the arrangements and the centre's colours (which the canvas does not draw) do not apply
  const asDrawn = drawn ? ui.input.bool(AS_DRAWN, {value: false, onValueChanged: (v) => {
    choice.asDrawn = !!v;
    for (const input of [width, fitSteps])
      input.enabled = !choice.asDrawn;
    colours.enabled = !choice.asDrawn && mapped;
    refresh();
  }}) : null;
  asDrawn?.setTooltip(AS_DRAWN_TOOLTIP);
  const isDrawn = () => !!choice.asDrawn && drawn !== null;
  /** The picture as drawn for `c`, made once for its options. */
  const drawnOf = (c: ExportChoice): Promise<Drawn | null> => {
    const o = drawnOptions(c);
    const key = JSON.stringify([o, !!c.transparent]);
    let job = drawnPictures.get(key);
    if (!job) {
      job = (async (): Promise<Drawn | null> => {
        const made = await drawn!.svg(o);
        if (!made)
          return null;
        const svg = c.transparent ? withoutBackground(made) : made;
        const size = svgPixels(svg);
        return size ? {svg, ...size, image: await svgImage(svg)} : null;
      })().catch((e) => {
        console.warn('Chem | the drawing was not pictured as drawn', e);
        return null;
      });
      drawnPictures.set(key, job);
    }
    return job;
  };
  /** The picture as drawn as a PNG at `dpi` (less when that is too large for a canvas: `dpi` says what it is): the
   * sketcher's own (white); a transparent one, the SVG drawn on a canvas with no background. */
  const drawnPng = async (dpi: number): Promise<{png: Blob, dpi: number} | null> => {
    const d = await drawnOf(choice);
    if (!d || !drawn)
      return null;
    const scale = Math.min(dpi / 96, CANVAS_MAX_SIDE / Math.max(d.width, d.height));
    const at = Math.round(scale * 96);
    if (!choice.transparent) {
      const png = await drawn.png({...drawnOptions(choice), scale});
      return png ? {png, dpi: at} : null;
    }
    const canvas = new OffscreenCanvas(Math.round(d.width * scale), Math.round(d.height * scale));
    canvas.getContext('2d')!.drawImage(d.image, 0, 0, canvas.width, canvas.height);
    return {png: await withResolution(await canvas.convertToBlob({type: 'image/png'}), scale * 96), dpi: at};
  };
  const formats = ui.input.choice('Format', {items: EXPORT_FORMATS, value: format, onValueChanged: (v) => {
    format = v ?? SVG;
    refresh();
  }});
  const transparent = ui.input.bool('Transparent background', {value: false, onValueChanged: (v) => {
    choice.transparent = !!v;
    refresh();
  }});
  const dialog = ui.dialog({title: EXPORT_TITLE, helpUrl: reactionsHelpUrl('export')});
  /** The picture drawn anew to fit the box as it is now; whether it was drawn. */
  const draw = (): boolean => {
    ui.empty(preview);
    const print = picture.printSizeOnScreen(choice);
    const zoom = print ? previewZoom(print.width, print.height, preview.clientWidth - 2 * PREVIEW_PAD,
      preview.clientHeight - 2 * PREVIEW_PAD) : 1;
    fittedTo = `${preview.clientWidth}x${preview.clientHeight}`;
    const canvas = print ? picture.preview(choice, zoom) : null;
    if (canvas)
      placed(canvas);
    preview.dataset.zoom = zoom.toFixed(3);
    return canvas !== null;
  };
  /** `canvas` into the box: centred in it; a transparent picture over a chequer. */
  const placed = (canvas: HTMLCanvasElement) => {
    Object.assign(canvas.style, {position: 'absolute', inset: '0', margin: 'auto',
      backgroundImage: choice.transparent ? 'repeating-conic-gradient(#e6e6e6 0% 25%, #ffffff 0% 50%)' : '',
      backgroundSize: choice.transparent ? '16px 16px' : ''});
    preview.append(canvas);
  };
  /** The picture as drawn fitted to the box as it is now, drawn anew at that size (its SVG's vectors); whether it was
   * drawn. */
  const fitDrawn = (): boolean => {
    ui.empty(preview);
    const d = drawnShown;
    if (!d)
      return false;
    const zoom = previewZoom(d.width, d.height, preview.clientWidth - 2 * PREVIEW_PAD,
      preview.clientHeight - 2 * PREVIEW_PAD);
    fittedTo = `${preview.clientWidth}x${preview.clientHeight}`;
    const ratio = Math.min(window.devicePixelRatio * zoom, CANVAS_MAX_SIDE / Math.max(d.width, d.height));
    const canvas = ui.canvas(Math.round(d.width * ratio), Math.round(d.height * ratio));
    Object.assign(canvas.style, {width: `${d.width * zoom}px`, height: `${d.height * zoom}px`, display: 'block'});
    const g = canvas.getContext('2d')!;
    if (!choice.transparent) {
      g.fillStyle = '#FFFFFF';
      g.fillRect(0, 0, canvas.width, canvas.height);
    }
    g.drawImage(d.image, 0, 0, canvas.width, canvas.height);
    placed(canvas);
    preview.dataset.zoom = zoom.toFixed(3);
    return true;
  };
  // the dialog resized: the picture fitted to the box again (the export itself is unchanged)
  let pending = 0;
  const resized = new ResizeObserver(() => {
    cancelAnimationFrame(pending);
    pending = requestAnimationFrame(() => {
      if (!preview.isConnected || fittedTo === `${preview.clientWidth}x${preview.clientHeight}`)
        return;
      if (isDrawn())
        fitDrawn();
      else if (!picture.error)
        draw();
    });
  });
  const enableSave = (on: boolean) => {
    const save = dialog.getButton('Save');
    if (save)
      ui.setClass(save, 'disabled', !on);
  };
  // CDXML: the drawing's page as Crux Sketch's reaction writer writes it, asked once, when CDXML is first chosen (the
  // format's own; of the picture's choices, only Map numbers applies to it)
  const startCdxml = () => {
    cdxml ??= options.cdxml ? options.cdxml() : reactionCdxml(reaction);
    cdxml.then((page) => {
      written = page;
      if (format === CDXML && dialog.root.isConnected)
        refreshCdxml();
    });
  };
  /** {@link AS_DRAWN}: the drawing's picture for the choice, once made (a later choice overtakes it), in the preview;
   * the note its size, the warning a picture wider than a double column. */
  const refreshDrawn = () => {
    const turn = ++drawnTurn;
    drawnShown = null;
    ui.empty(preview);
    Object.assign(preview.dataset, {width: '', height: '', cm: '', look: presetOf(choice.preset).look, percent: '',
      column: 'false', legend: 'false', rows: '0', drawn: 'true'});
    note.textContent = [DRAWN_DRAWING, ...(format === CDXML ? cdxmlLines() : [])].join('. ');
    warning.style.display = 'none';
    enableSave(format === CDXML && !!written);
    void drawnOf(choice).then((d) => {
      if (turn !== drawnTurn || !dialog.root.isConnected)
        return;
      drawnShown = d;
      const shown = fitDrawn();
      if (d) {
        Object.assign(preview.dataset, {width: `${d.width}`, height: `${d.height}`,
          cm: (d.width / 96 * 2.54).toFixed(1), percent: drawn!.takesOptions ? '100' : ''});
      }
      const lines = [d && shown ? drawnNote(d.width, d.height, drawn!.takesOptions) : DRAWN_NOT_MADE];
      if (format === CDXML)
        lines.push(...cdxmlLines());
      note.textContent = lines.join('. ');
      const warned = d && shown && format !== CDXML ? drawnWarning(d.width) : null;
      warningText.textContent = warned ?? '';
      warning.style.display = warned ? 'flex' : 'none';
      enableSave(format === CDXML ? !!written : shown);
    });
  };
  const refresh = () => {
    if (format === CDXML)
      startCdxml();
    if (isDrawn()) {
      refreshDrawn();
      return;
    }
    drawnTurn++; // a picture as drawn still being made is not shown
    preview.dataset.drawn = 'false';
    const l = picture.layout(choice);
    if (!l)
      ui.empty(preview);
    const canvas = l ? draw() : false;
    // what is previewed, for the tests: the size in units, the bonds against print size, a column's rows, the rows of
    // "Fit each step on one line" (and the zoom on the screen, set as it is drawn)
    Object.assign(preview.dataset, {width: `${l?.width ?? ''}`, height: `${l?.height ?? ''}`,
      cm: l ? (l.width / l.preset.unitsPerInch * 2.54).toFixed(1) : '', look: l?.preset.look ?? '',
      percent: `${l?.percent ?? ''}`, column: `${!!l?.column}`, legend: `${(l?.legend.length ?? 0) > 0}`,
      rows: `${l?.rows.length ?? 0}`});
    const lines = [l && canvas ? exportNote(l) : picture.error ?? 'The reaction could not be drawn'];
    if (options.textsLeftOut && !textsDrawn)
      lines.push(TEXTS_CDXML_ONLY);
    if ((options.moleculesLeftOut ?? 0) > 0 && format !== CDXML)
      lines.push(`${noStepWords(options.moleculesLeftOut!)} of the SVG and PNG; CDXML keeps them`);
    // CDXML: once written ({@link startCdxml})
    if (format === CDXML)
      lines.push(...cdxmlLines());
    // Map numbers on, drawn small (an engine that draws them at any size): said in the note
    const smallMaps = l && canvas && format !== CDXML ? mapNumbersNote(l, choice, reaction) : null;
    if (smallMaps)
      lines.push(smallMaps);
    note.textContent = lines.join('. ');
    const small = l && canvas && format !== CDXML ? exportWarning(l) : null;
    // Map numbers on, and the figure too small for them (an engine without `mapNumbersMinBond`): said, with the width
    // that shows them
    const doubleShows = () => {
      const double = picture.layout({...choice, width: DOUBLE_COLUMN});
      return !!double && !mapNumbersLeftOut(double);
    };
    const maps = l && canvas && format !== CDXML ? mapNumbersWarning(l, choice, reaction,
      choice.width === SINGLE_COLUMN && doubleShows()) : null;
    const warned = [small, maps].filter((w) => w !== null).join(' ');
    warningText.textContent = warned;
    warning.style.display = warned ? 'flex' : 'none';
    const save = dialog.getButton('Save');
    if (save)
      ui.setClass(save, 'disabled', format === CDXML ? !written : !canvas);
  };
  // what the note says of the CDXML page: being written, not written, written; and, a mapped page, whether it keeps
  // its atom mapping, as the Map numbers box says
  const cdxmlLines = (): string[] => {
    const lines = [written === undefined ? CDXML_WRITING : written === null ? CDXML_NOT_YET : CDXML_WRITTEN];
    if (written && cdxmlMapped(written))
      lines.push(choice.mapNumbers === 'show' ? CDXML_WITH_MAPPING : CDXML_WITHOUT_MAPPING);
    return lines;
  };
  // the CDXML lines and Save, once the page is written (the preview is the picture's, unchanged)
  const refreshCdxml = () => {
    const lines = note.textContent!.split('. ').filter((t) => t !== CDXML_WRITING);
    lines.push(...cdxmlLines());
    note.textContent = [...new Set(lines)].join('. ');
    const save = dialog.getButton('Save');
    if (save)
      ui.setClass(save, 'disabled', !written);
  };
  /** Why the picture was not made (Save's and Copy's warning). */
  const notDrawn = () => isDrawn() ? DRAWN_NOT_MADE.toLowerCase() : picture.error ?? 'the reaction could not be drawn';
  /** The PNG at `dpi`: the drawing's as drawn, or the reaction's as arranged. */
  const pngOf = (dpi: number): Promise<{png: Blob, dpi: number} | null> =>
    isDrawn() ? drawnPng(dpi) : picture.png(choice, dpi);
  const save = async () => {
    try {
      if (format === CDXML) {
        if (!written)
          return;
        DG.Utils.download(CDXML_FILE, cdxmlMapping(written, choice.mapNumbers === 'show'), 'chemical/x-cdxml');
        dialog.close();
        return;
      }
      if (format === SVG && isDrawn()) {
        const d = await drawnOf(choice);
        if (!d)
          throw new Error(notDrawn());
        DG.Utils.download(SVG_FILE, drawnSvgFile(d.svg, presetOf(choice.preset)), 'image/svg+xml');
      } else if (format === SVG) {
        const svg = picture.svg(choice);
        if (svg === null)
          throw new Error(notDrawn());
        DG.Utils.download(SVG_FILE, svg, 'image/svg+xml');
      } else {
        const png = await pngOf(format === PNG_600 ? 600 : 300);
        if (png === null)
          throw new Error(notDrawn());
        DG.Utils.download(PNG_FILE, png.png, 'image/png');
        if (png.dpi < (format === PNG_600 ? 600 : 300))
          grok.shell.info(`Saved at ${png.dpi} dpi, the most a picture this large allows`);
      }
      dialog.close();
    } catch (e) {
      grok.shell.warning(`Not saved: ${e instanceof Error ? e.message : e}`);
    }
  };
  const copy = async () => {
    try {
      const png = await pngOf(format === PNG_600 ? 600 : 300);
      if (png === null)
        throw new Error(notDrawn());
      await navigator.clipboard.write([new ClipboardItem({'image/png': png.png})]);
      grok.shell.info('Copied as a picture');
    } catch (e) {
      grok.shell.warning(`Not copied: ${e instanceof Error ? e.message : e}`);
    }
  };
  const choices = [preset, width, fitSteps, ...(asDrawn ? [asDrawn] : []), formats, transparent, colours, mapNumbers];
  choices.forEach((input) => dialog.add(input));
  dialog.add(preview).add(note).add(warning)
    .addButton('Save', () => save())
    .addButton(COPY_AS_PICTURE, () => copy(), 1, 'The PNG on the clipboard, at the chosen dpi (300 for SVG)');
  // the choices in as many columns as the dialog's width holds, one as it opens: a dialog made wider gives the preview
  // the rows it frees, so that a picture held by the box's height grows as the dialog widens too (widened, only
  // the grey grew)
  const columns = ui.div([], CHOICES_CLASS);
  Object.assign(columns.style, {display: 'grid', columnGap: `${CHOICE_COLUMN_GAP}px`,
    gridTemplateColumns: `repeat(auto-fill, minmax(${CHOICE_COLUMN_MIN}px, 1fr))`});
  preset.root.before(columns);
  columns.append(...choices.map((input) => input.root));
  dialog.onClose.subscribe(() => {
    resized.disconnect();
    cancelAnimationFrame(pending);
    picture.free();
  });
  // as wide as the preview's box at its opening size, its padding and the dialog's
  dialog.show({resizable: true, width: PREVIEW_OPEN_WIDTH + 4 * PREVIEW_PAD});
  // made taller than the window's 80% (the contents' own cap), the box follows still
  const contents = dialog.root.querySelector('.d4-dialog-contents') as HTMLElement | null;
  if (contents)
    contents.style.maxHeight = 'none';
  refresh();
  resized.observe(preview);
  const inputs = {style: preset, width, fitSteps, format: formats, transparent, colours, mapNumbers, asDrawn};
  panelInputs.set(dialog.root, inputs);
  return {dialog, inputs};
}

/** The reaction cell menu's items. */
export const EXPORT_ITEM = 'Export reaction...';
export const EDIT_ITEM = 'Edit reaction...';

/** A grid's context menu (`grok.events.onContextMenu`): on a reaction cell, in both modes, "Edit reaction..." opens
 * the editor a double-click opens (double-click had been the only way in), and
 * "Export reaction..." the Export panel on its value, when it has one; nothing else. */
export function addReactionExportItem(args: any): boolean {
  const cell = args?.args?.item;
  const menu = args?.args?.menu as DG.Menu | undefined;
  if (!menu || !cell?.isTableCell || cell.tableColumn?.semType !== DG.SEMTYPE.CHEMICAL_REACTION)
    return false;
  menu.item(EDIT_ITEM, () => editFromMenu(cell as DG.GridCell));
  const value: string = cell.cell?.valueString ?? '';
  if (value) {
    menu.item(EXPORT_ITEM, () => exportCellValue(value).catch((e) =>
      console.error('Chem | the Export panel did not open', e)));
  }
  return true;
}

/** "Export reaction..." on a cell's `value`: the Export panel on it; on a CDXML page's the reaction Crux Sketch reads
 * from the page, which the pictures draw, its CDXML the page itself (the panel once read the page as reaction text,
 * which Crux does not read: SVG and PNG were off, and CDXML a page written anew). A page that does not read is
 * given as it is: the panel says so, its CDXML still the page. */
export async function exportCellValue(value: string): Promise<DG.Dialog> {
  // the panel draws with Chem's Crux engine, which RDKit mode loads here; one that does not load: the panel says so
  await cruxForReactions('the Export panel');
  // "As drawn" for a value drawn with coordinates of its own (a CDXML page, an RXN block with coordinates): there is no
  // canvas, so the page's reader pictures it
  const asDrawn = cellPicture(value);
  if (!isCdxmlText(value))
    return openReactionExport(value, {asDrawn});
  const smiles = await cdxmlReactionSmiles(value);
  return openReactionExport(smiles ?? value, {cdxml: async () => value, moleculesLeftOut: moleculesOnNoStep(value),
    asDrawn});
}

/** "Edit reaction...": the cell's editor, as its double-click opens it. */
async function editFromMenu(cell: DG.GridCell): Promise<void> {
  try {
    const {openReactionCellEditor} = await import('../utils/reactions/reaction-cell-editor');
    await openReactionCellEditor(cell);
  } catch (e) {
    console.error('Chem | the reaction editor did not open', e);
  }
}

let menuInstalled = false;

/** Chem's init: the reaction cells' "Edit reaction..." and "Export reaction..." (once). */
export function installReactionExportMenu(): void {
  if (menuInstalled)
    return;
  menuInstalled = true;
  grok.events.onContextMenu.subscribe((args: any) => addReactionExportItem(args));
}
