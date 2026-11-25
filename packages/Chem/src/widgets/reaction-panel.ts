/**
 * The Reaction panel of the context panel: the current reaction drawn large, with a full-screen view, a route's steps
 * one by one (as Crux reads them: the steps the picture draws), and automatic atom mapping once Chem's Crux file has
 * the mapper (mapping happens here only, never over a column or into a cell, off the UI thread). Sections still to come
 * go between the steps and the mapping: the description (Crux's `describe`), the components with their roles, the exact
 * view of a template. It shows in Crux mode only (rendering/crux-reaction-drawing.ts) and draws as the cells do:
 * synchronously, through Crux, on canvases of its own.
 */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {parseMultiStepReaction} from '../rendering/rdkit-reaction-renderer';
import {CentreWidths, StepTexts, chemDrawOptions, cruxDrawsReactions, cruxDrawsSteps, cruxFittedPicture,
  cruxReactionNotes, cruxReactionSteps, cruxReadProblem, cruxReadsReaction, drawReactionByCrux, productAtomsMissing,
  stepTextsOfStep, tooLargeToDraw} from '../rendering/crux-reaction-drawing';
import {cruxMessage} from '../utils/reactions/reaction-document';
import {MappedReaction, mapReaction, mappedInTheValue, mapperExported, mapperReady}
  from '../utils/reactions/reaction-mapping';
import {Flag, checkMapping, checkReaction, checkTemplateMapping, checksModule} from '../utils/footguns/checks';
import {PICTURE_LEGEND_GAP, centreLegend, pictureWithLegend, setCentreLegend} from './reaction-centre-legend';
import {openReactionExport} from './reaction-export';
import {CDXML_DRAWING, CDXML_UNREAD, DrawnPicture, cdxmlReactionSmiles, cellPicture, isCdxmlText, moleculesOnNoStep,
  noStepWords} from '../rendering/cdxml-reaction-cells';
import {reactionsHelpUrl} from '../utils/reactions/reactions-help';

/** The large picture: its height as a share of the panel's width, and its least height (CSS px). */
const PICTURE_ASPECT = 0.4;
const PICTURE_MIN_HEIGHT = 120;
/** The large picture's bond length (px): larger than a cell's 22, so a small reaction or template draws larger than its
 * cell; a scheme too wide for the panel is shrunk to fit as before. */
const PICTURE_BOND_LENGTH = 32;
/** A step's picture in the Steps section: its height (CSS px); its width is the panel's. */
const STEP_HEIGHT = 110;
/** The full-screen view's steps' bond length (px; the cells' is 22, and Crux never enlarges a scheme to fill its box);
 * its picture's is fitted ({@link FULL_SCREEN_MAX_BOND}). */
const FULL_SCREEN_BOND_LENGTH = 44;
/** The full-screen dialog's room to the window's edges (CSS px). Datagrok's `fullScreen` show option made a dialog as
 * narrow as its contents (128 px): the dialog is placed and sized instead. */
const FULL_SCREEN_EDGE = 24;

/** A canvas as wide as its parent that draws `reaction`, `height(width)` high (at `bondLength` px, else Crux's own;
 * the molecules with Chem's draw options, as molecule cells), again whenever its width changes. */
/** Each picture's size subscription, which {@link letGo} ends. */
const sizeSubs = new WeakMap<HTMLCanvasElement, Subscription>();

/** Ends `canvas`'s size subscription: it is drawn no more, and nothing holds it. */
function letGo(canvas: HTMLCanvasElement): void {
  sizeSubs.get(canvas)?.unsubscribe();
  sizeSubs.delete(canvas);
}

function reactionCanvas(reaction: string, height: (width: number) => number, bondLength?: number,
  mapNumbers?: 'show' | 'hide', subs?: Subscription[], stepTexts?: StepTexts | null, step?: RouteStep):
  HTMLCanvasElement {
  // a route's step: drawn from the route where Crux draws a step alone (one layout, one scale for every card), else as
  // a reaction of its own with its texts alone
  const own = step !== undefined && cruxDrawsSteps();
  const drawn = step === undefined || own ? reaction : step.text;
  const texts = step === undefined || own ? stepTexts : stepTextsOfStep(stepTexts, step.index);
  const shownTexts = step === undefined ? stepTexts : stepTextsOfStep(stepTexts, step.index);
  const canvas = ui.canvas();
  canvas.style.width = '100%';
  // what this canvas draws (the tests read it), and how many step texts with it
  canvas.dataset.reaction = step?.text ?? reaction;
  canvas.dataset.stepTexts = `${(shownTexts ?? []).reduce((n, e) => n + (e ? e.above.length + e.below.length : 0), 0)}`;
  if (own)
    canvas.dataset.step = `${step!.index}`;
  let drawnWidth = 0;
  const draw = () => {
    const w = canvas.clientWidth;
    if (w < 1 || w === drawnWidth)
      return;
    drawnWidth = w;
    const h = height(w);
    const r = window.devicePixelRatio;
    canvas.width = Math.round(w * r);
    canvas.height = Math.round(h * r);
    canvas.style.height = `${h}px`;
    drawReactionByCrux(canvas.getContext('2d')!, drawn, 0, 0, w, h,
      {drawOptions: chemDrawOptions(), bondLength, mapNumbers, stepTexts: texts, step: own ? step!.index : undefined});
  };
  const sub = ui.onSizeChanged(canvas).subscribe(draw);
  sizeSubs.set(canvas, sub);
  subs?.push(new Subscription(() => letGo(canvas)));
  return canvas;
}

/** A route's step a card draws: its index (from 0) and its text, as the route's steps read (Crux's get_steps). */
interface RouteStep {
  index: number;
  text: string;
}

/** The route whose steps are `steps`, in Chem's route notation (steps joined by --**--): what the step cards draw step
 * by step, mapped steps among them. */
const routeOf = (steps: readonly string[]) => steps.join('--**--');

/** A step's picture in the full-screen view: its height (CSS px). */
const FULL_SCREEN_STEP_HEIGHT = 220;
/** The full-screen view's picture (at 2549 x 1304 it was drawn at 44 px, small in the middle of the window, and
 * still small at a 96 px cap): fitted to fill this share of its room's limiting
 * dimension, at a bond of at most {@link FULL_SCREEN_MAX_BOND} px, so a small reaction in a large window is drawn
 * large, never huge; its lines a share of its bond ({@link FULL_SCREEN_LINE_SHARE}), as a print-size drawing's are
 * (ACS 1996: 0.6 pt to a 14.4 pt bond, 4%; less reads right on screen), never under Chem's own width. */
export const FULL_SCREEN_FILL = 0.8;
export const FULL_SCREEN_MAX_BOND = 160;
export const FULL_SCREEN_LINE_SHARE = 0.028;
/** The full-screen picture's reaction centre widths (multiples of its line): its line already a share of its bond,
 * thicker than the cells' 0.7 px, Chem's 4x and 2x would be bands (11 px at a 100 px bond); 3x and 1.5x keep the
 * centre in about the panel's proportion to the bond (crux-core decision centre-width-option). */
export const FULL_SCREEN_CENTRE: CentreWidths = {centreWidth: 3, centreLinesWidth: 1.5};
/** The room between the full-screen picture and its legend, right under it (CSS px): the panel's. */
export const FULL_SCREEN_LEGEND_GAP = PICTURE_LEGEND_GAP;
/** A route's full-screen view: the share of the height its picture takes, the rest its steps' list (which scrolls by
 * itself: the dialog never does). */
const FULL_SCREEN_ROUTE_SHARE = 0.55;
/** The full-screen view's class (its flex column under the title: the tests find it by it). */
export const FULL_SCREEN_CLASS = 'chem-reaction-full-screen';

/** The full-screen view's picture and its legend, one group centred in the room a flex column gives `stage`: the
 * canvas the size of the drawing, fitted to the room the legend leaves ({@link cruxFittedPicture}), the legend right
 * under it, centred on it; drawn again whenever the room changes (`data-bond-length`, `data-line-width`: what is
 * drawn, which the tests read). The legend names what the picture and a route's `steps` colour. */
function fittedPicture(reaction: string, mapNumbers?: 'show' | 'hide', subs?: Subscription[],
  stepTexts?: StepTexts | null, steps: string[] = []): HTMLElement {
  const canvas = ui.canvas();
  Object.assign(canvas.style, {display: 'block', flex: 'none'});
  canvas.dataset.reaction = reaction;
  canvas.dataset.stepTexts = `${(stepTexts ?? []).reduce((n, e) => n + (e ? e.above.length + e.below.length : 0), 0)}`;
  const legend = centreLegend([reaction, ...steps]);
  // the panel's group (the legend right under the picture, centred on it), filling the room the column gives it
  const stage = pictureWithLegend(canvas, legend);
  Object.assign(stage.style, {flex: '1 1 0', minHeight: '0', justifyContent: 'center', overflow: 'hidden'});
  let drawnRoom = '';
  const draw = () => {
    const below = legend.style.display === 'none' ? 0 : legend.offsetHeight + FULL_SCREEN_LEGEND_GAP;
    const [room, roomHeight] = [stage.clientWidth, stage.clientHeight - below];
    if (room < 1 || roomHeight < 1 || `${room}x${roomHeight}` === drawnRoom)
      return;
    drawnRoom = `${room}x${roomHeight}`;
    const own = chemDrawOptions();
    const started = performance.now();
    const fit = cruxFittedPicture(reaction, room, roomHeight,
      {drawOptions: own, mapNumbers, stepTexts, centre: FULL_SCREEN_CENTRE},
      FULL_SCREEN_MAX_BOND, FULL_SCREEN_FILL);
    canvas.dataset.fitMs = `${Math.round(performance.now() - started)}`;
    // a route too long for its room at a readable bond: drawn taller, its room scrolling (from its top)
    const overflows = !!fit?.overflows;
    Object.assign(stage.style, {overflowY: overflows ? 'auto' : 'hidden',
      justifyContent: overflows ? 'flex-start' : 'center'});
    const [w, h] = fit ? [Math.min(room, Math.ceil(fit.width)),
      overflows ? Math.ceil(fit.height) : Math.min(roomHeight, Math.ceil(fit.height))] : [room, roomHeight];
    const lineWidth = fit ? Math.max(Number(own.bondLineWidth ?? 1),
      Math.round(fit.bondLength * FULL_SCREEN_LINE_SHARE * 10) / 10) : undefined;
    const drawOptions = lineWidth === undefined ? own : {...own, bondLineWidth: lineWidth};
    canvas.dataset.bondLength = fit ? `${Math.round(fit.bondLength * 10) / 10}` : '';
    canvas.dataset.column = fit?.columnWidth === undefined ? '' : `${fit.columnWidth}`;
    canvas.dataset.lineWidth = lineWidth === undefined ? '' : `${lineWidth}`;
    const r = window.devicePixelRatio;
    Object.assign(canvas.style, {width: `${w}px`, height: `${h}px`});
    canvas.width = Math.round(w * r);
    canvas.height = Math.round(h * r);
    drawReactionByCrux(canvas.getContext('2d')!, reaction, 0, 0, w, h,
      {drawOptions, mapNumbers, stepTexts, bondLength: fit?.bondLength, columnWidth: fit?.columnWidth,
        centre: FULL_SCREEN_CENTRE});
  };
  const sub = ui.onSizeChanged(stage).subscribe(draw);
  sizeSubs.set(canvas, sub);
  subs?.push(new Subscription(() => letGo(canvas)));
  return stage;
}

/** The reaction as large as the window allows, in a dialog that fills the window but for a margin (not a right-shifted
 * dialog whose only button was CANCEL): a view, closed by its x or Escape, with no buttons (as Datagrok's own picture
 * views, PowerGrid's image cell's). Its body never scrolls (the legend was under the fold of a dialog that scrolled, a
 * band of nothing over oversized buttons): a flex column, the picture and its legend one group centred in it, the
 * picture fitted to the room the legend leaves; a route's `steps` (as the panel shows them, mapped once mapped) in a
 * list of their own under it, which scrolls by itself. */
export function openReactionFullScreen(reaction: string, steps: string[] = [],
  mapNumbers?: 'show' | 'hide', stepTexts?: StepTexts | null): DG.Dialog {
  const subs: Subscription[] = [];
  const picture = fittedPicture(reaction, mapNumbers, subs, stepTexts, steps);
  const host = ui.div([picture], FULL_SCREEN_CLASS);
  Object.assign(host.style, {display: 'flex', flexDirection: 'column', flex: '1 1 auto', minHeight: '0',
    width: '100%', overflow: 'hidden'});
  if (steps.length > 1) {
    picture.style.flex = `${FULL_SCREEN_ROUTE_SHARE} 1 0`;
    const list = ui.div([ui.h3('Steps')]);
    Object.assign(list.style, {flex: `${1 - FULL_SCREEN_ROUTE_SHARE} 1 0`, minHeight: '0', overflowY: 'auto'});
    steps.forEach((step, i) => list.append(ui.label(`Step ${i + 1}`),
      reactionCanvas(routeOf(steps), () => FULL_SCREEN_STEP_HEIGHT, FULL_SCREEN_BOND_LENGTH, mapNumbers, subs,
        stepTexts, {index: i, text: step})));
    // the steps' colours are named by the route's legend, right under its picture
    host.append(list);
  }
  const dialog = ui.dialog({title: 'Reaction', helpUrl: reactionsHelpUrl('the-reaction-panel'), showFooter: false})
    .add(host)
    .show({resizable: true, x: FULL_SCREEN_EDGE, y: FULL_SCREEN_EDGE, width: window.innerWidth - 2 * FULL_SCREEN_EDGE,
      height: window.innerHeight - 2 * FULL_SCREEN_EDGE});
  // the contents fill the dialog under its title (Datagrok caps them at 80% of the window, d4.css), with no scroll
  const contents = dialog.root.querySelector('.d4-dialog-contents') as HTMLElement | null;
  if (contents)
    Object.assign(contents.style, {maxHeight: 'none', minHeight: '0', overflow: 'hidden'});
  const closed = dialog.onClose.subscribe(() => {
    subs.forEach((s) => s.unsubscribe());
    closed.unsubscribe();
  });
  return dialog;
}

/** The mapping's flags as lines of the Mapping section (`data-level`: what the tests read), one per flag; `step`
 * names the step they are of. */
function flagLines(flags: Flag[], step?: number): HTMLElement[] {
  return flags.map((f) => {
    const line = ui.divText(step === undefined ? f.text : `Step ${step}: ${f.text}`);
    line.dataset.level = f.level;
    line.style.color = f.level === 'warn' ? 'var(--orange-2)' : f.level === 'block' ? 'var(--red-3)' : '';
    // as wide as the panel, wrapped, never widening it: a long line widened a dialog sized to its content, and with it
    // the panel's picture (in "Reactions: Crux": a template's MAP-02 warn shown as the panel opens)
    line.style.contain = 'inline-size';
    return line;
  });
}

/** The limit on one mapping (5 s): a reaction, or a route's step, the mapper has not answered within it
 * is given up (its late answer ignored; the worker job's own limit is the same, so a job still running past it has its
 * worker replaced). It counts from the mapper being loaded (the model's one-time load is no reaction's time). */
export const MAP_TIMEOUT_MS = 5000;
/** The cap on the mapper's one-time load (its weights read, Chem's Crux workers started: {@link mapperReady}): past
 * it the panel says the mapper could not be loaded, and nothing waits for ever. */
export const MAP_LOAD_TIMEOUT_MS = 30000;
/** How long a panel waits before it maps (ms): a row a user passes by quicker than this is never mapped. */
const MAP_DELAY_MS = 300;
/** How long a route's step cards wait, after a step is mapped, before they are drawn again (ms), where every card is
 * drawn from the mapped route (Crux draws a step alone): steps mapped quicker than this draw the cards once, and the
 * last step draws them at once (5 redraws of a 5-step route's cards, ~85 ms, made one). */
const CARDS_REDRAW_MS = 200;

/** What the Mapping section is doing (`data-state` of `.chem-reaction-mapping`, what tests and bindings wait on):
 * waiting for its turn, loading the mapper (once a page), mapping; then done, given up (timeout), failed (the mapper's
 * error, or the mapper not loaded within its cap), or refused (a template, a value mapped already, an unreadable
 * component, a value too large) without asking the mapper. */
export type MappingState = 'waiting' | 'loading' | 'mapping' | 'done' | 'timeout' | 'failed' | 'refused';
/** The Mapping section's class. */
export const MAPPING_CLASS = 'chem-reaction-mapping';

/** What the Mapping section says when the mapper failed on a reaction (`Not mapped: …`) or a route's step (`Step N:
 * not mapped: …`): in words, the mapper's own text to the console (the panel once said `Step 1: not mapped: "" is
 * not a valid SMILES string`, the engine's text). */
export const MAPPER_FAILED_WORDS = 'Not mapped: the mapper could not map this reaction';
const MAPPER_FAILED_STEP = 'not mapped: the mapper could not map this step';
const MAPPER_NOT_LOADED = 'The mapper could not be loaded';
/** What the Mapping section says of a value Crux reads no reaction from (its |f:| groups out of range, say): the mapper
 * is not asked. */
export const NOT_READ_WORDS = 'Not mapped: the reaction could not be read';

/** {@link NOT_READ_WORDS} and why, as the app's status says it (Crux's words; its |f:| groups' in words), when
 * known. */
export function notReadWords(problem: string | null): string {
  const why = problem ? cruxMessage(problem) : '';
  return why ? `${NOT_READ_WORDS}: ${why.charAt(0).toLowerCase()}${why.slice(1)}` : NOT_READ_WORDS;
}

/** The side a reaction (or a route's step) has nothing on, its products first: a side the mapper cannot map, and is
 * never asked for (it said `"" is not a valid SMILES string`). Null when both sides have a molecule, and for text that
 * is not one `reactants>agents>products` (an RXN block, a route: its steps are judged one by one). */
export function emptySide(reaction: string): 'product' | 'reactant' | null {
  if (reaction.includes('\n'))
    return null;
  const parts = reaction.split(' |')[0].split('>');
  if (parts.length !== 3)
    return null;
  return parts[2].trim() === '' ? 'product' : parts[0].trim() === '' ? 'reactant' : null;
}

/** What the Mapping section says of a reaction, or of a route's step `step` (from 1), with nothing on `side`. */
export function noSideWords(side: 'product' | 'reactant', step?: number): string {
  return `${step === undefined ? 'The reaction' : `Step ${step}`} has no ${side}, so it can't be mapped`;
}

type Mapper = (reaction: string, timeoutMs: number) => Promise<MappedReaction>;
let mapper: Mapper | null = null;
let loader: (() => Promise<void>) | null = null;
let autoMapping = true;
/** The mapper's one-time load, once asked (again after it failed); `loaded` once it has come. */
let loading: Promise<void> | null = null;
let loaded = false;

/** Tests: the mapper the panels call instead of Chem's (RXNMapper on a Crux worker); null: Chem's again. A mapper set
 * here needs no load unless a loader is set too. */
export function setReactionMapperForTests(f: Mapper | null): void {
  mapper = f;
}

/** Tests: the mapper's load (instead of {@link mapperReady}); null: Chem's again. The load is asked again. */
export function setMapperLoaderForTests(f: (() => Promise<void>) | null): void {
  loader = f;
  loading = null;
  loaded = false;
}

/** Tests and measurements: whether the panels map at all (false: the Mapping section stays waiting). */
export function setAutoMappingForTests(on: boolean): void {
  autoMapping = on;
}

/** The panels' mappings, one at a time in the order asked: each waits for the one before, and a panel gone by its turn
 * maps nothing, so passing over rows queues no work for the worker. */
let mappingQueue: Promise<unknown> = Promise.resolve();

/** `p`'s value, or 'late' when it has not come within `ms`; a late value, or a late error, is ignored. */
async function within<T>(p: Promise<T>, ms: number): Promise<T | 'late'> {
  p.catch(() => {}); // an answer given up is nobody's: its error too
  let timer: ReturnType<typeof setTimeout> | undefined;
  const late = new Promise<'late'>((resolve) => timer = setTimeout(() => resolve('late'), ms));
  try {
    return await Promise.race([p, late]);
  } finally {
    clearTimeout(timer);
  }
}

/** Whether the panels must load the mapper before they map (Chem's mapper and not loaded yet, or a test's loader). */
const mapperToLoad = () => !loaded && (!mapper || !!loader);

/** The mapper loaded, or 'late' when its load has not finished within {@link MAP_LOAD_TIMEOUT_MS} (the load goes on:
 * the next panel waits for it again, with a cap of its own); the load's error thrown (and the load asked again next
 * time). */
async function mapperLoaded(): Promise<'ok' | 'late'> {
  if (!mapperToLoad())
    return 'ok';
  loading ??= (loader ?? mapperReady)().then(() => {
    loaded = true;
  }, (e) => {
    loading = null;
    throw e;
  });
  return await within(loading, MAP_LOAD_TIMEOUT_MS) === 'late' ? 'late' : 'ok';
}

/** The mapper's answer for `reaction`, or 'late' when it has not come within {@link MAP_TIMEOUT_MS} (the mapper loaded
 * before); a late answer, or a late error, is ignored. */
function mapWithin(reaction: string): Promise<MappedReaction | 'late'> {
  return within((mapper ?? mapReaction)(reaction, MAP_TIMEOUT_MS), MAP_TIMEOUT_MS);
}

/** A mapping done, in a user's words (not "confidence"): how sure the mapper is, as a percentage. */
/** The Export icon's class (the tests find it by it). */
export const EXPORT_ICON_CLASS = 'chem-reaction-export-icon';

export const mappedWords = (confidence: number, first = false): string =>
  `${first ? 'Mapped' : 'mapped'} automatically (the mapper is ${Math.round(confidence * 100)}% sure)`;

/** The class of the line over a CDXML page's panel saying how many of its molecules are on no step. */
export const NO_STEP_CLASS = 'chem-reaction-no-step';

/** The panel of a CDXML page (rendering/cdxml-reaction-cells.ts): the panel of the reaction Crux Sketch reads
 * from it (its reaction SMILES), once read; a line while it is read, or when it could not be. */
function cdxmlReactionPanel(cdxml: string, options: {mapNumbers?: 'show' | 'hide'}): DG.Widget {
  const host = ui.div([ui.divText(CDXML_DRAWING)]);
  const widget = new DG.Widget(host);
  let inner: DG.Widget | null = null;
  let detached = false;
  const detach = widget.detach.bind(widget);
  widget.detach = () => { // the inner panel's pictures let go with it
    detached = true;
    inner?.detach();
    detach();
  };
  cdxmlReactionSmiles(cdxml).then((smiles) => {
    if (detached)
      return;
    ui.empty(host);
    if (smiles === null)
      host.append(ui.divText(CDXML_UNREAD));
    else {
      // its Export's CDXML the page itself, its texts and arrow styles kept (not a page written anew from the
      // reaction's SMILES)
      inner = reactionPanel(smiles, {...options, cdxmlPage: cdxml});
      // the page's molecules on no step are in no picture of it, and the panel says so
      const left = moleculesOnNoStep(cdxml);
      if (left > 0) {
        const line = ui.divText(`${noStepWords(left)} of the reaction`, NO_STEP_CLASS);
        line.dataset.level = 'warn';
        Object.assign(line.style, {color: 'var(--orange-2)', contain: 'inline-size', margin: '0 0 4px 0'});
        host.append(line);
      }
      host.append(inner.root);
    }
  });
  return widget;
}

/** The panel for one reaction value; `options.mapNumbers` shows or hides map numbers in every picture of it (the
 * Reactions app's switch; left out, Crux's own choice, as the cells draw them); `options.stepTexts`, the texts over
 * and under its steps' arrows as a drawing holds them (the Reactions app's sketcher), drawn in every picture of it and
 * exported with it; `options.cdxmlPage`, the CDXML page `reaction` was read from (a CDXML cell's), which the Export
 * panel saves as its CDXML; `options.asDrawn`, the drawing as drawn when Export is clicked (the Reactions app's
 * sketcher), which the Export panel offers as "As drawn" (left out: the page, or a value with coordinates of its own,
 * as drawn). */
export function reactionPanel(reaction: string, options: {mapNumbers?: 'show' | 'hide', stepTexts?: StepTexts | null,
  cdxmlPage?: string, asDrawn?: () => DrawnPicture | null} = {}): DG.Widget {
  if (isCdxmlText(reaction))
    return cdxmlReactionPanel(reaction, options);
  const {mapNumbers, stepTexts, cdxmlPage} = options;
  const asDrawn = () => options.asDrawn ? options.asDrawn() : cellPicture(cdxmlPage ?? reaction);
  // what the panel holds, let go of by its widget's detach()
  const subs: Subscription[] = [];
  const pictureHeight = (w: number) => Math.max(PICTURE_MIN_HEIGHT, Math.round(w * PICTURE_ASPECT));
  let picture = reactionCanvas(reaction, pictureHeight, PICTURE_BOND_LENGTH, mapNumbers, subs, stepTexts);
  // the steps the picture draws, as Crux reads them (a restated intermediate is no step of its own); a value Crux reads
  // no reaction from, cut as RDKit's renderer cuts it
  const steps = cruxReactionSteps(reaction) ?? parseMultiStepReaction(reaction).flat();
  // what the panel shows now, mapped once mapped: full screen opens the same (full screen once drew the reaction
  // unmapped after mapping)
  let shown = reaction;
  const shownSteps = [...steps];
  const fullScreen = ui.iconFA('expand-alt', () => openReactionFullScreen(shown, shownSteps, mapNumbers, stepTexts),
    'Full screen');
  // the Export panel: what the panel shows, at print size
  const exportIcon = ui.iconFA('file-export', () => openReactionExport(shown, {mapNumbers, stepTexts,
    asDrawn: asDrawn(),
    ...(cdxmlPage ? {cdxml: async () => cdxmlPage, moleculesLeftOut: moleculesOnNoStep(cdxmlPage)} : {})}),
  'Export...');
  exportIcon.classList.add(EXPORT_ICON_CLASS);
  // the reaction centre's colours named right under the large picture, centred on it, as in full screen; a route's
  // steps, coloured as the route is, share it
  const pictureLegend = centreLegend([reaction]);
  const root = ui.divV([ui.divH([exportIcon, fullScreen], {style: {justifyContent: 'flex-end', gap: '8px'}}),
    pictureWithLegend(picture, pictureLegend)]);
  const stepCanvases: HTMLCanvasElement[] = [];
  if (steps.length > 1) {
    root.append(ui.h3('Steps'));
    steps.forEach((step, i) => {
      stepCanvases.push(reactionCanvas(routeOf(steps), () => STEP_HEIGHT, undefined, mapNumbers, subs, stepTexts,
        {index: i, text: step}));
      root.append(ui.label(`Step ${i + 1}`), stepCanvases[i]);
    });
  }
  // a template's notes, which no picture draws (in the context panel, not in the cells)
  const notes = cruxReactionNotes(reaction);
  if (notes.length) {
    const list = ui.divV(notes.map((n) => ui.divText(n)));
    list.classList.add('chem-reaction-notes');
    root.append(ui.h3('Template notes'), list);
  }
  // atom mapping: automatic (5 s a reaction, no button; errors caught; never for a template), on a Crux worker, after
  // MAP_DELAY_MS and one panel at a time; the picture is drawn unmapped at once and swapped for the mapped one when it
  // comes. A single reaction: the picture shows it mapped. A route: RXNMapper maps one reaction at a time, so each step
  // is mapped on its own, one after another (mapping a route whole said "not a valid reaction string"), drawn mapped
  // under Steps, and the route's picture is drawn again from its steps as mapped: what full screen and Export are given
  // (the route stayed unmapped there: Export's Map numbers and centre colours did nothing, and the legend sat under the
  // steps, not under the route's picture).
  if (mapperExported() || mapper) {
    const status = ui.divV([], MAPPING_CLASS);
    // its lines as wide as the panel, wrapped, never widening it (as the flag lines): "Not mapped: a template keeps…"
    // widened a test's dialog, and with it the panel's picture
    status.style.contain = 'inline-size';
    const say =(state: MappingState, lines: (string | HTMLElement)[]) => {
      status.dataset.state = state;
      ui.empty(status);
      lines.forEach((l) => status.append(typeof l === 'string' ? ui.divText(l) : l));
    };
    let gone = false; // the panel's widget detached (another row, another object): its mapping updates nothing
    subs.push(new Subscription(() => gone = true));
    // a value Crux reads no reaction from is said so first, and why (an unreadable typed text was called a
    // template); the mapper, asked for one, said "Cannot read properties of undefined"
    const unread = cruxDrawsReactions() && !tooLargeToDraw(reaction) && !cruxReadsReaction(reaction);
    // MAP-02 (Crux mode): a template's map numbers are its author's, what a run follows; RXNMapper reads SMILES only
    const templateFlags = unread ? [] : checkTemplateMapping(reaction);
    const template = templateFlags.length > 0;
    // a value with a mapping of its own keeps it (automatic mapping must not replace a mapping the value already has;
    // any mapped atom counts): drawn as written, the mapper never asked
    const ownMapping = !unread && !template && mappedInTheValue(reaction);
    // MAP-03 (Crux mode): a component no molecule reader reads is named, before the mapper's own words would
    const unreadable = checksModule() && !unread && !template && !ownMapping ?
      checkReaction(reaction, {template: false}).filter((f) => f.rule === 'componentUnreadable') : [];
    root.append(ui.h3('Mapping'), status);
    if (unread)
      say('refused', [notReadWords(cruxReadProblem(reaction))]);
    else if (template)
      say('refused', [...flagLines(templateFlags), 'Not mapped: a template keeps the map numbers written in it']);
    else if (ownMapping)
      say('refused', ['Mapped in the value itself']);
    else if (unreadable.length)
      say('refused', [...flagLines(unreadable), 'Not mapped: the mapper reads molecules only']);
    else if (tooLargeToDraw(reaction))
      say('refused', ['Not mapped: too large to map here']);
    else if (steps.length <= 1 && emptySide(reaction) !== null)
      say('refused', [noSideWords(emptySide(reaction)!)]);
    else {
      say('waiting', ['Mapping...']);
      const mapAll = async () => {
        // the mapper's one-time load, capped on its own; each reaction's 5 s counts from after it
        if (mapperToLoad())
          say('loading', ['Loading the mapper…']);
        try {
          if (await mapperLoaded() === 'late') {
            if (!gone)
              say('failed', ['The mapper could not be loaded']);
            return;
          }
        } catch (e) {
          console.warn('Chem | the mapper could not be loaded', e);
          if (!gone)
            say('failed', [MAPPER_NOT_LOADED]);
          return;
        }
        if (gone)
          return;
        say('mapping', ['Mapping...']);
        if (stepCanvases.length) {
          const lines: (string | HTMLElement)[] = [];
          let late = false;
          let failed = false;
          // the cards drawn again as the steps of the route mapped so far (`from`: the card to start at), where Crux
          // draws a step alone (one layout and one scale for all); else only the step mapped, by itself
          const drawCards = (from: number, to: number) => {
            for (let j = from; j < to; j++) {
              const next = reactionCanvas(routeOf(shownSteps), () => STEP_HEIGHT, undefined, mapNumbers, subs,
                stepTexts, {index: j, text: shownSteps[j]});
              letGo(stepCanvases[j]);
              stepCanvases[j].replaceWith(next);
              stepCanvases[j] = next;
            }
          };
          // the route's picture drawn again from its steps as mapped so far, coloured as the cards are; it is what the
          // panel shows now (full screen, Export), its legend right under it
          const drawRoute = () => {
            shown = routeOf(shownSteps);
            const next = reactionCanvas(shown, pictureHeight, PICTURE_BOND_LENGTH, mapNumbers, subs, stepTexts);
            letGo(picture);
            picture.replaceWith(next);
            picture = next;
            setCentreLegend(pictureLegend, [shown]);
          };
          let cardsTimer: ReturnType<typeof setTimeout> | null = null;
          const drawCardsNow = () => {
            if (cardsTimer !== null)
              clearTimeout(cardsTimer);
            cardsTimer = null;
            if (gone)
              return;
            drawCards(0, stepCanvases.length);
            drawRoute();
          };
          for (let i = 0; i < steps.length; i++) {
            const side = emptySide(steps[i]);
            if (side !== null) // skipped, said in words: no mapper asked, no failure
              lines.push(noSideWords(side, i + 1));
            else {
              try {
                const answer = await mapWithin(steps[i]);
                if (gone)
                  return;
                if (answer === 'late') {
                  late = true;
                  lines.push(`Step ${i + 1}: mapping took longer than ${MAP_TIMEOUT_MS / 1000} s`);
                } else {
                  shownSteps[i] = answer.mapped;
                  if (!cruxDrawsSteps()) {
                    drawCards(i, i + 1);
                    drawRoute();
                  } else {
                    if (cardsTimer !== null)
                      clearTimeout(cardsTimer);
                    cardsTimer = setTimeout(drawCardsNow, CARDS_REDRAW_MS);
                  }
                  lines.push(`Step ${i + 1}: ${mappedWords(answer.confidence)}`,
                    ...flagLines(checkMapping(answer.confidence, productAtomsMissing(steps[i])), i + 1));
                }
              } catch (e) {
                if (gone)
                  return;
                failed = true;
                console.warn(`Chem | the mapper did not map step ${i + 1}`, e);
                lines.push(`Step ${i + 1}: ${MAPPER_FAILED_STEP}`);
              }
            }
            if (i + 1 < steps.length)
              say('mapping', [...lines, 'Mapping...']);
            else {
              if (cardsTimer !== null)
                drawCardsNow();
              say(late ? 'timeout' : failed ? 'failed' : 'done', lines);
            }
          }
          return;
        }
        try {
          const answer = await mapWithin(reaction);
          if (gone)
            return;
          if (answer === 'late') {
            say('timeout', [`Mapping took longer than ${MAP_TIMEOUT_MS / 1000} s`]);
            return;
          }
          const next = reactionCanvas(answer.mapped, pictureHeight, PICTURE_BOND_LENGTH, mapNumbers, subs, stepTexts);
          letGo(picture);
          picture.replaceWith(next);
          picture = next;
          shown = answer.mapped;
          setCentreLegend(pictureLegend, [shown]);
          // (MAP-05's note, a value's own map numbers replaced, cannot arise: such a value is never mapped)
          say('done', [mappedWords(answer.confidence, true),
            ...flagLines(checkMapping(answer.confidence, productAtomsMissing(reaction)))]);
        } catch (e) {
          console.warn('Chem | the mapper did not map the reaction', e);
          if (!gone)
            say('failed', [MAPPER_FAILED_WORDS]);
        }
      };
      if (autoMapping) {
        const timer = setTimeout(() => {
          if (!gone)
            mappingQueue = mappingQueue.then(() => gone ? undefined : mapAll()).catch(() => {});
        }, MAP_DELAY_MS);
        subs.push(new Subscription(() => clearTimeout(timer)));
      }
    }
  }
  const widget = new DG.Widget(root);
  // every picture, those made later (a mapped one) too
  widget.subs.push(new Subscription(() => subs.forEach((s) => s.unsubscribe())));
  return widget;
}
