/**
 * Reaction drawing by Crux. In Crux mode (the "Chemical engine" setting, engine/chem-engine.ts), the reaction cell
 * renderer (rdkit-reaction-renderer.ts) and the Reaction panel (widgets/reaction-panel.ts) draw a reaction through the
 * scheme painter in Chem's Crux file (`ChemReactionScheme`), synchronously, as canvas commands on the canvas they are
 * given: no SVG, no image, no promise. The molecules inside a reaction are drawn with Chem's molecule draw options, the
 * ones molecule cells use (`drawOptions`). RDKit mode keeps RDKit's drawing as it was; until the Crux module is ready
 * RDKit draws, as for molecules (engine/README.md).
 */
import * as DG from 'datagrok-api/dg';
import {ChemicalEngine, engineInEffect} from '../engine/chem-engine';
import {ChemWasm, cruxLoadFailed, cruxModuleSync, loadCrux, noteCruxError} from '../engine/crux-module';
import {createRenderingOpts} from '../utils/chem-common-rdkit';
import {renderingOptions} from './rendering-settings';

/** How many read reactions are kept, one per distinct value drawn. A drawn one keeps its layout for the next paint of
 * the same box (0.3 ms instead of 2.5 ms), and with it ~5.5 KB a molecule plus ~1.2 KB an atom of the Crux heap, which
 * never shrinks: 1,000 hold ~100 MB. 200 would hold ~20 MB but made a scroll back up 68% slower (measured in node, the
 * scroll of 2,000 distinct reactions, and by leak-tests.ts's scroll benchmark in a browser), so the size
 * stays and the engine is asked to keep less a drawn scheme (a display list instead of the laid-out molecules). */
const SCHEME_CACHE_SIZE = 1000;
/** The largest reaction drawn, the cells and the panel drawing synchronously: every molecule's text at most
 * MAX_DRAWN_MOLECULE_LENGTH characters, the whole at most MAX_DRAWN_REACTION_LENGTH (a reaction SMILES, SMARTS, route
 * or step list; an RXN block is not judged, as Chem judges a molecule's SMILES only). A molecule's layout costs more
 * than its size: in node, a peptide's acetylation drew in 0.47 s at 250 residues (2,508 characters a
 * molecule), 2.1 s at 500, 6.8 s at 1,000; a 120-step route of drug-like molecules (5,300 characters) in 0.44 s.
 * A larger one is said to be too large (catalogue RXN-08, PERF-03: answer or refuse within 1 s). */
const MAX_DRAWN_MOLECULE_LENGTH = 2500;
const MAX_DRAWN_REACTION_LENGTH = 10000;
/** The reaction cells' muted text colour (a reaction written as text, too large to draw; the CDXML cells' messages). */
export const REACTION_CELL_TEXT_COLOUR = '#7f7f7f';
/** A value Crux reads no reaction from is written as its text, in this colour and font (CSS px). */
const TEXT_COLOUR = REACTION_CELL_TEXT_COLOUR;
const TEXT_FONT = '11px Roboto, sans-serif';
/** Agent names Chem adds to Crux's own table or replaces in it, by the canonical SMILES of the
 * whole agent: `{"<SMILES>": {"name": "...", "place": "above" | "below"}}`. Chem's own agent files fill it; there are
 * none yet. Crux reads it as the drawing option `agentNames` once its agent table exists, and ignores it until then. */
const CHEM_AGENT_NAMES = {};

type Scheme = InstanceType<ChemWasm['ChemReactionScheme']>;

/** A cache of `size` read reactions that frees the ones it evicts. */
function schemeCache(size: number): DG.LruCache<string, Scheme | null> {
  const cache = new DG.LruCache<string, Scheme | null>(size);
  cache.onItemEvicted = (scheme: Scheme | null) => scheme?.free();
  return cache;
}

/** The reactions read so far, by their text; null for a text Crux reads no reaction from. */
let schemes = schemeCache(SCHEME_CACHE_SIZE);

/** Measurements and tests: the reactions read so far freed, and the cache made again with room for `size`
 * (by default its own), so that a scroll can be timed at another cache size in the same session. */
export function resetCruxSchemeCache(size: number = SCHEME_CACHE_SIZE): void {
  const kept = (schemes as unknown as {V: (Scheme | null | undefined)[]}).V; // the cache's values (DG.LruCache)
  schemes = schemeCache(size);
  kept.forEach((scheme) => scheme?.free());
}

/** How many reactions Crux has drawn, how many values it drew as text, and how many it said were too large (tests). */
export const cruxReactionCounts = {drawn: 0, asText: 0, tooLarge: 0};

/** Whether `text` is too large for Crux to draw synchronously (MAX_DRAWN_MOLECULE_LENGTH): a molecule of it over
 * 2,500 characters, or the whole over 10,000, its CX part left out; an RXN block (a text with line breaks) never is. */
export function tooLargeToDraw(text: string): boolean {
  if (text.includes('\n'))
    return false;
  const body = text.split(' |')[0];
  return body.length > MAX_DRAWN_REACTION_LENGTH ||
    body.split(/[.>]/).some((molecule) => molecule.length > MAX_DRAWN_MOLECULE_LENGTH);
}

/** Whether Crux draws reactions now: Crux is the engine in effect (in Crux mode, once its module is ready). */
export function cruxDrawsReactions(): boolean {
  return engineInEffect() === ChemicalEngine.Crux;
}

/** Chem's Crux engine for the reaction features that are Crux's in both modes (the reaction editor and its menu's
 * Export work whatever the Chemical Engine): loaded once if RDKit mode has not loaded it, which
 * changes nothing else RDKit mode does (cells, panels and checks follow the engine in effect). False when it cannot be
 * loaded, or has failed in this session. */
export async function cruxForReactions(reason: string): Promise<boolean> {
  if (cruxLoadFailed())
    return false;
  if (cruxModuleSync())
    return true;
  try {
    await loadCrux(reason);
    return !cruxLoadFailed() && !!cruxModuleSync();
  } catch (e) {
    console.error(`Chem | Crux could not be loaded for ${reason}`, e);
    return false;
  }
}

/** The main thread's Crux module, or null until it is loaded or once it failed: after a trap (crux-module.ts
 * noteCruxError) the instance may be broken, and the reaction calls call it no more. */
function cruxNow(): ChemWasm | null {
  return cruxLoadFailed() ? null : cruxModuleSync() as unknown as ChemWasm | null;
}

/** The reaction Crux reads from `text` (in the notation it detects), kept for the next draw; null if it reads none. */
function schemeOf(crux: ChemWasm, text: string): Scheme | null {
  return schemes.getOrCreate(text, () => {
    try {
      return new crux.ChemReactionScheme(text, '');
    } catch (e) {
      noteCruxError(e); // a trap the module could not recover from: Crux is off for the session
      return null;
    }
  });
}

/** What a caller chooses about a picture besides its size: Chem's molecule draw options ({@link chemDrawOptions}), and
 * the bond length (px), Crux's own (22) when unset. */
export interface CruxPictureOptions {
  drawOptions: {[key: string]: any};
  bondLength?: number;
  /** Map numbers beside their atoms: `show` or `hide` (the Reactions app's switch); left out, Crux's own choice, as the
   * cells draw them. */
  mapNumbers?: 'show' | 'hide';
  /** The texts over and under each step's arrow, as the drawing holds them ({@link StepTexts}); left out, none. */
  stepTexts?: StepTexts | null;
  /** A route's step k alone (from 0), drawn with the route's own coordinates, every step of it at one scale in a box
   * (chem-wasm's `step`; {@link cruxDrawsSteps}); left out, the whole reaction. */
  step?: number;
  /** The scheme broken into rows as a page's column of this width is (px; chem-wasm's `column`, the Export panel's
   * Width): a long route wraps instead of shrinking into one row; left out, the box's own layout. */
  columnWidth?: number;
  /** The reaction centre's widths; left out, {@link CHEM_CENTRE} (full screen's lines are thicker, its centre less so
   * against them). */
  centre?: CentreWidths;
}

/** A text over or under a step's arrow as a drawing holds it (crux-sketch's `ReactionStepText`): its place,
 * its plain text, and the engine's drawing-JSON text object (its runs, faces, size and soft wrap), which Crux draws as
 * it is. */
export interface DrawnStepText {
  readonly place: 'above' | 'below';
  readonly text: string;
  readonly drawing?: Readonly<Record<string, unknown>>;
}

/** A step's texts as a drawing holds them (crux-sketch's `ReactionStepTexts`): `step` from 1, the scheme's k-th
 * step. */
export interface DrawnStepTexts {
  readonly step: number;
  readonly texts: readonly DrawnStepText[];
}

/** Crux's `stepTexts` option (chem-wasm's scheme painter): one entry a step, by the scheme's step index, its texts
 * over and under the arrow (each a drawing-JSON text object, or plain text); null for a step with none. */
export type StepTexts = ({above: unknown[], below: unknown[]} | null)[];

/** The `stepTexts` option of a drawing's step texts (each text's drawing object passed on untouched, its plain text
 * where it has none); null when no step has a text. */
export function stepTextsOf(steps: readonly DrawnStepTexts[] | null | undefined): StepTexts | null {
  const out: StepTexts = [];
  for (const s of steps ?? []) {
    if (!Number.isInteger(s?.step) || s.step < 1 || !s.texts?.length)
      continue;
    while (out.length < s.step)
      out.push(null);
    const entry = out[s.step - 1] ?? {above: [], below: []};
    for (const t of s.texts)
      (t.place === 'below' ? entry.below : entry.above).push(t.drawing ?? t.text);
    out[s.step - 1] = entry;
  }
  return out.some((e) => e !== null) ? out : null;
}

/** The `stepTexts` option for step `k` (from 0) drawn as a reaction of its own: its entry alone; undefined (none) when
 * there are no texts. */
export function stepTextsOfStep(texts: StepTexts | null | undefined, k: number): StepTexts | undefined {
  return texts ? [texts[k] ?? null] : undefined;
}

/** Whether Chem's Crux file draws a host's step texts (chem-wasm's `stepTexts`, which came with its
 * `get_option_notes`; an engine before it draws none, and the texts stay CDXML's alone). */
export function cruxDrawsStepTexts(): boolean {
  const crux = cruxNow();
  return !!crux && typeof (crux.ChemReactionScheme.prototype as {get_option_notes?: unknown}).get_option_notes ===
    'function';
}

/** Chem's molecule draw options for the cells of `col` (the package's without one): the options molecule cells draw
 * with (rdkit-cell-renderer.ts), the common ones and the Rendering settings. Crux draws a reaction's molecules with
 * them, and its arrow, plus signs and step numbers at their line width and label font. */
export function chemDrawOptions(col?: DG.Column | null): {[key: string]: any} {
  return createRenderingOpts(renderingOptions(col));
}

/** The reaction centre's widths Chem draws with, as multiples of its bond line (chem-wasm's `centreWidth`, a bond of one
 * line, and `centreLinesWidth`, each line of a double, triple or aromatic one). Crux's own, 2x and 1x, were chosen in
 * its screen look, whose line is 2 px; at Chem's 0.7 px they were hairlines (1.4 px, 0.7 px: the owner, "still too
 * thin"). 4x and 2x: 2.8 px and 1.4 px, and a double bond's two lines stay apart even in a cell (crux-core decision
 * centre-width-option: they touch from 5x and 2.5x). An engine before the keys ignores them. */
export const CHEM_CENTRE: CentreWidths = {centreWidth: 4, centreLinesWidth: 2};

/** A reaction centre's widths (multiples of the bond line): {@link CHEM_CENTRE}. */
export interface CentreWidths {
  centreWidth: number;
  centreLinesWidth: number;
}

/** The drawing options that change the picture, besides its size: the caller's, Chem's agent names, its reaction
 * centre's widths, and no template footnotes: Chem draws a template's notes in no picture, the Reaction panel lists
 * them as text (they belong in the context panel when the reaction is clicked, not in the cells;
 * {@link cruxReactionNotes}). */
function pictureOptions(options: CruxPictureOptions): object {
  return {bondLength: options.bondLength, agentNames: CHEM_AGENT_NAMES, drawOptions: options.drawOptions,
    footnotes: false, mapNumbers: options.mapNumbers, stepTexts: options.stepTexts ?? undefined, step: options.step,
    ...(options.centre ?? CHEM_CENTRE),
    ...(options.columnWidth === undefined ? {} : {column: true, width: options.columnWidth})};
}

/** Whether Chem's Crux file draws a route's step alone from the route's own layout (chem-wasm's `step` option,
 * so that an intermediate is drawn the same way in both its steps): a step past the last throws there, and an
 * engine before it ignores the key. Asked once a module. */
export function cruxDrawsSteps(): boolean {
  const crux = cruxNow();
  if (!crux)
    return false;
  if (drawsSteps === null || drawsStepsOf !== crux) {
    const scheme = schemeOf(crux, 'C>>CC');
    drawsStepsOf = crux;
    try {
      scheme?.get_size(JSON.stringify({step: STEP_PAST_THE_LAST}));
      drawsSteps = false;
    } catch {
      drawsSteps = scheme !== null;
    }
  }
  return drawsSteps;
}
/** {@link cruxDrawsSteps}'s answer and the module it is of; the step index it asks for. */
let drawsSteps: boolean | null = null;
let drawsStepsOf: ChemWasm | null = null;
const STEP_PAST_THE_LAST = 1000000;

/** The bond length (px) at which `reaction`, drawn with `options`, fills `fill` of the box `w` x `h` (CSS px) in its
 * limiting dimension, at most `maxBond`, and the drawing's size at it: Crux shrinks a scheme too large for its box but
 * never enlarges one, so a picture meant to fill a large box (the Reaction panel's full-screen view) asks for the bond
 * that fits, worked out from the scheme's size at Crux's own bond. Null where Crux reads no scheme or cannot size it
 * (Crux's own bond then). */
export function cruxFittedPicture(reaction: string, w: number, h: number, options: CruxPictureOptions,
  maxBond: number, fill: number = FIT_MARGIN): FittedPicture | null {
  const crux = cruxNow();
  if (!crux || w <= 0 || h <= 0 || tooLargeToDraw(reaction))
    return null;
  const scheme = schemeOf(crux, reaction);
  if (!scheme)
    return null;
  try {
    if (scheme.get_steps().length > 1)
      return fittedRoute(scheme, w, h, options, maxBond, fill);
    const [sw, sh] = Array.from(scheme.get_size(JSON.stringify({width: w, height: h,
      ...pictureOptions({...options, bondLength: FIT_BOND})}))) as number[];
    if (!(sw > 0 && sh > 0))
      return null;
    const bondLength = Math.max(1, Math.min(maxBond, FIT_BOND * Math.min(fill, FIT_MARGIN) * Math.min(w / sw, h / sh)));
    // a little room besides: a scheme's texts are not scaled with its bond
    const scale = bondLength / FIT_BOND;
    return {bondLength, width: sw * scale + FIT_PAD, height: sh * scale + FIT_PAD};
  } catch {
    return null;
  }
}
/** A picture fitted to a box ({@link cruxFittedPicture}): its bond and its size (CSS px); a route's column (its
 * width, which the drawing is given as `columnWidth`), and whether it is taller than the box. */
export interface FittedPicture {
  bondLength: number;
  width: number;
  height: number;
  columnWidth?: number;
  overflows?: boolean;
}

/** A route fitted to the box `w` x `h` (a long route fitted as one row was drawn at a 4-11 px bond): broken into rows
 * as a page's column of the box's width (`fill` of it) is, at the largest bond from {@link ROUTE_FLOOR_BOND} up to
 * `maxBond` whose rows fit {@link ROUTE_HEIGHT_FILL} of its height (8 sizes asked at most). Where even the floor's rows
 * are too tall, or too wide, it is drawn at {@link ROUTE_MIN_BOND} (the full-screen view's bond before it was fitted),
 * or the largest bond whose rows fit the column, taller than the box, whose room then scrolls: never one unreadable
 * row. */
function fittedRoute(scheme: Scheme, w: number, h: number, options: CruxPictureOptions, maxBond: number,
  fill: number): FittedPicture {
  const columnWidth = Math.floor(w * Math.min(fill, FIT_MARGIN));
  const room = h * ROUTE_HEIGHT_FILL;
  const sizeAt = (bondLength: number) => {
    const [cw, ch, percent] = Array.from(scheme.get_size(JSON.stringify(
      pictureOptions({...options, bondLength, columnWidth})))) as number[];
    return {width: cw, height: ch, fitsWidth: percent === undefined || percent >= 99.5};
  };
  const picture = (bondLength: number, at: {height: number}) => ({bondLength, width: columnWidth + FIT_PAD,
    height: at.height + FIT_PAD, columnWidth, overflows: at.height > room});
  const fits = (at: {height: number, fitsWidth: boolean}) => at.fitsWidth && at.height <= room;
  const floor = Math.min(ROUTE_FLOOR_BOND, maxBond);
  const floorSize = sizeAt(floor);
  if (!fits(floorSize)) { // too long for the room even at the floor: the old bond, or the largest whose rows fit
    let bond = Math.min(ROUTE_MIN_BOND, maxBond);
    let at = sizeAt(bond);
    if (!at.fitsWidth) {
      let [lo, hi] = [ROUTE_LEAST_BOND, bond];
      for (let i = 0; i < ROUTE_SEARCH_STEPS - 1 && hi - lo > 1; i++) {
        const mid = (lo + hi) / 2;
        if (sizeAt(mid).fitsWidth)
          lo = mid;
        else
          hi = mid;
      }
      bond = lo;
      at = sizeAt(bond);
    }
    return picture(bond, at);
  }
  const top = sizeAt(maxBond);
  if (fits(top))
    return picture(maxBond, top);
  let [lo, hi, best] = [floor, maxBond, floorSize];
  for (let i = 0; i < ROUTE_SEARCH_STEPS - 1 && hi - lo > 1; i++) {
    const mid = (lo + hi) / 2;
    const at = sizeAt(mid);
    if (fits(at)) {
      lo = mid;
      best = at;
    } else
      hi = mid;
  }
  return picture(lo, best);
}
/** {@link fittedRoute}'s least bond where the room allows it (px; at 22 px the default route was
 * too tall for a 1440 x 900 window's room, which scrolled, its legend under the fold), its bond where the room does not
 * (the full-screen view's before it was fitted), its least at all (a row's piece wider than the column even so), the
 * sizes it searches, and the share of the room's height its rows fill (a single reaction's {@link cruxFittedPicture}
 * `fill` is a look; a route's rows take what room there is). */
const ROUTE_FLOOR_BOND = 16;
const ROUTE_HEIGHT_FILL = 0.95;
const ROUTE_MIN_BOND = 44;
const ROUTE_LEAST_BOND = 8;
const ROUTE_SEARCH_STEPS = 7;

/** The bond {@link cruxFittedPicture} sizes a scheme at (Crux's own, px), the largest share of the box it fills (a
 * scheme's size at one bond is not exactly in proportion at another: its texts are not scaled), and the room it adds
 * to the drawing's size (CSS px). */
const FIT_BOND = 22;
const FIT_MARGIN = 0.95;
const FIT_PAD = 8;

/** Whether Crux reads a reaction from `text` (a route as one chain too); false until Chem's Crux file is loaded. */
export function cruxReadsReaction(text: string): boolean {
  const crux = cruxNow();
  return !!crux && schemeOf(crux, text) !== null;
}

/** Element symbols by atomic number - 1, as far as a reaction's product atoms are named
 * ({@link productAtomsMissing}). */
const SYMBOLS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ' +
  'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba').split(' ');

/** What {@link productAtomsMissing} adds to its atoms when they are one reactant's or agent's atoms once more, or more
 * times over: the step uses it more than once, which reaction SMILES does not write. */
export const REACTANT_REPEATED = '|repeated';

/** The heavy atoms one reaction step's products hold beyond its reactants and agents together, by element ("C 2,
 * O 1"; an N-acetylation with no acetyl source drawn was mapped onto a ring carbon, 63% sure), followed
 * by {@link REACTANT_REPEATED} when they are a whole number of times one reactant's or agent's heavy atoms, its other
 * elements those the products hold fewer of (else a self-condensation `CC=O>>CC(O)CC=O`, a homocoupling
 * `Brc1ccccc1>>c1ccc(-c2ccccc2)cc1`, a double acylation, each drawn once, were told to add a reagent); null when they
 * hold none more, when `step` is no `reactants>agents>products` (an RXN block, a route: its steps are judged one by
 * one), when a molecule does not read (a template), or when Chem's Crux file is not loaded. */
export function productAtomsMissing(step: string): string | null {
  const crux = tooLargeToDraw(step) || step.includes('\n') ? null : cruxNow();
  const sides = step.split(' |')[0].trim().split('>');
  if (!crux || sides.length !== 3)
    return null;
  const balance = new Map<number, number>();
  /** Each reactant's and agent's heavy atoms, by element. */
  const drawn: Map<number, number>[] = [];
  const count = (side: string, sign: number): boolean => {
    for (const smiles of side.split('.').filter((s) => s.trim())) {
      let mol: ReturnType<ChemWasm['get_mol']> = null;
      const own = new Map<number, number>();
      try {
        mol = crux.get_mol(smiles, '{}');
        if (!mol)
          return false;
        const json = JSON.parse(mol.get_json());
        const z0 = json.defaults?.atom?.z ?? 6;
        for (const m of json.molecules ?? []) {
          for (const a of m.atoms ?? []) {
            const z = a.z ?? z0;
            if (z > 1) {
              balance.set(z, (balance.get(z) ?? 0) + sign);
              own.set(z, (own.get(z) ?? 0) + 1);
            }
          }
        }
      } catch (e) {
        noteCruxError(e);
        return false;
      } finally {
        if (!cruxLoadFailed())
          mol?.delete();
      }
      if (sign < 0 && own.size)
        drawn.push(own);
    }
    return true;
  };
  if (!count(sides[2], 1) || !count(sides[0], -1) || !count(sides[1], -1))
    return null;
  const more = [...balance.entries()].filter(([, n]) => n > 0).sort((a, b) => a[0] - b[0]);
  if (!more.length)
    return null;
  // the atoms more are k times one molecule drawn: each of its elements k times over, but one the products hold fewer
  // of (a leaving group's: bromobenzene's Br in a homocoupling); and they hold no element it has not
  const repeated = drawn.some((own) => {
    const [z1, n1] = more[0];
    const k = (own.get(z1) ?? 0) > 0 ? n1 / own.get(z1)! : 0;
    return Number.isInteger(k) && k >= 1 && more.every(([z, n]) => own.get(z) === n / k) &&
      [...own.entries()].every(([z, n]) => (balance.get(z) ?? 0) === k * n || (balance.get(z) ?? 0) < 0);
  });
  return more.map(([z, n]) => `${SYMBOLS[z - 1] ?? `element ${z}`} ${n}`).join(', ') +
    (repeated ? REACTANT_REPEATED : '');
}

/** Why Crux reads no reaction from `text`, in its own words ("not a reaction"), or null when it reads one, or when
 * Chem's Crux file is not loaded (without it the Reaction panel called an unreadable text a template). */
export function cruxReadProblem(text: string): string | null {
  const crux = tooLargeToDraw(text) ? null : cruxNow();
  if (!crux || schemeOf(crux, text) !== null)
    return null;
  try {
    new crux.ChemReactionScheme(text, '').free();
    return null;
  } catch (e) {
    noteCruxError(e);
    return String(e instanceof Error ? e.message : e).replace(/^ChemReactionScheme: /, '').trim() || null;
  }
}

/**
 * The steps of `text` as Crux reads them, in order, each a reaction of its own (`ChemReactionScheme.get_steps`): a
 * chain cut by the restated-intermediate rule (`A.B>>C>>C.D>>E` is two steps, as its picture draws it), a step list at
 * its `--**--`, one reaction as one step. Null when Chem's Crux file is not loaded or reads no reaction from `text`.
 */
export function cruxReactionSteps(text: string): string[] | null {
  // a value too large to draw is not read for its steps either: get_steps grows with the square of a molecule's size
  // (a 1,000-residue acetylation: 0.4 s in node, 0.73 s in a browser, the context panel frozen)
  const crux = tooLargeToDraw(text) ? null : cruxNow();
  const scheme = crux ? schemeOf(crux, text) : null;
  return scheme ? scheme.get_steps() : null;
}

/**
 * A template's notes, as Crux words them (`ChemReactionScheme.get_notes`): what a generic atom or a query feature
 * stands for ("X: one of Cl, Br, I", "the amine's N takes 1 or 2 H"), one string each, which Chem's pictures leave out
 * and the Reaction panel lists. Empty for a reaction with none, or when Chem's Crux file is not loaded or reads no
 * reaction from `text`.
 */
export function cruxReactionNotes(text: string): string[] {
  const crux = tooLargeToDraw(text) ? null : cruxNow(); // as cruxReactionSteps: not read for its notes either
  const scheme: any = crux ? schemeOf(crux, text) : null;
  return scheme && typeof scheme.get_notes === 'function' ? scheme.get_notes() : [];
}

/**
 * The bonds Crux's picture of `text` colours as the reaction centre, by colour: `[made, broken, changed]`
 * (`ChemReactionScheme.get_centre_counts`; all 0 for an unmapped reaction), for the legend under a picture
 * (widgets/reaction-centre-legend.ts). Null when Chem's Crux file is not loaded, reads no reaction from `text`, or does
 * not count them yet.
 */
export function cruxCentreCounts(text: string): number[] | null {
  const crux = tooLargeToDraw(text) ? null : cruxNow(); // as cruxReactionSteps: not read for its counts either
  const scheme: any = crux ? schemeOf(crux, text) : null;
  if (!scheme || typeof scheme.get_centre_counts !== 'function')
    return null;
  try {
    return Array.from(scheme.get_centre_counts() as Int32Array);
  } catch (e) {
    noteCruxError(e);
    return null;
  }
}

/**
 * Draws `reaction` into the box `x, y, w, h` (CSS px) of `g`'s canvas through Crux, clipped to the box; the canvas
 * holds `devicePixelRatio` device pixels per CSS px, as the grid's and the panel's do. A value Crux reads no reaction
 * from, or cannot draw, is written as its text in grey (where RDKit's renderer drew a grey cross). Crux draws at its
 * own bond length (22 px) and shrinks a scheme that does not fit, never enlarges it; `options.bondLength` (px) asks
 * for another.
 */
export function drawReactionByCrux(g: CanvasRenderingContext2D, reaction: string,
  x: number, y: number, w: number, h: number, options: CruxPictureOptions): void {
  draw(g, reaction, x, y, w, h, window.devicePixelRatio, options);
}

/**
 * The key of `reaction`'s picture of `w` x `h` device px at `ratio` device px per CSS px: the text, the size, the
 * ratio and the options that change the picture, Chem's draw options among them, so that a change of the Rendering
 * settings draws anew (the reaction cell renderer keeps pictures by it, as it does RDKit's).
 */
export function cruxPictureKey(reaction: string, w: number, h: number, ratio: number,
  options: CruxPictureOptions): string {
  return `${w} || ${h} || ${ratio} || ${JSON.stringify(pictureOptions(options))} || ${reaction}`;
}

/**
 * `reaction` drawn by Crux on `canvas` (cleared first), `w` x `h` device px at `ratio` device px per CSS px, and
 * returned as its pixels: the reaction cell renderer keeps them and copies them on the next paint of the cell
 * (rdkit-reaction-renderer.ts), as it does RDKit's pictures.
 */
export function cruxReactionPixels(canvas: OffscreenCanvas, reaction: string, w: number, h: number,
  ratio: number, options: CruxPictureOptions): ImageData {
  const g = canvas.getContext('2d', {willReadFrequently: true})! as unknown as CanvasRenderingContext2D;
  g.save();
  g.setTransform(1, 0, 0, 1, 0, 0);
  g.clearRect(0, 0, w, h);
  g.restore();
  draw(g, reaction, 0, 0, w / ratio, h / ratio, ratio, options);
  return g.getImageData(0, 0, w, h);
}

/** drawReactionByCrux's drawing, at `ratio` device px per CSS px. */
function draw(g: CanvasRenderingContext2D, reaction: string,
  x: number, y: number, w: number, h: number, ratio: number, options: CruxPictureOptions): void {
  const crux = cruxNow();
  if (!crux || w <= 0 || h <= 0)
    return;
  g.save();
  try {
    g.setTransform(ratio, 0, 0, ratio, 0, 0);
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    if (tooLargeToDraw(reaction)) {
      g.font = TEXT_FONT;
      g.fillStyle = TEXT_COLOUR;
      g.textBaseline = 'middle';
      g.fillText(`Too large to draw here: ${reaction.length.toLocaleString('en-US')} characters`, x + 4, y + h / 2);
      cruxReactionCounts.tooLarge++;
      return;
    }
    const scheme = schemeOf(crux, reaction);
    try {
      if (!scheme)
        throw new Error('no reaction read');
      scheme.draw_to_canvas(g.canvas,
        JSON.stringify({width: w, height: h, offsetx: x, offsety: y, ...pictureOptions(options)}));
      cruxReactionCounts.drawn++;
    } catch (e) {
      noteCruxError(e); // a trap the module could not recover from: Crux off, RDKit draws the next paint
      g.font = TEXT_FONT;
      g.fillStyle = TEXT_COLOUR;
      g.textBaseline = 'middle';
      g.fillText(reaction, x + 4, y + h / 2);
      cruxReactionCounts.asText++;
    }
  } finally {
    g.restore();
  }
}

/** Draws a small reaction once, so that the first cell does not pay for the font and the first compilation (decision
 * 3's note: about 30 ms); Chem calls it when Crux is ready in Crux mode. */
export function warmCruxReactionDrawing(): void {
  cruxReactionPixels(new OffscreenCanvas(40, 20), 'C>>C', 40, 20, 1, {drawOptions: chemDrawOptions()});
}
