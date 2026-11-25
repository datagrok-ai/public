/**
 * Reaction cells that hold CDXML (a cell keeps what its column holds; Chem's demo table of such a column is
 * files/reactions-cdxml.csv, its pages' sources and licences in NOTICE.reactions-cdxml beside this file). Crux mode
 * only (RDKit mode never detects such a column, so it shows its text as before).
 * - A column whose values are CDXML pages is a reaction column of units `cdxml` ({@link detectCdxmlColumn}, called by
 *   detectors.js through Chem:detectCdxmlReactions, so that it runs where the engine setting is known).
 * - Its cells are drawn by Crux Sketch: one offscreen reaction sketcher per page reads each value (its CDXML module
 *   loads with the first) and gives the page's picture, its `svg` (the canvas's own renderer: arrows with their styles,
 *   texts, brackets, plus signs, cycles; framed to every object drawn), kept per value ({@link PICTURES_KEPT}); a page
 *   with objects away from its reaction (a key box) is pictured cut to its reaction ({@link reactionPage}). A cell
 *   whose picture is still being made is painted again when it is (its grid invalidated). One value is read at a time.
 * - What the page is as a reaction ({@link cdxmlReactionSmiles}: the sketcher's reaction SMILES), for the Reaction
 *   panel and the checks, which read reaction text.
 * The cell editor writes such a cell back as CDXML, the sketcher's `cdxml` getter (after an edit, `value` is reaction
 * SMILES, which drops the texts and the arrows' styles).
 */
import * as DG from 'datagrok-api/dg';
import type {SketcherHandle} from '../../vendor/crux-sketch/crux-sketch';
import {loadCruxSketch} from '../crux/crux-sketcher';
import {isCruxRequested} from '../engine/chem-engine';
import {REACTION_CELL_TEXT_COLOUR} from './crux-reaction-drawing';

/** A reaction column's units for CDXML pages (detectors.js's other units: smiles, smarts, rxn, route). */
export const CDXML_UNITS = 'cdxml';
/** How many values' pictures are kept (the oldest dropped first). */
export const PICTURES_KEPT = 200;
/** The offscreen sketcher's box (CSS px): `svg` is drawn at natural size whatever the box. */
const READER_BOX = {width: 800, height: 600};
/** Room kept free round a picture in its cell (CSS px). */
const CELL_PADDING = 4;
/** What a cell says while its picture is made, and when its page could not be read. */
export const CDXML_DRAWING = 'Drawing…';
export const CDXML_UNREAD = 'This CDXML page could not be read';
/** What a cell says of a page that fits it only with bonds shorter than {@link CDXML_MIN_BOND} CSS px: how to see it
 * (the mechanisms page's cell once said only "Too large to show at this size"; the cell editor opens it). */
export const CDXML_TOO_SMALL = 'Too large for this cell: double-click to open';
/** The cells' messages: their font size (CSS px), and the colour of a page not read (the others the reaction cells'
 * muted grey). */
const MESSAGE_FONT_PX = 12;
const MESSAGE_LINE = 1.3;
const UNREAD_COLOUR = '#a00';
export const CDXML_MIN_BOND = 4;

/** Whether `s` is a CDXML page: its first element CDXML, after an XML prolog and a doctype. */
export function isCdxmlText(s: string): boolean {
  return typeof s === 'string' && /^\s*(<\?xml[^>]*>\s*)?(<!DOCTYPE[^>]*>\s*)?<CDXML[\s>]/.test(s.slice(0, 4096));
}

/** A column of CDXML pages made a reaction column of units `cdxml`, in Crux mode; nothing in RDKit mode. */
export function detectCdxmlColumn(col: DG.Column): void {
  if (!isCruxRequested())
    return;
  col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
  col.meta.units = CDXML_UNITS;
  col.meta.cellRenderer = 'ChemicalReaction';
}

/** Whether a cell is drawn here: its column's units are `cdxml`, or its value is a CDXML page. */
export function isCdxmlCell(cell: DG.Cell): boolean {
  return cell.column?.meta.units === CDXML_UNITS || isCdxmlText(cell.value);
}

/** A value read: its picture (an image of its `svg`) and its size, its reaction SMILES; or why it was not read. */
interface Read {
  image: HTMLImageElement | null;
  width: number;
  height: number;
  /** A bond's length in the picture at its natural size (CSS px; the median of its bonds), 0 when it has none. */
  bond: number;
  smiles: string;
  error: string | null;
}

const reads = new Map<string, Read | Promise<Read | null>>();
/** A cell that asked for a value's read while it waits: whether it is on screen still (its grid shown, its row among
 * the rows the grid shows). */
export interface CdxmlAsker {
  onScreen: () => boolean;
}
/** The cells waiting for each value's read; a value someone needs whatever the screen shows (the Reaction panel's
 * reaction SMILES, the tests) is in {@link needed}. A read whose every cell has left the screen when its turn comes is
 * skipped, and dropped: a cell that comes back on screen asks again (a fast scroll through a 1,500-page table queued a
 * read per row passed, ~40 s of CPU after the scroll stopped). */
const askers = new Map<string, CdxmlAsker[]>();
const needed = new Set<string>();
/** The reads made and the reads skipped since the page loaded (the tests read them). */
export const cdxmlReadCounts = {read: 0, skipped: 0};
let reader: Promise<SketcherHandle> | null = null;
let readerHost: HTMLElement | null = null;
let queue: Promise<unknown> = Promise.resolve();

/** The page's one offscreen reaction sketcher, made on first use (again after it failed to start). */
function theReader(): Promise<SketcherHandle> {
  reader ??= (async () => {
    const host = document.createElement('div');
    Object.assign(host.style, {position: 'absolute', left: '-10000px', top: '0', width: `${READER_BOX.width}px`,
      height: `${READER_BOX.height}px`});
    document.body.append(host);
    readerHost = host;
    try {
      const sketch = await (await loadCruxSketch()).createSketcher(host,
        {document: 'reaction', toolset: 'publication', toolbar: {}});
      await sketch.loadFormat('cdxml');
      return sketch;
    } catch (e) {
      host.remove();
      throw e;
    }
  })().catch((e) => {
    reader = null;
    throw e;
  });
  return reader;
}

/** `svg` as an image, once it is decoded. */
function imageOf(svg: string): Promise<HTMLImageElement> {
  return new Promise((resolve, reject) => {
    const image = new Image();
    image.onload = () => resolve(image);
    image.onerror = () => reject(new Error('the picture could not be decoded'));
    image.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
  });
}

/** The width and height an SVG says it is (CSS px), or 0 x 0. */
function svgSize(svg: string): [number, number] {
  const head = svg.slice(0, svg.indexOf('>') + 1);
  const n = (name: string) => Number(new RegExp(`\\b${name}="([\\d.]+)`).exec(head)?.[1] ?? 0);
  return [n('width'), n('height')];
}

/** A bond's length in Crux Sketch's `svg` (its own units, CSS px at natural size): the median, over its bond groups
 * (ids `<prefix>bond-<i>`), of each one's length from atom centre to atom centre (its `data-atoms`' `data-x` and
 * `data-y`), else of its first line; 0 when it draws no bond. A line stops short of an atom's label: the pages keep
 * ChemDraw's labels (Ph, SmI2) contracted, and the lines' median put the demo's mechanisms page under
 * {@link CDXML_MIN_BOND} in a 400 x 200 cell (3.94 px) where its bonds are 4.20 px long. */
export function svgBondLength(svg: string): number {
  const doc = new DOMParser().parseFromString(svg, 'image/svg+xml');
  const lengths: number[] = [];
  const centres = new Map<string, number[]>(); // atom group id -> its centre
  for (const atom of Array.from(doc.querySelectorAll('[data-x][data-y]'))) {
    const xy = [Number(atom.getAttribute('data-x')), Number(atom.getAttribute('data-y'))];
    if (atom.id && xy.every(Number.isFinite))
      centres.set(atom.id, xy);
  }
  const centre = (prefix: string, atom: string) => centres.get(`${prefix}atom-${atom}`) ?? null;
  for (const group of Array.from(doc.querySelectorAll('[id*="bond-"]'))) {
    if (!/bond-\d+$/.test(group.id))
      continue;
    const prefix = group.id.replace(/bond-\d+$/, '');
    const [a, b] = (group.getAttribute('data-atoms') ?? '').trim().split(/\s+/).map((i) => centre(prefix, i));
    if (a && b) {
      lengths.push(Math.hypot(b[0] - a[0], b[1] - a[1]));
      continue;
    }
    const line = group.matches('line, path') ? group : group.querySelector('line, path');
    let ends: number[] = [];
    if (line?.tagName === 'line')
      ends = ['x1', 'y1', 'x2', 'y2'].map((a) => Number(line.getAttribute(a)));
    else if (line) // a path: its first two points
      ends = Array.from((line.getAttribute('d') ?? '').matchAll(/-?\d*\.?\d+(?:e-?\d+)?/gi), (m) => Number(m[0]));
    if (ends.length >= 4 && ends.slice(0, 4).every(Number.isFinite))
      lengths.push(Math.hypot(ends[2] - ends[0], ends[3] - ends[1]));
  }
  lengths.sort((a, b) => a - b);
  return lengths.length ? lengths[lengths.length >> 1] : 0;
}

/** Whether `cdxml`'s read is wanted when its turn comes: someone needs it, or a cell that asked for it is on screen. */
function wanted(cdxml: string): boolean {
  if (needed.has(cdxml))
    return true;
  return (askers.get(cdxml) ?? []).some((a) => {
    try {
      return a.onScreen();
    } catch {
      return false;
    }
  });
}

/** A grid cell's {@link CdxmlAsker}: its grid shown, its row among the rows the grid shows (a row's margin each
 * side). */
function cellAsker(gridCell: DG.GridCell): CdxmlAsker {
  const grid = gridCell.grid;
  const row = gridCell.gridRow;
  return {onScreen: () => {
    if (!grid || !grid.root.isConnected)
      return false;
    const shown = grid.vertScroll;
    return row >= Math.floor(shown.min) - 1 && row <= Math.ceil(shown.max) + 1;
  }};
}

/** `cdxml` read by the page's sketcher, one value at a time; null, read not, when its turn comes and it is no longer
 * {@link wanted}. */
function readValue(cdxml: string): Promise<Read | null> {
  const job = queue.then(async (): Promise<Read | null> => {
    if (!wanted(cdxml)) {
      cdxmlReadCounts.skipped++;
      return null;
    }
    cdxmlReadCounts.read++;
    try {
      const sketch = await theReader();
      sketch.setValue(cdxml);
      if (sketch.lastError !== null)
        return {image: null, width: 0, height: 0, bond: 0, smiles: '', error: sketch.lastError.message};
      let svg = sketch.svg;
      let smiles = '';
      try {
        smiles = sketch.smiles;
      } catch {
        // a page with no reaction text (a molecule, a query): pictured all the same
      }
      // the picture cut to the page's reaction, its key box and the like left out; the reaction SMILES is
      // the whole page's, the same
      const cut = reactionPage(cdxml);
      if (cut !== null) {
        sketch.setValue(cut);
        if (sketch.lastError === null && !sketch.isEmpty)
          svg = sketch.svg;
        else
          console.warn('Chem | a CDXML page cut to its reaction did not read: the page is pictured', sketch.lastError);
      }
      const image = svg ? await imageOf(svg) : null;
      const [width, height] = svgSize(svg);
      return {image, width: width || image?.naturalWidth || 0, height: height || image?.naturalHeight || 0,
        bond: svgBondLength(svg), smiles, error: null};
    } catch (e) {
      console.error('Chem | a CDXML cell could not be read', e);
      return {image: null, width: 0, height: 0, bond: 0, smiles: '',
        error: e instanceof Error ? e.message : String(e)};
    }
  });
  queue = job;
  return job;
}

/** `cdxml`'s read: at once when made, else null, the read started (once) and `done` called when it is made. `asker`:
 * the cell that asks, whose read is skipped if it has left the screen when its turn comes; none, the read is needed
 * whatever the screen shows. */
function readNow(cdxml: string, done: () => void, asker?: CdxmlAsker): Read | null {
  const kept = reads.get(cdxml);
  if (kept && !(kept instanceof Promise)) {
    reads.delete(cdxml); // the newest last: the oldest is dropped first
    reads.set(cdxml, kept);
    return kept;
  }
  if (asker)
    askers.set(cdxml, [...(askers.get(cdxml) ?? []).slice(-MAX_ASKERS + 1), asker]);
  else
    needed.add(cdxml);
  if (!kept) {
    const p = readValue(cdxml).then((r) => {
      askers.delete(cdxml);
      needed.delete(cdxml);
      if (r === null) { // skipped: asked again when a cell shows it
        reads.delete(cdxml);
        return null;
      }
      reads.set(cdxml, r);
      while (reads.size > PICTURES_KEPT)
        reads.delete(reads.keys().next().value!);
      return r;
    });
    reads.set(cdxml, p);
    p.then((r) => r !== null && done());
  } else
    kept.then((r) => r !== null && done());
  return null;
}
/** The cells kept for a value waiting to be read (the latest: a cell painted again asks again). */
const MAX_ASKERS = 8;

/** Tests: `cdxml`'s read queued as a cell's that is on screen while `onScreen()` says so; resolves once its turn has
 * come: true when it was read, false when it was skipped. */
export async function queueCdxmlReadForTests(cdxml: string, onScreen: () => boolean): Promise<boolean> {
  readNow(cdxml, () => {}, {onScreen});
  const r = await reads.get(cdxml);
  return r !== null && r !== undefined;
}

/** The reaction SMILES of CDXML page `cdxml` (the steps it draws), or null when it was not read. */
export async function cdxmlReactionSmiles(cdxml: string): Promise<string | null> {
  readNow(cdxml, () => {});
  const read = await reads.get(cdxml);
  return read && read.error === null && read.smiles ? read.smiles : null;
}

/** A step's lists of the page's objects, by their ids (ChemDraw's `<step>`): its reactants, its products and what is
 * over and under its arrow. */
const STEP_LISTS = ['ReactionStepReactants', 'ReactionStepProducts', 'ReactionStepObjectsAboveArrow',
  'ReactionStepObjectsBelowArrow'];

/**
 * How many molecules of CDXML `page` are on no reaction step: its structures (a `<fragment>` of its own, not one inside
 * another or inside an atom node, as an abbreviation's) whose id no step lists, by itself or by a group that holds
 * it. Crux Sketch's reaction text holds the steps' molecules only, so a reaction SMILES or RXN made from the page
 * leaves these out, and says so ({@link noStepWords}). 0 for a page with no step (no reaction text
 * is made from it) or one that does not parse.
 */
export function moleculesOnNoStep(page: string): number {
  let doc: Document;
  try {
    doc = new DOMParser().parseFromString(page, 'text/xml');
  } catch {
    return 0;
  }
  const steps = Array.from(doc.getElementsByTagName('step'));
  if (doc.getElementsByTagName('parsererror').length || !steps.length)
    return 0;
  const listed = new Set<string>();
  for (const step of steps) {
    for (const name of STEP_LISTS)
      (step.getAttribute(name) ?? '').split(/\s+/).filter((id) => id).forEach((id) => listed.add(id));
  }
  const onStep = (fragment: Element): boolean => {
    for (let e: Element | null = fragment; e; e = e.parentElement) {
      if (listed.has(e.getAttribute('id') ?? ''))
        return true;
    }
    return false;
  };
  return Array.from(doc.getElementsByTagName('fragment'))
    .filter((f) => !f.parentElement?.closest('fragment, n') && !onStep(f)).length;
}

/** A step's lists of the page's objects besides {@link STEP_LISTS}: its arrows and its plus signs. */
const STEP_DRAWN_LISTS = ['ReactionStepArrows', 'ReactionStepPlusses'];
/** A box on a CDXML page: left, top, right, bottom (pt, y down). */
type PageBox = [number, number, number, number];

/** The smallest box holding `a` and `b` (either may be null). */
function unionBox(a: PageBox | null, b: PageBox | null): PageBox | null {
  return !a ? b : !b ? a :
    [Math.min(a[0], b[0]), Math.min(a[1], b[1]), Math.max(a[2], b[2]), Math.max(a[3], b[3])];
}

/** `e`'s box on its CDXML page: its BoundingBox, else the box of its points (p, Head3D, Tail3D, Center3D) and of what
 * it holds; null when it says none. */
function pageBoxOf(e: Element): PageBox | null {
  const numbers = (name: string) => (e.getAttribute(name) ?? '').trim().split(/[\s,]+/).filter((v) => v).map(Number);
  const b = numbers('BoundingBox');
  if (b.length === 4 && b.every(Number.isFinite))
    return [Math.min(b[0], b[2]), Math.min(b[1], b[3]), Math.max(b[0], b[2]), Math.max(b[1], b[3])];
  let box: PageBox | null = null;
  for (const name of ['p', 'Head3D', 'Tail3D', 'Center3D']) {
    const p = numbers(name);
    if (p.length >= 2 && Number.isFinite(p[0]) && Number.isFinite(p[1]))
      box = unionBox(box, [p[0], p[1], p[0], p[1]]);
  }
  for (const child of Array.from(e.children))
    box = unionBox(box, pageBoxOf(child));
  return box;
}

/**
 * CDXML `page` cut to its reaction, for its cell's picture (the Suzuki and N-sulfonylation pages of the demo table
 * drew their key box too, the reaction a few pixels high at the opening row height and the key's labels on one
 * another), or null when nothing would be cut (or the page has no step, or does not parse). What is kept of the page's
 * objects: its schemes; those its steps list (molecules, arrows, plus signs, the texts over and under the arrows),
 * or that hold one (a group); an object beside them, its centre within two bonds of their box (a plus sign, an agent
 * not listed); a text within two bonds of what is kept, one after another (a title, a name, its amount, its role); an
 * object with no box. The cell's value stays the whole page: the editor opens it and Export writes it.
 */
export function reactionPage(page: string): string | null {
  let doc: Document;
  try {
    doc = new DOMParser().parseFromString(page, 'text/xml');
  } catch {
    return null;
  }
  const sheet = doc.getElementsByTagName('page')[0];
  const steps = Array.from(doc.getElementsByTagName('step'));
  if (doc.getElementsByTagName('parsererror').length || !sheet || !steps.length)
    return null;
  const listed = new Set<string>();
  for (const step of steps) {
    for (const name of [...STEP_LISTS, ...STEP_DRAWN_LISTS])
      (step.getAttribute(name) ?? '').split(/\s+/).filter((id) => id).forEach((id) => listed.add(id));
  }
  const onStep = (e: Element) => listed.has(e.getAttribute('id') ?? '') ||
    listed.has(e.getAttribute('SupersededBy') ?? '') || // an arrow's older graphic
    Array.from(e.querySelectorAll('[id]')).some((d) => listed.has(d.getAttribute('id') ?? ''));
  const objects = Array.from(sheet.children);
  const kept = new Set(objects.filter((e) => e.tagName === 'scheme' || onStep(e)));
  let region: PageBox | null = null;
  for (const e of kept)
    region = unionBox(region, pageBoxOf(e));
  if (!region)
    return null;
  const room = 2 * (Number(doc.documentElement.getAttribute('BondLength')) || 14.4);
  const stepsBox = region;
  const centreBeside = (b: PageBox) => {
    const [x, y] = [(b[0] + b[2]) / 2, (b[1] + b[3]) / 2];
    return x >= stepsBox[0] - room && x <= stepsBox[2] + room && y >= stepsBox[1] - room && y <= stepsBox[3] + room;
  };
  for (const e of objects) {
    const b = kept.has(e) ? null : pageBoxOf(e);
    if (!kept.has(e) && (!b || (e.tagName !== 't' && centreBeside(b))))
      kept.add(e);
  }
  const gap = (a: PageBox, b: PageBox) => Math.max(0, a[0] - b[2], b[0] - a[2], a[1] - b[3], b[1] - a[3]);
  for (let grown = true; grown;) {
    grown = false;
    for (const e of objects) {
      const b = kept.has(e) || e.tagName !== 't' ? null : pageBoxOf(e);
      if (b && gap(b, region!) <= room) {
        kept.add(e);
        region = unionBox(region, b);
        grown = true;
      }
    }
  }
  const dropped = objects.filter((e) => !kept.has(e));
  if (!dropped.length)
    return null;
  dropped.forEach((e) => e.remove());
  return new XMLSerializer().serializeToString(doc);
}

/** What a reaction text made from a page with `n` molecules on no step says of them ({@link moleculesOnNoStep}). */
export function noStepWords(n: number): string {
  return `${n} ${n === 1 ? 'molecule' : 'molecules'} on no step ${n === 1 ? 'was' : 'were'} left out`;
}

/** Reaction `text` (any notation Crux Sketch reads) as a CDXML page, written by Crux Sketch's reaction writer on the
 * page's reader (the Export panel's CDXML), one at a time with the cells' reads; null while the vendored
 * sketcher writes no reaction as CDXML (its getter then throws: "a reaction is not written as CDXML yet"), or for a
 * text it does not read. */
export function reactionCdxml(text: string): Promise<string | null> {
  const job = queue.then(async (): Promise<string | null> => {
    try {
      const sketch = await theReader();
      sketch.setValue(text);
      return sketch.lastError === null ? sketch.cdxml || null : null;
    } catch (e) {
      console.warn('Chem | the reaction was not written as CDXML', e);
      return null;
    }
  });
  queue = job;
  return job;
}

/** What a drawing's picture takes: crux-sketch's picture API (`svgWith` and `png`), whose defaults are the
 * sketcher's own (its map numbers as shown, a 20 pt bond, its atom colours). */
export interface DrawnPictureOptions {
  /** Map numbers drawn beside their atoms, or not. */
  mapNumbers?: boolean;
  /** The mean bond, in points (1 pt = 4/3 CSS px): ACS 1996's is 14.4. */
  bondLength?: number;
  /** Atoms in element colours, or all black. */
  colour?: boolean;
}

/** A drawing as it is drawn (the Export panel's "As drawn"): its picture as an SVG, in CSS px at its natural size, and
 * as a PNG at `scale` times that (crux-sketch's `png`: a white background, its pHYs 96 x scale dpi). */
export interface DrawnPicture {
  /** Whether the pictures follow the options; false: the sketcher's own size, colours and map numbers (a crux-sketch
   * before its picture API). */
  readonly takesOptions: boolean;
  svg(options: DrawnPictureOptions): Promise<string | null>;
  png(options: DrawnPictureOptions & {scale: number}): Promise<Blob | null>;
}

/** The picture getters of a Crux sketcher, the element's or its handle's; `svgWith` from crux-sketch's picture API. */
interface PictureGetters {
  readonly svg: string;
  png(options?: {scale?: number}): Promise<Blob | null>;
  svgWith?: (options?: DrawnPictureOptions) => string;
}

/** `sketch`'s drawing as a {@link DrawnPicture}: crux-sketch's `svgWith` and `png(options)` where it has them, else
 * its `svg` and `png({scale})`, the options left out. */
export function sketcherPicture(sketch: PictureGetters): DrawnPicture {
  const takesOptions = typeof sketch.svgWith === 'function';
  return {
    takesOptions,
    svg: async (options) => (takesOptions ? sketch.svgWith!(options) : sketch.svg) || null,
    png: (options) => sketch.png(takesOptions ? options : {scale: options.scale}),
  };
}

/** Whether reaction `value` has coordinates of its own to draw as drawn: a CDXML page, or an RXN block (V2000 or
 * V3000) with an atom off the origin. A reaction SMILES or SMARTS has none: the sketcher would lay it out. */
export function drawnWithCoordinates(value: string): boolean {
  if (isCdxmlText(value))
    return true;
  if (!/^\s*\$RXN/.test(value))
    return false;
  const coordinate = /^-?\d+\.\d+$/;
  return value.split('\n').some((line) => {
    const f = line.trim().split(/\s+/);
    const [x, y] = f[0] === 'M' && f[1] === 'V30' ? [f[4], f[5]] : [f[0], f[1]];
    return f.length > 3 && coordinate.test(x ?? '') && coordinate.test(y ?? '') &&
      (Number(x) !== 0 || Number(y) !== 0);
  });
}

/** A cell's value as drawn ({@link drawnWithCoordinates}): its picture made on the page's reader, one at a time with
 * the cells' reads (the page as a whole for a CDXML page, as ChemDraw draws it); null for a value with no
 * coordinates. */
export function cellPicture(value: string): DrawnPicture | null {
  if (!drawnWithCoordinates(value))
    return null;
  let takesOptions = true;
  const onReader = <T>(make: (sketch: SketcherHandle) => Promise<T | null>): Promise<T | null> => {
    const job = queue.then(async (): Promise<T | null> => {
      try {
        const sketch = await theReader();
        sketch.setValue(value);
        return sketch.lastError === null ? await make(sketch) : null;
      } catch (e) {
        console.warn('Chem | the cell\'s drawing was not pictured', e);
        return null;
      }
    });
    queue = job;
    return job;
  };
  return {
    get takesOptions() {
      return takesOptions;
    },
    svg: (options) => onReader(async (sketch) => {
      const picture = sketcherPicture(sketch as unknown as PictureGetters);
      takesOptions = picture.takesOptions;
      return picture.svg(options);
    }),
    png: (options) => onReader((sketch) => sketcherPicture(sketch as unknown as PictureGetters).png(options)),
  };
}

/** Paints CDXML cell `gridCell` into `g` at (x, y, w, h), in the canvas's own pixels (`ratio` of them a CSS px, as
 * the grid's reaction renderer paints): its picture fitted to the cell, never larger than its natural size; a line
 * while it is made (the grid painted again when it is) or when its page was not read. */
export function drawCdxmlCell(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
  gridCell: DG.GridCell, ratio = 1): void {
  const value: string = gridCell.cell.value ?? '';
  if (!value.trim()) // an empty cell: nothing to read or draw
    return;
  const read = readNow(value, () => gridCell.grid?.invalidate(), cellAsker(gridCell));
  const pad = CELL_PADDING * ratio;
  g.save();
  try {
    // in the canvas's own pixels whatever transform the grid's context carries (at a device
    // pixel ratio over 1 the grid's context is scaled, and the cells' pictures, scaled twice, were cut off and slid
    // against a sideways scroll), and kept to the cell, as the other reaction cells are drawn
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    const fit = read?.image ? cdxmlFit(read, w, h, ratio) : null;
    if (read?.image && fit !== null && !fit.tooSmall) {
      const [dw, dh] = [read.width * fit.scale, read.height * fit.scale];
      g.drawImage(read.image, x + (w - dw) / 2, y + (h - dh) / 2, dw, dh);
      return;
    }
    if (fit !== null && fit.scale <= 0) // a cell with no room
      return;
    g.font = `${MESSAGE_FONT_PX * ratio}px Roboto, Roboto Local, sans-serif`;
    g.fillStyle = read?.error ? UNREAD_COLOUR : REACTION_CELL_TEXT_COLOUR;
    g.textBaseline = 'middle';
    g.textAlign = 'center';
    const message = read?.error ? CDXML_UNREAD : fit?.tooSmall ? CDXML_TOO_SMALL : CDXML_DRAWING;
    // centred, never squeezed to the cell's width
    const line = MESSAGE_FONT_PX * MESSAGE_LINE * ratio;
    const lines = legibleLines((t) => g.measureText(t).width, message, Math.max(0, w - 2 * pad), h - 2 * pad,
      MESSAGE_FONT_PX * ratio, line);
    lines.forEach((t, i) => g.fillText(t, x + w / 2, y + h / 2 + (i - (lines.length - 1) / 2) * line));
  } finally {
    g.restore();
  }
}

/** The lines of a cell's `message` that are legible in `room` x `height` px (`width(t)`: a text's width; `font`: the
 * font's size, `line`: a line's height), never squeezed (in a box one line high, such as Datagrok's view header, the
 * message was drawn squashed to a smear): the message on one line; else broken after its colon, when two lines fit;
 * else its first words that fit, with an ellipsis; else an ellipsis alone; else nothing. */
export function legibleLines(width: (t: string) => number, message: string, room: number, height: number,
  font: number, line: number): string[] {
  const fits = (t: string) => width(t) <= room;
  if (height < font)
    return [];
  if (fits(message))
    return [message];
  const colon = message.indexOf(': ');
  const two = [message.slice(0, colon + 1), message.slice(colon + 2)];
  if (colon > 0 && height >= 2 * line && two.every(fits))
    return two;
  const words = message.split(' ');
  for (let k = words.length - 1; k > 0; k--) {
    const shortened = `${words.slice(0, k).join(' ').replace(/[:,]$/, '')}…`;
    if (fits(shortened))
      return [shortened];
  }
  return fits('…') ? ['…'] : [];
}

/** A page's picture in a cell `w` x `h` canvas px (`ratio` of them a CSS px): the scale that fits it whole, the cell's
 * padding kept, never larger than its natural size; and whether its bonds would then be shorter than
 * {@link CDXML_MIN_BOND} CSS px, too small to read (the cell then says so instead, as Chem's cells say a molecule too
 * large to draw). A page whose bonds are not known is never too small. */
export function cdxmlFit(read: {width: number, height: number, bond: number}, w: number, h: number, ratio = 1):
  {scale: number, tooSmall: boolean} {
  const pad = CELL_PADDING * ratio;
  if (read.width <= 0 || read.height <= 0)
    return {scale: 0, tooSmall: false};
  const scale = Math.min(ratio, (w - 2 * pad) / read.width, (h - 2 * pad) / read.height);
  return {scale, tooSmall: scale > 0 && read.bond > 0 && read.bond * scale / ratio < CDXML_MIN_BOND};
}

/** Tests: the picture CDXML cell value `value` is drawn from, once read (its size and its bond, CSS px at natural
 * size), or null when it was not read. */
export async function cdxmlPictureForTests(value: string): Promise<{width: number, height: number, bond: number} |
  null> {
  readNow(value, () => {});
  const read = await reads.get(value);
  return read && read.error === null && read.image ? {width: read.width, height: read.height, bond: read.bond} : null;
}

/** Tests: every read dropped and the reader let go, so that the next read starts afresh. */
export async function resetCdxmlReadsForTests(): Promise<void> {
  reads.clear();
  askers.clear();
  needed.clear();
  const r = reader;
  reader = null;
  try {
    (await r)?.destroy();
  } catch {
    // a reader that did not start: nothing to let go
  }
  readerHost?.remove();
  readerHost = null;
}
