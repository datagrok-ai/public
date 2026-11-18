/**
 * Category `Leaks: Crux`: an action repeated N times leaves nothing behind. Each test does the action once (the
 * module, the caches and the platform warm), takes a baseline, repeats the action N times, and checks that
 * - the live Crux engine objects of Chem's main-thread module (its wasm-bindgen classes: ChemReactionScheme, JSMol,
 *   Reaction, ...) are back at the baseline. An object counts from the moment the module registers it for finalization
 *   until its `free()`, whatever GC does: an object left to the garbage collector holds its wasm memory until a GC that
 *   may never come, and wasm memory never shrinks;
 * - the event listeners the action added to what outlives it (window, document, elements still in the page) are
 *   removed again.
 * The DOM of a closed view is not checked: the platform keeps every closed view's root (a plain `DG.View` with a
 * canvas, opened and closed 10 times: none of the 10 roots collected after GC, +114 nodes and +71 listeners a view),
 * so it would measure the platform, not Chem.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {ChemicalEngine, isCruxRequested} from '../engine/chem-engine';
import {ChemWasm, armCruxTrapForTests, cruxHeapBytes, loadCrux} from '../engine/crux-module';
import {chemDrawOptions, cruxReactionCounts, cruxReactionPixels, cruxReactionSteps, cruxReadsReaction,
  drawReactionByCrux, resetCruxSchemeCache} from '../rendering/crux-reaction-drawing';
import {inEngine} from './engine-utils';
import {MAPPING_CLASS, reactionPanel} from '../widgets/reaction-panel';
import {DOCUMENT_KEY, forgetDocument} from '../utils/reactions/reaction-document';
import {PackageFunctions, getMoleculeRenderer} from '../package';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';

/** How many times each action is repeated after its warm-up. */
const ROUNDS = 10;
/** The scheme cache's size in the cache test (resetCruxSchemeCache): small, so that the test draws well past it. Chem's
 * own size stays 1,000 (crux-reaction-drawing.ts SCHEME_CACHE_SIZE): each drawn scheme keeps its layout, ~110 KB of the
 * Crux heap, so 1,000 hold ~107 MB (measured in a browser), but 200 made a scroll back up 68% slower (the scroll
 * benchmark below); the engine's slim-down of a drawn scheme is the fix proposed instead. */
const TEST_CACHE_SIZE = 50;
/** Distinct reactions the cache test draws: four times its cache. */
const DISTINCT_REACTIONS = 4 * TEST_CACHE_SIZE;
/** The longest a reaction cell may hold the main thread for one value (ms): a 1,000-residue reaction drew in 7 s
 * (catalogue RXN-08, PERF-03), the tab frozen meanwhile. */
const HUGE_REACTION_BUDGET_MS = 1000;
/** The longest the Reaction panel may take to build for a reaction too large to draw (ms): it shows words, reads
 * nothing (the catalogue's PERF-03 gives the context panel 1 s; the panel is one of its sections). */
const HUGE_PANEL_BUDGET_MS = 100;
/** The scroll benchmark: rows, rows a page (a 900 px grid of 120 px rows), the cell's size (CSS px), the scheme cache
 * sizes compared (1,000, its first size, and Chem's), and the reps of each (interleaved). */
const SCROLL_ROWS = 2000;
const SCROLL_PAGE = 7;
const SCROLL_CELL_W = 600;
const SCROLL_CELL_H = 120;
const SCROLL_CACHE_SIZES = [1000, 200];
const SCROLL_REPS = 4;
/** Listeners left on long-lived targets after ROUNDS rounds, over the baseline, still counted as no leak: a leak of one
 * a round leaves ROUNDS; one-off lazy listeners (a first tooltip, a first menu) leave a few. */
const LISTENER_SLACK = 3;

/** The picture byte budget test: cells drawn, and the budget (Chem's CRUX_PICTURE_BYTES, 64 MB). */
const PICTURES_DRAWN = 100;
const PICTURE_BUDGET_BYTES = 64 * 1048576;

/** The bytes of the pictures a renderer's picture cache keeps: its `bytes` (a PictureCache), else the ImageData a
 * DG.LruCache holds (its values, `V`). */
function keptPictureBytes(cache: unknown): number {
  const c = cache as {bytes?: number, V?: (ImageData | undefined)[]};
  return typeof c.bytes === 'number' ? c.bytes : (c.V ?? []).reduce((n, v) => n + (v?.data?.byteLength ?? 0), 0);
}

type AnyClass = abstract new (...args: any[]) => any;

/**
 * The live objects of the module's classes: `FinalizationRegistry.prototype.register` and `unregister` watched while a
 * test runs (wasm-bindgen registers each object it makes or wraps, with the object as its unregister token, and
 * unregisters it in `free()`). Objects made before the watch are not counted.
 */
class LiveObjects {
  private readonly live = new Map<string, Set<object>>();
  private readonly register = FinalizationRegistry.prototype.register;
  private readonly unregister = FinalizationRegistry.prototype.unregister;

  constructor(private readonly classes: [string, AnyClass][]) {
    const self = this; // eslint-disable-line @typescript-eslint/no-this-alias
    const proto = FinalizationRegistry.prototype as any;
    proto.register = function(target: object, held: unknown, token?: object) {
      const entry = self.classes.find(([, c]) => target instanceof c);
      if (entry) {
        if (!self.live.has(entry[0]))
          self.live.set(entry[0], new Set());
        self.live.get(entry[0])!.add(target);
      }
      return self.register.call(this, target, held, token);
    };
    proto.unregister = function(token: object) {
      for (const set of self.live.values())
        set.delete(token);
      return self.unregister.call(this, token);
    };
  }

  /** The live count by class, the classes with none left out. */
  counts(): {[cls: string]: number} {
    const out: {[cls: string]: number} = {};
    for (const [cls, set] of this.live) {
      if (set.size)
        out[cls] = set.size;
    }
    return out;
  }

  stop(): void {
    (FinalizationRegistry.prototype as any).register = this.register;
    (FinalizationRegistry.prototype as any).unregister = this.unregister;
  }
}

/** The module's classes that make engine objects (those with `free`), by their names in the module. */
function engineClasses(crux: ChemWasm): [string, AnyClass][] {
  return Object.entries(crux).filter(([, v]) => typeof v === 'function' && typeof v.prototype?.free === 'function')
    .map(([name, v]) => [name, v as AnyClass]);
}

interface ListenerEntry {target: EventTarget; type: string; listener: unknown; capture: boolean; chem: boolean;
  origin: string}

/** The package a stack frame runs in: Chem for its bundle's files (at Chem's web root), another package for its files
 * (".../files/<Package>/..."), null for the platform's own code. */
function packageOfFrame(frame: string): string | null {
  const root = _package?.webRoot;
  if (root ? frame.includes(root) : /\/files\/Chem\//.test(frame))
    return 'Chem';
  return /\/files\/([^/]+)\//.exec(frame)?.[1] ?? null;
}

/** Who added a listener, from the stack taken in the watch's own addEventListener: the innermost package on it (the
 * nearest frame to the call), else the platform. The stack's first frame is the watch itself, in Chem's test bundle,
 * so it is left out: counted, it made every listener Chem's. And the innermost package decides, not any package on the
 * stack: a listener another package adds while handling an event a Chem call fires has Chem's frames further down
 * (Grokky's keyup listeners on the "Search everywhere" input, added as a closed app view brings the Home view back,
 * were counted as Chem's on cold runs: 33 -> 39). */
function originOf(stack: string): {chem: boolean, origin: string} {
  for (const frame of stack.split('\n').slice(2)) { // "Error", then the watch's own frame
    const pkg = packageOfFrame(frame);
    if (pkg !== null)
      return {chem: pkg === 'Chem', origin: pkg};
  }
  return {chem: false, origin: 'platform'};
}

/** Frames kept in the watch's stacks: deep enough to reach the package below the platform's frames of a listener the
 * platform adds for a package's call (V8 keeps 10 by default). */
const WATCH_STACK_FRAMES = 200;
/** V8's Error, whose `stackTraceLimit` the watch raises while it takes a stack. */
const v8Error = Error as unknown as {stackTraceLimit: number};

/** The event listeners added and not removed while a test runs (`EventTarget.prototype` watched), each marked by who
 * added it ({@link originOf}). */
class Listeners {
  private readonly entries: ListenerEntry[] = [];
  private readonly add = EventTarget.prototype.addEventListener;
  private readonly remove = EventTarget.prototype.removeEventListener;

  constructor() {
    const self = this; // eslint-disable-line @typescript-eslint/no-this-alias
    const capture = (options: unknown) => typeof options === 'boolean' ? options : !!(options as any)?.capture;
    EventTarget.prototype.addEventListener = function(type: string, listener: any, options?: any) {
      if (listener && !(typeof options === 'object' && options?.once)) {
        // the stack taken here, so that its first frame is this function (originOf leaves it out)
        const limit = v8Error.stackTraceLimit;
        v8Error.stackTraceLimit = WATCH_STACK_FRAMES;
        const stack = new Error().stack ?? '';
        v8Error.stackTraceLimit = limit;
        self.entries.push({target: this, type, listener, capture: capture(options), ...originOf(stack)});
      }
      return self.add.call(this, type, listener, options);
    };
    EventTarget.prototype.removeEventListener = function(type: string, listener: any, options?: any) {
      const c = capture(options);
      const i = self.entries.findIndex((e) => e.target === this && e.type === type && e.listener === listener &&
        e.capture === c);
      if (i >= 0)
        self.entries.splice(i, 1);
      return self.remove.call(this, type, listener, options);
    };
  }

  /** The listeners still on a target that outlives a closed view: window, document, an element in the page. */
  onLongLived(): ListenerEntry[] {
    return this.entries.filter(({target}) => target === window || target === document ||
      (target instanceof Node && target.isConnected));
  }

  stop(): void {
    EventTarget.prototype.addEventListener = this.add;
    EventTarget.prototype.removeEventListener = this.remove;
  }
}

/** What stays after an action: the live engine objects by class, the listeners on long-lived targets. */
interface Residue {objects: {[cls: string]: number}; listeners: number; where: string[]; others: string[]}

/**
 * Runs `action` once (the warm-up), then ROUNDS times, watching the module's objects and the page's listeners; returns
 * the residue after the warm-up and after the rounds.
 */
async function repeated(action: () => Promise<void>): Promise<{base: Residue; end: Residue}> {
  const crux = await loadCrux('tests');
  const objects = new LiveObjects(engineClasses(crux));
  const listeners = new Listeners();
  const residue = (): Residue => {
    const left = listeners.onLongLived();
    const chem = left.filter((e) => e.chem);
    return {objects: objects.counts(), listeners: chem.length,
      where: chem.map((e) => `${(e.target as any)?.constructor?.name}.${e.type}`),
      others: left.filter((e) => !e.chem)
        .map((e) => `${(e.target as any)?.constructor?.name}.${e.type} by ${e.origin}`)};
  };
  try {
    await action();
    await DG.delay(300);
    const base = residue();
    for (let i = 0; i < ROUNDS; i++)
      await action();
    await DG.delay(300);
    return {base, end: residue()};
  } finally {
    listeners.stop();
    objects.stop();
  }
}

/** The residue's growth over the baseline, worded for a failure. */
function growth(base: Residue, end: Residue): string {
  const classes = new Set([...Object.keys(base.objects), ...Object.keys(end.objects)]);
  const objects = [...classes].map((c) => `${c} ${base.objects[c] ?? 0} -> ${end.objects[c] ?? 0}`).join(', ');
  const extra = end.where.filter((w) => !base.where.includes(w));
  const othersExtra = [...end.others];
  for (const o of base.others) {
    const i = othersExtra.indexOf(o);
    if (i >= 0)
      othersExtra.splice(i, 1);
  }
  return `live objects: ${objects || 'none'}; Chem's listeners on long-lived targets ${base.listeners} -> ` +
    `${end.listeners}` + (extra.length ? ` (${extra.slice(0, 8).join(', ')})` : '') +
    `; others' (the platform's, other packages') ${base.others.length} -> ${end.others.length}` +
    (othersExtra.length ? ` (${othersExtra.slice(0, 8).join(', ')})` : '');
}

/** Checks that `end` is back at `base`: the same live objects, at most LISTENER_SLACK listeners more. */
function expectBaseline(base: Residue, end: Residue, what: string): void {
  const said = growth(base, end);
  const classes = new Set([...Object.keys(base.objects), ...Object.keys(end.objects)]);
  for (const c of classes)
    expect(end.objects[c] ?? 0, base.objects[c] ?? 0, `${what}: live ${c} after ${ROUNDS} rounds (${said})`);
  expect(end.listeners - base.listeners <= LISTENER_SLACK, true, `${what}: listeners left (${said})`);
}

category('Leaks: Crux', () => {
  let keptDocument: string | null = null;
  // the app shuts the context and help panels as it opens (reactions-app.ts shutSidePanels): put back after
  let keptPanels: [boolean, boolean] | null = null;

  before(async () => {
    keptDocument = window.localStorage.getItem(DOCUMENT_KEY);
    keptPanels = [grok.shell.windows.showContextPanel, grok.shell.windows.showHelp];
    if (!chemCommonRdKit.moduleInitialized) { // the reaction renderer is made with RDKit's module (its RDKit mode)
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  after(async () => {
    if (keptDocument === null)
      window.localStorage.removeItem(DOCUMENT_KEY);
    else
      window.localStorage.setItem(DOCUMENT_KEY, keptDocument);
    if (keptPanels)
      [grok.shell.windows.showContextPanel, grok.shell.windows.showHelp] = keptPanels;
  });

  test('the Reactions app opened and closed 10 times: no engine object or listener left behind', async () => {
    forgetDocument(); // the default example each time: what one open draws, every open draws
    const openAndClose = async () => {
      const view = await PackageFunctions.reactionsApp() as DG.View;
      grok.shell.addView(view);
      try {
        await awaitCheck(() => [...view.root.querySelectorAll('canvas')].some((c) => c.width > 0 && c.height > 0),
          'the app drew nothing', 10000);
        // the app's Reaction pane maps its steps by itself: every round waits for the answer, so that each
        // round makes what the warm-up made (a mapping landing before the close in some rounds only would add a mapped
        // step's scheme the baseline lacks); the first one loads the mapper (capped at 30 s)
        await awaitCheck(() => [...view.root.querySelectorAll(`.${MAPPING_CLASS}`)]
          .every((s) => !['waiting', 'loading', 'mapping'].includes((s as HTMLElement).dataset.state ?? '')),
        'the app\'s Reaction pane is still mapping', 60000);
        await DG.delay(300); // the pane's own drawings after the first
      } finally {
        view.close();
      }
    };
    const {base, end} = await repeated(openAndClose);
    console.log(`Leaks: Crux, the Reactions app: ${growth(base, end)}`); // kept in the run's logs, passed or not
    expectBaseline(base, end, 'the Reactions app');
  }, {timeout: 120000});

  test('reaction cells: the scheme cache keeps at most its size and frees every scheme it drops', async () => {
    const crux = await loadCrux('tests');
    const mols = (await grok.dapi.files.readCsv('System:AppData/Chem/tests/smi10K.csv')).col('smiles')!.toList();
    resetCruxSchemeCache(TEST_CACHE_SIZE);
    const objects = new LiveObjects(engineClasses(crux));
    const heapBefore = cruxHeapBytes();
    try {
      const canvas = new OffscreenCanvas(600, 120);
      const options = {drawOptions: chemDrawOptions()};
      for (let i = 0; i < DISTINCT_REACTIONS; i++) // as a reaction column's cells: a picture each, at one size
        cruxReactionPixels(canvas, `${mols[8000 + i]}.CC(=O)Cl>CCN(CC)CC>${mols[8001 + i]}`, 600, 120, 1, options);
      const live = objects.counts()['ChemReactionScheme'] ?? 0;
      const grew = ((cruxHeapBytes() - heapBefore) / 1048576).toFixed(1);
      expect(live <= TEST_CACHE_SIZE, true, `${live} drawn schemes not freed after ${DISTINCT_REACTIONS} distinct ` +
        `reactions through a cache of ${TEST_CACHE_SIZE} (the heap grew ${grew} MB)`);
    } finally {
      objects.stop();
      resetCruxSchemeCache();
    }
  }, {timeout: 60000});

  // Crux's cell pictures were kept by count alone (100, DG.LruCache's default): a picture is its cell's w x h x DPR^2 x
  // 4 bytes, so 100 tall cells at DPR 2 held 384 MB (300 x 400 CSS px). 100 distinct reactions drawn into 600 x 400
  // device px cells (960,000 B a picture, 96 MB unbounded): the pictures kept stay within 64 MB.
  test('reaction cells: Crux\'s pictures kept within their byte budget', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      await loadCrux('tests');
      const mols = (await grok.dapi.files.readCsv('System:AppData/Chem/tests/smi10K.csv')).col('smiles')!.toList();
      const renderer = await PackageFunctions.rdKitReactionRenderer();
      const canvas = ui.canvas(600, 400);
      const options = {drawOptions: chemDrawOptions()};
      for (let i = 0; i < PICTURES_DRAWN; i++) {
        const reaction = `${mols[6000 + i]}.CC(=O)Cl>>${mols[6001 + i]}`;
        renderer._drawReactionByCrux(0, 0, 600, 400, canvas, reaction, 1, options);
      }
      const kept = keptPictureBytes(renderer.cruxRendersCache);
      expect(kept <= PICTURE_BUDGET_BYTES, true, `${(kept / 1048576).toFixed(1)} MB of pictures kept after ` +
        `${PICTURES_DRAWN} cells of 600 x 400 device px (budget ${PICTURE_BUDGET_BYTES / 1048576} MB)`);
    });
  }, {timeout: 120000});

  // RDKit's cell pictures have the same byte budget as Crux's: the same 100 distinct reactions drawn by RDKit into
  // 600 x 400 px cells (960,000 B a picture, 96 MB unbounded) keep within 64 MB.
  test('reaction cells: RDKit\'s pictures kept within the same byte budget', async () => {
    const mols = (await grok.dapi.files.readCsv('System:AppData/Chem/tests/smi10K.csv')).col('smiles')!.toList();
    const renderer = await PackageFunctions.rdKitReactionRenderer();
    for (let i = 0; i < PICTURES_DRAWN; i++)
      renderer._fetchRender(600, 400, `${mols[6200 + i]}.CC(=O)Cl>>${mols[6201 + i]}`);
    const kept = keptPictureBytes(renderer.reactionRendersCache);
    expect(kept <= PICTURE_BUDGET_BYTES, true, `${(kept / 1048576).toFixed(1)} MB of RDKit pictures kept after ` +
      `${PICTURES_DRAWN} cells of 600 x 400 px (budget ${PICTURE_BUDGET_BYTES / 1048576} MB)`);
  }, {timeout: 120000});

  test('a 1,000-residue reaction (20,033 characters): its cell answers within 1 s, drawn or refused', async () => {
    await loadCrux('tests');
    const peptide = (n: number) => 'NCC(=O)' + 'NC(C)C(=O)'.repeat(n) + 'O';
    const text = `${peptide(1000)}.CC(=O)Cl>>${peptide(1000).replace(/^N/, 'CC(=O)N')}`;
    const shown = () => cruxReactionCounts.drawn + cruxReactionCounts.asText + cruxReactionCounts.tooLarge;
    const before = shown();
    const start = performance.now();
    cruxReactionPixels(new OffscreenCanvas(600, 120), text, 600, 120, 1, {drawOptions: chemDrawOptions()});
    const ms = performance.now() - start;
    expect(ms < HUGE_REACTION_BUDGET_MS, true, `the cell took ${ms.toFixed(0)} ms (budget ${HUGE_REACTION_BUDGET_MS})`);
    expect(shown(), before + 1, 'the cell shows the reaction or words');
  }, {timeout: 60000});

  test('a 1,000-residue reaction: the Reaction panel answers within 100 ms and reads no scheme of it', async () => {
    const crux = await loadCrux('tests');
    const peptide = (n: number) => 'NCC(=O)' + 'NC(C)C(=O)'.repeat(n) + 'O';
    const text = `${peptide(1000)}.CC(=O)Cl>>${peptide(1000).replace(/^N/, 'CC(=O)N')}`;
    const objects = new LiveObjects(engineClasses(crux));
    let widget: DG.Widget | null = null;
    try {
      const start = performance.now();
      widget = reactionPanel(text);
      const ms = performance.now() - start;
      // get_steps alone took 0.4-0.73 s on this value (grows with the square of a molecule's size): the panel froze
      expect(ms < HUGE_PANEL_BUDGET_MS, true, `the panel took ${ms.toFixed(0)} ms (budget ${HUGE_PANEL_BUDGET_MS})`);
      expect(objects.counts()['ChemReactionScheme'] ?? 0, 0, 'schemes read for a reaction too large to draw');
    } finally {
      objects.stop();
      widget?.detach();
    }
  }, {timeout: 60000});

  // a trap (a Rust panic, an out-of-bounds access) may leave the instance broken: a free() after one once read out of
  // bounds. The module recovers inside the call, so that a user never has to reload (crux/crux-recover.ts): the
  // trapped instance is never called again (each object rebuilt on a fresh one at its next call, a free of an old
  // one a no-op), and Crux stays on. crux-core's real test trap (the next engine call panics).
  test('a trap in a reaction drawing is recovered from: the cell is drawn, Crux stays on, the reaction calls go on',
    async () => {
      await inEngine(ChemicalEngine.Crux, async () => {
        await loadCrux('tests');
        const g = ui.canvas(300, 100).getContext('2d')!;
        const draw = (text: string) =>
          drawReactionByCrux(g, text, 0, 0, 300, 100, {drawOptions: chemDrawOptions()});
        const {drawn, asText} = cruxReactionCounts;
        expect(armCruxTrapForTests(1), true, 'the main thread\'s module is loaded');
        draw('CCCO>>CCCOC(C)=O');
        expect(cruxReactionCounts.drawn, drawn + 1, 'drawn, a trap in it');
        expect(isCruxRequested(), true, 'Crux stays on after a trap in a reaction drawing');
        expect(armCruxTrapForTests(1, 2), true);
        draw('CCCCO>>CCCCOC(C)=O');
        expect(cruxReactionCounts.asText, asText + 1, 'trapped on a fresh instance too: that cell written as text');
        expect(isCruxRequested(), true, 'Crux stays on');
        expect(cruxReadsReaction('CCCN>>CCCNC(C)=O'), true, 'a reaction read after the traps');
        expect(cruxReactionSteps('CCCS>>CCCSC')?.length, 1, 'steps after the traps');
      });
    }, {timeout: 60000});

  // `grok test --benchmark` only: the cell drawing of a 2,000-row reaction column scrolled down a page at a time and
  // back up, the scheme cache at 1,000 (its first size) and at Chem's size, interleaved; each page's draw time (ms),
  // medians over the reps of each pass's median, 95th percentile and worst page.
  test('reaction cells scrolled down and back up: page draw times by scheme cache size', async () => {
    if (!DG.Test.isInBenchmark)
      return;
    const mols = (await grok.dapi.files.readCsv('System:AppData/Chem/tests/smi10K.csv')).col('smiles')!.toList();
    const values = Array.from({length: SCROLL_ROWS}, (_, i) => `${mols[i]}.CC(=O)Cl>CCN(CC)CC>${mols[i + 1]}`);
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('reaction', values),
      DG.Column.fromStrings('molecule', mols.slice(0, SCROLL_ROWS))]);
    df.col('reaction')!.semType = 'ChemicalReaction';
    df.col('molecule')!.semType = DG.SEMTYPE.MOLECULE;
    const grid = df.plot.grid();
    const reactions = await PackageFunctions.rdKitReactionRenderer();
    const molecules = getMoleculeRenderer();
    const g = ui.canvas(SCROLL_CELL_W, SCROLL_CELL_H).getContext('2d')!;
    const stats = (xs: number[]) => {
      const s = [...xs].sort((a, b) => a - b);
      return {median: s[s.length >> 1], p95: s[Math.floor(0.95 * (s.length - 1))], max: s[s.length - 1]};
    };
    const out: {[leg: string]: {[pass: string]: number[]}} = {};
    /** Cells whose drawing threw, by leg: how many, and the first (a run once stopped on an uncaught
     * WebAssembly.Exception, leg unknown): the benchmark goes on and says so. */
    const thrown: {[leg: string]: {count: number, first: string}} = {};
    /** One leg: the column's cells drawn a page at a time down, then back up; each page's time. */
    const scroll = async (leg: string, column: string, renderer: DG.GridCellRenderer) => {
      const pages: number[] = [];
      for (let first = 0; first < SCROLL_ROWS; first += SCROLL_PAGE)
        pages.push(first);
      const page = (first: number): number => {
        const start = performance.now();
        for (let i = first; i < Math.min(SCROLL_ROWS, first + SCROLL_PAGE); i++) {
          const cell = grid.cell(column, i);
          try {
            renderer.render(g, 0, 0, SCROLL_CELL_W, SCROLL_CELL_H, cell, cell.style);
          } catch (e) {
            const t = (thrown[leg] ??= {count: 0, first: ''});
            if (t.count++ === 0) {
              t.first = `row ${i} (${String(cell.cell.value).substring(0, 160)}): ` +
                `${(e as object)?.constructor?.name} ${e instanceof Error ? e.message : String(e)}`;
            }
          }
        }
        return performance.now() - start;
      };
      const down = pages.map(page);
      await DG.delay(0);
      const up = [...pages].reverse().map(page);
      await DG.delay(0);
      const m = (out[leg] ??= {});
      for (const [pass, times] of [['down', down], ['up', up]] as [string, number[]][]) {
        for (const [k, v] of Object.entries(stats(times)))
          (m[`${pass} ${k}`] ??= []).push(v);
        (m[`${pass} total`] ??= []).push(times.reduce((a, b) => a + b, 0));
      }
    };
    const emptyMolecules = () => {
      (molecules as any).cacheGeneration = -1;
      (molecules as any)._syncEngine();
    };
    await inEngine(ChemicalEngine.Crux, async () => {
      for (let rep = 0; rep < SCROLL_REPS; rep++) {
        for (const size of rep % 2 ? [...SCROLL_CACHE_SIZES].reverse() : SCROLL_CACHE_SIZES) {
          resetCruxSchemeCache(size); // both caches empty: each leg starts cold
          reactions.cruxRendersCache = new DG.LruCache<string, ImageData>();
          await scroll(`cache ${size}`, 'reaction', reactions);
        }
      }
      resetCruxSchemeCache();
      // for reference: molecule cells of the same grid, the reactions' first molecules (Crux mode)
      emptyMolecules();
      await scroll('molecules, Crux', 'molecule', molecules);
    });
    // for reference only (RDKit mode is not changed): the same reactions drawn by RDKit
    await inEngine(ChemicalEngine.RDKit, async () => {
      reactions.reactionRendersCache = new DG.LruCache<string, ImageData>();
      await scroll('reactions, RDKit', 'reaction', reactions);
    });
    const medians = Object.fromEntries(Object.entries(out).map(([leg, m]) => [leg, Object.fromEntries(
      Object.entries(m).map(([k, v]) => [k, +[...v].sort((a, b) => a - b)[v.length >> 1].toFixed(2)]))]));
    console.log(`Leaks: Crux scroll: ${JSON.stringify(medians)} (all: ${JSON.stringify(out)}; thrown: ` +
      `${JSON.stringify(thrown)})`);
    return JSON.stringify({...medians, thrown});
  }, {benchmark: true, timeout: 900000});
});
