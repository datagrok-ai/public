// Leaks and speed (features/leaks/leaks-perf.feature): what
// repeated use and pathological input must never do in Crux mode. Read through the page and Chrome's DevTools protocol:
// - engine objects: FinalizationRegistry.prototype's register/unregister watched from a "notes" step (wasm-bindgen
// registers each object it makes and unregisters it in free()), by kind: a reaction scheme (it has get_steps), a
// Crux sketcher molecule (it has toMolblock), else its class name; held weakly, counted after a forced GC (CDP);
// - wasm memory: every WebAssembly memory an instance of the page exports, recorded as made by an init script the
// Background's first step adds (Runtime.queryObjects did not answer on Datagrok's page), read after a forced GC, each
// compared with itself;
// - detached sketcher roots: the crux-sketch elements ever connected (recorded by the same init script, weakly) alive
// after a forced GC and not in the document;
// - main-thread tasks: a PerformanceObserver of long tasks (over 50 ms) from a mark;
// - frames: the time from a scroll to the second animation frame after it.
// What Chem's package tests check of the same (category "Leaks: Crux", src/tests/leak-tests.ts) is not repeated.
import {CDPSession, Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {el, expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';
import {doubleClickArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {runTransformation} from './footguns.js';

declare const grok: any;
declare const DG: any;

/** Chem's scheme cache (src/rendering/crux-reaction-drawing.ts SCHEME_CACHE_SIZE): the most reaction schemes the
 * cells keep alive. */
const SCHEME_CACHE_BOUND = 1000;
/** The scroll steps: the frames measured, a page of rows at a time. */
const FRAME_SETTLE_FRAMES = 2;

/** The longest a DevTools call or a liveness check may take before the step says the page does not answer (ms). */
const CDP_ANSWER_MS = 30000;

/** `p`, or an error saying what did not answer within `ms`. */
function within<T>(p: Promise<T>, ms: number, what: string): Promise<T> {
  let timer: ReturnType<typeof setTimeout> | undefined;
  const late = new Promise<never>((_, reject) => {
    timer = setTimeout(() => reject(new Error(`${what} did not answer within ${ms / 1000} s`)), ms);
  });
  return Promise.race([p, late]).finally(() => clearTimeout(timer));
}

/** The page's DevTools session, after checking that the page's main thread answers (a script and a task): a call
 * into a page whose main thread is busy waits for it, and the step should say so rather than run out of time. */
const sessions = new WeakMap<Page, CDPSession>();
async function cdp(page: Page): Promise<CDPSession> {
  try {
    await within(page.evaluate(() => new Promise((r) => setTimeout(r, 0))), CDP_ANSWER_MS,
      'the main thread of the page (a script and a task)');
  } catch (e) {
    throw new Error(`${(e as Error).message}; where it was: ${await whereBusy(page)}`);
  }
  let s = sessions.get(page);
  if (!s) {
    s = await page.context().newCDPSession(page);
    sessions.set(page, s);
  }
  return s;
}

/** Where a busy main thread is: the debugger paused into it (a running script is interrupted; a thread blocked
 * outside JavaScript does not pause, and the answer says so), its top frames, then resumed. */
async function whereBusy(page: Page): Promise<string> {
  const s = await page.context().newCDPSession(page);
  try {
    const paused = new Promise<string>((resolve) => s.once('Debugger.paused', (e: any) => resolve(e.callFrames
      .slice(0, 10).map((f: any) => `${f.functionName || '(anonymous)'} ${String(f.url).split('/').pop()}:` +
        `${f.location.lineNumber + 1}`).join(' < '))));
    await within(s.send('Debugger.enable'), 10000, 'DevTools Debugger.enable');
    await within(s.send('Debugger.pause'), 10000, 'DevTools Debugger.pause');
    return await within(paused, 10000, 'a debugger pause (no script running: blocked outside JavaScript?)');
  } catch (e) {
    return `not paused: ${(e as Error).message}`;
  } finally {
    await s.send('Debugger.resume').catch(() => undefined);
    await s.detach().catch(() => undefined);
  }
}

/** A DevTools call that fails saying which one did not answer. */
function send<T>(s: CDPSession, method: string, params?: object): Promise<T> {
  return within(s.send(method as any, params as any) as Promise<T>, CDP_ANSWER_MS, `DevTools ${method}`);
}

/** A forced GC through the session; if it does not answer, the error says where the main thread was (a GC right
 * after the cell editor rounds once waited 110 s, until the test's time ran out). */
async function gcOn(page: Page, s: CDPSession): Promise<void> {
  try {
    await send(s, 'HeapProfiler.collectGarbage');
  } catch (e) {
    throw new Error(`${(e as Error).message}; where the main thread was: ${await whereBusy(page)}`);
  }
}

/** In the page, from its first script on (an init script, before the login's navigation): every WebAssembly memory an
 * instance exports (window.__leakMemsFromStart, with the module it came from: Chem's chem-wasm, crux-sketch's, RDKit's
 * or another), and every crux-sketch element ever connected, held weakly (window.__leakSketchers). DevTools'
 * Runtime.queryObjects found them after the fact on a small page, but on Datagrok's it did not answer in 30 s, so they
 * are recorded as they are made. */
function watchFromStart(): void {
  const w = window as any;
  if (w.__leakMemsFromStart)
    return;
  const mems: {tag: string, memory: WebAssembly.Memory}[] = [];
  w.__leakMemsFromStart = mems;
  const tag = (exported: object) => {
    const names = Object.keys(exported ?? {});
    return names.some((k) => k.startsWith('chemreactionscheme_')) ? 'chem-wasm' :
      names.some((k) => /sketch/i.test(k)) ? 'crux-sketch' :
        names.some((k) => /get_rxn|_malloc|emscripten/i.test(k)) ? 'rdkit' : `other (${names.length} exports)`;
  };
  const note = (instance: any) => {
    const memory = instance?.exports?.memory;
    if (memory instanceof WebAssembly.Memory && !mems.some((m) => m.memory === memory))
      mems.push({tag: tag(instance.exports), memory});
  };
  const instantiate = WebAssembly.instantiate;
  (WebAssembly as any).instantiate = async function(...args: any[]) {
    const r = await (instantiate as any).apply(this, args);
    note(r?.instance ?? r);
    return r;
  };
  const streaming = WebAssembly.instantiateStreaming;
  if (streaming) {
    (WebAssembly as any).instantiateStreaming = async function(...args: any[]) {
      const r = await (streaming as any).apply(this, args);
      note(r?.instance);
      return r;
    };
  }
  const Instance = WebAssembly.Instance;
  const Counted = function(module: WebAssembly.Module, imports?: any) {
    const i = new Instance(module, imports);
    note(i);
    return i;
  } as any;
  Counted.prototype = Instance.prototype;
  (WebAssembly as any).Instance = Counted;
  // crux-sketch elements: their class's connectedCallback records each one, weakly
  const sketchers: WeakRef<Element>[] = [];
  w.__leakSketchers = sketchers;
  const define = CustomElementRegistry.prototype.define;
  CustomElementRegistry.prototype.define = function(name: string, ctor: CustomElementConstructor, options?: any) {
    if (name === 'crux-sketch') {
      const connected = ctor.prototype.connectedCallback;
      ctor.prototype.connectedCallback = function(this: Element, ...args: any[]) {
        if (!sketchers.some((r) => r.deref() === this))
          sketchers.push(new WeakRef(this));
        return connected?.apply(this, args);
      };
    }
    return define.call(this, name, ctor, options);
  };
}

export const watchesFromStart = Given('the page records the WebAssembly memories and Crux sketchers it makes',
  async (page: Page) => {
    await page.addInitScript(watchFromStart);
    // a page an earlier feature left in the shell is not navigated again by the login step: reloaded here, so that
    // the script runs before the shell's first script
    const inShell = await page.evaluate(() => typeof (window as any).grok?.shell?.closeAll === 'function')
      .catch(() => false);
    if (inShell) {
      await page.reload({waitUntil: 'domcontentloaded', timeout: 180000});
      await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
    }
  }, {tier: 'api', description: 'an init script, before the login\'s navigation (a page already in the shell is reloaded): every WebAssembly memory an instance exports and every crux-sketch element ever connected, recorded as made (weakly for the elements); the steps that read memories and detached sketchers need it'});

/** After a forced GC: the byte length of every WebAssembly memory the page made, in the order made. */
async function wasmMemories(page: Page): Promise<number[]> {
  await gcOn(page, await cdp(page));
  const mems = await page.evaluate(() => ((window as any).__leakMemsFromStart as {tag: string, memory: any}[] | undefined)
    ?.map((m) => m.memory.buffer.byteLength));
  if (!mems)
    throw new Error('no memories recorded: the Background must begin with "the page records the WebAssembly ' +
      'memories and Crux sketchers it makes"');
  return mems;
}

/** The page's memories by module (for messages). */
const memoryTags = (page: Page) => page.evaluate(() =>
  ((window as any).__leakMemsFromStart as {tag: string}[] | undefined)?.map((m) => m.tag) ?? []);

/** After a forced GC: the crux-sketch elements alive and not in the document. */
async function detachedSketchers(page: Page): Promise<number> {
  await gcOn(page, await cdp(page));
  const n = await page.evaluate(() => ((window as any).__leakSketchers as WeakRef<Element>[] | undefined)
    ?.filter((r) => r.deref() && !r.deref()!.isConnected).length);
  if (n === undefined)
    throw new Error('no sketchers recorded: the Background must begin with "the page records the WebAssembly ' +
      'memories and Crux sketchers it makes"');
  return n;
}

/** In the page: engine objects counted from now on, by kind, until freed (window.__leakObjects); size observations
 * of the page's elements counted the same way (window.__leakObserved). Installed once a page. Everything is held
 * weakly, so the watch keeps nothing alive: an object the page dropped without free() is freed by wasm-bindgen's
 * finalizer once collected, and is no longer counted after a forced GC. */
async function watchObjects(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__leakObjects)
      return;
    const live: {[kind: string]: Set<WeakRef<object>>} = {};
    const refOf = new WeakMap<object, WeakRef<object>>();
    w.__leakObjects = live;
    const kindOf = (t: any) => typeof t?.get_steps === 'function' && typeof t?.draw_to_canvas === 'function' ?
      'reaction scheme' : typeof t?.toMolblock === 'function' ? 'Crux sketcher molecule' :
        `other ${t?.constructor?.name ?? '?'}`;
    const proto = FinalizationRegistry.prototype as any;
    const register = proto.register;
    const unregister = proto.unregister;
    proto.register = function(target: object, held: unknown, token?: object) {
      const ref = new WeakRef(target);
      refOf.set(token ?? target, ref); // wasm-bindgen's token is the object itself
      (live[kindOf(target)] ??= new Set()).add(ref);
      return register.call(this, target, held, token);
    };
    proto.unregister = function(token: object) {
      const ref = refOf.get(token);
      if (ref) {
        for (const set of Object.values(live))
          set.delete(ref);
      }
      return unregister.call(this, token);
    };
    // ResizeObserver observations (ui.onSizeChanged): each target observed and not unobserved, by observer
    const observed = new Map<WeakRef<ResizeObserver>, Set<WeakRef<Element>>>();
    const setOf = new WeakMap<ResizeObserver, Set<WeakRef<Element>>>();
    w.__leakObserved = observed;
    const ro = ResizeObserver.prototype as any;
    const observe = ro.observe;
    const unobserve = ro.unobserve;
    const disconnect = ro.disconnect;
    ro.observe = function(target: Element, options?: any) {
      let set = setOf.get(this);
      if (!set) {
        set = new Set();
        setOf.set(this, set);
        observed.set(new WeakRef(this), set);
      }
      if (!Array.from(set).some((r) => r.deref() === target))
        set.add(new WeakRef(target));
      return observe.call(this, target, options);
    };
    ro.unobserve = function(target: Element) {
      const set = setOf.get(this);
      set?.forEach((r) => r.deref() === target && set.delete(r));
      return unobserve.call(this, target);
    };
    ro.disconnect = function() {
      setOf.get(this)?.clear();
      return disconnect.call(this);
    };
  });
}

/** A forced garbage collection in the page (DevTools), after checking that its main thread answers. */
async function collectGarbage(page: Page): Promise<void> {
  await gcOn(page, await cdp(page));
}

/** In the page, after a forced GC: the engine objects made since the watch began, neither freed nor collected, by
 * kind. */
async function liveObjects(page: Page): Promise<{[kind: string]: number}> {
  await collectGarbage(page);
  return page.evaluate(() => Object.fromEntries(Object.entries(
    (window as any).__leakObjects as {[kind: string]: Set<WeakRef<object>>})
    .map(([k, s]) => [k, Array.from(s).filter((r) => r.deref() !== undefined).length])));
}

/** What a "notes" step keeps for its Then steps (per page). */
interface Noted {objects?: {[kind: string]: number}; mems?: number[]; detached?: number; canvases?: number;
  observed?: number}
const noted = new WeakMap<Page, Noted>();
const note = (page: Page, more: Noted) => noted.set(page, {...(noted.get(page) ?? {}), ...more});
function notedOf(page: Page): Noted {
  const n = noted.get(page);
  if (!n)
    throw new Error('nothing was noted: a "user notes ..." step comes first');
  return n;
}

// ---- the molecule cell editor ----

/** The cell editor of `area` opened (a double-click), Crux in it ready; with `edit`, a bond drawn (the molecule set
 * with one more carbon on its last atom, through the sketcher's host, as footguns.ts's "puts ... into the Crux
 * sketcher" does); then Cancel, and the dialog gone. */
async function editorRound(page: Page, area: string, edit: boolean): Promise<void> {
  await doubleClickArea(page, area, el('grid'));
  await expect.poll(() => page.evaluate(() => {
    const dialogs = Array.from(document.querySelectorAll('.d4-dialog')).filter((d) => d.querySelector('.crux-sketcher'));
    const root = dialogs[dialogs.length - 1]?.querySelector('.crux-sketcher');
    const impl = root ? DG.Widget.find(root) : null;
    return !!impl?.isInitialized && !impl.isRenderPending;
  }), {message: 'the cell editor with Crux ready', timeout: pollMs(30000)}).toBe(true);
  if (edit) {
    await page.evaluate(() => {
      const dialogs = Array.from(document.querySelectorAll('.d4-dialog')).filter((d) => d.querySelector('.crux-sketcher'));
      const root = dialogs[dialogs.length - 1]!.querySelector('.crux-sketcher')!;
      let host: any = null;
      for (let e = root.parentElement; e !== null && host === null; e = e.parentElement) {
        const w = DG.Widget.find(e);
        if (w && w.onSketcherReady)
          host = w;
      }
      if (!host)
        throw new Error('no sketcher host round Crux');
      host.setSmiles(`${host.getSmiles()}C`);
    });
  }
  await page.evaluate(() => {
    const d = Array.from(document.querySelectorAll('.d4-dialog')).filter((x) => x.querySelector('.crux-sketcher')).pop()!;
    const cancel = Array.from(d.querySelectorAll('.d4-dialog-footer button'))
      .find((b) => /^cancel$/i.test((b.textContent ?? '').trim())) as HTMLElement | undefined;
    if (!cancel)
      throw new Error('no CANCEL button in the cell editor');
    cancel.click();
  });
  await expect.poll(() => page.evaluate(() => document.querySelectorAll('.d4-dialog .crux-sketcher').length),
    {message: 'the cell editor closed', timeout: pollMs(15000)}).toBe(0);
}

export const editorWarmed = Given('user has opened and closed the cell editor of {string} {int} times',
  async (page: Page, area: string, times: number) => {
    await sketcherIs(page, 'Crux'); // the scenario's sketcher is Crux (its objects are what it counts)
    await watchObjects(page);
    for (let i = 0; i < times; i++)
      await editorRound(page, area, true);
  }, {tier: 'ui', description: 'the session\'s molecule sketcher pinned to Crux (the account\'s choice back at feature end), then the cell editor opened by a double-click on the cell, Crux ready, a bond drawn, CANCEL, the dialog gone, as many times as said'});

export const editorRounds = When('user opens and closes the cell editor of {string} {int} times, drawing a bond each time',
  async (page: Page, area: string, times: number) => {
    for (let i = 0; i < times; i++)
      await editorRound(page, area, true);
  }, {tier: 'ui', description: 'the cell editor opened by a double-click, Crux ready, a bond drawn (one more carbon on the last atom, set through the sketcher\'s host), CANCEL, the dialog gone; as many times as said'});

export const noteSketcher = When('user notes the live Crux sketcher objects and the wasm memory of every Crux instance',
  async (page: Page) => {
    await watchObjects(page);
    note(page, {objects: await liveObjects(page), mems: await wasmMemories(page), detached: await detachedSketchers(page)});
  }, {tier: 'api', description: 'after a forced GC: the engine objects made since the watch began and not freed, by kind; every WebAssembly memory\'s size; the detached crux-sketch elements'});

export const sketcherObjectsAsNoted = Then('the live Crux sketcher objects are as many as noted', async (page: Page) => {
  const before = notedOf(page).objects!;
  const now = await liveObjects(page);
  for (const kind of new Set([...Object.keys(before), ...Object.keys(now)])) {
    if (kind === 'reaction scheme')
      continue;
    expect(now[kind] ?? 0, `live ${kind} objects (noted ${before[kind] ?? 0}; now ${JSON.stringify(now)})`)
      .toBe(before[kind] ?? 0);
  }
}, {description: 'every kind of engine object but reaction schemes (the cells\' cache): as many made and not freed as when noted, exactly'});

export const memoryGrewAtMost = Then('the wasm memory of every Crux instance grew by at most {int} MiB',
  async (page: Page, mib: number) => {
    const before = notedOf(page).mems!;
    const now = await wasmMemories(page);
    const tags = await memoryTags(page);
    for (let i = 0; i < before.length; i++) {
      expect(now[i] - before[i], `wasm memory ${i + 1} of ${before.length} (${tags[i]}): ` +
        `${(before[i] / 1048576).toFixed(1)} → ${(now[i] / 1048576).toFixed(1)} MiB`)
        .toBeLessThanOrEqual(mib * 1048576);
    }
  }, {description: 'each WebAssembly memory noted (the Crux instances and RDKit\'s among them), compared with itself after a forced GC'});

export const noMoreDetached = Then('the page holds no more detached sketcher roots than before', async (page: Page) => {
  const before = notedOf(page).detached!;
  const now = await detachedSketchers(page);
  expect(now, `detached crux-sketch elements after a forced GC (noted ${before})`).toBeLessThanOrEqual(before);
}, {description: 'the crux-sketch elements alive after a forced GC and not in the document: no more than noted'});

// ---- the Reaction panel ----

export const openCorpus = Given('user opens the reaction corpus {string} of Chem', async (page: Page, file: string) => {
  const path = `System:AppData/Chem/${file.replace(/^files\//, '')}`;
  await page.evaluate(async (p) => {
    const df = await grok.dapi.files.readCsv(p);
    df.name = 'reactions-crux';
    await grok.data.detectSemanticTypes(df);
    grok.shell.addTableView(df);
    grok.shell.windows.showContextPanel = true;
  }, path);
  await expect.poll(() => page.evaluate(() => grok.shell.tv?.dataFrame?.name === 'reactions-crux' &&
    grok.shell.tv.dataFrame.col('reaction')?.semType === 'ChemicalReaction' && !grok.shell.tv.grid.isRenderPending),
  {message: 'the corpus open, its reaction column typed', timeout: pollMs(60000)}).toBe(true);
}, {tier: 'api', description: 'Chem\'s own file (files/... is System:AppData/Chem/...) read and opened as a table view, its semantic types detected, the context panel shown'});

/** Row `row` (1-based) of the reaction column made current by a click on its cell (scrolled into view first); the
 * Reaction pane's picture of its value awaited. (Setting the current cell and grok.shell.o from a script left the
 * panel on row 1 in a BDD run, though it followed in a probe of a running Datagrok: a click is what a user does.) */
async function makeCurrent(page: Page, row: number): Promise<number> {
  await page.evaluate((r) => grok.shell.tv.grid.scrollToCell('reaction', r - 1), row);
  const c = viewers.centerOf(await viewers.hitArea(page, el('grid'), `cell ${row} of reaction`, true));
  await page.mouse.click(c.x, c.y);
  return page.evaluate(async (r) => {
    const df = grok.shell.tv.dataFrame;
    const want = df.get('reaction', r - 1);
    const t0 = performance.now();
    for (let k = 0; k < 250; k++) {
      const pane = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((h) => (h.textContent ?? '').trim() === 'Reaction');
      if (pane && !pane.classList.contains('expanded'))
        (pane as HTMLElement).click();
      if (Array.from(document.querySelectorAll('canvas[data-reaction]')).some((c: any) => c.dataset.reaction === want &&
        c.width > 0))
        return performance.now() - t0;
      await new Promise((res) => setTimeout(res, 20));
    }
    throw new Error(`the Reaction pane drew no picture of row ${r}`);
  }, row);
}

export const madeCurrentOnce = Given('user has made rows {int} to {int} current in turn', async (page: Page,
  from: number, to: number) => {
  await watchObjects(page);
  for (let r = from; r <= to; r++) {
    await makeCurrent(page, r);
    // the pane maps a reaction by itself: its mapped picture's scheme is made here too, so that a row the
    // rounds leave current long enough to be mapped adds no scheme the count did not see (RXNMapper answers the same)
    await mappingSettled(page);
  }
  // then once as the rounds switch them, without waiting: chem-wasm's memory reaches its high-water mark in the first
  // quick round and stays there (measured: 5.50 | 5.50 | 5.50 | 5.50 | 5.50 MiB after rounds 1-5, from 4.3 noted after
  // the settled round), so the note is taken past it and the rounds measure growth per switch, which a leak shows
  for (let r = from; r <= to; r++)
    await makeCurrent(page, r);
  await mappingSettled(page);
}, {tier: 'api', description: 'each row\'s reaction cell made the current cell and the context panel\'s object, its Reaction pane opened and its picture drawn, and its automatic mapping answered, in turn; then the rows once more in quick turn, as the rounds go'});

export const madeCurrentRounds = When('user makes rows {int} to {int} current in turn, {int} times', async (page: Page,
  from: number, to: number, times: number) => {
  // each round's wasm memories said (a rise of 1.19 MiB over 5 rounds once needed telling apart): a rise in the first
  // round alone is a high-water mark, the same rise every round a leak
  const rounds: string[] = [];
  for (let t = 0; t < times; t++) {
    for (let r = from; r <= to; r++)
      await makeCurrent(page, r);
    rounds.push((await wasmMemories(page)).map((b) => (b / 1048576).toFixed(2)).join('/'));
  }
  console.log(`LEAKS rounds ${from}-${to} x${times}: wasm memories (MiB, ${(await memoryTags(page)).join('/')}) ` +
    `after each round: ${rounds.join(' | ')}`);
}, {tier: 'api', description: 'as "has made rows ... current in turn", the rows over and over'});

/** In the page: the Reaction pane's canvases in the document, and the size observations of canvases. */
const paneCounts = async (page: Page) => {
  await collectGarbage(page);
  return page.evaluate(() => {
    const observed = (window as any).__leakObserved as Map<WeakRef<ResizeObserver>, Set<WeakRef<Element>>>;
    let canvases = 0;
    for (const [observer, set] of observed) {
      if (observer.deref())
        set.forEach((r) => canvases += r.deref() instanceof HTMLCanvasElement ? 1 : 0);
    }
    return {canvases: document.querySelectorAll('canvas[data-reaction]').length, observed: canvases};
  });
};

/** The Reaction pane's automatic mapping answered wherever one is under way: a count taken while a
 * mapping is still to come would see its mapped picture's scheme appear after it. At most 15 s. */
const mappingSettled = (page: Page) => page.evaluate(async () => {
  const going = () => Array.from(document.querySelectorAll('.chem-reaction-mapping'))
    .some((s) => ['waiting', 'loading', 'mapping'].includes((s as HTMLElement).dataset.state ?? ''));
  for (let k = 0; k < 300 && going(); k++)
    await new Promise((r) => setTimeout(r, 50));
});

export const notePanel = When('user notes the live reaction schemes, the Reaction panel\'s canvases and their size listeners',
  async (page: Page) => {
    await watchObjects(page);
    await mappingSettled(page);
    const {canvases, observed} = await paneCounts(page);
    note(page, {objects: await liveObjects(page), mems: await wasmMemories(page), canvases, observed});
  }, {tier: 'api', description: 'the reaction schemes made since the watch began and not freed; the Reaction pane\'s canvases in the page and the canvases observed for their size (ui.onSizeChanged); every WebAssembly memory\'s size'});

export const schemesBounded = Then('the live reaction schemes are at most the cell cache\'s bound', async (page: Page) => {
  const live = (await liveObjects(page))['reaction scheme'] ?? 0;
  expect(live, `reaction schemes made and not freed (Chem's scheme cache keeps at most ${SCHEME_CACHE_BOUND})`)
    .toBeLessThanOrEqual(SCHEME_CACHE_BOUND);
}, {description: `the reaction schemes made since the watch began and not freed: at most Chem's scheme cache's ${SCHEME_CACHE_BOUND}`});

export const schemesAsNoted = Then('the live reaction schemes are as many as noted', async (page: Page) => {
  const before = notedOf(page).objects?.['reaction scheme'] ?? 0;
  await mappingSettled(page);
  const live = (await liveObjects(page))['reaction scheme'] ?? 0;
  expect(live, `reaction schemes made since the watch began, neither freed nor collected (noted ${before})`)
    .toBe(before);
}, {description: 'D-LEAK-3: the reaction schemes made since the watch began and alive after a forced GC, exactly as many as noted (the rows\' schemes were made before the note, so showing them again makes none)'});

export const paneAsNoted = Then('the Reaction panel\'s canvases and size listeners are as many as noted', async (page: Page) => {
  const before = notedOf(page);
  await mappingSettled(page);
  const now = await paneCounts(page);
  expect(now.canvases, 'the Reaction pane\'s canvases in the page').toBe(before.canvases);
  expect(now.observed, 'canvases observed for their size').toBeLessThanOrEqual(before.observed!);
}, {description: 'the Reaction pane\'s canvases in the page as many as noted, and no more canvases observed for their size'});

export const mainMemoryGrew = Then('the wasm memory of the main thread\'s Crux grew by at most {int} MiB',
  (page: Page, mib: number) => memoryGrewAtMost(page, mib),
{description: 'every WebAssembly memory of the page (the main thread\'s Crux among them) compared with itself after a forced GC'});

// ---- scrolling ----

/** A table of `n` distinct values: `kind` reaction (`A.CC(=O)Cl>CCN(CC)CC>B`), route (four molecules) or molecule,
 * from Chem's smi10K, typed before its view opens. */
async function openGenerated(page: Page, name: string, n: number, kinds: string[]): Promise<void> {
  await page.evaluate(async ([nm, rows, ks]) => {
    const mols = (await grok.dapi.files.readCsv('System:AppData/Chem/tests/smi10K.csv')).col('smiles').toList();
    const values: string[] = [];
    for (let i = 0; i < rows; i++) {
      const k = ks[i % ks.length];
      const m = (j: number) => mols[(4 * i + j) % mols.length];
      values.push(k === 'molecule' ? m(0) : k === 'route' ? `${m(0)}>>${m(1)}>>${m(2)}>>${m(3)}` :
        `${m(0)}.CC(=O)Cl>CCN(CC)CC>${m(1)}`);
    }
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings(ks[0] === 'molecule' ? 'molecule' : 'reaction', values)]);
    df.name = nm;
    await grok.data.detectSemanticTypes(df);
    grok.shell.addTableView(df);
  }, [name, n, kinds] as [string, number, string[]]);
  await expect.poll(() => page.evaluate((nm) => grok.shell.tv?.dataFrame?.name === nm && !grok.shell.tv.grid.isRenderPending,
    name), {message: `the table "${name}" open and drawn`, timeout: pollMs(60000)}).toBe(true);
}

export const scrollTable = Given('user opens a table of 2,000 reactions and routes made from the reaction corpus',
  async (page: Page) => {
    await watchObjects(page);
    await openGenerated(page, 'reactions-2000', 2000, ['reaction', 'reaction', 'reaction', 'route']);
  }, {tier: 'api', description: '2,000 distinct values, three reactions to a route, made of Chem\'s smi10K molecules (the corpus has 41 rows: repeated, every page after the first would hit the caches; dispute D-LEAK-1)'});

export const sizeColumn = Given('the reaction column of the grid is {int} px wide and its rows {int} px tall',
  async (page: Page, width: number, height: number) => {
    await page.evaluate(([w, h]) => {
      const grid = grok.shell.tv.grid;
      grid.props.rowHeight = h;
      grid.columns.byName(grid.dataFrame.columns.byIndex(0).name).width = w;
    }, [width, height] as [number, number]);
    await page.waitForFunction(() => !grok.shell.tv.grid.isRenderPending, null, {timeout: pollMs(30000)});
  }, {tier: 'api', description: 'the grid\'s row height and its first data column\'s width set'});

/** In the page: the grid scrolled a page at a time from the top to the last row; each frame's time (ms). */
const scrollFrames = (page: Page) => page.evaluate(async (settle) => {
  const grid = grok.shell.tv.grid;
  const df = grid.dataFrame;
  const col = df.columns.byIndex(0).name;
  const perPage = Math.max(1, Math.floor(grid.root.clientHeight / grid.props.rowHeight) - 1);
  grid.scrollToCell(col, 0);
  await new Promise((r) => setTimeout(r, 300));
  const frames: number[] = [];
  for (let top = perPage; top < df.rowCount + perPage; top += perPage) {
    const t0 = performance.now();
    grid.scrollToCell(col, Math.min(df.rowCount - 1, top));
    for (let f = 0; f < settle; f++)
      await new Promise((r) => requestAnimationFrame(() => r(null)));
    frames.push(performance.now() - t0);
  }
  return frames;
}, FRAME_SETTLE_FRAMES);

const scrolled = new WeakMap<Page, number[]>();

export const scrollDown = When('user scrolls the grid from the first row to the last, one page at a time',
  async (page: Page) => {
    scrolled.set(page, await scrollFrames(page));
  }, {tier: 'api', description: 'the grid scrolled to each next page of rows in turn; each page\'s time from the scroll to the second animation frame after it'});

export const noFrameOver = Then('no frame took more than {int} ms', async (page: Page, ms: number) => {
  const frames = scrolled.get(page)!;
  const worst = Math.max(...frames);
  expect(worst, `the worst of ${frames.length} frames (median ${[...frames].sort((a, b) => a - b)[frames.length >> 1]
    .toFixed(1)} ms)`).toBeLessThanOrEqual(ms);
}, {description: 'the worst page\'s time of the scroll'});

export const twiceMolecules = Then('the worst frame is at most twice the worst frame of the same grid scrolled over a molecule column',
  async (page: Page) => {
    const frames = scrolled.get(page)!;
    const size = await page.evaluate(() => {
      const grid = grok.shell.tv.grid;
      return {w: grid.columns.byIndex(1).width, h: grid.props.rowHeight, n: grid.dataFrame.rowCount, name: grok.shell.tv.dataFrame.name};
    });
    await openGenerated(page, 'molecules-same', size.n, ['molecule']);
    await page.evaluate(([w, h]) => {
      const grid = grok.shell.tv.grid;
      grid.props.rowHeight = h;
      grid.columns.byIndex(1).width = w;
    }, [size.w, size.h] as [number, number]);
    await page.waitForFunction(() => !grok.shell.tv.grid.isRenderPending, null, {timeout: pollMs(30000)});
    const molecules = await scrollFrames(page);
    await page.evaluate((nm) => {
      const view = grok.shell.tableView('molecules-same');
      view?.close();
      const back = grok.shell.tableView(nm);
      if (back)
        grok.shell.v = back;
    }, size.name);
    expect(Math.max(...frames), `the worst reaction frame against twice the worst molecule frame (${Math.max(...molecules)
      .toFixed(1)} ms)`).toBeLessThanOrEqual(2 * Math.max(...molecules));
  }, {tier: 'api', description: 'a table of as many distinct molecules (smi10K) in a grid of the same column width and row height, scrolled the same way; its view closed after'});

// ---- huge values ----

/** In the page: long tasks recorded from now on (window.__leakTasks: [start, duration] in ms). */
const watchTasks = (page: Page) => page.evaluate(() => {
  const tasks: [number, number][] = [];
  (window as any).__leakTasks = tasks;
  (window as any).__leakTaskObserver?.disconnect();
  const o = new PerformanceObserver((list) => list.getEntries().forEach((e) => tasks.push([e.startTime, e.duration])));
  o.observe({type: 'longtask', buffered: false});
  (window as any).__leakTaskObserver = o;
});
const longestTask = (page: Page) => page.evaluate(async () => {
  await new Promise((r) => setTimeout(r, 200)); // the observer's callback of the last task
  return Math.max(0, ...((window as any).__leakTasks as [number, number][]).map(([, d]) => d));
});

/** In the page: the time of the next click on a dialog's OK (window.__leakOkAt), from now on. */
const watchOk = (page: Page) => page.evaluate(() => {
  const w = window as any;
  w.__leakOkAt = undefined;
  if (w.__leakOkWatched)
    return;
  w.__leakOkWatched = true;
  document.addEventListener('click', (e) => {
    if ((e.target as Element | null)?.closest?.('.d4-dialog-footer [name="button-OK"]'))
      w.__leakOkAt = performance.now();
  }, true);
});

/** In the page: the time a column is next added to the current table (window.__leakColumnAddedAt): Chem's dialog adds
 * the products column empty and fills it in the same task, the platform's. */
const watchColumnAdded = (page: Page) => page.evaluate(() => {
  const w = window as any;
  w.__leakColumnAddedAt = undefined;
  w.__leakColumnSub?.unsubscribe();
  w.__leakColumnSub = grok.shell.tv.dataFrame.onColumnsAdded.subscribe(() => {
    w.__leakColumnAddedAt ??= performance.now();
  });
});

/** The long tasks seen since watchTasks, split at the OK click: from it on (a task running into it counts as after),
 * the products column's creation (the task holding the moment the column was added) set aside; every task as [start
 * relative to OK, duration] ms for the message. */
const tasksAroundOk = (page: Page) => page.evaluate(async () => {
  await new Promise((r) => setTimeout(r, 200));
  const w = window as any;
  w.__leakColumnSub?.unsubscribe();
  const ok = (w.__leakOkAt as number | undefined) ?? 0;
  const added = w.__leakColumnAddedAt as number | undefined;
  const tasks = w.__leakTasks as [number, number][];
  const isCreation = ([s, d]: [number, number]) => added !== undefined && s <= added && added <= s + d;
  const after = tasks.filter(([s, d]) => s + d >= ok && !isCreation([s, d])).map(([, d]) => d);
  const creation = tasks.filter(isCreation).map(([, d]) => d);
  return {longest: Math.max(0, ...after), creation: creation.length ? Math.max(...creation) : null,
    addedAt: added === undefined ? null : Math.round(added - ok),
    all: tasks.map(([s, d]) => [Math.round(s - ok), Math.round(d)])};
});

/** A table of one reaction value, its long tasks watched from before it opens. */
async function openOne(page: Page, name: string, column: string, value: string): Promise<void> {
  await watchTasks(page);
  await page.evaluate(async ([nm, c, v]) => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings(c, [v])]);
    df.name = nm;
    df.col(c).semType = 'ChemicalReaction';
    grok.shell.addTableView(df);
  }, [name, column, value] as [string, string, string]);
}

export const chainTable = Given('user opens a table {string} with a reaction column {string} holding a {int}-step chain {string}',
  (page: Page, name: string, column: string, steps: number, start: string) => {
    const chain = Array.from({length: steps + 1}, (_, i) => `${'C'.repeat(i + 1)}O`).join('>>');
    if (!chain.startsWith(start.replace(/…$/, '')))
      throw new Error(`the chain CO>>CCO>>… does not start as "${start}"`);
    return openOne(page, name, column, chain);
  }, {tier: 'api', description: 'a one-row table, its column typed ChemicalReaction, holding the chain CO>>CCO>>CCCO>>… of as many steps as said (the text given is its start); the page\'s long tasks recorded from before it opens'});

export const gridDraws = When('the grid draws {string}', async (page: Page, area: string) => {
  const [, row, column] = /^cell (\d+) of (.+)$/.exec(area) ?? [];
  if (!row)
    throw new Error(`"${area}" is no "cell N of <column>"`);
  await page.evaluate(([c, r]) => {
    grok.shell.tv.grid.scrollToCell(c, r - 1);
    grok.shell.tv.grid.invalidate();
  }, [column, +row] as [string, number]);
  await page.waitForFunction(() => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(() => r(true)))),
    null, {timeout: pollMs(120000)});
  await page.waitForFunction(() => !grok.shell.tv.grid.isRenderPending, null, {timeout: pollMs(120000)});
}, {tier: 'api', description: 'the cell scrolled into view and the grid repainted: two animation frames passed and no render is pending'});

export const firstDrawAtMost = Then('the cell\'s first draw took at most {int} ms', async (page: Page, ms: number) => {
  expect(await longestTask(page), 'the longest main-thread task since the table opened (ms)').toBeLessThanOrEqual(ms);
}, {description: 'the longest task since the table step (the cell\'s first draw is one task of the grid\'s first paint)'});

// ---- mapping ----

// The page checks each mapping's time; the mapping worker's memory is checked at engine level, in crux-core (memcheck's
// reaction-map section).

/** Per page: the next corpus row to try, the last row mapped, and the times of the mappings the When step made. */
const mapping = new WeakMap<Page, {next: number; last: number; times: [number, number][]}>();

/** In the Reaction pane of the current row: its automatic mapping's answer awaited (the panel maps by itself, with no
 * button): the Mapping section's `data-state` no longer waiting, loading or mapping. The time from
 * the pane's picture being drawn (makeCurrent) to the answer, which takes in the panel's 300 ms wait before it maps;
 * the state and the section's words. */
const mappingAnswer = (page: Page, value: string) => page.evaluate(async (want) => {
  const picture = Array.from(document.querySelectorAll('canvas[data-reaction]'))
    .find((c) => (c as HTMLCanvasElement).dataset.reaction === want);
  const root = picture?.closest('.d4-accordion-pane-content') ?? document.body;
  const status = root.querySelector('.chem-reaction-mapping') as HTMLElement | null;
  if (!status)
    return {ms: -1, state: 'none', text: 'no Mapping section in the Reaction pane'};
  const t0 = performance.now();
  for (let k = 0; k < 1600; k++) { // 40 s: the mapper's one-time load is capped at 30 s
    const state = status.dataset.state ?? '';
    if (!['waiting', 'loading', 'mapping'].includes(state))
      return {ms: performance.now() - t0, state, text: status.textContent ?? ''};
    await new Promise((r) => setTimeout(r, 25));
  }
  return {ms: performance.now() - t0, state: status.dataset.state ?? '', text: 'no answer within 40 s'};
}, value);

/** `n` mappings in the Reaction pane: the corpus's single reactions in turn (a route maps step by step, so it is left
 * out), each row made current by a click and its automatic mapping awaited; a row the pane refuses to map is passed over. The
 * rows go round again when they run out; the same row is never mapped twice in a row (the pane would not be rebuilt). */
async function mapReactions(page: Page, n: number, keep: boolean): Promise<void> {
  const state = mapping.get(page) ?? {next: 1, last: 0, times: []};
  mapping.set(page, state);
  const kinds: string[] = await page.evaluate(() => grok.shell.tv.dataFrame.col('kind').toList());
  let tried = 0;
  for (let done = 0; done < n;) {
    if (tried++ > 2 * kinds.length)
      throw new Error(`no row of the corpus was mapped after ${tried - 1} tries (${done} of ${n} mapped)`);
    const row = state.next;
    state.next = row % kinds.length + 1;
    if (!kinds[row - 1].startsWith('single'))
      continue;
    if (row === state.last)
      await makeCurrent(page, row === 1 ? 2 : 1);
    await makeCurrent(page, row);
    const value: string = await page.evaluate((r) => grok.shell.tv.dataFrame.get('reaction', r - 1), row);
    const {ms, state: answered, text} = await mappingAnswer(page, value);
    if (ms < 0)
      throw new Error(`row ${row}: ${text}`);
    if (answered === 'refused')
      continue;
    state.last = row;
    tried = 0;
    done++;
    if (keep)
      state.times.push([row, Math.round(ms)]);
  }
}

export const mappedSome = Given('user has mapped {int} reactions in the Reaction panel', (page: Page, n: number) =>
  mapReactions(page, n, false),
{tier: 'ui', description: 'the corpus\'s single reactions in turn: the row\'s cell clicked, the Reaction pane\'s automatic mapping awaited (its data-state no longer waiting, loading or mapping; a value the pane refuses without asking the mapper, a template or an unreadable one, is passed over); the first mapping waits for the mapper\'s one-time load, so these are not timed'});

export const mapMore = When('user maps {int} more reactions in the Reaction panel, one after another', (page: Page,
  n: number) => mapReactions(page, n, true),
{tier: 'ui', description: 'as "has mapped ... reactions", going on through the rows (round again when they run out), each mapping timed from the pane\'s picture drawn to its answer (the Mapping section: "Mapped (confidence ...)", "Not mapped: ..." or "Mapping took longer than 5 s"; the panel waits 300 ms before it maps)'});

export const eachMapping = Then('each mapping answered within {int} s', async (page: Page, s: number) => {
  const times = mapping.get(page)?.times ?? [];
  expect(times.length, 'mappings timed').toBeGreaterThan(0);
  expect(Math.max(...times.map(([, ms]) => ms)), `the slowest mapping (ms); every one as [row, ms]: ` +
    JSON.stringify(times)).toBeLessThanOrEqual(s * 1000);
}, {description: 'every mapping the When step made, from the pane\'s picture drawn to the Mapping section\'s answer'});

// ---- Transform ----

export const bigMolecules = Given('user opens a table of 50,000 molecules {string} cut to 50,000 rows', async (page: Page,
  file: string) => {
  await page.evaluate(async (f) => {
    const all = (await grok.dapi.files.readBinaryDataFrames(`System:AppData/Chem/tests/${f}.zip`))[0];
    const keep = DG.BitSet.create(all.rowCount, (i: number) => i < 50000);
    const df = all.clone(keep);
    df.name = f;
    await grok.data.detectSemanticTypes(df);
    grok.shell.addTableView(df);
  }, file);
  await expect.poll(() => page.evaluate((f) => grok.shell.tv?.dataFrame?.name === f && grok.shell.tv.dataFrame.rowCount,
    file), {message: `the table ${file} of 50,000 rows`, timeout: pollMs(120000)}).toBe(50000);
}, {tier: 'api', description: 'Chem\'s files/tests/<name>.zip read, its first 50,000 rows opened as a table view, typed'});

export const runOnMolecules = When('user runs the transformation {string} on the molecule column', async (page: Page,
  name: string) => {
  const column = await page.evaluate(() => grok.shell.tv.dataFrame.columns.toList()
    .find((c: any) => c.semType === 'Molecule')?.name ?? '');
  if (!column)
    throw new Error('the table has no molecule column');
  await watchTasks(page);
  await watchOk(page);
  await watchColumnAdded(page);
  await runTransformation(page, name, column);
  // done when Chem says so ("Reaction complete: N products, ...", reactions.ts); timed from the dialog's OK
  await expect.poll(() => page.evaluate(() => ((window as any).__guardBalloons as string[])
    .some((b) => /Reaction complete:/.test(b))), {message: 'the transformation\'s "Reaction complete" message',
    timeout: pollMs(300000)}).toBe(true);
  transformSeconds.set(page, await page.evaluate(() => (performance.now() - (window as any).__leakOkAt) / 1000));
}, {tier: 'ui', description: 'GUARD\'s "user runs the transformation ... on the ... column" on the table\'s molecule column, the page\'s long tasks recorded from the menu on and the dialog\'s OK click timed, until Chem\'s "Reaction complete: ..." message shows; its time from OK kept'});

const transformSeconds = new WeakMap<Page, number>();

export const noTaskOverMs = Then('no main-thread task took more than {int} ms while it ran, the products column\'s creation aside',
  async (page: Page, ms: number) => {
    const {longest, addedAt, all} = await tasksAroundOk(page);
    expect(longest, 'the longest main-thread task from the dialog\'s OK on, the products column\'s creation aside ' +
      `(ms); the column added ${addedAt} ms after OK; every long task as [start from OK, ms]: ${JSON.stringify(all)}`)
      .toBeLessThanOrEqual(ms);
  }, {description: 'D-LEAK-4: the longest task the long-task observer saw from the dialog\'s OK click on (a task running into the click counts), leaving out the one task in which the products column was added (the platform builds its category dictionary there); the dialog\'s own opening comes before OK and is listed in the message, not judged'});

export const columnCreatedWithin = Then('the platform created the products column in one task of at most {int} ms',
  async (page: Page, ms: number) => {
    const {creation, addedAt, all} = await tasksAroundOk(page);
    expect(addedAt, 'the products column was added (ms after OK)').not.toBeNull();
    // a creation under 50 ms is no long task at all: the observer sees none, and that is within any budget
    expect(creation ?? 0, `the task in which the products column was added, ${addedAt} ms after OK (ms; the ` +
      `platform's: a 50,000-value string column's category dictionary built and sorted); every long task as ` +
      `[start from OK, ms]: ${JSON.stringify(all)}`).toBeLessThanOrEqual(ms);
  }, {description: 'D-LEAK-4: the long task holding the moment the products column was added to the table (Chem adds it empty and fills it in that task; the cost is the platform\'s string column build), at most this long; none (under 50 ms) passes'});

export const finishedWithin = Then('it finished within {int} s with a message counting products, errors and rows with no reaction',
  async (page: Page, s: number) => {
    const message = await page.evaluate(() => ((window as any).__guardBalloons as string[])
      .find((b) => /Reaction complete:/.test(b)) ?? '');
    expect(message, 'Chem\'s message counts products, errors and rows with no reaction')
      .toMatch(/\d+ products, \d+ errors, \d+ no reaction/);
    expect(transformSeconds.get(page)!, 'the transformation\'s time to its message (s)').toBeLessThanOrEqual(s);
  }, {description: 'Chem\'s "Reaction complete: N products, N errors, N no reaction, ..." message, and the time from the start to it'});
