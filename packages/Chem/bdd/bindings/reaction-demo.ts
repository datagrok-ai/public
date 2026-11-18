// The reaction showcase clips (features/demos/reaction-sketcher.feature and reaction-tables.feature: a demo of
// Chem's reaction features, filmed with the BDD demo recorder). What no other binding has: the
// Reaction Sketcher opened on a reaction of the clip's own, and the gestures on Crux's canvas that act on what Crux
// draws rather than on an atom or a toolbar control: an arrow drawn after a molecule, an arrow's "+" handle and its
// menu, a condition written over an arrow and dragged, a molecule pasted where the pointer is, a drag from atom to
// atom. They read Crux's canvas inside its open shadow root (its arrows `g[data-arrow]`, its texts
// `g[data-object][data-kind="text"]`, its overlays and menus by test id: crux-sketch docs/conventions/test-ids.md),
// in whichever sketcher the page shows: the reaction editor's dialog over the app's. The checks read the app's Text,
// a route step by step, each molecule canonical by RDKit. What the canvas does with each gesture is crux-sketch's own
// BDD (features/rxn); what the app and the cells keep is Chem's package tests ('Reactions: app', 'Reactions: CDXML').
import {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, expect, gestures, locate, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

/** Crux in the page: the reaction editor's (a dialog over the app) first, else the Reaction Sketcher's. */
const SKETCHES = ['.d4-dialog:has(.chem-reaction-cell-editor) crux-sketch', '.chem-reactions-app-sketcher crux-sketch'];
/** reaction-document.ts DOCUMENT_KEY and its version: what the app reopens on. */
const DOCUMENT_KEY = 'Chem.ReactionsApp.document';
const DOCUMENT_VERSION = 2;

// ---------------------------------------------------------------- the app on a reaction of the clip's own

/** What the browser kept as the app's document before the feature's first clip, put back at feature end (once). */
const keeping = new WeakSet<Page>();
async function keepDocument(page: Page): Promise<void> {
  if (keeping.has(page))
    return;
  keeping.add(page);
  const kept = await page.evaluate((k) => window.localStorage.getItem(k), DOCUMENT_KEY);
  atFeatureEnd(page, async () => {
    keeping.delete(page);
    await page.evaluate(([k, v]) => v === null ? window.localStorage.removeItem(k) : window.localStorage.setItem(k, v),
      [DOCUMENT_KEY, kept] as [string, string | null]);
  });
}

/** Done when the Reaction Sketcher is the current view and its sketcher draws its reaction. */
async function sketcherShown(page: Page): Promise<void> {
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
    {message: 'the current view', timeout: pollMs(30000)}).toBe('Reaction Sketcher');
  await expect.poll(() => page.evaluate((sel) => {
    const el = document.querySelector(sel) as any;
    try {
      return !!el && el.isEmpty === false && el.isPending === false;
    } catch {
      return false;
    }
  }, SKETCHES[1]), {message: 'the app\'s sketcher, drawn', timeout: pollMs(60000)}).toBe(true);
}

/** Opens the Reaction Sketcher on `text` (kept in the browser as the app keeps a document). */
async function openSketcher(page: Page, title: string, text: string, mapNumbers: boolean): Promise<void> {
  await keepDocument(page);
  await page.evaluate(async ([k, doc]) => {
    DG.Dialog.getOpenDialogs().forEach((d: any) => d.close());
    for (const v of [...grok.shell.views].filter((v: any) => v.name === 'Reaction Sketcher'))
      v.close();
    window.localStorage.setItem(k, doc);
    grok.shell.addView(await grok.functions.call('Chem:reactionsApp'));
  }, [DOCUMENT_KEY, JSON.stringify({version: DOCUMENT_VERSION, title, text, mapNumbers})] as [string, string]);
  await sketcherShown(page);
  // the drawing's first fit and the preview beside it
  await page.waitForTimeout(500);
  await parkPointer(page);
}

/** The pointer rested in the canvas's empty top-right corner, where it shows no tooltip. */
async function parkPointer(page: Page): Promise<void> {
  const {box} = await canvasOf(page);
  await page.mouse.move(box.x + box.width * 0.9, box.y + box.height * 0.08);
}

export const opensOnExample = Given('the Reaction Sketcher opens on its example', async (page: Page) => {
  await keepDocument(page);
  await page.evaluate((k) => {
    for (const v of [...grok.shell.views].filter((v: any) => v.name === 'Reaction Sketcher'))
      v.close();
    window.localStorage.removeItem(k);
  }, DOCUMENT_KEY);
}, {tier: 'api', description: 'nothing kept in the browser as the app\'s document, so the app opens on its default ' +
  'example; what was kept comes back at feature end'});

export const showsExample = Then('the Reaction Sketcher should show its example route', async (page: Page) => {
  await sketcherShown(page);
  await expect.poll(async () => (await arrows(page)).length, {message: 'the example\'s five arrows',
    timeout: pollMs(10000)}).toBe(5);
}, {description: 'the app is the current view, its sketcher drawing the default example: five step arrows'});

export const sketcherOnExample = Given('the Reaction Sketcher is open on its example', async (page: Page) => {
  await keepDocument(page);
  await page.evaluate(async (k) => {
    DG.Dialog.getOpenDialogs().forEach((d: any) => d.close());
    for (const v of Array.from(grok.shell.views as any[]).filter((v: any) => v.name === 'Reaction Sketcher'))
      v.close();
    window.localStorage.removeItem(k);
    grok.shell.addView(await grok.functions.call('Chem:reactionsApp'));
  }, DOCUMENT_KEY);
  await sketcherShown(page);
  await page.waitForTimeout(500);
  await parkPointer(page);
}, {tier: 'api', description: 'the app run with nothing kept in the browser: its default example (five mapped steps), ' +
  'drawn; what the browser kept comes back at feature end'});

export const sketcherOn = Given('the Reaction Sketcher is open on {string}', (page: Page, text: string) =>
  openSketcher(page, '', text, false),
{tier: 'api', description: 'Chem\'s Reaction Sketcher (Chem:reactionsApp) on this reaction or route as the document the ' +
  'app reopens on, untitled, map numbers off; drawn; what the browser kept comes back at feature end'});

export const sketcherOnTitled = Given('the Reaction Sketcher is open on {string}, titled {string}',
  (page: Page, text: string, title: string) => openSketcher(page, title, text, false),
  {tier: 'api', description: 'the same, with the document\'s title'});

// ---------------------------------------------------------------- Crux's canvas, in viewport px

type Point = {x: number, y: number};
type Box = {x: number, y: number, width: number, height: number};

/** In the page: the selector of the sketcher shown (the editor's dialog over the app). */
async function sketchSelector(page: Page): Promise<string> {
  const sel = await page.evaluate((sels) => sels.find((s) => {
    const el = document.querySelector(s) as HTMLElement | null;
    return !!el && el.getBoundingClientRect().width > 0;
  }) ?? null, SKETCHES);
  if (sel === null)
    throw new Error('no Crux sketcher shown: neither the reaction editor nor the Reaction Sketcher');
  return sel;
}

/** The arrows Crux draws, in drawing order: each its tail and its head's tip (viewport px). */
async function arrows(page: Page): Promise<{tail: Point, head: Point}[]> {
  const sel = await sketchSelector(page);
  return page.evaluate((s) => {
    const svg = document.querySelector(s)!.shadowRoot!.querySelector('[data-testid="canvas"]')!;
    const nums = (n: Element | null) => (n?.getAttribute('d') ?? '').match(/-?\d+(?:\.\d+)?/g)?.map(Number) ?? [];
    const at = (n: SVGGraphicsElement, x: number, y: number) => {
      const p = new DOMPoint(x, y).matrixTransform(n.getScreenCTM()!);
      return {x: p.x, y: p.y};
    };
    return Array.from(svg.querySelectorAll('g[data-arrow]')).map((g) => {
      const shaft = g.querySelector('path.arrow-shaft') as SVGGraphicsElement;
      const head = g.querySelector('path.arrow-head') as SVGGraphicsElement;
      const [a, b] = [nums(shaft), nums(head)];
      return {tail: at(shaft, a[0], a[1]), head: at(head, b[0], b[1])};
    });
  }, sel);
}

async function arrowOf(page: Page, n: number): Promise<{tail: Point, head: Point}> {
  let all: {tail: Point, head: Point}[] = [];
  await expect.poll(async () => (all = await arrows(page)).length, {message: `arrow ${n} drawn`,
    timeout: pollMs(10000)}).toBeGreaterThanOrEqual(n);
  return all[n - 1];
}

const middleOf = (a: {tail: Point, head: Point}): Point => ({x: (a.tail.x + a.head.x) / 2, y: (a.tail.y + a.head.y) / 2});

/** Crux's drawn bond at the current zoom (px), and its canvas box (viewport px). */
async function canvasOf(page: Page): Promise<{bond: number, box: Box}> {
  const sel = await sketchSelector(page);
  return page.evaluate((s) => {
    const el = document.querySelector(s) as any;
    const r = el.shadowRoot.querySelector('[data-testid="canvas"]').getBoundingClientRect();
    return {bond: el.positions.bondLength as number, box: {x: r.x, y: r.y, width: r.width, height: r.height}};
  }, sel);
}

/** The box of the molecules' atoms (viewport px): every atom, or those of the atoms numbered `from` to `to`. */
async function atomsBox(page: Page, from?: number, to?: number): Promise<Box> {
  const sel = await sketchSelector(page);
  return page.evaluate(([s, f, t]) => {
    const el = document.querySelector(s as string) as any;
    const r = el.getBoundingClientRect();
    const ps = (el.positions.atoms as ({x: number, y: number} | null)[])
      .filter((p, i) => p !== null && (f === null || (i >= (f as number) && i <= (t as number)))) as {x: number, y: number}[];
    const xs = ps.map((p) => p.x + r.x);
    const ys = ps.map((p) => p.y + r.y);
    return {x: Math.min(...xs), y: Math.min(...ys), width: Math.max(...xs) - Math.min(...xs),
      height: Math.max(...ys) - Math.min(...ys)};
  }, [sel, from ?? null, to ?? null] as [string, number | null, number | null]);
}

/** A drag from `a` to `b`, as a person draws one: pressed, moved in steps, released. */
async function drag(page: Page, a: Point, b: Point, steps = 12): Promise<void> {
  await page.mouse.move(a.x, a.y);
  await page.mouse.down();
  await page.mouse.move((a.x + b.x) / 2, (a.y + b.y) / 2, {steps: Math.ceil(steps / 2)});
  await page.mouse.move(b.x, b.y, {steps: Math.ceil(steps / 2)});
  await page.mouse.up();
}

// ---------------------------------------------------------------- gestures

export const arrowUnder = When('user draws an arrow {int} bond lengths long under the drawing',
  async (page: Page, bonds: number) => {
    const m = await atomsBox(page);
    const {bond} = await canvasOf(page);
    // four bonds down: room above the arrow for a reagent, and its product clear of the row above
    const y = m.y + m.height + 4 * bond;
    await drag(page, {x: m.x, y}, {x: m.x + bonds * bond, y});
  }, {tier: 'ui', description: 'with the tool chosen (the reaction arrow tool), a drag four bonds under the drawing\'s ' +
    'atoms, from their left, that many bonds across: Crux continues the route from the row\'s last product'});

export const clickTool = When('user clicks the {string} tool of the sketcher', async (page: Page, name: string) => {
  const sel = await sketchSelector(page);
  const at = await page.evaluate(([s, n]) => {
    const el = document.querySelector(s) as any;
    const c = el.layout.controls[n];
    if (!c)
      return null;
    const r = el.getBoundingClientRect();
    return {x: r.x + c.x + c.width / 2, y: r.y + c.y + c.height / 2};
  }, [sel, name] as [string, string]);
  expect(at, `Crux shows the control "${name}"`).not.toBeNull();
  await page.mouse.move(at!.x, at!.y, {steps: 8});
  await page.mouse.down();
  await page.mouse.up();
}, {tier: 'ui', description: 'a click on a control of Crux\'s toolbars by its name in Crux\'s toolbar configuration ' +
  '(its layout: reaction.text, reaction.map, reaction.map.open, reaction.cycles, select, copy, paste)'});

export const hoverArrow = When('user hovers over the middle of arrow {int}', async (page: Page, n: number) => {
  const m = middleOf(await arrowOf(page, n));
  await page.mouse.move(m.x - 4, m.y - 2);
  await page.mouse.move(m.x, m.y, {steps: 3});
}, {tier: 'ui', description: 'the pointer on the middle of the arrow (arrows counted from 1, in drawing order): a step ' +
  'arrow shows its "+" handle'});

/** The arrow's "+" handle as the canvas shows it (overlay.step-add): its middle, null while hidden. */
async function addHandle(page: Page): Promise<Point | null> {
  const sel = await sketchSelector(page);
  return page.evaluate((s) => {
    const g = document.querySelector(s)!.shadowRoot!.querySelector('[data-testid="overlay.step-add"]');
    if (!g || g.getAttribute('display') === 'none')
      return null;
    const r = g.getBoundingClientRect();
    return r.width > 0 ? {x: r.x + r.width / 2, y: r.y + r.height / 2} : null;
  }, sel);
}

export const clickAddHandle = When('user clicks the "+" handle of the arrow', async (page: Page) => {
  let at: Point | null = null;
  await expect.poll(async () => (at = await addHandle(page)) !== null, {message: 'the arrow\'s "+" handle shown',
    timeout: pollMs(5000)}).toBe(true);
  await page.mouse.move(at!.x, at!.y, {steps: 3});
  await page.mouse.down();
  await page.mouse.up();
}, {tier: 'ui', description: 'a click on the "+" handle a hovered step arrow shows (overlay.step-add): its menu opens'});

export const chooseInMenu = When('user chooses {string} in the sketcher\'s menu', async (page: Page, item: string) => {
  const sel = await sketchSelector(page);
  const entry = page.locator(`${sel} [role^="menuitem"]`).filter({visible: true}).filter({hasText: item}).first();
  await expect(entry, `"${item}" in Crux's menu`).toBeVisible({timeout: pollMs(5000)});
  await entry.hover();
  await entry.click();
}, {tier: 'ui', description: 'an item of the menu Crux shows (a handle\'s, an atom\'s, a text\'s), by its name'});

export const pasteAt = When('user pastes {string} at {int}% across and {int}% down the canvas',
  async (page: Page, text: string, across: number, down: number) => {
    const {box} = await canvasOf(page);
    await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);
    await page.evaluate((t) => navigator.clipboard.writeText(t), text);
    await page.mouse.move(box.x + box.width * across / 100, box.y + box.height * down / 100, {steps: 6});
    await gestures.press(page, 'Control+V');
  }, {tier: 'ui', description: 'the text on the clipboard, the pointer on that point of Crux\'s canvas, Ctrl+V: Crux ' +
    'pastes where the pointer is (or where a mode it is in puts it)'});

export const pasteHere = When('user pastes {string}', async (page: Page, text: string) => {
  await page.context().grantPermissions(['clipboard-read', 'clipboard-write']);
  await page.evaluate((t) => navigator.clipboard.writeText(t), text);
  await gestures.press(page, 'Control+V');
}, {tier: 'ui', description: 'the text on the clipboard, Ctrl+V with the pointer where it is'});

export const writeOver = When('user writes {string} over arrow {int}', async (page: Page, text: string, n: number) => {
  const m = middleOf(await arrowOf(page, n));
  const {bond} = await canvasOf(page);
  await page.mouse.move(m.x, m.y - bond / 3, {steps: 6});
  await page.mouse.down();
  await page.mouse.up();
  const sel = await sketchSelector(page);
  const editor = page.locator(`${sel} [data-testid="canvas.text-editor"]`).first();
  await expect(editor, 'the text editor over the arrow').toBeVisible({timeout: pollMs(5000)});
  await page.keyboard.type(text, {delay: 90});
  await gestures.press(page, 'Escape');
  await expect(editor, 'the text editor closed').toBeHidden({timeout: pollMs(5000)});
}, {tier: 'ui', description: 'with the Text tool chosen, a click just above the middle of the arrow, the text typed, ' +
  'Escape: one edit, the step\'s conditions'});

/** The texts Crux draws: their words and boxes (viewport px). */
async function texts(page: Page): Promise<{words: string, box: Box}[]> {
  const sel = await sketchSelector(page);
  return page.evaluate((s) => Array.from(document.querySelector(s)!.shadowRoot!.querySelectorAll('g[data-object][data-kind="text"]'))
    .map((g) => {
      const r = g.getBoundingClientRect();
      return {words: Array.from(g.querySelectorAll('text')).map((t) => t.textContent ?? '').join(' '),
        box: {x: r.x, y: r.y, width: r.width, height: r.height}};
    }), sel);
}

async function textOf(page: Page, words: string): Promise<{words: string, box: Box}> {
  let all: {words: string, box: Box}[] = [];
  const plain = (s: string) => s.replace(/[₀-₉]/g, (c) => String(c.charCodeAt(0) - 0x2080)).replace(/\s+/g, '');
  await expect.poll(async () => (all = await texts(page)).some((t) => plain(t.words) === plain(words)),
    {message: `the text "${words}" drawn`, timeout: pollMs(5000)}).toBe(true);
  return all.find((t) => plain(t.words) === plain(words))!;
}

export const dragTextUnder = When('user drags the text {string} under arrow {int}', async (page: Page, words: string,
  n: number) => {
  const t = await textOf(page, words);
  const a = await arrowOf(page, n);
  const {bond} = await canvasOf(page);
  const from = {x: t.box.x + t.box.width / 2, y: t.box.y + t.box.height / 2};
  const m = middleOf(a);
  await drag(page, from, {x: m.x, y: m.y + t.box.height / 2 + bond / 3}, 16);
}, {tier: 'ui', description: 'with a selection tool, the text pressed in its middle and dragged just under the arrow\'s ' +
  'middle, in steps: Crux shows where it will land while it moves'});

export const rightClickText = When('user right-clicks the text {string}', async (page: Page, words: string) => {
  const t = await textOf(page, words);
  await page.mouse.move(t.box.x + t.box.width / 2, t.box.y + t.box.height / 2, {steps: 6});
  await page.mouse.down({button: 'right'});
  await page.mouse.up({button: 'right'});
}, {tier: 'ui', description: 'a right-click in the middle of the text Crux draws: its menu opens'});

export const dragAtoms = When('user drags from the {string} area to the {string} area of the sketcher',
  async (page: Page, from: string, to: string) => {
    const sel = await sketchSelector(page);
    const at = (name: string) => page.evaluate(([s, i]) => {
      const el = document.querySelector(s as string) as any;
      const p = el.positions.atoms[i as number];
      const r = el.getBoundingClientRect();
      return p ? {x: p.x + r.x, y: p.y + r.y} : null;
    }, [sel, Number(name.replace(/^atom /, ''))] as [string, number]);
    const [a, b] = [await at(from), await at(to)];
    expect(a, `${from} drawn`).not.toBeNull();
    expect(b, `${to} drawn`).not.toBeNull();
    await drag(page, a!, b!, 16);
  }, {tier: 'ui', description: 'a drag from one atom of Crux to another ("atom N", by the molfile\'s numbering), in steps ' +
    '(the Map tool pairs them)'});

export const clickCanvasBelow = When('user clicks the canvas {int} bond lengths under the drawing', async (page: Page,
  bonds: number) => {
  const m = await atomsBox(page);
  const {bond} = await canvasOf(page);
  const p = {x: m.x + m.width / 2, y: m.y + m.height + bonds * bond};
  await page.mouse.move(p.x, p.y, {steps: 6});
  await page.mouse.down();
  await page.mouse.up();
}, {tier: 'ui', description: 'a click on the empty canvas that many bonds under the drawing\'s atoms, at their middle'});

// ---------------------------------------------------------------- checks: the app's Text, step by step

/** In the page: the app's Text as steps, each [reactants, agents, products], every molecule canonical by RDKit and each
 * side sorted (a CX part left out). */
function stepsOf(page: Page, text?: string): Promise<string[][][]> {
  return page.evaluate(async (given) => {
    const t = given ?? (document.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement | null)?.value ?? '';
    const rdkit = await grok.functions.call('Chem:getRdKitModule');
    const canonical = (smiles: string) => {
      let mol: any = null;
      try {
        mol = rdkit.get_mol(smiles);
        return mol?.is_valid() ? mol.get_smiles() : `(unread) ${smiles}`;
      } catch {
        return `(unread) ${smiles}`;
      } finally {
        mol?.delete();
      }
    };
    return t.split(' |')[0].trim().split('--**--').map((step: string) => step.split('>')
      .map((side: string) => side === '' ? [] : side.split('.').map(canonical).sort()));
  }, text ?? null);
}

export const routeSteps = Then('the Reaction Sketcher\'s Text should be the route {string}', async (page: Page,
  route: string) => {
  const want = JSON.stringify(await stepsOf(page, route));
  await expect.poll(async () => JSON.stringify(await stepsOf(page)), {message: 'the app\'s Text, step by step',
    timeout: pollMs(10000)}).toBe(want);
}, {tier: 'api', description: 'the app\'s Text holds this route (steps joined by --**--): each step\'s reactants, agents ' +
  'and products, each molecule read by RDKit in any SMILES, the map numbers kept as written'});

export const textShown = Then('the sketcher should draw the text {string}', async (page: Page, words: string) => {
  await textOf(page, words);
}, {tier: 'ui', description: 'Crux draws a text of these words (subscript digits read as digits, spaces aside)'});

export const textSized = Then('the sketcher should draw the text {string} in {int} pt', async (page: Page,
  words: string, pt: number) => {
  const sel = await sketchSelector(page);
  const plain = (s: string) => s.replace(/[₀-₉]/g, (c) => String(c.charCodeAt(0) - 0x2080)).replace(/\s+/g, '');
  // the size as the drawing keeps it: its CDXML's text runs (ChemDraw's size, in points)
  await expect.poll(() => page.evaluate(([s, w]) => {
    const el = document.querySelector(s) as any;
    const doc = new DOMParser().parseFromString(String(el.cdxml), 'text/xml');
    const norm = (x: string) => x.replace(/[₀-₉]/g, (c) => String(c.charCodeAt(0) - 0x2080)).replace(/\s+/g, '');
    const t = Array.from(doc.getElementsByTagName('t')).find((n) => n.parentElement?.tagName !== 'n' &&
      norm(Array.from(n.getElementsByTagName('s')).map((r) => r.textContent ?? '').join('')) === w);
    return t ? Array.from(t.getElementsByTagName('s')).map((r) => Number(r.getAttribute('size'))) : null;
  }, [sel, plain(words)] as [string, string]), {message: `the text "${words}"'s size`, timeout: pollMs(5000)})
    .toEqual(expect.arrayContaining([pt]));
}, {tier: 'api', description: 'the text\'s runs are that size in the drawing\'s CDXML (ChemDraw\'s points)'});

/** The Reaction Sketcher's right pane: its scrolling box, and the heading of one of its sections. */
const PREVIEW = '.chem-reactions-app-drawing';

export const scrollPreview = When('user scrolls the Reaction Sketcher\'s preview to its {word}', async (page: Page,
  part: string) => {
  const heading = ({steps: 'Steps', mapping: 'Mapping'} as {[k: string]: string})[part];
  expect(heading, `a section of the preview: steps or mapping, not "${part}"`).toBeTruthy();
  const pane = page.locator(PREVIEW).first();
  const box = (await pane.boundingBox())!;
  const h3 = pane.locator('h3', {hasText: heading}).first();
  await expect(h3, `the preview's ${heading} section`).toBeAttached({timeout: pollMs(15000)});
  await page.mouse.move(box.x + box.width / 2, box.y + Math.min(box.height, 700) / 2, {steps: 6});
  // notch by notch, as a person turns the wheel, until the section's heading is in the pane's upper third
  for (let i = 0; i < 40; i++) {
    const top = await h3.evaluate((e) => e.getBoundingClientRect().top);
    const view = await pane.evaluate((e) => {
      const scroller = e.parentElement!;
      const r = scroller.getBoundingClientRect();
      return {top: r.top, bottom: r.bottom, end: scroller.scrollTop + scroller.clientHeight >= scroller.scrollHeight - 1};
    });
    if (top <= view.top + (view.bottom - view.top) / 3 || view.end)
      break;
    await page.mouse.wheel(0, 100);
    await page.waitForTimeout(60);
  }
}, {tier: 'ui', description: 'the mouse wheel turned over the right pane, notch by notch, until the section (Steps or ' +
  'Mapping) heads its upper third or the pane can scroll no further'});

/** In the page: the right pane's Mapping section: its state (reaction-panel.ts MappingState) and its lines. */
function previewMapping(page: Page): Promise<{state: string | null, lines: {level: string | null, text: string}[]}> {
  return page.evaluate((p) => {
    const status = document.querySelector(`${p} .chem-reaction-mapping`) as HTMLElement | null;
    return {state: status?.dataset.state ?? null, lines: status ? Array.from(status.children).map((c) => ({
      level: (c as HTMLElement).dataset.level ?? null, text: (c.textContent ?? '').trim()})) : []};
  }, PREVIEW);
}

async function mappedLine(page: Page, test: (l: {level: string | null, text: string}) => boolean, what: string) {
  let s: Awaited<ReturnType<typeof previewMapping>> | null = null;
  try {
    await expect.poll(async () => (s = await previewMapping(page)).lines.some(test),
      {message: what, timeout: pollMs(90000)}).toBe(true);
  } catch (e) {
    throw new Error(`${what}: the Mapping section is ${JSON.stringify(s)}\n${e instanceof Error ? e.message : e}`);
  }
}

export const previewSure = Then('the Reaction Sketcher\'s preview should say how sure the mapper is of step {int}',
  (page: Page, step: number) => mappedLine(page, (l) => new RegExp(`^Step ${step}:.*the mapper is \\d+% sure`).test(l.text),
    `step ${step} mapped, how sure the mapper is said`),
  {description: 'the right pane\'s Mapping section says step N was mapped and how sure the mapper is (polled while ' +
    'RXNMapper maps each step on a Crux worker)'});

export const previewWarns = Then('the Reaction Sketcher\'s preview should warn of step {int} {string}', (page: Page,
  step: number, words: string) => mappedLine(page, (l) => l.level === 'warn' && l.text.startsWith(`Step ${step}:`) &&
    l.text.includes(words), `step ${step}'s warning "${words}"`),
{description: 'the right pane\'s Mapping section has an orange line for step N holding these words'});

export const everyStepMapped = Then('every step of the Reaction Sketcher\'s Text should carry map numbers',
  async (page: Page) => {
    const read = () => page.evaluate(() => {
      const t = (document.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement | null)?.value ?? '';
      return t.split(' |')[0].split('--**--').map((s) => (s.match(/:\d+\]/g) ?? []).length);
    });
    await expect.poll(async () => (await read()).every((n) => n > 0), {message: 'map numbers in every step of the Text',
      timeout: pollMs(60000)}).toBe(true);
  }, {tier: 'api', description: 'each step of the app\'s Text (--**--) has atoms with map numbers (Auto-map answers a ' +
    'step at a time)'});

export const rightClickAtom = When('user right-clicks atom {int} of the sketcher', async (page: Page, i: number) => {
  const sel = await sketchSelector(page);
  const at = await page.evaluate(([s, n]) => {
    const el = document.querySelector(s as string) as any;
    const p = el.positions.atoms[n as number];
    const r = el.getBoundingClientRect();
    return p ? {x: p.x + r.x, y: p.y + r.y} : null;
  }, [sel, i] as [string, number]);
  expect(at, `atom ${i} drawn`).not.toBeNull();
  await page.mouse.move(at!.x, at!.y, {steps: 8});
  await page.mouse.down({button: 'right'});
  await page.mouse.up({button: 'right'});
  await expect(page.locator(`${sel} [data-testid="menu.context"]`), 'Crux\'s menu').toBeVisible({timeout: pollMs(5000)});
}, {tier: 'ui', description: 'a right-click on the atom (by the molfile\'s numbering): Crux\'s own menu of the atom'});

export const clickLastAtom = When('user clicks on the last atom of the sketcher', async (page: Page) => {
  const sel = await sketchSelector(page);
  const at = await page.evaluate((s) => {
    const el = document.querySelector(s) as any;
    const r = el.getBoundingClientRect();
    const ps = el.positions.atoms as ({x: number, y: number} | null)[];
    const p = [...ps].reverse().find((q) => q !== null);
    return p ? {x: p.x + r.x, y: p.y + r.y} : null;
  }, sel);
  expect(at, 'an atom drawn').not.toBeNull();
  await page.mouse.move(at!.x, at!.y, {steps: 8});
  await page.mouse.down();
  await page.mouse.up();
}, {tier: 'ui', description: 'a click on the atom Crux numbers last (the molecule added last)'});

export const resizeDialog = When('user resizes {element} by {int} and {int} pixels from its corner',
  async (page: Page, target: ElementRef, dx: number, dy: number) => {
    const dialog = (await locate(page, target)).filter({visible: true}).last();
    const b = (await dialog.boundingBox())!;
    const from = {x: b.x + b.width - 3, y: b.y + b.height - 3};
    await drag(page, from, {x: from.x + dx, y: from.y + dy}, 16);
    await page.waitForTimeout(300);
  }, {tier: 'ui', description: 'the dialog dragged by its bottom-right corner, in steps'});

export const paneSays = Then('the Reaction pane should say {string}', async (page: Page, words: string) => {
  const pane = page.locator('.grok-prop-panel .d4-accordion-pane').filter({has: page.locator('.d4-accordion-pane-header',
    {hasText: /^Reaction$/})}).first();
  await expect(pane, `the Reaction pane saying "${words}"`).toContainText(words, {timeout: pollMs(15000)});
}, {description: 'the context panel\'s Reaction pane holds these words'});

export const pointAtLabel = When('user points at the label {string} in the sketcher', async (page: Page,
  label: string) => {
  const sel = await sketchSelector(page);
  const found: {at: Point | null} = {at: null};
  await expect.poll(async () => {
    found.at = await page.evaluate(([s, l]) => {
      const e = document.querySelector(s)!.shadowRoot!.querySelector(`[data-label="${l}"]`);
      const r = e?.getBoundingClientRect();
      return r && r.width > 0 ? {x: r.x + r.width / 2, y: r.y + r.height / 2} : null;
    }, [sel, label] as [string, string]);
    return found.at !== null;
  }, {message: `the label ${label} drawn`, timeout: pollMs(10000)}).toBe(true);
  await page.mouse.move(found.at!.x, found.at!.y, {steps: 10});
  await page.waitForTimeout(600);
}, {tier: 'ui', description: 'the pointer moved onto the atom Crux draws with this label, and rested there'});

export const labelsShown = Then('the sketcher should draw the labels {string} and {string}', async (page: Page,
  a: string, b: string) => {
  const sel = await sketchSelector(page);
  await expect.poll(() => page.evaluate(([s, want]) => {
    const labels = Array.from(document.querySelector(s as string)!.shadowRoot!.querySelectorAll('[data-label]'))
      .map((e) => e.getAttribute('data-label'));
    return (want as string[]).filter((w) => labels.includes(w));
  }, [sel, [a, b]] as [string, string[]]), {message: `the labels ${a} and ${b} drawn`, timeout: pollMs(10000)})
    .toEqual([a, b]);
}, {tier: 'ui', description: 'Crux draws atoms labelled so (their data-label: a contracted label kept whole)'});
