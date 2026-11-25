/**
 * Category `Reactions: app`: the Reaction Sketcher (utils/reactions/reactions-app.ts) on a real view, in its layout:
 * ONE reaction, the editing area on the left, the reaction drawn and told of on the right.
 * It opens cold on the default example (a five-step mapped route), loads an example from the ribbon's Examples menu and
 * keeps it, keeps an edit in the browser and reopens with it, resets to the default example, and exports the reaction
 * as SVG and PNG. The document this browser keeps is set aside before the category and put back after it, so a run
 * leaves a person's document as it was.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {ChemWasm, armCruxTrapForTests, cruxLoadFailed, cruxModuleSync, loadCrux} from '../engine/crux-module';
import {REAGENT_HINT_NAMES, agentNames} from '../utils/reactions/reaction-add-panel';
import {REFUSED_MESSAGE} from '../crux/crux-recover';
import {DOCUMENT_KEY, DocumentReaction, FRAGMENT_GROUPS_UNREAD, NO_ARROW, PAGE_PADDING, ReactionPicture, TOO_LARGE,
  defaultExample,
  documentPng, documentSvg, examples, forgetDocument, loadDocument, sameReactionText, saveDocument}
  from '../utils/reactions/reaction-document';
import {DRAWING_LOSS_WORDS} from '../utils/reactions/reaction-losses';
import {ADD_TO_TABLE_BUTTON, DRAWING_UNREAD, EDITOR_FAILED, EDITOR_MIN_WIDTH, EDITOR_SHARE, EXAMPLES_MENU,
  EXPORT_BUTTON, LOSSES_CLASS, NEW_TABLE_NAME, NOTHING_TO_ADD, NO_ARROW_CLASS, PAGE_READ_MS, PREVIEW_MIN_WIDTH,
  REACTION_COLUMN, SKETCHER_FULL_HEIGHT, SKETCHER_REFUSED_CLASS, SVG_FILE, TITLE_COLUMN, UNTITLED, reactionUnitsOf,
  readsAs, sketcherBoxHeight} from '../utils/reactions/reactions-app';
import {ACS_1996, CDXML, CDXML_NOT_YET, CDXML_WRITTEN, EXPORT_TITLE, NOTE_CLASS, PREVIEW_CLASS, SINGLE_COLUMN,
  TEXTS_CDXML_ONLY, drawnOptions, exportInputsOf} from '../widgets/reaction-export';
import {cdxmlReactionSmiles, sketcherPicture} from '../rendering/cdxml-reaction-cells';
import {DrawnStepTexts, cruxDrawsStepTexts, stepTextsOf} from '../rendering/crux-reaction-drawing';
import {CDXML_DEMO_FILE, EDITED_ROW, captions} from './reaction-cdxml-tests';
import {REACTIONS_HELP_FILE, reactionsHelpUrl} from '../utils/reactions/reactions-help';
import {SET_VALUE_ARRANGES, arrangesOnWrite} from '../utils/reactions/reaction-arrange';
import {arrowRows} from '../utils/reactions/reaction-cell-editor';
import {PackageFunctions} from '../package';
import {loadFileAsText} from './utils';
import {ReactionSketchArea, SKETCH_LOADING, SKETCH_LOADING_CLASS} from '../utils/reactions/reaction-cell-editor';
import {CruxSketchLoaderForTests, setCruxSketchLoaderForTests} from '../crux/crux-sketcher';

/** The width the export tests lay the reaction out at (CSS px). */
const WIDTH = 900;
/** Files > App Data > Chem, where `grok publish` puts Chem's files/ (what the help page tells a person to open). */
const APP_DATA_CHEM = 'System:AppData/Chem';

/** Pixels with any ink (alpha above 0) in the canvas, or with a colour other than white (`onWhite`). */
function ink(canvas: HTMLCanvasElement, onWhite = false): number {
  const data = canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height).data;
  let n = 0;
  for (let i = 0; i < data.length; i += 4) {
    if (onWhite ? data[i] < 250 || data[i + 1] < 250 || data[i + 2] < 250 : data[i + 3] > 0)
      n++;
  }
  return n;
}

/** The right pane: the reaction drawn and told of. */
const drawing = (view: DG.View) => view.root.querySelector('.chem-reactions-app-drawing') as HTMLElement;
/** The large picture of the right pane (the Reaction pane's first canvas), or null. */
const picture = (view: DG.View) => drawing(view)?.querySelector('canvas') as HTMLCanvasElement | null;
/** The left pane's editing area: the text box. */
const textBox = (view: DG.View) => view.root.querySelector('.chem-reactions-app-text') as HTMLTextAreaElement;
/** The Crux sketcher of the editing area (a crux-sketch element), or null. */
const sketchEl = (view: DG.View) => view.root.querySelector('.chem-reactions-app-sketcher crux-sketch') as any;
/** Crux for the tests that compare reactions (loaded by {@link keepingTheDocument}). */
let crux: ChemWasm | null = null;
/** A drawn bond's least length, px, for a reaction opened arranged (readable, not 20% in a row). */
// with the context and help panels shut by the app: 15 px measured at 1440x900, so the floor is measured - 1 (to be 15
// once the measure clears 16.5 with margin, e.g. after a box-shaped Arrange)
const MIN_ARRANGED_BOND_PX = 14;
/** Whether the editing area's sketcher has read the document (its reaction document, its value the document's text:
 * opened arranged by an Arrange edit, the sketcher's text of the same reaction). */
const sketcherShows = (view: DG.View, text: string) => sketchEl(view)?.document === 'reaction' && (() => {
  try {
    const value = sketchEl(view).value;
    return value === text || (!SET_VALUE_ARRANGES && arrangesOnWrite(text) && crux !== null &&
      sameReactionText(crux, value, text));
  } catch {
    return false;
  }
})();
/** The title the right pane shows. */
const shownTitle = (view: DG.View) => drawing(view)?.querySelector('.chem-reactions-app-title')?.textContent ?? '';
/** What the status line under the editing area says. */
const status = (view: DG.View) => view.root.querySelector('.chem-reactions-app-status')!.textContent ?? '';

/** Whether Crux Sketch's left bar in `el` shows every group, the reaction arrow tool inside it, with no scroll button
 * (crux-sketch's `toolbar.left.scroll-forward`, shown while more of the bar lies below). */
function leftBarWhole(el: HTMLElement | null): boolean {
  const root = el?.shadowRoot;
  const bar = root?.querySelector('[data-testid="toolbar.left"]') as HTMLElement | null;
  const arrow = root?.querySelector('[data-testid="toolbar.reaction.arrows"]') as HTMLElement | null;
  if (!bar || !arrow)
    return false;
  const forward = root!.querySelector('[data-testid="toolbar.left.scroll-forward"]') as HTMLElement | null;
  const scrolls = !!forward && forward.getBoundingClientRect().height > 0 &&
    getComputedStyle(forward).display !== 'none' && getComputedStyle(forward).visibility !== 'hidden';
  const b = bar.getBoundingClientRect();
  const a = arrow.getBoundingClientRect();
  return !scrolls && bar.scrollHeight <= bar.clientHeight + 1 && a.height > 0 && a.top >= b.top - 0.5 &&
    a.bottom <= b.bottom + 0.5;
}

/** The app opened on the shell, once its reaction is drawn. */
async function openApp(): Promise<DG.View> {
  const view = await PackageFunctions.reactionsApp() as DG.View;
  grok.shell.addView(view);
  await awaitCheck(() => !!picture(view) && ink(picture(view)!) > 500, 'the app\'s reaction was not drawn', 8000);
  return view;
}

/** Types `text` into the editing area, as a person's keystrokes do. */
function type(view: DG.View, text: string): void {
  const box = textBox(view);
  box.value = text;
  box.dispatchEvent(new Event('input', {bubbles: true}));
}

/** The ribbon's button `caption` (the ribbon holds each in a box of its own). */
function ribbonButton(view: DG.View, caption: string): HTMLElement {
  const item = view.getRibbonPanels().flat().find((b) => b.textContent === caption)!;
  return item.tagName === 'BUTTON' ? item : item.querySelector('button')!;
}

/** The documents as pictures read by Crux. */
async function picturesOf(docs: DocumentReaction[]): Promise<ReactionPicture[]> {
  const crux: ChemWasm = await loadCrux('tests');
  return docs.map((r) => new ReactionPicture(crux, r));
}

/** The document the browser keeps, given back after the category's tests (they open the app on documents of their
 * own). */
function keepingTheDocument(): void {
  let kept: string | null = null;
  let putPanelsBack: (() => void) | null = null;

  before(async () => {
    kept = window.localStorage.getItem(DOCUMENT_KEY);
    putPanelsBack = shellPanelsKept();
    crux ??= await loadCrux('tests');
  });

  after(async () => {
    if (kept === null)
      window.localStorage.removeItem(DOCUMENT_KEY);
    else
      window.localStorage.setItem(DOCUMENT_KEY, kept);
    putPanelsBack?.();
  });
}

/** The shell's context and help panels as they are now, and what puts them back so (the app shuts both as it opens,
 * {@link shutSidePanels}: a category that opens it leaves the shell as it found it). */
export function shellPanelsKept(): () => void {
  const windows = grok.shell.windows;
  const [context, help] = [windows.showContextPanel, windows.showHelp];
  return () => {
    windows.showContextPanel = context;
    windows.showHelp = help;
  };
}

/** The app opens with the context and help panels shut, as other Datagrok apps do, its sketcher given the room; a
 * person may open them again (at 1440x900 with both open, the sketcher was 359 px wide and the opened route's bonds
 * 8 px). */
async function panelsShut(): Promise<void> {
  forgetDocument();
  grok.shell.windows.showContextPanel = true;
  grok.shell.windows.showHelp = true;
  const view = await openApp();
  try {
    expect(grok.shell.windows.showContextPanel, false, 'the context panel shut');
    expect(grok.shell.windows.showHelp, false, 'the help panel shut');
    // the room was the context panel's: opened again by a person, it takes it back, and the app leaves it open
    const split = (view.root.querySelector('.chem-reactions-app-editor') as HTMLElement).parentElement!.parentElement!
      .parentElement!;
    const shut = split.getBoundingClientRect().width;
    grok.shell.windows.showContextPanel = true;
    await awaitCheck(() => split.getBoundingClientRect().width < shut - 100,
      `the view as wide with the context panel open: ${shut} px`, 3000);
    await delay(300);
    expect(grok.shell.windows.showContextPanel, true, 'the context panel opened again stays open');
  } finally {
    view.close();
  }
}

/** The box a reaction sketcher is made in ({@link ReactionSketchArea}) is under a cover until the sketcher has read the
 * reaction, so that a click then is the cover's; uncovered, a click is the box's. */
async function boxCovered(): Promise<void> {
  const box = document.createElement('div');
  Object.assign(box.style, {position: 'fixed', left: '20px', top: '20px', width: '400px', height: '300px',
    zIndex: '100000', background: 'white'});
  document.body.append(box);
  const area = new ReactionSketchArea(box);
  try {
    const cover = box.querySelector(`:scope > .${SKETCH_LOADING_CLASS}`) as HTMLElement | null;
    expect(cover !== null, true, 'a cover in the box');
    expect(cover!.textContent, SKETCH_LOADING);
    expect(area.loading, true, 'loading as it is made');
    const hit = () => document.elementFromPoint(220, 40); // where the sketcher's top bar is drawn
    expect(hit() === cover, true, `a click as the box opens reaches ${hit()?.className}`);
    area.loading = false;
    expect(getComputedStyle(cover!).display, 'none', 'uncovered once the sketcher has read the reaction');
    expect(hit() === box, true, `a click then reaches ${hit()?.className}`);
  } finally {
    area.detach();
    box.remove();
  }
}

/** Crux Sketch, for the tests' loader, as a module whose sketchers do not start (createSketcher rejects). */
const NOT_STARTING = {
  createSketcher: async () => {
    throw new Error('test: the sketcher does not start');
  },
} as unknown as Awaited<ReturnType<CruxSketchLoaderForTests>>;

/** The app opened while Crux Sketch does not start: the box uncovered and saying so, a click on it reaching it (not
 * the cover), the Text still the document's. */
async function notStarted(): Promise<void> {
  forgetDocument();
  setCruxSketchLoaderForTests(async () => NOT_STARTING);
  let view: DG.View | null = null;
  try {
    view = await openApp();
    const host = view.root.querySelector('.chem-reactions-app-sketcher') as HTMLElement;
    await awaitCheck(() => (host.textContent ?? '').includes(EDITOR_FAILED), 'the box did not say it failed', 10000);
    const cover = host.querySelector(`:scope > .${SKETCH_LOADING_CLASS}`) as HTMLElement;
    expect(getComputedStyle(cover).display, 'none', 'the cover gone');
    const r = host.getBoundingClientRect();
    const hit = document.elementFromPoint(r.left + r.width / 2, r.top + r.height / 2);
    expect(hit !== cover && host.contains(hit), true, `a click on the box reaches ${hit?.className}`);
    expect(textBox(view).value, defaultExample().text, 'the Text');
    expect(status(view), readsAs(5), 'the status');
  } finally {
    setCruxSketchLoaderForTests(null);
    view?.close();
  }
}

/** A document reopened on its kept drawing is read under the sketcher's cover, which lifted before the drawing was read
 * (a click on the empty sketcher then became the document). The cover stays at most PAGE_READ_MS: past it (here a CDXML
 * module that never loads), or when the drawing does not read (an empty page), the text is drawn and the status says
 * the kept drawing could not be read. */
async function keptDrawingCovered(): Promise<void> {
  const text = 'CC(=O)O.OCC>>CC(=O)OCC';
  const cover = (view: DG.View) => view.root.querySelector(
    `.chem-reactions-app-sketcher > .${SKETCH_LOADING_CLASS}`) as HTMLElement | null;
  const covered = (view: DG.View) => !!cover(view) && getComputedStyle(cover(view)!).display !== 'none';
  // a drawing that does not read: the text drawn, said in the status
  saveDocument({title: '', text, mapNumbers: false, drawing: '<CDXML><page id="1"/></CDXML>'});
  let view = await openApp();
  try {
    await awaitCheck(() => sketcherShows(view, text) && !covered(view), 'the text was not drawn', 15000);
    await awaitCheck(() => status(view).includes(DRAWING_UNREAD), `the status: ${status(view)}`, 3000);
    expect(textBox(view).value, text, 'the Text');
  } finally {
    view.close();
  }
  // a drawing whose CDXML module never loads: covered meanwhile, then the text drawn, said in the status
  const page = (DG.DataFrame.fromCsv(await loadFileAsText(CDXML_DEMO_FILE)).col('cdxml')!
    .toList() as string[])[EDITED_ROW];
  setCruxSketchLoaderForTests(async (load) => {
    const module = await load();
    const createSketcher: typeof module.createSketcher = async (...args) => {
      const sketch = await module.createSketcher(...args);
      return new Proxy(sketch, {get: (target, key) => {
        if (key === 'loadFormat')
          return () => new Promise<void>(() => {});
        const v = (target as any)[key];
        return typeof v === 'function' ? v.bind(target) : v;
      }});
    };
    return {...module, createSketcher} as typeof module;
  });
  saveDocument({title: '', text, mapNumbers: false, drawing: page});
  view = await openApp();
  try {
    const started = Date.now();
    await awaitCheck(() => !!sketchEl(view), 'the sketcher did not start', 15000);
    await delay(1000);
    expect(covered(view), true, 'covered while the kept drawing is read');
    await awaitCheck(() => !covered(view), `the cover still up past ${PAGE_READ_MS} ms`, PAGE_READ_MS + 5000);
    expect(Date.now() - started >= PAGE_READ_MS - 1500, true, `the cover kept ${Date.now() - started} ms`);
    await awaitCheck(() => sketcherShows(view, text), 'the text was not drawn', 5000);
    await awaitCheck(() => status(view).includes(DRAWING_UNREAD), `the status: ${status(view)}`, 3000);
  } finally {
    setCruxSketchLoaderForTests(null);
    view.close();
  }
}

/** The default route opened: arranged in rows (reaction-arrange.ts), its bonds long enough to read,
 * the Text as written, and no edit to undo when the write itself arranges. */
async function opensArranged(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  try {
    const text = defaultExample().text;
    await awaitCheck(() => sketcherShows(view, text), 'the sketcher did not show the document', 15000);
    const el = sketchEl(view);
    expect(typeof el.arrange, 'function', 'the vendored sketcher arranges');
    await awaitCheck(() => arrowRows(el, el.positions.bondLength) >= 2,
      `the route in ${arrowRows(el, el.positions.bondLength)} row(s)`, 5000);
    expect(el.positions.bondLength >= MIN_ARRANGED_BOND_PX, true,
      `bonds drawn ${el.positions.bondLength.toFixed(1)} px long`);
    expect(textBox(view).value, text, 'the Text as the example writes it');
    expect(status(view), readsAs(5));
    if (SET_VALUE_ARRANGES)
      expect(el.canUndo, false, 'opening arranged is no edit to undo');
  } finally {
    view.close();
  }
}

/** A text typed is one edit of the drawing (setValue's `history: 'edit'`): Undo brings the document back. */
async function textTypedUndone(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  const text = 'CC(=O)O.OCC>>CC(=O)OCC';
  try {
    const before = defaultExample().text;
    await awaitCheck(() => sketcherShows(view, before), 'the sketcher did not show the document', 15000);
    const el = sketchEl(view);
    type(view, text);
    await awaitCheck(() => sketcherShows(view, text), 'the sketcher did not take the text', 5000);
    expect(el.canUndo, true, 'the text typed is an edit to undo');
    expect(el.undo(), true, 'Undo');
    await awaitCheck(() => sketcherShows(view, before), 'Undo did not bring the document back', 5000);
  } finally {
    view.close();
  }
}

/** A route typed (drawn as read, not arranged), then laid out by a person's Arrange from the Clean Up and Layout menu:
 * the Text stays as typed. */
async function layoutOnly(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  try {
    await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document', 15000);
    const typed = examples()[0].text; // the paracetamol route, three steps
    type(view, typed);
    textBox(view).dispatchEvent(new Event('blur')); // drawn now, not after the pause
    await awaitCheck(() => sketcherShows(view, typed), 'the sketcher did not take the text', 5000);
    const host = view.root.querySelector('.chem-reactions-app-sketcher') as HTMLElement;
    const changes = () => Number(DG.Widget.find(host)?.getWidgetStatus()?.values?.changes ?? -1);
    const before = changes();
    const root = sketchEl(view).shadowRoot as ShadowRoot;
    (root.querySelector('[data-testid="toolbar.clean.open"]') as HTMLElement).click();
    const arrange = () => root.querySelector('[data-testid="menu.clean-up.arrange"]') as HTMLElement | null;
    await awaitCheck(() => arrange() !== null, 'no Arrange in the menu', 3000);
    arrange()!.click();
    await awaitCheck(() => changes() === before + 1, `Arrange made no edit (changes ${before} -> ${changes()})`, 5000);
    expect(textBox(view).value, typed, 'the Text as typed');
    expect(status(view), readsAs(3));
  } finally {
    view.close();
  }
}

/** The editing pane (the split's box round the editor) is {@link EDITOR_SHARE}% of the view, at least
 * {@link EDITOR_MIN_WIDTH} px, and the sketcher's box fills its width. */
async function paneShare(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  try {
    const box = (view.root.querySelector('.chem-reactions-app-editor') as HTMLElement).parentElement!.parentElement!;
    const split = box.parentElement!;
    const preview = split.lastElementChild as HTMLElement;
    const widths = () => ({split: split.getBoundingClientRect().width, editor: box.getBoundingClientRect().width,
      preview: preview.getBoundingClientRect().width});
    const shared = (w = widths()) => Math.abs(w.editor - Math.min(Math.max(EDITOR_MIN_WIDTH,
      w.split * EDITOR_SHARE / 100), w.split - PREVIEW_MIN_WIDTH - 4)) <= 2 && w.preview >= PREVIEW_MIN_WIDTH - 1;
    expect(shared(), true, `the panes: ${JSON.stringify(widths())}`);
    const sketcher = view.root.querySelector('.chem-reactions-app-sketcher') as HTMLElement;
    expect(sketcher.getBoundingClientRect().width > widths().editor - 40, true,
      `the sketcher ${sketcher.getBoundingClientRect().width} px in a pane of ${widths().editor}`);
    // the view narrowed as the help pane narrows it (where the preview was once ~75 px): the editing area gives way
    split.style.width = '560px';
    await awaitCheck(() => shared(), `narrowed: ${JSON.stringify(widths())}`, 3000);
    split.style.width = '100%';
    await awaitCheck(() => shared(), `widened again: ${JSON.stringify(widths())}`, 3000);
    // a drag on the divider: the person's widths stay, the preview's least width still kept
    const divider = split.querySelector('.ui-split-h-divider') as HTMLElement;
    divider.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
    document.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
    split.style.width = '500px';
    await delay(500);
    expect(widths().preview >= PREVIEW_MIN_WIDTH - 1, true, `after a drag, narrowed: ${JSON.stringify(widths())}`);
  } finally {
    view.close();
  }
}

/** The Text emptied, a name refused by the Reagent slot, then CCO and CC=O added as reactant and product: the Text,
 * the status and the preview have the reaction; no refusal line left. */
async function addedOnly(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  try {
    await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document', 15000);
    type(view, '');
    textBox(view).dispatchEvent(new Event('blur'));
    await awaitCheck(() => sketchEl(view).isEmpty === true, 'the drawing was not emptied', 5000);
    const slot = (role: string) =>
      view.root.querySelector(`.chem-reactions-add-slot[data-role="${role}"] input`) as HTMLInputElement;
    const enter = (role: string, text: string) => {
      slot(role).value = text;
      slot(role).dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true}));
    };
    const line = view.root.querySelector('.chem-reactions-add-refused') as HTMLElement;
    enter('agent', 'Xylophone'); // a name no agent table holds (not HATU, which the table holds)
    expect(line.style.display !== 'none', true, 'the name refused');
    enter('reactant', 'CCO');
    enter('product', 'CC=O');
    expect(slot('reactant').value + slot('product').value, '', 'the slots emptied once added');
    expect(line.style.display, 'none', 'the refusal line cleared by the next addition');
    const crux = await loadCrux('tests');
    await awaitCheck(() => sameReactionText(crux, textBox(view).value, 'CCO>>CC=O'),
      `the Text: "${textBox(view).value}"`, 5000);
    expect(status(view), readsAs(1));
    await awaitCheck(() => !!picture(view) && ink(picture(view)!) > 200, 'the preview did not draw it', 5000);
  } finally {
    view.close();
  }
}

/** The default example's title dropped once its reaction is replaced (the Title box empty, "Untitled" in it, the right
 * pane untitled, the browser's copy too); a layout-only change keeps it; a title typed stays through later edits. */
async function titleFollows(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  const title = () => view.root.querySelector('.chem-reactions-app-editor input[type="text"]') as HTMLInputElement;
  try {
    expect(title().value, defaultExample().title, 'the example\'s title as it opens');
    type(view, defaultExample().text + ' '); // the same reaction: the title stays
    expect(title().value, defaultExample().title, 'the same reaction');
    type(view, 'CC(=O)O.OCC>>CC(=O)OCC');
    expect(title().value, '', 'the example\'s title dropped');
    expect(title().placeholder, UNTITLED);
    expect(loadDocument()?.title, '', 'kept so');
    await awaitCheck(() => shownTitle(view) === '', `the right pane's title: ${shownTitle(view)}`, 3000);
    title().value = 'My ester';
    title().dispatchEvent(new Event('input', {bubbles: true}));
    await awaitCheck(() => loadDocument()?.title === 'My ester', 'the title typed was not kept', 3000);
    type(view, 'CC(=O)O.OCCC>>CC(=O)OCCC');
    expect(title().value, 'My ester', 'a title of one\'s own stays');
    expect(loadDocument()?.title, 'My ester');
  } finally {
    view.close();
  }
}

/** A step added by the Add panel's Next step wiped the example's title to "Untitled" (so did the drawing pasted back as
 * CDXML). The example's title stays through edits of its reaction: a step added, a step deleted; it goes when the Text
 * is emptied (Clear), comes back with the reaction (Clear undone), and goes for a reaction with none of the example's
 * reactants and products. */
async function titleKept(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  const example = defaultExample();
  const title = () => view.root.querySelector('.chem-reactions-app-editor input[type="text"]') as HTMLInputElement;
  const stepsOf = (text: string) => {
    const scheme = new crux!.ChemReactionScheme(text, '');
    try {
      return scheme.get_steps();
    } finally {
      scheme.free();
    }
  };
  const kept = async (what: string) => {
    expect(title().value, example.title, `${what}: the Title box`);
    expect(loadDocument()?.title, example.title, `${what}: the document kept`);
    await awaitCheck(() => shownTitle(view) === example.title, `${what}: the right pane's title: ${shownTitle(view)}`,
      3000);
  };
  try {
    await awaitCheck(() => sketchEl(view)?.document === 'reaction' && sketchEl(view).isEmpty === false,
      'the sketcher did not show the example', 15000);
    const next = view.root.querySelector('.chem-reactions-add-slot[data-role="next"] input') as HTMLInputElement;
    next.value = 'CC(=O)O';
    next.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true, cancelable: true}));
    await awaitCheck(() => stepsOf(textBox(view).value).length === stepsOf(example.text).length + 1,
      `a step added: "${textBox(view).value.slice(0, 80)}..."`, 5000);
    await kept('a step added (Next step)');
    type(view, stepsOf(example.text).slice(1).join('--**--'));
    await kept('a step deleted');
    type(view, '');
    expect(title().value, '', 'the Text emptied: the title dropped');
    expect(title().placeholder, UNTITLED);
    await awaitCheck(() => shownTitle(view) === '', `emptied: the right pane's title: ${shownTitle(view)}`, 3000);
    type(view, example.text); // as Ctrl+Z brings it back
    await kept('the reaction back');
    type(view, 'CC(=O)O.OCC>>CC(=O)OCC'); // acetic acid is an agent of the example's, no reactant or product
    expect(title().value, '', 'a reaction with none of the example\'s molecules: the title dropped');
  } finally {
    view.close();
  }
}

/** The right pane's pictures (the reaction's, then each step's), drawn. */
const pictures = (view: DG.View) => [...drawing(view).querySelectorAll('canvas')] as HTMLCanvasElement[];
const pixelsOf = (c: HTMLCanvasElement) => c.getContext('2d')!.getImageData(0, 0, c.width, c.height).data;
const samePixels = (a: Uint8ClampedArray, b: Uint8ClampedArray) =>
  a.length === b.length && a.every((v, i) => v === b[i]);
/** The Map numbers switch of the app's preview (in its row of controls at its top). */
const mapSwitch = (view: DG.View) =>
  view.root.querySelector('.chem-reactions-app-preview-controls input[type="checkbox"]') as HTMLInputElement;

/** The Map numbers switch, on the default route and on a reaction the right pane draws at its own bond length. On the
 * route every picture is drawn again with the numbers; whether that adds ink there depends on the pane's width (the
 * route's numbers widen it, so a 430-px picture drew it smaller, 5796 -> 4997 pixels, and a step 3132 -> 2298): the ink
 * of the numbers is read where the reaction is not made smaller to fit (the Diels-Alder example, 2667 -> 2982 at
 * 430 px). */
async function mapNumbersDrawn(): Promise<void> {
  forgetDocument();
  let view = await openApp();
  try {
    const route = pictures(view);
    expect(route.length, 6, 'the route\'s picture and its five steps\'');
    await awaitCheck(() => route.every((c) => ink(c) > 0), 'the steps were not drawn', 5000);
    const before = route.map(pixelsOf);
    mapSwitch(view).click();
    // the input says its change in a task of its own
    await awaitCheck(() => loadDocument()?.mapNumbers === true, 'the switch was not kept', 3000);
    await awaitCheck(() => pictures(view).length === route.length && pictures(view).every((c, i) => c !== route[i] &&
      ink(c) > 0), 'the pictures were not made again', 5000);
    pictures(view).forEach((c, i) =>
      expect(samePixels(pixelsOf(c), before[i]), false, `picture ${i}: drawn as before, no numbers`));
  } finally {
    view.close();
  }
  const example = examples().find((e) => e.title.startsWith('Atom-mapped: a Diels-Alder'))!;
  saveDocument({...example, mapNumbers: false});
  view = await openApp();
  try {
    const hidden = picture(view)!;
    const before = ink(hidden);
    mapSwitch(view).click();
    await awaitCheck(() => loadDocument()?.mapNumbers === true, 'the switch was not kept', 3000);
    await awaitCheck(() => picture(view) !== hidden && ink(picture(view)!) > 0, 'the drawing was not made again',
      3000);
    expect(ink(picture(view)!) > before, true, `the map numbers drawn: ${before} -> ${ink(picture(view)!)} pixels`);
  } finally {
    view.close();
  }
}

/** The demo page {@link EDITED_ROW} (texts over its arrows) drawn in the app's sketcher, its reaction in the Text: the
 * Export panel says the texts are in CDXML only, and its CDXML is the drawing's, every caption of the page kept. */
async function exportsTexts(): Promise<void> {
  forgetDocument();
  // read as Chem's tests read files (package-test's package: package.ts's is never given the platform in this bundle)
  const page = (DG.DataFrame.fromCsv(await loadFileAsText(CDXML_DEMO_FILE)).col('cdxml')!
    .toList() as string[])[EDITED_ROW];
  const smiles = await cdxmlReactionSmiles(page);
  expect(!!smiles, true, 'the demo page read');
  const view = await openApp();
  const download = DG.Utils.download;
  const got: [string, BlobPart][] = [];
  DG.Utils.download = (name: string, content: BlobPart) => {
    got.push([name, content]);
  };
  const panel = () => DG.Dialog.getOpenDialogs().find((d) => d.title === EXPORT_TITLE);
  const note = () => panel()?.root.querySelector(`.${NOTE_CLASS}`)?.textContent ?? '';
  const hasTexts = () => (sketchEl(view).warnings as {kind: string, code?: string}[])
    .some((w) => w.kind === 'format' && w.code === 'texts');
  try {
    await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document', 15000);
    type(view, smiles!);
    textBox(view).dispatchEvent(new Event('blur'));
    await awaitCheck(() => [1, 2, 3, 4].some((n) => status(view) === readsAs(n)), `the status: ${status(view)}`, 5000);
    // the page as the sketcher reads it, its texts over the arrow (what a person draws or pastes); the Text holds its
    // reaction already
    sketchEl(view).setValue(page);
    await awaitCheck(() => !sketchEl(view).isPending && hasTexts(), 'the page\'s texts were not drawn', 15000);
    // the texts over the arrows are the drawing's step texts (crux-sketch's reactionSteps),
    // which the pictures draw: the panel no longer says they are CDXML's alone
    const steps = sketchEl(view).reactionSteps as DrawnStepTexts[] | null | undefined;
    const texts = stepTextsOf(steps);
    expect(texts !== null, true, `the drawing's step texts: ${JSON.stringify(steps ?? null).slice(0, 200)}`);
    ribbonButton(view, EXPORT_BUTTON).click();
    await awaitCheck(() => !!panel(), 'the Export panel opens', 5000);
    expect(cruxDrawsStepTexts(), true, 'Chem\'s Crux file draws a host\'s step texts');
    expect(note().includes(TEXTS_CDXML_ONLY), false, `drawn, not said to be CDXML's alone: ${note()}`);
    // the right pane draws them too
    await awaitCheck(() => Number(picture(view)?.dataset.stepTexts ?? 0) > 0,
      `the right pane drew no step texts (${picture(view)?.dataset.stepTexts})`, 5000);
    panel()!.input('Format').value = CDXML;
    await awaitCheck(() => note().includes(CDXML_WRITTEN) || note().includes(CDXML_NOT_YET), note(), 15000);
    panel()!.getButton('Save').click();
    await awaitCheck(() => got.length === 1, `one download: ${note()}`, 5000);
    const written = String(got[0][1]);
    expect(captions(page).length > 0, true, 'the demo page has texts');
    const lost = captions(page).filter((t) => !captions(written).includes(t));
    expect(lost.length, 0, `texts not kept: ${lost.join(' | ')}`);
  } finally {
    DG.Utils.download = download;
    panel()?.close();
    view.close();
  }
}

/** Chem's Crux file has the agent table: the Reagent slot's hint names names, and HATU typed
 * there adds its structure over the arrow, the Text getting its SMILES; the slot emptied, no refusal said. */
async function namesAdded(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  try {
    const names = agentNames((cruxModuleSync() as unknown as {chemAgentTable?: () => string} | null)
      ?.chemAgentTable?.() ?? null);
    const hatu = names.smiles.get('hatu');
    expect(!!hatu, true, `Chem's Crux file has the agent table, HATU in it (${names.names.length} names)`);
    await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document', 15000);
    type(view, '');
    textBox(view).dispatchEvent(new Event('blur'));
    await awaitCheck(() => sketchEl(view).isEmpty === true, 'the drawing was not emptied', 5000);
    const slot = (role: string) =>
      view.root.querySelector(`.chem-reactions-add-slot[data-role="${role}"] input`) as HTMLInputElement;
    const enter = (role: string, text: string) => {
      slot(role).value = text;
      slot(role).dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true}));
    };
    slot('agent').dispatchEvent(new Event('focus'));
    expect(slot('agent').placeholder, REAGENT_HINT_NAMES);
    enter('reactant', 'CC(=O)O');
    enter('reactant', 'NCC');
    enter('agent', 'HATU');
    enter('product', 'CC(=O)NCC');
    const line = view.root.querySelector('.chem-reactions-add-refused') as HTMLElement;
    expect(line.style.display, 'none', `nothing refused: ${line.textContent}`);
    expect(slot('agent').value, '', 'the Reagent slot emptied once added');
    const crux = await loadCrux('tests');
    await awaitCheck(() => sameReactionText(crux, textBox(view).value, `CC(=O)O.NCC>${hatu}>CC(=O)NCC`),
      `the Text: "${textBox(view).value}"`, 5000);
    expect(textBox(view).value.includes('*'), false, 'HATU as its structure, never a "*"');
  } finally {
    view.close();
  }
}

/** A text Crux cannot read: a message under the editing area, gone once the text reads. */
async function unreadText(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  try {
    type(view, 'CC(=O)O.OCC>>CC(=O)OCC(');
    expect(status(view).startsWith('One reaction') || status(view).startsWith('A route'), false, status(view));
    expect(status(view).length > 0, true, 'a message');
    type(view, 'CC(=O)O.OCC>>CC(=O)OCC');
    expect(status(view), 'One reaction.');
  } finally {
    view.close();
  }
}

/** A route Crux Sketch does not open (step 2 has no product, so step 3 does not start from it: the text a Next step bug
 * once left): "a route that is neither a chain of steps nor a tree of chains opens in a later version". The editor's
 * refusal test uses it too: the vendored crux-sketch opens convergent routes. */
export const REFUSED_ROUTE ='c1ccccc1>>Oc1ccccc1O--**--Oc1ccccc1O>>--**--C1COc2ccccc2O1>>O=[N+]([O-])c1ccc2OCCOc2c1';

/** The app reopened on a stored route its drawing refused showed an empty canvas ("Draw reactants, a Reaction arrow,
 * then products") while the Text, the right pane and Export held the route. The drawing area says it did not open the
 * text, and why; the text stays the reaction, drawn on the right, and Export takes the text (no CDXML of the empty
 * drawing); a text that opens takes the line away. */
async function refusedOnOpen(): Promise<void> {
  saveDocument({title: '', text: REFUSED_ROUTE, mapNumbers: false});
  const view = await openApp();
  try {
    const said = () => view.root.querySelector(`.${SKETCHER_REFUSED_CLASS}`) as HTMLElement | null;
    await awaitCheck(() => !!said() && said()!.style.display !== 'none',
      'the drawing area does not say it did not open the text', 10000);
    const words = said()!.textContent ?? '';
    expect(words.startsWith('The drawing cannot open this text: ') && words.includes('chain'), true, words);
    expect(textBox(view).value, REFUSED_ROUTE, 'the text kept');
    expect(picture(view)?.dataset.reaction, REFUSED_ROUTE, 'the right pane draws the text');
    const dialog = (ribbonButton(view, EXPORT_BUTTON).click(), DG.Dialog.getOpenDialogs().find((d) =>
      d.title === EXPORT_TITLE));
    try {
      expect(dialog?.root.querySelector(`.${PREVIEW_CLASS}`) !== null, true, 'Export opens on the text');
    } finally {
      dialog?.close();
    }
    type(view, 'CC(=O)O.OCC>>CC(=O)OCC');
    await awaitCheck(() => said() === null || said()!.style.display === 'none',
      'a text the drawing opens: the line stays', 5000);
  } finally {
    view.close();
  }
}

/** Whether the right pane's `canvas` draws reaction `text`: as written, or as the Reaction panel draws it once it has
 * mapped it by itself (a single reaction, a moment after it is shown: reaction-panel.ts swaps in the mapped picture,
 * its molecules in the mapper's order), compared with the map numbers left out. */
function drawsReaction(crux: ChemWasm, canvas: HTMLCanvasElement | null, text: string): boolean {
  const drawn = canvas?.dataset.reaction;
  return drawn !== undefined && (drawn === text ||
    sameReactionText(crux, drawn.replace(/(\[[^\]:]+):\d+\]/g, '$1]'), text));
}

/** The right pane went blank during every drag in the drawing. While a button is held in the editing area the pane
 * keeps its picture, drawn again once the button is let go; and as it is drawn again a picture shows all along (the new
 * page drawn hidden over the old, which stays until it is drawn). The pane's own mapping is its picture's alone: the
 * Text and the document keep the reaction as stored (a run once read the pane's mapped picture,
 * "O[C:4]([CH3:5])=[O:6].[CH3:1][CH2:2][OH:3]>>...", where the test expected the text). */
async function keptWhileDragged(): Promise<void> {
  const [before, after] = ['CC(=O)O.OCC>>CC(=O)OCC', 'CC(=O)Cl.N>>CC(N)=O'];
  saveDocument({title: '', text: before, mapNumbers: false});
  const crux = await loadCrux('tests');
  const view = await openApp();
  try {
    await awaitCheck(() => sketcherShows(view, before), 'the sketcher did not open the document', 8000);
    // the pane's mapping done (or given up: its 5 s cap) before the watch, so that its swap is not taken for a drag's
    for (let t = 0; t < 60 && !/:\d+\]/.test(picture(view)?.dataset.reaction ?? ''); t++)
      await delay(100);
    expect(textBox(view).value, before, 'the Text as stored: no map numbers added');
    expect(loadDocument()?.text, before, 'the document as stored: no map numbers added');
    const shown = () => ([...drawing(view).querySelectorAll('canvas')] as HTMLCanvasElement[]).filter((c) =>
      c.style.height !== '' && getComputedStyle(c).visibility !== 'hidden');
    let blank = 0;
    const watch = new MutationObserver(() => {
      if (shown().length === 0)
        blank++;
    });
    watch.observe(drawing(view), {childList: true, subtree: true, attributes: true});
    try {
      sketchEl(view).dispatchEvent(new PointerEvent('pointerdown', {bubbles: true, composed: true}));
      type(view, after);
      await delay(1500); // past the Text's pause and the pane's
      expect(drawsReaction(crux, picture(view), before), true,
        `a button held: the pane keeps its picture: ${picture(view)?.dataset.reaction}`);
      window.dispatchEvent(new PointerEvent('pointerup', {bubbles: true}));
      await awaitCheck(() => drawsReaction(crux, picture(view), after) && shown().length > 0,
        'let go: the pane drawn again', 5000);
      expect(blank, 0, 'the pane never blank');
      expect(drawing(view).querySelectorAll('.chem-reactions-app-title').length, 1, 'one page left');
    } finally {
      watch.disconnect();
      window.dispatchEvent(new PointerEvent('pointerup', {bubbles: true}));
    }
  } finally {
    view.close();
  }
}

/** A press and release on the editing area's Crux canvas at `p` (the element's `positions` frame), as a person's. */
function pressAt(el: any, p: {x: number, y: number}): void {
  const canvas = el.shadowRoot.querySelector('[data-testid="canvas"]') as HTMLElement;
  const r = (el as HTMLElement).getBoundingClientRect();
  const init = (buttons: number): PointerEventInit => ({bubbles: true, cancelable: true, composed: true,
    clientX: r.left + p.x, clientY: r.top + p.y, pointerId: 1, pointerType: 'mouse', isPrimary: true, button: 0,
    buttons});
  canvas.dispatchEvent(new PointerEvent('pointerdown', init(1)));
  canvas.dispatchEvent(new PointerEvent('pointerup', init(0)));
}

/** Whether the editing area's drawing holds a catalytic cycle (Crux Sketch says reaction text drops one). */
const drawsCycle = (view: DG.View) => {
  try {
    return (sketchEl(view)?.warnings ?? []).some((w: {code?: string}) => w.code === 'cycles');
  } catch {
    return false;
  }
};

/** The demo page {@link EDITED_ROW} (texts over its arrows and captions) and its reaction as the reader writes it. */
async function demoPage(): Promise<{page: string, smiles: string}> {
  const page = (DG.DataFrame.fromCsv(await loadFileAsText(CDXML_DEMO_FILE)).col('cdxml')!
    .toList() as string[])[EDITED_ROW];
  const smiles = await cdxmlReactionSmiles(page);
  expect(!!smiles, true, 'the demo page read');
  return {page, smiles: smiles!};
}

/** CDXML `page` with its texts left out (every text that is no atom's label, and the steps' objects over and under
 * their arrows that were texts): the same structures, so the same reaction, and nothing a reaction text drops for
 * texts. */
function withoutTexts(page: string): string {
  const doc = new DOMParser().parseFromString(page, 'text/xml');
  const removed = new Set<string>();
  for (const t of Array.from(doc.getElementsByTagName('t'))) {
    if (t.parentElement?.tagName !== 'n') {
      removed.add(t.getAttribute('id') ?? '');
      t.remove();
    }
  }
  for (const step of Array.from(doc.getElementsByTagName('step'))) {
    for (const name of ['ReactionStepObjectsAboveArrow', 'ReactionStepObjectsBelowArrow']) {
      const ids = (step.getAttribute(name) ?? '').split(/\s+/).filter((id) => id && !removed.has(id));
      if (ids.length)
        step.setAttribute(name, ids.join(' '));
      else
        step.removeAttribute(name);
    }
  }
  return new XMLSerializer().serializeToString(doc);
}

/** Whether the editing area's drawing has texts that reaction text drops (Crux Sketch's `texts` loss). */
const drawsTexts = (view: DG.View) => ((sketchEl(view)?.warnings ?? []) as {kind: string, code?: string}[])
  .some((w) => w.kind === 'format' && w.code === 'texts');

/** Crux Sketch's own control `id` (its Undo, its Redo: a person's edit, `source: 'user'`) in the editing area. */
function sketchControl(view: DG.View, id: string): HTMLElement {
  const control = sketchEl(view)?.shadowRoot?.querySelector(`[data-testid="${id}"]`) as HTMLElement | null;
  expect(control !== null, true, `Crux Sketch's ${id}`);
  return control!;
}

/**
 * The line under the status naming what the Text does not keep of the drawing follows the drawing. A document reopened
 * on its kept drawing says it as the drawing opens (it said nothing until an edit); and a person's edit that leaves the
 * Text as it is says it at once (it waited for the next change of the Text): here Undo and Redo between the demo page
 * without its texts and with them, the same reaction. The page is written into the drawing by the host (no change the
 * app hears) and the Text is never typed: the demo route's reaction SMILES is not read back by the sketcher as written
 * ("the sketcher did not take the text").
 */
async function lossesFollowDrawing(): Promise<void> {
  forgetDocument();
  const {page, smiles} = await demoPage();
  const bare = withoutTexts(page);
  const said = (view: DG.View) => view.root.querySelector(`.${LOSSES_CLASS}`) as HTMLElement;
  const saysTexts = (view: DG.View) => said(view).style.display !== 'none' &&
    (said(view).textContent ?? '').includes(DRAWING_LOSS_WORDS['texts']);
  saveDocument({title: '', text: smiles, mapNumbers: false, drawing: page});
  const view = await openApp();
  try {
    // reopened on the kept drawing: said as it opens
    await awaitCheck(() => drawsTexts(view), 'reopened: the kept drawing was not opened', 15000);
    expect(textBox(view).value, smiles, 'reopened: the Text as stored');
    await awaitCheck(() => saysTexts(view),
      `reopened: the line does not name the texts: "${said(view).textContent}"`, 3000);
    // the page without its texts, then with them, written by the host; a person's Undo goes back to the one without
    sketchEl(view).setValue(bare, {history: 'edit'});
    await awaitCheck(() => !sketchEl(view).isPending && !drawsTexts(view),
      `the page without its texts was not drawn: ${sketchEl(view).lastError?.message ?? ''}`, 15000);
    sketchEl(view).setValue(page, {history: 'edit'});
    await awaitCheck(() => !sketchEl(view).isPending && drawsTexts(view), 'the page was not drawn again', 15000);
    sketchControl(view, 'toolbar.undo').click();
    await awaitCheck(() => !drawsTexts(view) && !saysTexts(view),
      `Undo, no texts drawn: the line still names them: "${said(view).textContent}"`, 5000);
    const text = textBox(view).value;
    // a person's Redo: the texts again, the reaction (so the Text) as it was
    sketchControl(view, 'toolbar.redo').click();
    await awaitCheck(() => drawsTexts(view), 'Redo did not bring the texts back', 5000);
    if (textBox(view).value !== text) // the drawing's reaction written otherwise: the Text's own change said it
      console.warn(`Chem | lossesFollowDrawing: the Redo changed the Text (${text} -> ${textBox(view).value})`);
    await awaitCheck(() => saysTexts(view), `Redo: the line does not name the texts: "${said(view).textContent}"`,
      3000);
  } finally {
    view.close();
  }
}

/** The right pane draws the drawing's texts over its arrows only while the drawing holds the document. A document
 * reopened on its kept drawing drew none of them until an edit. A Text typed reached the drawing 500 ms after the
 * keystroke and the pane 400 ms after: the pane drew the previous drawing's texts over the new text's arrows. */
async function stepTextsFollowDrawing(): Promise<void> {
  forgetDocument();
  const {page, smiles} = await demoPage();
  const crux = await loadCrux('tests');
  const texts = (view: DG.View) => Number(picture(view)?.dataset.stepTexts ?? 0);
  saveDocument({title: '', text: smiles, mapNumbers: false, drawing: page});
  const view = await openApp();
  try {
    // reopened on the kept drawing: the pane draws its texts once the drawing is read
    await awaitCheck(() => drawsTexts(view), 'reopened: the kept drawing was not opened', 15000);
    await awaitCheck(() => texts(view) > 0, `reopened: the pane drew no texts over its arrows (${texts(view)})`, 5000);
    // a Text typed: the drawing is the text's, with no texts, and so is the pane, from its first picture of it
    const typed = 'CC(=O)O.OCC>>CC(=O)OCC';
    type(view, typed);
    await awaitCheck(() => drawsReaction(crux, picture(view), typed), 'the pane did not draw the text typed', 5000);
    expect(texts(view), 0, 'the pane drew the previous drawing\'s texts over the text typed');
    await awaitCheck(() => !drawsTexts(view) && sketcherShows(view, typed), 'the drawing did not take the text', 5000);
    await delay(800); // past any picture drawn again
    expect(texts(view), 0, 'texts drawn over the text typed, later');
  } finally {
    view.close();
  }
}

/** A cycle drawn with the Cycle tool (Pd -> PhBr -> toluene -> Pd) came back as a straight route of 3 steps when
 * Reactions was opened again; two cycles as a refusal. The app kept only the text, which holds no cycle. The drawing is
 * kept beside it while it holds what the text drops, and the app reopens on it, the Text as written; a text typed drops
 * it (the text is the document again). */
async function cycleKept(): Promise<void> {
  const molecules = '[Pd].Brc1ccccc1.Cc1ccccc1';
  forgetDocument();
  let view = await openApp(); // on the default example (a text with no arrow draws no picture for openApp to wait on)
  let text = '';
  try {
    await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not open the example', 15000);
    type(view, molecules);
    textBox(view).dispatchEvent(new Event('blur')); // drawn at once
    await awaitCheck(() => sketchEl(view)?.document === 'reaction' && (sketchEl(view).positions?.atoms ?? [])
      .filter((a: unknown) => a).length >= 15, 'the sketcher did not draw the three molecules', 10000);
    const el = sketchEl(view);
    const tool = (el.shadowRoot.querySelector('[data-testid="toolbar.reaction.cycle.3"]') ??
      el.shadowRoot.querySelector('[data-testid="toolbar.reaction.cycles"]')) as HTMLElement | null;
    expect(tool !== null, true, 'the Cycle of 3 tool');
    tool!.click();
    for (const atom of [0, 1, 8]) { // Pd, bromobenzene's Br, toluene's methyl: each molecule in turn
      pressAt(el, el.positions.atoms[atom]);
      await delay(150);
    }
    await awaitCheck(() => drawsCycle(view) && textBox(view).value.split('--**--').length === 3,
      `no cycle of 3 steps drawn: ${textBox(view).value}`, 5000);
    text = textBox(view).value;
    await awaitCheck(() => (loadDocument()?.drawing ?? '').includes('crux-cycle') && loadDocument()?.text === text,
      'the drawing was not kept with the text', 10000);
  } finally {
    view.close();
  }
  view = await openApp();
  try {
    await awaitCheck(() => drawsCycle(view), 'reopened: the cycle is not drawn (a straight route?)', 10000);
    expect(textBox(view).value, text, 'reopened: the Text as written');
    expect(view.root.querySelector(`.${SKETCHER_REFUSED_CLASS}`)?.getAttribute('style')?.includes('none') ?? true,
      true, 'nothing refused');
    type(view, 'CC(=O)O.OCC>>CC(=O)OCC');
    await awaitCheck(() => loadDocument()?.text === 'CC(=O)O.OCC>>CC(=O)OCC' && loadDocument()?.drawing === undefined,
      'a text typed: the kept drawing dropped', 5000);
  } finally {
    view.close();
  }
}

/** There was no way to put a route into a table (the app, the pane's icons, the Reactions menu). The ribbon's Add to
 * table puts the reaction there: with no open table that has a reaction column, a new table (its reaction as written, a
 * route with its steps separated by --**--, read as a reaction column of units route; its title), opened; with one
 * open, a row added to it, the app kept in view. A text that reads as no reaction is not added, and it says so. */
async function routeIntoTable(): Promise<void> {
  forgetDocument();
  grok.shell.closeAll();
  const view = await openApp();
  let table: DG.DataFrame | null = null;
  const warnings: string[] = [];
  const warning = grok.shell.warning;
  grok.shell.warning = (s: string | HTMLElement) => {
    warnings.push(String(typeof s === 'string' ? s : s.textContent));
  };
  try {
    const example = defaultExample();
    expect(reactionUnitsOf(example.text), 'route', 'the default route\'s units');
    expect(`${reactionUnitsOf('CCO>>CC=O')} ${reactionUnitsOf('[Na+].[Cl-]>>[Na]Cl')} ` +
      `${reactionUnitsOf('[#6:1]>>[#6:1]O')} ${reactionUnitsOf('$RXN\n')}`, 'smiles smiles smarts rxn', 'units');
    const tables = () => Array.from(grok.shell.tableViews).map((v) => v.dataFrame);
    expect(tables().length, 0, 'no table open');
    ribbonButton(view, ADD_TO_TABLE_BUTTON).click();
    await awaitCheck(() => tables().some((t) => t.name === NEW_TABLE_NAME), 'no new table', 10000);
    table = tables().find((t) => t.name === NEW_TABLE_NAME)!;
    const column = table.col(REACTION_COLUMN)!;
    expect(`${table.rowCount} ${column.get(0) === example.text} ${table.col(TITLE_COLUMN)?.get(0) === example.title}`,
      '1 true true', 'one row: the route as written, its title');
    expect(`${column.semType} ${column.meta.units}`, `${DG.SEMTYPE.CHEMICAL_REACTION} route`, 'a reaction column');
    expect(grok.shell.v?.root !== view.root && grok.shell.tv?.dataFrame === table, true, 'the new table opened');
    // back in the app, the next reaction: a row of that table, the app kept in view
    grok.shell.v = view;
    type(view, 'CCO>>CC=O');
    textBox(view).dispatchEvent(new Event('blur'));
    await awaitCheck(() => status(view) === readsAs(1), `the status: ${status(view)}`, 5000);
    ribbonButton(view, ADD_TO_TABLE_BUTTON).click();
    await awaitCheck(() => table!.rowCount === 2, `no row added: ${table!.rowCount}`, 5000);
    expect(`${column.get(1)} ${tables().length} ${grok.shell.v?.root === view.root}`, 'CCO>>CC=O 1 true',
      'the reaction a row of the open table; no new table; the app in view');
    // a text that reads as no reaction
    type(view, 'CCO');
    textBox(view).dispatchEvent(new Event('blur'));
    await awaitCheck(() => status(view) !== readsAs(1), `the status: ${status(view)}`, 5000);
    ribbonButton(view, ADD_TO_TABLE_BUTTON).click();
    await awaitCheck(() => warnings.includes(NOTHING_TO_ADD), `not said: ${warnings.join(' | ')}`, 3000);
    expect(table.rowCount, 2, 'nothing added');
  } finally {
    grok.shell.warning = warning;
    view.close();
    if (table)
      grok.shell.closeTable(table);
  }
}

/** The Export panel's "As drawn" is the sketcher's own picture of the drawing as it is drawn (crux-sketch's svgWith, or
 * its svg before the picture API), not the reaction arranged: the default route's every atom, the sketcher's layers,
 * sized in points for ACS 1996. */
async function asDrawnFromTheApp(): Promise<void> {
  forgetDocument();
  const view = await openApp();
  const download = DG.Utils.download;
  const got: [string, BlobPart][] = [];
  DG.Utils.download = (name: string, content: BlobPart) => {
    got.push([name, content]);
  };
  const panel = () => DG.Dialog.getOpenDialogs().find((d) => d.title === EXPORT_TITLE);
  const note = () => panel()?.root.querySelector(`.${NOTE_CLASS}`)?.textContent ?? '';
  const atoms = (svg: string) => (svg.match(/data-atom=/g) ?? []).length;
  try {
    await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document', 15000);
    ribbonButton(view, EXPORT_BUTTON).click();
    await awaitCheck(() => !!panel(), 'the Export panel opens', 5000);
    const inputs = exportInputsOf(panel()!.root);
    expect(!!inputs?.asDrawn, true, 'As drawn offered');
    expect(inputs!.asDrawn!.value, false, 'off');
    inputs!.asDrawn!.value = true;
    await awaitCheck(() => note().includes(', as drawn, '), `previewed as drawn: ${note()}`, 10000);
    panel()!.getButton('Save').click();
    await awaitCheck(() => got.length === 1, 'one download', 5000);
    const saved = String(got[0][1]);
    const own = await sketcherPicture(sketchEl(view)).svg(drawnOptions({preset: ACS_1996, width: SINGLE_COLUMN,
      transparent: false, centre: false, mapNumbers: 'hide'})) ?? '';
    expect(saved.includes('data-layer=') && /width='[\d.]+pt'/.test(saved), true,
      `the sketcher's picture, in points: ${saved.slice(0, 300)}`);
    expect(atoms(saved) > 0 && atoms(saved) === atoms(own), true,
      `the drawing's every atom: ${atoms(saved)} saved, ${atoms(own)} drawn`);
  } finally {
    DG.Utils.download = download;
    panel()?.close();
    view.close();
  }
}

category('Reactions: app', () => {
  keepingTheDocument();

  test('the default example and the nine Examples: Crux reads and draws each, routes as their steps', async () => {
    const pictures = await picturesOf([defaultExample(), ...examples()]);
    try {
      expect(pictures.length, 10, 'the default example and nine examples');
      for (const p of pictures) {
        expect(p.frame(WIDTH).height > 40, true, `${p.reaction.title}: a frame ${p.frame(WIDTH).height} px high`);
        expect(p.error, null, p.reaction.title);
      }
      // the default example: five mapped steps, the map numbers hidden; then the paracetamol route (a chain), the
      // Diels-Alder RXN, ..., the imatinib route (a step list)
      expect(pictures.map((p) => p.steps().length).join(' '), '5 3 1 1 1 1 1 1 1 5', 'the steps Crux reads');
      expect(defaultExample().mapNumbers, false, 'the default example\'s map numbers hidden');
    } finally {
      pictures.forEach((p) => p.free());
    }
  });

  test('the app opens cold on the default example, both panes present: the text on the left, ' +
    'the route drawn on the right', async () => {
    forgetDocument();
    const view = await openApp();
    try {
      const example = defaultExample();
      expect(textBox(view).value, example.text, 'the editing area holds the default example');
      expect(shownTitle(view), example.title, 'the right pane\'s title');
      expect(status(view), 'A route of 5 steps.');
      // nothing on the view names the engine, dumps a file or says "confidence"
      for (const word of ['Crux', 'RDKit', 'M  END', 'confidence'])
        expect(view.root.textContent!.includes(word), false, `the view says "${word}"`);
      expect(ink(picture(view)!) > 5000, true, `the route drawn: ${ink(picture(view)!)} pixels`);
      expect(picture(view)!.dataset.reaction, example.text, 'the right pane draws the document');
      expect(drawing(view).textContent!.includes('Steps'), true, 'the steps told under the picture');
      expect(window.localStorage.getItem(DOCUMENT_KEY), null, 'opening keeps nothing');
    } finally {
      view.close();
    }
  });

  test('an Examples entry loads in place of the document and is kept; the app reopens with it', async () => {
    forgetDocument();
    const example = examples()[1]; // the Diels-Alder reaction, mapped
    let view = await openApp();
    try {
      const menu = view.ribbonMenu.find(EXAMPLES_MENU);
      expect(!!menu, true, 'the Examples menu');
      menu.find(example.title).click();
      await awaitCheck(() => textBox(view).value === example.text, 'the example was not loaded', 3000);
      expect(shownTitle(view), example.title, 'the right pane\'s title');
      await awaitCheck(() => picture(view)?.dataset.reaction === example.text, 'the example was not drawn', 3000);
      expect(status(view), 'One reaction.');
      expect(JSON.stringify(loadDocument()), JSON.stringify(example), 'kept at once');
    } finally {
      view.close();
    }
    view = await openApp();
    try {
      expect(textBox(view).value, example.text, 'the reopened app\'s reaction');
      expect(shownTitle(view), example.title);
    } finally {
      view.close();
    }
  });

  // "A route of five mapped steps..." stayed after Clear, a new drawing and a step deleted
  test('a reaction typed in place of the example drops its title; a title of one\'s own stays', titleFollows);
  test('an example\'s title stays through its edits; Clear drops it, its undo brings it back', titleKept);

  test('an edit is kept in the browser and drawn; the app reopens with it', async () => {
    forgetDocument();
    const text = 'CC(=O)Cl.OCC>>CCOC(C)=O';
    let view = await openApp();
    try {
      type(view, text);
      expect(loadDocument()?.text, text, 'kept at once');
      expect(status(view), 'One reaction.');
      await awaitCheck(() => picture(view)?.dataset.reaction === text && ink(picture(view)!) > 500,
        'the edited reaction was not drawn', 3000);
    } finally {
      view.close();
    }
    view = await openApp();
    try {
      expect(textBox(view).value, text, 'the reopened app\'s reaction');
      expect(ink(picture(view)!) > 500, true, 'drawn');
    } finally {
      view.close();
    }
  });

  test('a kept document that cannot be read, of version 1, or with no reaction typed ' +
    'gives the default example', async () => {
    const page = JSON.stringify({version: 1, reactions: [{title: 'Mine', text: 'CC>>CO', mapNumbers: true}]});
    for (const stored of ['{"version": 2, "title": ',
      '{"version": 99, "title": "a", "text": "CC>>CO", "mapNumbers": true}',
      '{"version": 2, "title": 1, "text": "CC>>CO", "mapNumbers": true}', page,
      '{"version": 2, "title": "Empty", "text": "  ", "mapNumbers": true}']) {
      window.localStorage.setItem(DOCUMENT_KEY, stored);
      expect(loadDocument(), null, stored);
      const view = await openApp();
      try {
        expect(textBox(view).value, defaultExample().text, stored);
      } finally {
        view.close();
      }
    }
  });

  test('Map numbers: the switch is kept, and the right pane draws as it says', mapNumbersDrawn);

  // the switch sat at the foot of the editing area, under the sketcher, the Add panel
  // and the Text, out of sight at 1440 x 900
  test('Map numbers: the switch is at the top of the preview, in sight without scrolling, and off on the default ' +
    'example', async () => {
    forgetDocument();
    const view = await openApp();
    try {
      const input = mapSwitch(view);
      expect(input !== null, true, 'the switch in the preview\'s controls');
      expect(drawing(view).contains(input), true, 'in the preview');
      expect(input.checked, false, 'off on the default example');
      const pane = drawing(view).parentElement!.getBoundingClientRect();
      const r = input.getBoundingClientRect();
      expect(r.height > 0 && r.top >= pane.top && r.bottom <= Math.min(pane.bottom, window.innerHeight), true,
        `in sight: ${Math.round(r.top)}-${Math.round(r.bottom)} in ${Math.round(pane.top)}-${Math.round(pane.bottom)}`);
      const title = drawing(view).querySelector('.chem-reactions-app-title')!.getBoundingClientRect();
      expect(r.top <= title.top, true, 'over the preview\'s title, the first of its controls');
    } finally {
      view.close();
    }
  });

  test('Reset to example brings the default example back and forgets the document', async () => {
    saveDocument({title: 'Mine', text: 'CC>>CO', mapNumbers: true});
    const view = await openApp();
    try {
      expect(textBox(view).value, 'CC>>CO', 'the kept document');
      ribbonButton(view, 'Reset to example').click();
      const dialog = DG.Dialog.getOpenDialogs().find((d) => d.title === 'Reset to example')!;
      expect(!!dialog, true, 'asked first');
      dialog.getButton('OK').click();
      await awaitCheck(() => textBox(view).value === defaultExample().text, 'the example did not come back', 3000);
      expect(shownTitle(view), defaultExample().title);
      expect(window.localStorage.getItem(DOCUMENT_KEY), null, 'nothing kept: the app opens on the example again');
    } finally {
      view.close();
    }
  });

  test('a text that reads as no reaction: its message under the editing area', unreadText);

  test('Export SVG: the reaction as one SVG, titled and drawn', async () => {
    const pictures = await picturesOf([defaultExample()]);
    try {
      const svg = documentSvg(pictures, WIDTH);
      const doc = new DOMParser().parseFromString(svg, 'image/svg+xml');
      expect(doc.querySelector('parsererror'), null, 'well formed');
      const root = doc.documentElement;
      expect(root.getAttribute('width'), `${WIDTH + 2 * PAGE_PADDING}px`);
      const texts = [...root.children].filter((e) => e.tagName === 'text').map((e) => e.textContent);
      expect(JSON.stringify(texts), JSON.stringify([defaultExample().title]), 'the title, as text');
      const drawn = [...root.children].filter((e) => e.tagName === 'svg');
      expect(drawn.length, 1, 'the picture');
      const image = new Image();
      image.src = URL.createObjectURL(new Blob([svg], {type: 'image/svg+xml'}));
      await image.decode();
      const p = pictures[0];
      const top = +drawn[0].getAttribute('y')! - p.frame(WIDTH).boxTop;
      const canvas = document.createElement('canvas');
      canvas.width = WIDTH;
      canvas.height = p.frame(WIDTH).height;
      canvas.getContext('2d')!.drawImage(image, -PAGE_PADDING, -top);
      expect(ink(canvas, true) > 2000, true, `${ink(canvas, true)} pixels`);
      URL.revokeObjectURL(image.src);
    } finally {
      pictures.forEach((p) => p.free());
    }
  });

  test('Export PNG: the reaction at 300 dpi, drawn', async () => {
    const pictures = await picturesOf([defaultExample()]);
    try {
      const png = await documentPng(pictures, WIDTH);
      expect(png.type, 'image/png');
      const bytes = new Uint8Array(await png.arrayBuffer());
      const at = bytes.findIndex((_, i) => String.fromCharCode(...bytes.subarray(i, i + 4)) === 'pHYs');
      expect(at, 8 + 25 + 4, 'a pHYs chunk after the header');
      expect(new DataView(bytes.buffer).getUint32(at + 4), 11811, 'pixels per metre: 300 dpi');
      const bitmap = await createImageBitmap(png);
      expect(bitmap.width, Math.round((WIDTH + 2 * PAGE_PADDING) * 300 / 96), 'the width at 300 dpi');
      const canvas = document.createElement('canvas');
      [canvas.width, canvas.height] = [bitmap.width, bitmap.height];
      canvas.getContext('2d')!.drawImage(bitmap, 0, 0);
      expect(ink(canvas, true) > 50000, true, `${ink(canvas, true)} pixels drawn`);
    } finally {
      pictures.forEach((p) => p.free());
    }
  });

  // one button, the Export panel on the reaction, in place of two one-click buttons and no choice; Export dropped the
  // text over the arrow, and the pictures drew none of it
  test('Export: CDXML is the drawing as drawn, its texts kept; the pictures draw its step texts',
    exportsTexts);

  test('the ribbon\'s Export... opens the Export panel on the reaction; Save downloads what it previews', async () => {
    forgetDocument();
    const view = await openApp();
    const download = DG.Utils.download;
    const got: [string, BlobPart, string | undefined][] = [];
    DG.Utils.download = (name: string, content: BlobPart, type?: string) => {
      got.push([name, content, type]);
    };
    const exportPanel = () => DG.Dialog.getOpenDialogs().find((d) => d.title === EXPORT_TITLE);
    try {
      expect(view.getRibbonPanels().flat().some((b) => /Export (SVG|PNG)/.test(b.textContent ?? '')), false,
        'the one-click buttons are gone');
      ribbonButton(view, EXPORT_BUTTON).click();
      await awaitCheck(() => !!exportPanel(), 'the Export panel opens', 5000);
      const preview = exportPanel()!.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
      await awaitCheck(() => !!preview.querySelector('canvas'), 'its preview is drawn', 5000);
      exportPanel()!.getButton('Save').click();
      await awaitCheck(() => got.length === 1, 'one download', 5000);
    } finally {
      DG.Utils.download = download;
      exportPanel()?.close();
      view.close();
    }
    expect(JSON.stringify(got.map(([name, , type]) => [name, type])), JSON.stringify([[SVG_FILE, 'image/svg+xml']]),
      'ACS 1996, a single column, SVG: the panel\'s first choices');
    expect(typeof got[0][1] === 'string' && (got[0][1] as string).includes(`width='235pt'`), true,
      'the SVG, 8.3 cm wide');
  });

  // there was no way to put a route into a table
  test('Add to table: the route into a new table, opened, its reaction column read as a route; the next reaction a ' +
    'row of that table; nothing to add, said', routeIntoTable);

  // the Export panel's "As drawn" exports the canvas's drawing, not an arranged scheme
  test('the ribbon\'s Export... offers As drawn: the sketcher\'s own picture of the drawing, as SVG in points',
    asDrawnFromTheApp);

  // the editing area is the Crux sketcher, the text box under it, each following the other; the app opens with the
  // reaction tools in query mode, as a cell's editor (the owner, 2026-11-17: it opened with the publication tools, its
  // Atom Properties without the query features), the text as it was
  test('the editing area is the Crux sketcher: a reaction document with the reaction tools in query mode, ' +
    'the text under it', async () => {
    forgetDocument();
    const view = await openApp();
    try {
      await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document',
        15000);
      expect(sketchEl(view).toolset, 'reaction');
      expect(sketchEl(view).mode, 'query', 'the query tools offered');
      expect(textBox(view).value, defaultExample().text, 'the text under the sketcher');
    } finally {
      view.close();
    }
  });

  // on a page's first open the sketcher's tools showed about half a second before it had read
  // the route, and Clear canvas clicked then did nothing (no change, the route drawn after it)
  test('the sketcher\'s box is covered until it has read the reaction: no click reaches it before', boxCovered);

  // at a fixed 480 px an arranged route's bonds were ~11 px; the sketcher is where the work happens
  // with the help pane open the preview was ~75 px
  test('the sketcher\'s pane takes 58% of the view, never under 480 px; the preview never under 280 px',
    paneShare);

  test('the app opens with the context and help panels shut, as other apps do; opened again, they stay', panelsShut);

  // after a name the drawing refused, CCO and CC=O added as reactant and product were drawn, but the
  // Text stayed empty, the preview asked for a reaction and Export had none, until a canvas edit
  // reagent names work in the Add panel (Chem's Crux file has the agent table)
  test('with Chem\'s agent table, HATU in the Reagent slot adds its structure over the arrow', namesAdded);

  test('a reaction made only with the Add panel, after a refused name, reaches the Text, the status and the preview',
    addedOnly);

  // a sketcher that does not start must not leave its box under the cover, taking every click
  test('a sketcher that does not start: the cover goes, the box says so, and the Text still works', notStarted);

  test('a kept drawing is read under the cover, a few seconds at most; one not read: its text drawn, said',
    keptDrawingCovered);

  test('an edit in the sketcher writes the text and the drawing: Clear canvas empties both', async () => {
    forgetDocument();
    const view = await openApp();
    try {
      await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document',
        15000);
      const cover = view.root.querySelector(`.chem-reactions-app-sketcher > .${SKETCH_LOADING_CLASS}`) as HTMLElement;
      expect(!!cover && getComputedStyle(cover).display === 'none', true, 'uncovered once the document is read');
      (sketchEl(view).shadowRoot.querySelector('[data-testid="toolbar.clear"]') as HTMLElement).click();
      await awaitCheck(() => textBox(view).value === '', 'the text did not follow the sketcher', 5000);
      await awaitCheck(() => picture(view) === null && drawing(view).textContent!.includes('Type or paste a reaction'),
        'the drawing did not follow the sketcher', 5000);
    } finally {
      view.close();
    }
  });

  // a molecule drawn alone turned the Text into a MOL file and the preview into an
  // error dump. Crux Sketch's value for a drawing with no arrow is a MOL block: the Text gets its SMILES, and the
  // status and the preview say in one line that a reaction needs an arrow.
  test('a molecule drawn alone is no error: the Text gets its SMILES, never a MOL file; one line asks for an arrow',
    async () => {
      forgetDocument();
      const view = await openApp();
      try {
        await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document',
          15000);
        const el = sketchEl(view);
        const press = (tool: string) =>
          (el.shadowRoot.querySelector(`[data-testid="toolbar.${tool}"]`) as HTMLElement).click();
        press('clear');
        await awaitCheck(() => textBox(view).value === '', 'Clear canvas did not empty the Text', 5000);
        press('bond.single'); // a click on the empty canvas draws a bond: ethane
        const canvas = el.shadowRoot.querySelector('[data-testid="canvas"]') as HTMLElement;
        const r = canvas.getBoundingClientRect();
        const init = (buttons: number): PointerEventInit => ({bubbles: true, cancelable: true, composed: true,
          clientX: r.left + r.width / 2, clientY: r.top + r.height / 2, pointerId: 1, pointerType: 'mouse',
          isPrimary: true, button: 0, buttons});
        canvas.dispatchEvent(new PointerEvent('pointerdown', init(1)));
        canvas.dispatchEvent(new PointerEvent('pointerup', init(0)));
        await awaitCheck(() => textBox(view).value !== '', 'the bond drawn did not reach the Text', 5000);
        expect(String(el.value).includes('M  END'), true, 'Crux Sketch gives a MOL block for a drawing with no arrow');
        expect(textBox(view).value, 'CC', 'the Text: ethane\'s SMILES');
        expect(status(view), NO_ARROW, 'the status');
        await awaitCheck(() => drawing(view).querySelector(`.${NO_ARROW_CLASS}`)?.textContent === NO_ARROW,
          'the preview did not say a reaction needs an arrow', 3000);
        expect((view.root.textContent ?? '').includes('M  END'), false, 'no MOL file on the view');
      } finally {
        view.close();
      }
    });

  // Arrange in the app must leave the Text exactly as it was. Crux Sketch gives an edited drawing's
  // reaction in its own writing, so the app keeps the person's Text while the reaction is the same one.
  test('one reaction however it is written: the molecules of each side in any order and spelling, map numbers kept',
    async () => {
      const crux = await loadCrux('tests');
      const same = (a: string, b: string) => sameReactionText(crux, a, b);
      expect(same('CC(=O)O.OCC>>CC(=O)OCC', 'C(C)O.OC(C)=O>>CCOC(C)=O'), true, 'reordered, respelt');
      expect(same('CC(=O)O.OCC>CN(C)C=O.O>CC(=O)OCC', 'OC(C)=O.CCO>O.CN(C)C=O>CCOC(C)=O'), true, 'the agents too');
      expect(same(defaultExample().text, defaultExample().text), true, 'the default route');
      expect(same('CC(=O)O.OCC>C[N+](C)=C(C)C.F[P-](F)(F)(F)(F)F>CC(=O)OCC',
        'CC(=O)O.OCC>C[N+](C)=C(C)C.F[P-](F)(F)(F)(F)F>CC(=O)OCC |f:2.3|'), true, 'ions grouped by a CX part');
      expect(same('CC(=O)O.OCC>>CC(=O)OCC', 'CC(=O)O.OCCC>>CC(=O)OCCC'), false, 'a molecule changed');
      expect(same('[CH3:1][OH:2]>>[CH2:1]=[O:2]', '[CH3:2][OH:1]>>[CH2:2]=[O:1]'), false, 'renumbered: a Map edit');
      expect(same('CCO>>CC=O', 'CCO>>CC=O.O'), false, 'a product added');
    });

  // a document opens arranged, which Clean Up leaves as it is; a route typed is drawn as read, and
  // Arrange (the Clean Up and Layout menu) lays it out anew: a person's edit of its layout only
  test('a drawing changed only in its layout (Arrange) keeps the Text as written', layoutOnly);

  // the default route opened at 20% as one long row; a reaction or route with no
  // coordinates of its own opens arranged in rows (reaction-arrange.ts)
  test('the default route opens arranged: in rows, its bonds long enough to read, the Text as written', opensArranged);

  // crux-sketch's setValue(text, {history: 'edit'}): a text typed is one edit of the drawing, its Undo taking it back
  test('a text typed is one edit of the drawing: Undo takes the drawing back to the document before it',
    textTypedUndone);

  test('an edit in the text is drawn in the sketcher', async () => {
    forgetDocument();
    const view = await openApp();
    const text = 'CC(=O)O.OCC>>CC(=O)OCC';
    try {
      await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document',
        15000);
      type(view, text);
      await awaitCheck(() => sketcherShows(view, text), 'the sketcher did not take the text', 5000);
    } finally {
      view.close();
    }
  });

  // the app's sketcher is the "crux sketcher widget" the cell editor's is, for the demo clips' steps
  test('the sketcher is a widget the platform finds: an area per atom, and the atoms drawn read', async () => {
    forgetDocument();
    const view = await openApp();
    const text = 'CC(=O)Cl.NCc1ccccc1>>CC(=O)NCc1ccccc1';
    try {
      await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document',
        15000);
      type(view, text);
      await awaitCheck(() => sketcherShows(view, text), 'the sketcher did not take the text', 5000);
      expect(view.root.querySelector('[data-u2-name="cruxSketch"] crux-sketch') === sketchEl(view), true,
        'the sketcher inside a widget named cruxSketch');
      const widget = DG.Widget.find(view.root.querySelector('.chem-reactions-app-sketcher') as HTMLElement);
      expect(widget !== null, true, 'a widget on the sketcher\'s host');
      await awaitCheck(() => widget!.getWidgetStatus().values?.atoms === 23, 'the status did not read 23 atoms', 5000);
      const status = widget!.getWidgetStatus();
      expect(['atom 0', 'atom 22', 'bond 0', 'canvas'].every((a) => !!status.hitAreas?.[a]), true,
        'an area per atom and bond, and the canvas');
      expect(status.values?.arrows, 1, 'the arrow drawn');
    } finally {
      view.close();
    }
  });

  // at 440 px the left bar scrolled and hid the arrow tool
  test('the sketcher is tall enough for its left bar: in a 1080-px window every group shows, the arrow tool too, ' +
    'no scroll button', async () => {
    forgetDocument();
    const view = await openApp();
    try {
      await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document',
        15000);
      const host = view.root.querySelector('.chem-reactions-app-sketcher') as HTMLElement;
      expect(Math.round(host.getBoundingClientRect().height), sketcherBoxHeight(), 'the box\'s height for this window');
      if (window.innerHeight < 1080) // a shorter window than the package tests' 1920x1080: the box as tall as it allows
        return;
      expect(sketcherBoxHeight(), SKETCHER_FULL_HEIGHT, 'the full height in a 1080-px window');
      for (const toolset of ['reaction', 'publication'] as const) {
        sketchEl(view).toolset = toolset; // the gear's choice; no edit
        await awaitCheck(() => leftBarWhole(sketchEl(view)),
          `${toolset} tools: the left bar scrolls or hides the arrow tool`, 3000);
      }
      sketchEl(view).toolset = 'reaction';
    } finally {
      view.close();
    }
  });

  test('closing the app destroys its sketcher', async () => {
    forgetDocument();
    const view = await openApp();
    await awaitCheck(() => sketcherShows(view, defaultExample().text), 'the sketcher did not show the document', 15000);
    const el = sketchEl(view);
    view.close();
    // destroyed: empty and inert, a write throws (crux-sketch API-037)
    await awaitCheck(() => {
      try {
        el.setValue('CC>>CO');
        return false;
      } catch {
        return true;
      }
    }, 'the sketcher still takes writes after the app closed', 5000);
  });

  // one scheme was left per open; a plain view's subs are not ended by View.close()
  test('closing the app frees the reaction it read', async () => {
    forgetDocument();
    const view = await openApp();
    const crux: ChemWasm = await loadCrux('tests');
    const proto = crux.ChemReactionScheme.prototype as any;
    const free = proto.free;
    let frees = 0;
    proto.free = function(this: unknown) {
      frees++;
      return free.call(this);
    };
    try {
      view.close();
      await awaitCheck(() => frees > 0, 'the app freed nothing when its view closed', 5000);
    } finally {
      proto.free = free;
    }
  });

  // a 1,000-residue reaction took Crux 7 s to draw, the tab frozen; the cells' rule applies
  test('a reaction too large to draw is not drawn: the status says so at once, and no picture is made', async () => {
    forgetDocument();
    const view = await openApp();
    try {
      const big = `${'C'.repeat(3000)}O>>${'C'.repeat(3000)}N`;
      const t0 = performance.now();
      type(view, big);
      await awaitCheck(() => !!drawing(view).querySelector('.chem-reactions-app-too-large'),
        'the right pane does not say the reaction is too large', 3000);
      expect(status(view), TOO_LARGE);
      expect(picture(view), null, 'no picture drawn');
      expect(performance.now() - t0 < 2000, true, `answered in ${Math.round(performance.now() - t0)} ms`);
      const [p] = await picturesOf([{title: 'Big', text: big, mapNumbers: false}]);
      try {
        expect(p.error, TOO_LARGE);
        expect(p.steps().length, 0);
        expect(p.svg(0, 0, WIDTH).startsWith('<text'), true, 'the export writes it as text');
      } finally {
        p.free();
      }
    } finally {
      view.close();
    }
  });

  // the right pane, its status and the Export's "Not saved" say why a value is not drawn in words, never in Crux's
  // notation (a help screenshot once showed "the fragment group f:[0, 1, 2] names components out of range or of two
  // roles" there)
  test('a value whose |f:| groups Crux does not read is said in words, not in Crux\'s notation', async () => {
    const crux: ChemWasm = await loadCrux('tests');
    for (const text of ['CC>>CC |f:0.1.2|', 'CC.O>>CC |f:0.2|']) {
      const p = new ReactionPicture(crux, {title: 'Groups', text, mapNumbers: false});
      try {
        expect(p.error, FRAGMENT_GROUPS_UNREAD, text);
        expect(p.svg(0, 0, WIDTH).includes('fragment group'), false, `${text}: the export's text`);
      } finally {
        p.free();
      }
    }
  });

  // a wasm trap may leave Crux's instance broken (a double free, for one). The module recovers inside the call
  // (crux/crux-recover.ts) and Crux stays on; a call that traps on a fresh instance too fails alone. crux-core's real
  // test trap (the next engine call panics), not a stub.
  test('a trap while the picture is drawn is recovered from; trapped twice, that picture alone is written as text',
    async () => {
      const crux: ChemWasm = await loadCrux('tests');
      const doc: DocumentReaction = {title: 'Ester', text: 'CC(=O)O.OCC>>CC(=O)OCC', mapNumbers: false};
      const pictures = [new ReactionPicture(crux, doc), new ReactionPicture(crux, doc)];
      const [once, twice] = pictures;
      try {
        expect(once.steps().length, 1, 'read before the trap');
        expect(armCruxTrapForTests(1), true, 'the main thread\'s module is loaded');
        expect(once.svg(0, 0, WIDTH).startsWith('<svg'), true, 'drawn, a trap in it');
        expect(once.error, null);
        expect(armCruxTrapForTests(1, 2), true);
        expect(twice.svg(0, 0, WIDTH).startsWith('<text'), true, 'trapped twice: written as text');
        expect(twice.error, REFUSED_MESSAGE, 'its words name no engine');
        expect(cruxLoadFailed(), false, 'Crux stays on');
        const next = new ReactionPicture(crux, doc);
        pictures.push(next);
        expect(next.svg(0, 0, WIDTH).startsWith('<svg'), true, 'a new picture is drawn');
        expect(once.steps().length, 1, 'the first picture still reads');
      } finally {
        for (const p of pictures)
          p.free();
      }
    });

  test('a stored text the drawing does not open: it says so and why; the text is drawn and exported', refusedOnOpen);

  test('a drag in the drawing: the right pane keeps its picture, never blank, drawn again once let go',
    keptWhileDragged);

  test('a catalytic cycle drawn is kept: the app reopens on the drawing, the Text as written; a text typed drops it',
    cycleKept);

  test('the line of what the Text does not keep follows the drawing: an edit that leaves the Text, a reopened drawing',
    lossesFollowDrawing);

  test('the right pane draws the drawing\'s texts over its arrows only while the drawing holds the document',
    stepTextsFollowDrawing);

  // "Reaction Sketcher", not "Reactions" (Apps showed it as Reactions | Reactions)
  test('the app is called Reaction Sketcher: its entry in Apps and its view', async () => {
    const apps = DG.Func.find({package: 'Chem'}).filter((f) => f.options['role'] === 'app');
    expect(apps.some((f) => f.name === 'Reaction Sketcher' || f.friendlyName === 'Reaction Sketcher'), true,
      `Chem's apps: ${apps.map((f) => f.friendlyName || f.name).join(', ')}`);
    expect(apps.some((f) => f.name === 'Reactions' || f.friendlyName === 'Reactions'), false,
      'no app called Reactions');
    forgetDocument();
    const view = await openApp();
    try {
      expect(view.name, 'Reaction Sketcher', 'the view\'s name');
    } finally {
      view.close();
    }
  });
});

// the help checks, a category of their own (the app's category is long enough)
category('Reactions: help', () => {
  keepingTheDocument();

  // the help pane said "Help unavailable" for Reset to example and the Reaction dialog
  test('the app, Reset to example, full screen and the Export panel open Chem\'s Reactions help page, which is served',
    async () => {
      forgetDocument();
      const view = await openApp();
      const dialog = (title: string) => DG.Dialog.getOpenDialogs().find((d) => d.title === title);
      const helpOf = (title: string) => dialog(title)?.helpUrl ?? '';
      const atSection = (url: string | null, section: string) => !!url?.endsWith(`${REACTIONS_HELP_FILE}#${section}`);
      try {
        expect(atSection(view.helpUrl, 'the-reaction-sketcher'), true, `the view's help: ${view.helpUrl}`);
        ribbonButton(view, 'Reset to example').click();
        expect(atSection(helpOf('Reset to example'), 'the-reaction-sketcher'), true, helpOf('Reset to example'));
        dialog('Reset to example')!.close();
        const fullScreen = drawing(view).querySelector('i.fa-expand-alt') as HTMLElement;
        expect(!!fullScreen, true, 'the full-screen button');
        fullScreen.click();
        await awaitCheck(() => !!dialog('Reaction'), 'full screen did not open', 5000);
        expect(atSection(helpOf('Reaction'), 'the-reaction-panel'), true, helpOf('Reaction'));
        dialog('Reaction')!.close();
        ribbonButton(view, EXPORT_BUTTON).click();
        await awaitCheck(() => !!dialog(EXPORT_TITLE), 'the Export panel did not open', 5000);
        expect(atSection(helpOf(EXPORT_TITLE), 'export'), true, helpOf(EXPORT_TITLE));
        dialog(EXPORT_TITLE)!.close();
        // the page the help pane loads: served from the package, with the sections the links point at
        const response = await fetch(reactionsHelpUrl());
        expect(response.ok, true, `the page is served: ${response.status}`);
        const page = (await response.text()).replace(/\r\n/g, '\n');
        expect(page.startsWith('# Reactions\n'), true, page.slice(0, 80));
        for (const heading of ['The Reaction Sketcher', 'The Reaction panel', 'Editing a reaction cell', 'Export'])
          expect(page.includes(`\n## ${heading}\n`), true, `the section "${heading}"`);
        // where the example tables are (a person once took 2:20 to find one), and they are there: the folder the page
        // names, by its path (package.ts's `_package` is never given the platform in the tests' bundle: its files were
        // System:AppData/ itself)
        expect(page.includes('**Files** > **App Data** > **Chem**'), true, 'the page says where the tables are');
        for (const file of ['test-reactions.csv', 'reactions-cdxml.csv']) {
          expect(page.includes(`\`${file}\``), true, `the page names ${file}`);
          expect(await grok.dapi.files.exists(`${APP_DATA_CHEM}/${file}`), true, `${file} in App Data/Chem`);
        }
      } finally {
        for (const title of ['Reset to example', 'Reaction', EXPORT_TITLE])
          dialog(title)?.close();
        view.close();
      }
    });
});
