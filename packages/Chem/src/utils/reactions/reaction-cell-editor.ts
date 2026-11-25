/**
 * The reaction cell editor: a reaction cell's menu "Edit reaction..." and a double-click on any reaction cell
 * (package.ts editReactionCell) open the Crux sketcher on its reaction in a dialog, in both modes (the reaction editor
 * is Crux's whatever the Chemical Engine, so a reaction is edited the same way in either), as a reaction document
 * (crux-sketch's reaction mode) with the reaction tools (toolset `reaction`: a cell may hold a template with real
 * SMARTS, so the query tools are on). OK, the dialog's own button, writes the sketcher's `value`: the reaction in the
 * notation the cell gave, the cell's own string while nothing was edited, so OK without an edit writes nothing; Cancel
 * writes nothing. The footer's middle, between the dialog's edge and OK and Cancel, is the foot-gun flags' fixed slot
 * (utils/footguns/footer-slot.ts), which {@link ReactionCellEditor.onChange} feeds; a block there disables OK, which
 * the dialog honours for a click and for Enter. The same in RDKit mode, flags included (the reaction editor is the Crux
 * editor whatever the engine).
 */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject} from 'rxjs';
import type {SketchWarning, SketcherHandle} from '../../../vendor/crux-sketch/crux-sketch';
import {focusDrawingSoon, keepSketcherKeys, loadCruxSketchForReactions} from '../../crux/crux-sketcher';
import {cruxForReactions} from '../../rendering/crux-reaction-drawing';
import {sketchReactionMapper} from './reaction-sketch-mapper';
import {CDXML_UNITS} from '../../rendering/cdxml-reaction-cells';
import {reactionsHelpUrl} from './reactions-help';
import {arrangeOpened} from './reaction-arrange';
import {reactionTextFormat} from './reaction-losses';
import type {FooterSlot} from '../footguns/footer-slot';

/** The dialog's title, and its content's class (what the tests find it by). */
export const REACTION_EDITOR_TITLE = 'Reaction';
export const REACTION_EDITOR_CLASS = 'chem-reaction-cell-editor';
/** The class of the line saying the cell's text was not read. */
export const REACTION_EDITOR_REFUSED_CLASS = 'chem-reaction-cell-editor-refused';
/** The sketcher's box as the dialog opens (CSS px); the dialog resizes, and the box fills it as it does (else the box
 * keeps its size in a dialog made larger, empty space below and right of it). */
const EDITOR_WIDTH = 860;
/** The dialog's contents' padding (Datagrok's .d4-dialog-contents: 12 px) and border, around the box (CSS px). */
const DIALOG_PADDING = 12;
const DIALOG_BORDER = 1;
/** The box's height in which Crux Sketch's left bar shows every group with no scroll button (measured with crux-sketch
 * d0b8bb9's one reaction bar with the reaction tools: 710 px; 580 with the bars before it). A 1440x900 window still
 * gives it 720. */
export const EDITOR_FULL_HEIGHT = 716;
/** What the dialog adds around the box (its title and footer, 105 px measured, and the line saying why a cell was not
 * read), and the room kept free above and below it in the window (CSS px). */
const DIALOG_CHROME = 132;
const WINDOW_MARGIN = 24;
// A length of the dialog's style in whole pixels: Datagrok reads a dialog's left and top back as integers
// (Modal.clipToScreen threw FormatException at a left of 622.4px, measured on a screen scaled 125%).
const px = (n: number): string => `${Math.round(n)}px`;
/** The least height the box is given in a window too short for {@link EDITOR_FULL_HEIGHT}: its bars then scroll;
 * and the least width the dialog is made in a narrow window (CSS px). */
const EDITOR_MIN_HEIGHT = 320;
const EDITOR_MIN_WIDTH = 480;
/** The header's full-screen switch (the dialog had none, in its title bar or its menus), and its tip in and out of
 * full screen (not "Full screen" while full screen). */
export const EDITOR_FULL_SCREEN_CLASS = 'chem-reaction-cell-editor-full-screen';
export const EDITOR_FULL_SCREEN = 'Full screen';
export const EDITOR_FULL_SCREEN_OFF = 'Exit full screen';

/** What the editor says of the sketcher's reason `message` for not opening a cell (a message a user sees names no
 * engine and asks for no reload): the reason itself, unless it names the engine or asks for a reload (the vendored
 * sketcher's engine stop does), then the plain {@link NOT_READ_NOW}, the reason going to the console. */
export const NOT_READ_NOW = 'it could not be read just now';
export function plainRefusal(message: string): string {
  if (!/crux|rdkit|reload/i.test(message))
    return message;
  console.warn(`Chem | the reaction editor did not open the cell: ${message}`);
  return NOT_READ_NOW;
}

/** What the editor says, once, of a cell it does not open, for `reason` (its sketcher's, in plain words). */
export function refusedWords(reason: string): string {
  return `The editor cannot open this reaction: ${reason.replace(/[.\s]+$/, '')}. Cancel closes it; the cell ` +
    'stays as it is.';
}

/** The box's height in a window `windowHeight` px tall: tall enough for the left bar's every group, unless the window
 * is too short for that and the dialog; then as tall as the dialog can be in it, the bar scrolling (its scroll
 * buttons). */
export function editorBoxHeight(windowHeight: number = window.innerHeight): number {
  return Math.max(EDITOR_MIN_HEIGHT, Math.min(EDITOR_FULL_HEIGHT, windowHeight - DIALOG_CHROME - 2 * WINDOW_MARGIN));
}

type Rect = {x: number, y: number, width: number, height: number};

/** The cover over a reaction sketcher's box until Crux Sketch has read and drawn the reaction, and what it says (on a
 * page's first open the sketcher's tools show about half a second before it has read the reaction, and a click on them
 * then is lost: Clear canvas did nothing, and the route was drawn after it). */
export const SKETCH_LOADING_CLASS = 'chem-reaction-sketch-loading';
export const SKETCH_LOADING = 'Drawing the reaction…';

/** The rows the reaction arrows Crux draws in `element` stand in (a route arranged breaks into rows at its arrows):
 * their middles' heights, a new row wherever two differ by more than three bond lengths (`bondLength`, px). */
export function arrowRows(element: HTMLElement | null, bondLength: number): number {
  const ys = [...element?.shadowRoot?.querySelectorAll('g[data-arrow]') ?? []].map((g) => {
    const r = g.getBoundingClientRect();
    return r.top + r.height / 2;
  }).sort((a, b) => a - b);
  return ys.filter((y, i) => i === 0 || y - ys[i - 1] > 3 * bondLength).length;
}

/** The box Crux draws a reaction in (the cell editor's, the Reactions app's), a widget named `cruxSketch`
 * as Chem's molecule sketcher is (HOST-043), so that the platform's tests find it as the BDD library's "crux sketcher
 * widget": its status (`getWidgetStatus()`) gives a hit area for every atom and bond Crux draws (`atom 0`, `bond 0`, by
 * index), the canvas, and the readings `ready`, `pending`, `atoms`, `bonds`, `arrows` (the reaction arrows drawn),
 * `rows` (the rows they stand in), `bondLength` (px, rounded) and
 * `changes` (a person's edits); and `isRenderPending` / `onRendered` to settle on. Its host sets `sketch` once Crux
 * Sketch has read the reaction, counts `changes`, clears `loading`, and completes `onRendered` and detaches it as it
 * lets go. While `loading`, the box is under a cover ({@link SKETCH_LOADING_CLASS}) that takes every click. */
export class ReactionSketchArea extends DG.Widget {
  sketch: SketcherHandle | null = null;
  changes = 0;
  /** Its host does not open the reaction: the box says so in place of a sketcher. */
  refused = false;
  readonly onRendered = new Subject<void>();
  private loadingNow = true;
  // Over the box until the sketcher has read the reaction (or could not start): a click on the tools before then is
  // the cover's, so none is lost to a sketcher not ready yet.
  private readonly cover = ui.divText(SKETCH_LOADING, SKETCH_LOADING_CLASS);

  /** `root` is the box, positioned (its host's): the cover fills it. */
  constructor(root: HTMLElement) {
    super(root);
    this.name = 'cruxSketch';
    this.addStatusProvider('crux-sketch', () => this.status());
    // the keys the sketcher takes stay in it (Ctrl+A in the cell editor's canvas also selected every row of the grid
    // behind the dialog); its undo and redo are its own, not the platform's
    keepSketcherKeys(root);
    Object.assign(this.cover.style, {position: 'absolute', inset: '0', zIndex: '2', display: 'flex',
      alignItems: 'center', justifyContent: 'center', background: 'rgba(255, 255, 255, 0.7)', color: 'var(--grey-5)',
      cursor: 'progress'});
    root.append(this.cover);
  }

  /** True until the sketcher has read the reaction or could not start; the box is covered meanwhile. */
  get loading(): boolean {
    return this.loadingNow;
  }

  set loading(on: boolean) {
    this.loadingNow = on;
    this.cover.style.display = on ? 'flex' : 'none';
  }

  get isRenderPending(): boolean {
    return this.loading || (this.sketch?.isPending ?? false);
  }

  private status(): {hitAreas: {[name: string]: Rect}, values: {[name: string]: string | number | boolean}} {
    const sketch = this.sketch;
    const element = this.root.querySelector<HTMLElement>(':scope > crux-sketch');
    // a reaction the editor does not open: settled, with nothing drawn and nothing to edit (the BDD's steps on its OK
    // and Cancel wait for a settled box; a demo scenario once waited for a sketcher that is no longer there)
    if (this.refused && !this.loading && this.root.isConnected) {
      return {hitAreas: {}, values: {ready: true, pending: false, refused: true, atoms: 0, bonds: 0, arrows: 0,
        changes: this.changes}};
    }
    if (sketch === null || element === null || !this.root.isConnected)
      return {hitAreas: {}, values: {ready: false, pending: this.isRenderPending, changes: this.changes}};
    const {atoms, bonds, bondLength} = sketch.positions;
    const at = element.getBoundingClientRect();
    const root = this.root.getBoundingClientRect();
    const [ox, oy] = [at.left - root.left, at.top - root.top];
    const half = Math.max(4, bondLength / 4);
    const hitAreas: {[name: string]: Rect} = {};
    const add = (name: string, p: {x: number, y: number} | null) => {
      if (p !== null)
        hitAreas[name] = {x: ox + p.x - half, y: oy + p.y - half, width: 2 * half, height: 2 * half};
    };
    atoms.forEach((p, i) => add(`atom ${i}`, p));
    bonds.forEach((p, i) => add(`bond ${i}`, p));
    const canvas = sketch.layout.canvas;
    if (canvas)
      hitAreas['canvas'] = {x: ox + canvas.x, y: oy + canvas.y, width: canvas.width, height: canvas.height};
    const arrows = element.shadowRoot?.querySelectorAll('g[data-arrow]').length ?? 0;
    return {hitAreas, values: {ready: true, pending: sketch.isPending, atoms: atoms.length, bonds: bonds.length, arrows,
      rows: arrowRows(element, bondLength), bondLength: Math.round(bondLength), changes: this.changes}};
  }
}

/** An open reaction cell editor: its dialog, its sketcher once loaded, and what a person changed. */
export class ReactionCellEditor {
  readonly dialog: DG.Dialog;
  /** The cell's text as the editor opened: what OK compares with (the TPL-12 check too). */
  readonly original: string;
  /** The sketcher, once Crux Sketch is loaded and has read the cell (null before, or when it could not start). */
  sketch: SketcherHandle | null = null;
  /** True once a person changed the reaction: only then does OK write. */
  edited = false;
  /** The text OK would write, each time a person changes the reaction (the footer's flags read it); completed when the
   * dialog closes, which ends every subscription to it. */
  readonly onChange = new Subject<string>();
  /** Resolves once the sketcher has read the cell (or could not start). */
  readonly ready: Promise<void>;
  /** The footer's foot-gun slot (footguns/editors.ts `watchReactionEditor`, set as the editor opens, in both modes): a
   * block
   * there disables OK, and OK writes nothing while it shows one (Enter included). */
  slot: Pick<FooterSlot, 'blocked'> | null = null;
  /** Why Crux Sketch did not read the cell (its `lastError`: a convergent route "opens in a later version", RXN-06), or
   * null: the dialog says it above the canvas, and OK then writes nothing, whatever is drawn. */
  refused: string | null = null;
  // The sketcher's value once the cell opened arranged by an Arrange edit (reaction-arrange.ts): a person's edits
  // undone come back to it, and OK then writes nothing. Null otherwise.
  private arrangedValue: string | null = null;
  // The line above the canvas saying why the cell was not read (shown only then).
  private readonly notice = ui.divText('', REACTION_EDITOR_REFUSED_CLASS);
  // The sketcher's box, which gives way when the window is too small for the dialog.
  private readonly box: HTMLElement;
  // The dialog's place and size before it went full screen (its style's), or null while it is not.
  private beforeFullScreen: {left: string, top: string, width: string, height: string, box: string} | null = null;

  constructor(readonly cell: DG.GridCell) {
    this.original = cell.cell.value ?? '';
    // the same in both modes: Crux Sketch, whatever the Chemical Engine
    const host = this.box = ui.div([], REACTION_EDITOR_CLASS);
    // the box fills the dialog's room (its contents are a flex column), from its opening height; never below the
    // least the bars scroll in
    Object.assign(host.style, {flex: '1 1 auto', height: `${editorBoxHeight()}px`,
      minHeight: `${EDITOR_MIN_HEIGHT}px`, position: 'relative'});
    Object.assign(this.notice.style,
      {display: 'none', color: 'var(--red-3)', margin: '0 0 6px 0', maxWidth: `${EDITOR_WIDTH}px`});
    const titled = {title: REACTION_EDITOR_TITLE, helpUrl: reactionsHelpUrl('editing-a-reaction-cell')};
    const body = ui.divV([this.notice, host]);
    Object.assign(body.style, {flex: '1 1 auto', minHeight: '0'});
    this.dialog = ui.dialog(titled).add(body)
      .onOK(() => this.write())
      .show({resizable: true, width: EDITOR_WIDTH + 2 * (DIALOG_PADDING + DIALOG_BORDER)});
    // the sketcher's keys pressed anywhere in the dialog (but a text field) never reach the grid behind it (Ctrl+A over
    // the canvas, the dialog focused as it opened, selected every row of the table)
    keepSketcherKeys(this.dialog.root, true);
    // kept in the window as it opens and as the window is resized
    const fit = () => this.fitToWindow();
    // made taller than the window's 80% (the contents' own cap), the box follows still
    const contents = this.dialog.root.querySelector('.d4-dialog-contents') as HTMLElement | null;
    if (contents)
      contents.style.maxHeight = 'none';
    this.fitToWindow();
    window.addEventListener('resize', fit);
    this.addFullScreenSwitch();
    const area = new ReactionSketchArea(host);
    const closed = this.dialog.onClose.subscribe(() => {
      closed.unsubscribe();
      window.removeEventListener('resize', fit);
      this.onChange.complete();
      this.sketch?.destroy();
      area.onRendered.complete();
      area.detach();
    });
    this.ready = this.start(area);
  }

  /**
   * The dialog kept inside a window `width` x `height` (at 1024 x 700 its OK, Cancel and x were off-screen,
   * Escape the only way out): no wider or taller than the window less its margins, the sketcher's box giving way (to
   * its least height, its bars then scrolling), and moved back in where it hangs out.
   */
  fitToWindow(width: number = window.innerWidth, height: number = window.innerHeight): void {
    const root = this.dialog.root;
    const box = this.box;
    if (!box || !root.isConnected)
      return;
    if (this.fullScreen) { // full screen: the window's size again
      this.fillWindow();
      return;
    }
    const [roomWidth, roomHeight] = [width - 2 * WINDOW_MARGIN, height - 2 * WINDOW_MARGIN];
    let r = root.getBoundingClientRect();
    if (r.width > roomWidth)
      root.style.width = px(Math.max(EDITOR_MIN_WIDTH, roomWidth));
    if (root.style.height) { // a height of the person's: the box, flexed, follows it
      if (r.height > roomHeight)
        root.style.height = px(Math.max(EDITOR_MIN_HEIGHT, roomHeight));
    } else { // the box as tall as the window lets it be, up to its full height: smaller, and taller again
      const now = box.getBoundingClientRect().height;
      const wanted = Math.max(EDITOR_MIN_HEIGHT, Math.min(EDITOR_FULL_HEIGHT, roomHeight - (r.height - now)));
      if (Math.abs(wanted - now) > 0.5)
        box.style.height = `${wanted}px`;
    }
    r = root.getBoundingClientRect();
    const left = Math.max(WINDOW_MARGIN, Math.min(r.left, width - WINDOW_MARGIN - r.width));
    const top = Math.max(WINDOW_MARGIN, Math.min(r.top, height - WINDOW_MARGIN - r.height));
    if (Math.abs(left - r.left) > 0.5 || Math.abs(top - r.top) > 0.5) {
      // the dialog is placed by its style's left and top: Datagrok's are fixed, in the window's own coordinates
      // (as measured); one placed in a parent, in its parent's
      const fixed = getComputedStyle(root).position === 'fixed';
      const parent = fixed || !root.offsetParent ? {left: 0, top: 0} : root.offsetParent.getBoundingClientRect();
      Object.assign(root.style, {left: px(left - parent.left), top: px(top - parent.top)});
    }
  }

  /** Whether the dialog fills the window ({@link toggleFullScreen}). */
  get fullScreen(): boolean {
    return this.beforeFullScreen !== null;
  }

  /** The dialog filling the window less its margins, the sketcher's box with it, or back to its place and size: the
   * header's switch, beside its help and close, or a double-click on its title. */
  toggleFullScreen(): void {
    const root = this.dialog.root;
    if (!this.box)
      return;
    const before = this.beforeFullScreen;
    if (before) {
      this.beforeFullScreen = null;
      Object.assign(root.style, {left: before.left, top: before.top, width: before.width, height: before.height});
      this.box.style.height = before.box;
      this.fitToWindow();
    } else {
      this.beforeFullScreen = {left: root.style.left, top: root.style.top, width: root.style.width,
        height: root.style.height, box: this.box.style.height};
      this.fillWindow();
    }
    const icon = root.querySelector(`.${EDITOR_FULL_SCREEN_CLASS}`) as HTMLElement | null;
    if (icon) {
      icon.classList.toggle('fa-expand', !this.fullScreen);
      icon.classList.toggle('fa-compress', this.fullScreen);
      icon.setAttribute('aria-label', this.fullScreenTip());
      ui.tooltip.hide(); // the tip shown at the click was the other state's: the next hover shows this one's
    }
  }

  /** What the full-screen switch says now: in or out of full screen. */
  fullScreenTip(): string {
    return this.fullScreen ? EDITOR_FULL_SCREEN_OFF : EDITOR_FULL_SCREEN;
  }

  // Full screen: the dialog the window's size less its margins, its box (flexed) taking the room.
  private fillWindow(): void {
    const root = this.dialog.root;
    const fixed = getComputedStyle(root).position === 'fixed';
    const parent = fixed || !root.offsetParent ? {left: 0, top: 0} : root.offsetParent.getBoundingClientRect();
    Object.assign(root.style, {left: px(WINDOW_MARGIN - parent.left), top: px(WINDOW_MARGIN - parent.top),
      width: px(Math.max(EDITOR_MIN_WIDTH, window.innerWidth - 2 * WINDOW_MARGIN)),
      height: px(Math.max(EDITOR_MIN_HEIGHT, window.innerHeight - 2 * WINDOW_MARGIN))});
    if (this.box)
      this.box.style.height = `${EDITOR_MIN_HEIGHT}px`; // its least: the flex gives it the rest
  }

  // The header's full-screen switch, before its help and close icons, and the title's double-click.
  private addFullScreenSwitch(): void {
    const header = this.dialog.root.querySelector('.d4-dialog-header') as HTMLElement | null;
    if (!header)
      return;
    const icon = ui.iconFA('expand', () => this.toggleFullScreen());
    ui.tooltip.bind(icon, () => this.fullScreenTip()); // read at each hover: a second bind left the first tip shown
    icon.setAttribute('aria-label', this.fullScreenTip());
    icon.classList.add(EDITOR_FULL_SCREEN_CLASS);
    Object.assign(icon.style, {fontSize: '14px', margin: '0 8px', color: 'var(--grey-5)', cursor: 'pointer',
      order: '99'});
    // the header moves the dialog on a press: the switch's own press is not a move
    icon.addEventListener('mousedown', (e) => e.stopPropagation());
    const title = header.querySelector('.d4-dialog-title');
    header.insertBefore(icon, title ? title.nextSibling : header.firstChild);
    header.addEventListener('dblclick', (e) => {
      if (e.target !== icon)
        this.toggleFullScreen();
    });
  }

  /** The text OK would write now: a CDXML cell's page as CDXML (its texts and arrow styles kept: after an edit Crux
   * Sketch's `value` is reaction SMILES; its `cdxml` the whole page), else the sketcher's
   * value. */
  get value(): string {
    if (!this.edited) // no person's edit (an Arrange on opening is none): OK writes nothing
      return this.original;
    if (this.sketch && this.holdsCdxml) {
      try {
        return this.sketch.cdxml;
      } catch (e) {
        console.error('Chem | the edited page could not be written as CDXML', e);
        return this.original; // nothing written but the cell's own page
      }
    }
    return this.sketch?.value ?? this.original;
  }

  /** Whether the cell's column holds CDXML pages (units `cdxml`, rendering/cdxml-reaction-cells.ts). */
  get holdsCdxml(): boolean {
    return this.cell.cell.column?.meta.units === CDXML_UNITS;
  }

  /** What the sketcher says its drawing loses in each notation (Crux Sketch's `warnings`; none before it is loaded
   * or when it cannot say): the footer names what OK's text does not keep (reaction-losses.ts). */
  get warnings(): readonly SketchWarning[] {
    try {
      return this.sketch?.warnings ?? [];
    } catch {
      return [];
    }
  }

  private async start(area: ReactionSketchArea): Promise<void> {
    try {
      // the reaction document and the reaction tools: crux-sketch's createSketcher options for reactions (passed as a
      // value, not a literal, so that this compiles against an older vendored crux-sketch too); Auto-map by Chem's
      // mapper (an option a sketcher without the Map tool ignores), written only on OK as any edit. The mapper runs on
      // Chem's Crux engine, which RDKit mode loads here (the editor works in both modes); one that does not load leaves
      // Auto-map mapping nothing
      const [module] = await Promise.all([loadCruxSketchForReactions(), cruxForReactions('the reaction editor')]);
      const options = {value: this.original, document: 'reaction', toolset: 'reaction',
        mapReaction: sketchReactionMapper() ?? undefined} as const;
      const sketch = await module.createSketcher(area.root, options);
      if (!this.dialog.root.isConnected) { // closed while Crux loaded: nothing to edit
        sketch.destroy();
        return;
      }
      // a reaction with no coordinates of its own opens arranged (reaction-arrange.ts): no person's edit
      if (this.original.trim() !== '' && !this.holdsCdxml)
        this.arrangedValue = arrangeOpened(sketch, this.original, cellReadFormat(this.cell));
      sketch.on('render', () => area.onRendered.next());
      sketch.on('change', ({source}) => {
        if (source !== 'user')
          return;
        area.changes++;
        this.edited = true;
        this.onChange.next(this.value);
      });
      this.sketch = area.sketch = sketch;
      const error = sketch.lastError;
      if (this.original.trim() !== '' && error !== null) {
        this.refused = plainRefusal(error.message);
        this.sayRefused(area, sketch);
      } else // the keyboard the drawing's as the dialog opens, as a click on it gives it (Ctrl+A included)
        focusDrawingSoon(area.root, this.dialog.root);
    } catch (e) {
      // a message a user sees names no engine; the error itself goes to the console
      console.error('Chem | the reaction editor did not start', e);
      area.root.append(ui.divText('The drawing area could not start. Cancel leaves the cell as it is.'));
    } finally {
      area.loading = false;
    }
  }

  // The cell refused ({@link refused}): said once, in the box in place of the sketcher, which is let go (not above
  // the drawing and again in the footer, over an empty sketcher that invited a drawing, "Draw reactants, a Reaction
  // arrow, then products"). OK is off; Cancel closes.
  private sayRefused(area: ReactionSketchArea, sketch: SketcherHandle): void {
    this.sketch = area.sketch = null;
    area.refused = true;
    sketch.destroy();
    area.root.querySelector(':scope > crux-sketch')?.remove();
    this.notice.textContent = refusedWords(this.refused!);
    Object.assign(this.notice.style, {display: 'flex', position: 'absolute', inset: '0', margin: '0',
      alignItems: 'center', justifyContent: 'center', textAlign: 'center', padding: '24px', maxWidth: 'none',
      color: 'var(--grey-6)', background: 'var(--grey-1)', fontSize: '14px'});
    area.root.append(this.notice);
    const ok = this.dialog.getButton('OK');
    if (ok)
      ok.disabled = true;
  }

  // OK: the reaction into the cell, once a person changed it into another text.
  private write(): void {
    if (!this.edited || this.slot?.blocked || this.refused !== null)
      return;
    const value = this.value;
    if (value === this.original || value === this.arrangedValue)
      return;
    try {
      this.cell.setValue(value, true);
    } catch {
      this.cell.cell.value = value;
    }
  }
}

/** How the editor reads `cell`'s value (a template that is valid SMILES too opened read as SMILES, its bracket atoms
 * without an H count drawn as radicals): as reaction SMARTS when its column holds templates (units `smarts`, as
 * detectors.js gives a column with SMARTS in it) or the value is written in SMARTS; else null, as Crux Sketch detects
 * it. */
export function cellReadFormat(cell: DG.GridCell): 'smarts' | null {
  const value: string = cell.cell.value ?? '';
  return cell.cell.column?.meta.units === 'smarts' || (value.trim() !== '' && reactionTextFormat(value) === 'smarts') ?
    'smarts' : null;
}

/** Opens the reaction cell editor on `cell`. */
export function editReactionCell(cell: DG.GridCell): ReactionCellEditor {
  return new ReactionCellEditor(cell);
}

/** The editor a reaction cell's double-click opens (Chem's cell editor), with the foot-gun footer slot (RXN-15/16/17,
 * TPL-12; in both modes), which OK honours too: its cell menu's "Edit reaction..." opens it as well. */
export async function openReactionCellEditor(cell: DG.GridCell): Promise<ReactionCellEditor> {
  const {watchReactionEditor} = await import('../footguns/editors');
  const editor = editReactionCell(cell);
  editor.slot = watchReactionEditor(editor);
  return editor;
}
