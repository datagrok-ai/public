/**
 * The Reaction Sketcher app (it was "Reactions", in Chem's Reactions folder of Apps; reactions in Crux): ONE reaction
 * or route (utils/reactions/reaction-document.ts). On the left the editing area, behind one small interface
 * ({@link EditingArea}): the Crux sketcher's reaction document (toolset `reaction`: every tool, the query tools among
 * them, as in a reaction cell's editor; the owner, 2026-11-17) with the reaction's text under
 * it, and what Crux reads under them ("Crux reads 5 steps.", or Crux's error), the title and the map numbers' switch.
 * On the right the reaction drawn large, and under it what the Reaction pane tells of it (widgets/reaction-panel.ts:
 * full screen, the steps, the template's notes, Map atoms). The ribbon's Examples menu loads one of nine titled
 * examples in place of the document. The app opens on the default example (a five-step mapped route) until a person
 * changes something, keeps every change in the browser and reopens with it; "Reset to example" brings the default
 * example back. The ribbon's Export... opens the Export panel on the reaction (widgets/reaction-export.ts: presets,
 * widths, SVG and PNG, previewed at print size). Under the status, a line names what the text a drawing wrote does not
 * keep of it (a label such as HATU written as "*", the drawing's texts; reactions/reaction-losses.ts).
 */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {ChemWasm, cruxModuleSync, loadCrux} from '../../engine/crux-module';
import {reactionPanel} from '../../widgets/reaction-panel';
import type {SetValueOptions, SketchWarning, SketcherHandle} from '../../../vendor/crux-sketch/crux-sketch';
import {loadCruxSketchForReactions} from '../../crux/crux-sketcher';
import {Flag, flagsText} from '../footguns/rules';
import {drawingLossFlags, drawingLosses} from './reaction-losses';
import {sketchReactionMapper} from './reaction-sketch-mapper';
import {DrawnStepTexts, StepTexts, stepTextsOf, tooLargeToDraw} from '../../rendering/crux-reaction-drawing';
import {ReactionSketchArea, plainRefusal} from './reaction-cell-editor';
import {ExportOptions, exportCellValue, openReactionExport} from '../../widgets/reaction-export';
import {reactionsHelpUrl} from './reactions-help';
import {AddPanel, AddRole, agentNames} from './reaction-add-panel';
import {DocumentReaction, NO_ARROW, PAGE_PADDING, ReactionPicture, TOO_LARGE, defaultExample, examples,
  forgetDocument, lastProduct, loadDocument, noReactionArrow, sameReactionText, saveDocument, shareAMolecule}
  from './reaction-document';
import {arrangeOpened, writeArranged} from './reaction-arrange';
import {isCdxmlText, sketcherPicture} from '../../rendering/cdxml-reaction-cells';

/** The app's name: its view's, and Apps' (package.ts reactionsApp). */
export const APP_NAME = 'Reaction Sketcher';

/** The editing area's pane, on the left: this share of the view, percent, never under {@link EDITOR_MIN_WIDTH} px
 * (the sketcher is where the work happens; at a fixed 480 px an arranged route's bonds were ~11 px).
 * The divider still moves. */
export const EDITOR_SHARE = 58;
export const EDITOR_MIN_WIDTH = 480;
/** The preview's least width, on the right (CSS px): below it the editing area gives way (with the help pane open the
 * preview was ~75 px, its title one word a line, the legend cut). */
export const PREVIEW_MIN_WIDTH = 280;
/** ui.splitH's divider (CSS px). */
const DIVIDER_WIDTH = 4;
const EDITOR_PANE_WIDTH = `max(${EDITOR_MIN_WIDTH}px, ${EDITOR_SHARE}%)`;

/** The split's two boxes held at their shares through any resize (ui.splitH turns its boxes' widths into pixels at
 * each one): the editing area {@link EDITOR_SHARE}% and never under {@link EDITOR_MIN_WIDTH} px, the preview never
 * under {@link PREVIEW_MIN_WIDTH} px; once the divider is dragged the person's widths stay, the preview's least width
 * kept. */
export function keepShares(split: HTMLElement, editor: HTMLElement, preview: HTMLElement): Subscription {
  const share = () => {
    editor.style.width = EDITOR_PANE_WIDTH;
    preview.style.width = `calc(100% - ${EDITOR_PANE_WIDTH} - ${DIVIDER_WIDTH}px)`;
  };
  editor.style.maxWidth = `calc(100% - ${PREVIEW_MIN_WIDTH + DIVIDER_WIDTH}px)`;
  preview.style.minWidth = `${PREVIEW_MIN_WIDTH}px`;
  share();
  let dragged = false;
  split.querySelector('.ui-split-h-divider')?.addEventListener('mousedown', () => dragged = true);
  return ui.onSizeChanged(split).subscribe(() => {
    if (!dragged)
      share();
  });
}
/** The drawing on the right is made again this long after the last keystroke (ms). */
const DRAWING_DELAY_MS = 400;
/** The longest the right pane keeps its picture while the next one is drawn (rather than going blank). */
const SWAP_MS = 1500;
/** A Text edit reaches the sketcher this long after the last keystroke, or at once when the field is left (ms): one
 * edit of the drawing a pause in typing (a Text replaced is one undo step). */
const TEXT_TO_DRAWING_MS = 500;
/** The notations the editing area reads, said under it. */
const NOTATIONS = 'Reaction SMILES or SMARTS (a CX part too), a route (A>>B>>C, steps separated by --**--, ' +
  'route JSON), an RXN block or an RD file.';
/** The exported files' names (the Export panel's, widgets/reaction-export.ts). */
export const SVG_FILE = 'reaction.svg';
export const PNG_FILE = 'reaction.png';
/** The ribbon's button that opens the Export panel. */
export const EXPORT_BUTTON = 'Export...';
/** The ribbon's button that puts the reaction into a table (there was no way to put a route into a table):
 * a row added to the open table that has a reaction column, else a new table of one, opened. A ribbon button beside
 * Export..., as the app's other commands are (Datagrok's apps keep their view's commands on the ribbon), named for
 * what it does; the table is Datagrok's own, the reaction in its column as text, drawn by its cells. */
export const ADD_TO_TABLE_BUTTON = 'Add to table';
const ADD_TO_TABLE_TOOLTIP = 'The reaction as a row of the open table with a reaction column, or of a new table';
/** The new table's name and its columns: the reaction (its text: a route's steps separated by --**--, Chem's route
 * notation) and its title. */
export const NEW_TABLE_NAME = 'Reactions';
export const REACTION_COLUMN = 'reaction';
export const TITLE_COLUMN = 'title';
/** What Add to table says with no reaction to add. */
export const NOTHING_TO_ADD = 'There is no reaction to add: draw or type one first';

/** The units Chem's detector gives a reaction column of `text` (detectors.js reactionUnits): a route (steps joined by
 * --**--, or a chain A>>B>>C), an RXN block, SMARTS or SMILES. */
export function reactionUnitsOf(text: string): string {
  if (text.includes('--**--') || />>[^>]*>>/.test(text) || /^\s*(\{|\[\s*\{)/.test(text)) // route JSON too
    return 'route';
  if (/^\s*\$(RXN|RDFILE)/.test(text))
    return 'rxn';
  return /[&;!$~]|\[#\d/.test(text) ? 'smarts' : 'smiles';
}

/** The open table a reaction is added to, and its reaction column (not a column of CDXML pages): the current table when
 * it has one, else the last table view's that has one; null when no open table has one. */
export function reactionTable(): {table: DG.DataFrame, column: DG.Column} | null {
  const columnOf = (table: DG.DataFrame | null | undefined) => table?.columns.bySemTypeAll(DG.SEMTYPE.CHEMICAL_REACTION)
    .find((c) => c.meta.units !== 'cdxml') ?? null;
  const tables = [grok.shell.t, ...Array.from(grok.shell.tableViews).map((v) => v.dataFrame).reverse()];
  for (const table of tables) {
    const column = columnOf(table);
    if (table && column)
      return {table, column};
  }
  return null;
}

/** Reaction `text`, titled `title`, put into a table: a row added to {@link reactionTable}'s (its title column set too
 * when it has one), or a new table ({@link NEW_TABLE_NAME}: the reaction and, when titled, the title), opened. What it
 * was added to, and the row (from 0). */
export async function addToTable(text: string, title: string): Promise<{table: DG.DataFrame, row: number}> {
  const open = reactionTable();
  if (open) {
    const {table, column} = open;
    table.rows.addNew();
    const row = table.rowCount - 1;
    column.set(row, text);
    const titles = table.col(TITLE_COLUMN);
    if (titles && titles.type === DG.TYPE.STRING)
      titles.set(row, title);
    grok.shell.info(`Added to ${table.name} as row ${row + 1}`);
    return {table, row};
  }
  const columns = [DG.Column.fromStrings(REACTION_COLUMN, [text])];
  if (title)
    columns.push(DG.Column.fromStrings(TITLE_COLUMN, [title]));
  const table = DG.DataFrame.fromColumns(columns);
  table.name = NEW_TABLE_NAME;
  await grok.data.detectSemanticTypes(table);
  const column = table.col(REACTION_COLUMN)!;
  if (column.semType !== DG.SEMTYPE.CHEMICAL_REACTION) { // as Chem's detector makes it
    column.semType = DG.SEMTYPE.CHEMICAL_REACTION;
    column.meta.units = reactionUnitsOf(text);
    column.meta.cellRenderer = 'ChemicalReaction';
  }
  grok.shell.addTableView(table);
  return {table, row: 0};
}

/** What the Title box says while the document has no title (an example's dropped once its reaction is changed). */
export const UNTITLED = 'Untitled';
/** The ribbon menu that lists the examples. */
export const EXAMPLES_MENU = 'Examples';
/** The class of the preview's row of its own controls, at its top: the Map numbers switch (off by default, and in
 * sight without scrolling; it sat under the sketcher, the Add panel and the Text). */
export const PREVIEW_CONTROLS_CLASS = 'chem-reactions-app-preview-controls';
/** The class of the line under the status saying what the text does not keep of the drawing (what tests find it by). */
export const LOSSES_CLASS = 'chem-reactions-app-losses';
/** The class of the line over the drawing area while it does not hold the document's text (a text it cannot open:
 * the app reopened on a stored route the drawing refused showed its canvas empty under the route drawn on the
 * right). */
export const SKETCHER_REFUSED_CLASS = 'chem-reactions-app-sketcher-refused';

/** What the drawing area says of a text it cannot open, for `reason` (the sketcher's, in plain words). */
export function sketcherRefusedWords(reason: string): string {
  return `The drawing cannot open this text: ${reason.replace(/[.\s]+$/, '')}. The text below is the reaction, ` +
    'drawn on the right and exported as it is: correct it there, or draw here to start a new one.';
}

/** The class of the preview's line for a text with no reaction arrow (a lone molecule). */
export const NO_ARROW_CLASS = 'chem-reactions-app-no-arrow';
/** The longest the drawing area keeps its cover up while a reopened document's kept drawing is read (ms; the cover
 * once lifted before it was read, so a click on the empty sketcher meanwhile became the document; a CDXML module that
 * never loads must not keep it covered): past it the text is drawn, as for a kept drawing that does not read. */
export const PAGE_READ_MS = 8000;
/** What the status adds when the document's kept drawing was not read (it did not read, or not within
 * {@link PAGE_READ_MS}): the sketcher draws the text, which holds less. */
export const DRAWING_UNREAD = 'Its kept drawing could not be read: it is drawn from its text.';

/** What the app says when its sketcher or its engine does not start (a message a user sees names no engine and dumps
 * no error text; the error goes to the console). */
/** What the Add panel says when the drawing has not started yet. */
export const ADD_NOT_READY = 'The drawing area is not ready yet';
export const EDITOR_FAILED = 'The drawing area could not start. The reaction text below still works.';
export const ENGINE_FAILED = 'Reactions cannot be drawn here: the drawing engine did not load.';

/** The status line under the editing area: what the text reads as (no engine named). */
export function readsAs(steps: number): string {
  return steps === 1 ? 'One reaction.' : `A route of ${steps} steps.`;
}

/** What the Export panel takes from the editing area's drawing ({@link EditingArea.drawingExport}). */
export type DrawingExport = Pick<ExportOptions, 'cdxml' | 'textsLeftOut' | 'stepTexts' | 'asDrawn'>;

/**
 * The app's editing area: what edits the reaction's text, and tells when a person changed it. The app's is the Crux
 * sketcher's reaction document, toolset `reaction` ({@link SketcherEditingArea}), with a text box under "Text"
 * ({@link TextEditingArea}, an editing area of its own).
 */
export interface EditingArea {
  readonly root: HTMLElement;
  /** The reaction's text; setting it is no person's change, and calls no listener. */
  value: string;
  /** `listener` is called with the text after each change a person makes. */
  onChange(listener: (text: string) => void): void;
  /** What the text, as a drawing changed by a person wrote it, does not keep of the drawing (a label written as "*",
   * the texts, the arrows' styles; reactions/reaction-losses.ts): none for a text typed or loaded. */
  readonly losses: readonly Flag[];
  focus(): void;
  /** What the Export panel takes from the drawing itself: its CDXML as drawn, whether it has texts that the pictures,
   * drawn from the text, leave out, its steps' texts, which the pictures draw, and its picture as drawn (the
   * panel's "As drawn"). None for a text box. */
  drawingExport(): DrawingExport;
  /** The document's `text` shown, as setting {@link value} does, its drawing opened from `drawing` (a CDXML page the
   * app kept: DocumentReaction.drawing) when there is one and it reads, else drawn from the text. */
  open(text: string, drawing?: string): void;
  /** `listener` is called after a person's change of the drawing with its CDXML page while the drawing holds what its
   * text does not keep (reaction-losses.ts), else with null. Never for a text box. */
  onDrawing(listener: (drawing: string | null) => void): void;
  /** `listener` is called once the drawing area has opened a document's kept drawing ({@link open}'s `drawing`), which
   * it reads after its CDXML module is in: what the app shows of the drawing ({@link losses}) follows it then (else a
   * document reopened on its kept drawing said nothing of what its text drops). Never for a text box. */
  onShown(listener: () => void): void;
  /** The document's kept drawing was not read (it did not read, or not within {@link PAGE_READ_MS}): its text is
   * drawn. False for a text box, and once a text is typed, a document loaded or the drawing edited. */
  readonly drawingUnread: boolean;
  /** Lets go of what it holds (its listeners; the sketcher's engine), as the app closes. */
  dispose(): void;
}

/** The editing area as a text box (monospace, no spelling check), the notations it reads said under it. */
class TextEditingArea implements EditingArea {
  private readonly input = ui.input.textArea('Reaction', {onValueChanged: (v) => this.changed(v ?? '')});
  private readonly listeners: ((text: string) => void)[] = [];
  // True while the text is set from outside, so that the input's change event is no person's.
  private setting = false;
  readonly root: HTMLElement;
  /** A text box drops nothing: what is typed is the reaction. */
  readonly losses: readonly Flag[] = [];

  constructor() {
    const box = this.input.input as HTMLTextAreaElement;
    box.classList.add('chem-reactions-app-text');
    Object.assign(box.style, {fontFamily: 'monospace', fontSize: '12px'});
    box.rows = 12;
    box.setAttribute('spellcheck', 'false');
    const help = ui.divText(NOTATIONS, {style: {color: 'var(--grey-5)', fontSize: '11px'}});
    this.root = ui.divV([ui.form([this.input]), help]);
  }

  get value(): string {
    return this.input.value ?? '';
  }

  set value(text: string) {
    this.setting = true;
    try {
      this.input.value = text;
    } finally {
      this.setting = false;
    }
  }

  onChange(listener: (text: string) => void): void {
    this.listeners.push(listener);
  }

  drawingExport(): DrawingExport {
    return {};
  }

  open(text: string): void {
    this.value = text;
  }

  onDrawing(): void {} // a text box has no drawing

  onShown(): void {} // nor a drawing of its own to open

  readonly drawingUnread = false;

  focus(): void {
    this.input.input.focus();
  }

  /** The text box itself (its focus leaving it draws what was typed at once). */
  get box(): HTMLTextAreaElement {
    return this.input.input as HTMLTextAreaElement;
  }

  dispose(): void {
    this.listeners.length = 0;
  }

  private changed(text: string): void {
    if (!this.setting)
      this.listeners.forEach((l) => l(text));
  }
}

/** The text of what `sketch` holds after a person's edit: its value, the reaction in the notation read; a drawing with
 * no arrow yet, which Crux Sketch's value gives as a MOL block, as its SMILES (else drawing one
 * molecule turned the Text into a MOL file and the preview into an error dump). */
function textOf(sketch: SketcherHandle): string {
  const value = sketch.value;
  // a drawing opened from the document's kept CDXML page (DocumentReaction.drawing) writes its reaction as text
  if (isCdxmlText(value)) {
    try {
      return sketch.smiles || sketch.smarts || '';
    } catch {
      return '';
    }
  }
  if (!value.includes('M  END') || value.trimStart().startsWith('$RXN'))
    return value;
  try {
    return sketch.smiles || sketch.smarts || value; // a query has no SMILES: its SMARTS
  } catch {
    return value; // neither written: the block, which the preview still says needs an arrow
  }
}

/** The drawing opened arranged (`value`, or null when it was not arranged by an edit) and the text it opened on. */
function withText(value: string | null, text: string): {value: string, text: string} | null {
  return value === null ? null : {value, text};
}

/** The main thread's agent table (chem-wasm `chemAgentTable()`, a JSON of SMILES to names), null before the module
 * is loaded or when its build has none. */
function agentTableJson(): string | null {
  const crux = cruxModuleSync() as unknown as {chemAgentTable?: () => string} | null;
  try {
    return crux?.chemAgentTable?.() ?? null;
  } catch (e) {
    console.warn('Chem | the agent table did not load', e);
    return null;
  }
}

/** The sketcher's drawing as a CDXML page (Crux Sketch's writer, its module loaded first), or null. */
async function drawingCdxml(sketch: SketcherHandle): Promise<string | null> {
  try {
    await sketch.loadFormat('cdxml');
    return sketch.cdxml || null;
  } catch (e) {
    console.warn('Chem | the drawing was not written as CDXML', e);
    return null;
  }
}

/** The texts over and under each step's arrow in `sketch`'s drawing, as Crux's `stepTexts` (crux-sketch's
 * `reactionSteps`: a sketcher without it gives none), or null when there are none. */
function stepTextsOfSketch(sketch: SketcherHandle): StepTexts | null {
  try {
    return stepTextsOf((sketch as unknown as {reactionSteps?: readonly DrawnStepTexts[] | null}).reactionSteps);
  } catch (e) {
    console.warn('Chem | the drawing\'s step texts were not read', e);
    return null;
  }
}

/** What `sketch` says its drawing loses in each notation (its `warnings`), or none when it cannot say. */
function warningsOf(sketch: SketcherHandle): readonly SketchWarning[] {
  try {
    return sketch.warnings;
  } catch {
    return [];
  }
}

/** The height of the Crux sketcher's box in the editing area in which its left bar shows every group (the text tool,
 * the arrow styles, the map tool among them) with no scroll button, in either toolset (CSS px; crux-sketch d0b8bb9:
 * the publication tools need 690, the gear's reaction tools, one tool more (the query bond), 712: at 690 their bar
 * scrolled, 640 px of tools in 618; the cell editor's box is 716 too. At 580, on an older sketcher's bars, Text and
 * Reaction map were hidden behind the bar's scroll button at 1440x900). */
export const SKETCHER_FULL_HEIGHT = 716;
/** How far down the window the box starts (the shell's bar, the view's ribbon, the Title row: 76 px measured), and the
 * room kept free under it (CSS px). */
const SKETCHER_TOP = 76;
const SKETCHER_MARGIN = 24;
/** The least height the box is given in a window too short for {@link SKETCHER_FULL_HEIGHT}: its bars then scroll. */
const SKETCHER_MIN_HEIGHT = 320;

/** The box's height in a window `windowHeight` px tall: tall enough for the left bar's every group, unless the view is
 * shorter than that below the box's top; then the view's height there, the bar scrolling (its scroll buttons). */
export function sketcherBoxHeight(windowHeight: number = window.innerHeight): number {
  return Math.max(SKETCHER_MIN_HEIGHT, Math.min(SKETCHER_FULL_HEIGHT, windowHeight - SKETCHER_TOP - SKETCHER_MARGIN));
}

/**
 * The editing area as the Crux sketcher: a reaction document with the publication tools (the standalone app's
 * default; the gear switches to the reaction tools, never stored), the text box under it, under
 * "Text", each following the other: a drawing changed writes its text into the box, a text typed is drawn.
 */
class SketcherEditingArea implements EditingArea {
  readonly root: HTMLElement;
  private readonly host = ui.div([], 'chem-reactions-app-sketcher');
  // The host as the "crux sketcher widget" the cell editor's box is: its atoms' areas and readings, for the tests.
  private readonly area: ReactionSketchArea;
  private readonly text = new TextEditingArea();
  private readonly listeners: ((text: string) => void)[] = [];
  private sketch: SketcherHandle | null = null;
  private current = '';
  private disposed = false;
  // The drawing as a document opened arranged by an Arrange edit (its value) and the text it opened on: a person's
  // edits undone back to it give that text again, as written (reaction-arrange.ts). Null otherwise.
  private opened: {value: string, text: string} | null = null;
  /** What the text the drawing last wrote does not keep of it; emptied when a text is typed or set. */
  losses: readonly Flag[] = [];
  // A Text edit waiting to be drawn ({@link TEXT_TO_DRAWING_MS}).
  private textTimer: ReturnType<typeof setTimeout> | undefined;
  // The text is too large for Crux to draw in time (the cells' rule, `tooLargeToDraw`): the sketcher is left empty,
  // says so, and a drawing made in it meanwhile is not taken as the text.
  private tooLarge = false;
  /** Why the sketcher did not open the text (its `lastError`, in plain words): it is left empty and says so, the
   * text staying the document; a drawing made in it then takes the text's place. Null when it holds the text. */
  refused: string | null = null;
  private readonly notice = ui.divText('', 'chem-reactions-app-sketcher-too-large');
  /** The Add panel: reactant, reagent and product slots, each addition one edit of the drawing. */
  readonly addPanel: AddPanel;
  // The kept drawing (a CDXML page) the document opens on in place of its text, or null: the text is drawn.
  private page: string | null = null;
  // Counts the drawing's and the Text's changes: a drawing's CDXML written for a change that a later one overtook is
  // not kept.
  private generation = 0;
  private readonly drawingListeners: ((drawing: string | null) => void)[] = [];
  private readonly shownListeners: (() => void)[] = [];
  // The kept drawing being read ({@link showPage}: its CDXML module loading), by its turn: the sketcher does not hold
  // the document meanwhile. 0 when none is.
  private opening = 0;
  private openings = 0;
  // The read of the document's kept drawing {@link show} started, or null: the cover stays up while it lasts.
  private pageRead: Promise<void> | null = null;
  drawingUnread = false;

  /** `same(a, b)`: whether two reaction texts are one reaction, so that a drawing changed only in its layout keeps the
   * Text as written; `lastOf(text)`: the SMILES of its route's last product, which the Add panel's Next step draws
   * from, or null. */
  constructor(private readonly same: (a: string, b: string) => boolean = (a, b) => a === b,
    private readonly lastOf: (text: string) => string | null = () => null) {
    Object.assign(this.host.style,
      {height: `${sketcherBoxHeight()}px`, position: 'relative', border: '1px solid var(--grey-2)'});
    this.area = new ReactionSketchArea(this.host);
    Object.assign(this.notice.style, {display: 'none', position: 'absolute', left: '48px', right: '48px', top: '30%',
      textAlign: 'center', color: 'var(--grey-6)', pointerEvents: 'none', zIndex: '1', padding: '8px 12px',
      background: 'rgba(255, 255, 255, 0.92)'});
    this.addPanel = new AddPanel({add: (t, role) => this.add(t, role),
      lastProduct: () => this.tooLarge ? null : this.lastOf(this.current)}, () => agentNames(agentTableJson()));
    this.root = ui.divV([ui.div([this.host, this.notice], {style: {position: 'relative'}}), this.addPanel.root,
      ui.h3('Text'), this.text.root]);
    this.text.onChange((t) => {
      this.current = t;
      this.losses = [];
      this.drawingUnread = false;
      this.page = null; // a text typed is the document: its drawing is the text's
      this.generation++;
      this.listeners.forEach((l) => l(t));
      // the drawing follows a pause in typing, as one edit (a host's change: no change event comes back as a person's)
      clearTimeout(this.textTimer);
      this.textTimer = setTimeout(() => this.showTyped(), TEXT_TO_DRAWING_MS);
    });
    this.text.box.addEventListener('blur', () => this.showTyped());
    this.start();
  }

  // The Text typed and not drawn yet: into the sketcher now, as one edit its undo takes back.
  private showTyped(): void {
    if (this.textTimer === undefined)
      return;
    clearTimeout(this.textTimer);
    this.textTimer = undefined;
    this.show(this.current, true);
  }

  // `text` into the sketcher, unless it is too large to draw in time: then the sketcher is emptied and says so. A
  // text typed (`edit`) is one edit of the drawing, Ctrl+Z taking it back (crux-sketch's `history: 'edit'`, which a
  // sketcher before it ignores, starting a new history as any value written does); a document loaded starts anew.
  private show(text: string, edit = false): void {
    this.tooLarge = tooLargeToDraw(text);
    this.refused = null;
    this.pageRead = null;
    this.say();
    if (!this.sketch)
      return;
    if (!edit && !this.tooLarge && this.page !== null) { // the document's kept drawing, read in place of its text
      this.pageRead = this.showPage(this.sketch, this.page, text);
      return;
    }
    const options: SetValueOptions & {history?: 'edit'} = edit && !this.tooLarge ? {history: 'edit'} : {};
    if (edit || this.tooLarge)
      this.sketch.setValue(this.tooLarge ? '' : text, options);
    else // a document loaded: a reaction with no coordinates of its own drawn arranged
      this.opened = withText(writeArranged(this.sketch, text, options), text);
    this.checkOpened(this.sketch, text);
  }

  // The document's kept drawing `page` (CDXML) into `sketch`, its text `text` staying the Text: a catalytic cycle, the
  // texts, the arrows' styles come back as drawn (not a cycle reopened as a straight route). A page that does not
  // read, or that the area left meanwhile (another document, a text typed), is not drawn: the text is, as before.
  private async showPage(sketch: SketcherHandle, page: string, text: string): Promise<void> {
    const turn = this.opening = ++this.openings;
    try {
      await this.readPage(sketch, page, text);
    } finally {
      if (this.opening === turn)
        this.opening = 0;
    }
    if (turn === this.openings && !this.disposed) // what the app shows of the drawing follows it
      this.shownListeners.forEach((l) => l());
  }

  private async readPage(sketch: SketcherHandle, page: string, text: string): Promise<void> {
    const current = () => !this.disposed && this.sketch === sketch && this.page === page && this.current === text;
    try {
      await sketch.loadFormat('cdxml');
    } catch (e) {
      console.warn('Chem | the kept drawing could not be read; the text is drawn', e);
    }
    if (!current())
      return;
    sketch.setValue(page, {format: 'cdxml'});
    if (sketch.lastError !== null || sketch.isEmpty) {
      console.warn('Chem | the kept drawing did not read; the text is drawn', sketch.lastError);
      this.page = null;
      this.drawingUnread = true;
      this.opened = withText(writeArranged(sketch, text, {}), text);
    } else {
      this.opened = {value: sketch.value, text};
      // what its text does not keep of the drawing reopened: the losses a person's edit of it would name
      const lost = text.trim() === '' ? null : drawingLosses(warningsOf(sketch), text);
      this.losses = lost ? drawingLossFlags(lost) : [];
    }
    this.checkOpened(sketch, text);
  }

  // Whether `sketch` opened `text`: one it refused (its lastError, its drawing empty) is said over it ({@link
  // refused}), as one too large is.
  private checkOpened(sketch: SketcherHandle, text: string): void {
    const error = this.tooLarge || text.trim() === '' ? null : sketch.lastError;
    this.refused = error !== null && (sketch.value ?? '').trim() === '' ? plainRefusal(error.message) : null;
    this.say();
  }

  // The line over the drawing area: the text too large to draw, or refused, or none.
  private say(): void {
    this.notice.textContent = this.tooLarge ? TOO_LARGE + ' here: edit its text below.' :
      this.refused !== null ? sketcherRefusedWords(this.refused) : '';
    this.notice.classList.toggle(SKETCHER_REFUSED_CLASS, !this.tooLarge && this.refused !== null);
    this.notice.style.display = this.tooLarge || this.refused !== null ? '' : 'none';
  }

  private async start(): Promise<void> {
    try {
      const module = await loadCruxSketchForReactions();
      const shown = this.current;
      // Auto-map by Chem's mapper (crux-sketch's reaction-map option, which a sketcher before it ignores: passed as
      // a value, not a literal, so that this compiles against it); a kept drawing is read once the sketcher is up. The
      // reaction tools, query tools included (query mode): a reaction drawn here may be a template to run, as a cell's
      // (the owner, 2026-11-17: the app opened in molecule mode, its Atom Properties without the query features); the
      // gear's choice still switches to the publication tools, and neither changes the reaction's text
      const options = {value: tooLargeToDraw(shown) || this.page !== null ? '' : shown, document: 'reaction',
        toolset: 'reaction', mapReaction: sketchReactionMapper() ?? undefined} as const;
      const sketch = await module.createSketcher(this.host, options);
      if (this.disposed) {
        sketch.destroy();
        return;
      }
      this.sketch = this.area.sketch = sketch;
      if (this.current !== shown || this.page !== null) // a text set while the sketcher loaded, or a kept drawing
        this.show(this.current);
      else if (!tooLargeToDraw(shown)) { // the document opened arranged (reaction-arrange.ts)
        this.opened = withText(arrangeOpened(sketch, shown), shown);
        this.checkOpened(sketch, shown);
      }
      sketch.on('render', () => this.area.onRendered.next());
      sketch.on('change', ({source}) => {
        if (source !== 'user')
          return;
        this.area.changes++;
        this.drawingChanged(sketch);
      });
    } catch (e) {
      // a message a user sees names no engine; the error itself goes to the console
      console.error('Chem | the Reaction Sketcher sketcher did not start', e);
      this.host.append(ui.divText(EDITOR_FAILED));
    } finally {
      const reading = this.pageRead;
      if (reading) // a reopened document's kept drawing: under the cover until it is read
        void this.coverWhile(reading);
      else
        this.area.loading = false;
    }
  }

  // The cover kept up while `reading` (a kept drawing's read) lasts, at most {@link PAGE_READ_MS}: past it the read
  // is let go (it draws nothing when it ends), the text is drawn, and the app hears the drawing went unread.
  private async coverWhile(reading: Promise<void>): Promise<void> {
    let timer: ReturnType<typeof setTimeout> | undefined;
    const late = await Promise.race([reading.then(() => false, () => false),
      new Promise<boolean>((resolve) => timer = setTimeout(() => resolve(true), PAGE_READ_MS))]);
    clearTimeout(timer);
    if (late && !this.disposed && this.sketch !== null && this.pageRead === reading) {
      console.warn(`Chem | the kept drawing was not read within ${PAGE_READ_MS} ms; the text is drawn`);
      this.openings++; // the late read's turn is over: it opens nothing and tells no one
      this.opening = 0;
      this.page = null;
      this.show(this.current);
      this.drawingUnread = true;
      this.shownListeners.forEach((l) => l());
    }
    if (!this.disposed)
      this.area.loading = false;
  }

  // The drawing changed (a person's edit, or an addition from the Add panel): its text into Text, unless it is the
  // same reaction as Text (a move, Arrange, Clean Up keep the Text as the person wrote it).
  private drawingChanged(sketch: SketcherHandle): void {
    if (this.tooLarge)
      return;
    if (this.refused !== null) { // a drawing in place of the text it did not open: the drawing is the reaction now
      this.refused = null;
      this.say();
    }
    this.drawingUnread = false; // the drawing is the person's now
    // a drawing edit is newer than a Text not drawn yet: the drawing wins
    clearTimeout(this.textTimer);
    this.textTimer = undefined;
    const next = this.opened !== null && sketch.value === this.opened.value ? this.opened.text : textOf(sketch);
    if (!this.same(this.current, next))
      this.current = next;
    this.text.value = this.current;
    const lost = this.current.trim() === '' ? null : drawingLosses(warningsOf(sketch), this.current);
    this.losses = lost ? drawingLossFlags(lost) : [];
    this.listeners.forEach((l) => l(this.current));
    this.keepDrawing(sketch, !!lost && (lost.objects.length > 0 || lost.labels > 0));
  }

  // The drawing kept beside the text while it holds what the text drops (`holds`: a cycle, texts, arrow styles, a
  // label), as a CDXML page written once its module is in; else none: the listeners told either way, a page overtaken
  // by a later change dropped.
  private keepDrawing(sketch: SketcherHandle, holds: boolean): void {
    const generation = ++this.generation;
    if (!holds) {
      this.page = null;
      this.drawingListeners.forEach((l) => l(null));
      return;
    }
    void drawingCdxml(sketch).then((page) => {
      if (generation !== this.generation || this.disposed || this.sketch !== sketch)
        return;
      this.page = page;
      this.drawingListeners.forEach((l) => l(page));
    });
  }

  // The Add panel's addition: `text` into the drawing as `role` (crux-sketch's addMolecule: one edit, one undo
  // step), then its text into Text; null when added, else why not (the sketcher's reason, in plain words).
  private add(text: string, role: AddRole): string | null {
    const sketch = this.sketch;
    if (!sketch || this.tooLarge)
      return sketch ? TOO_LARGE : ADD_NOT_READY;
    // crux-sketch's lastError keeps an earlier refusal after a later addition succeeds (every addition after a refused
    // HATU was taken as refused, the Text never told): refused only when it leaves a new error
    const before = sketch.lastError;
    sketch.addMolecule(text, {role});
    const error = sketch.lastError;
    if (error !== null && error !== before)
      return plainRefusal(error.message);
    this.area.changes++;
    this.drawingChanged(sketch);
    return null;
  }

  get value(): string {
    return this.current;
  }

  set value(text: string) {
    this.open(text);
  }

  open(text: string, drawing?: string): void {
    clearTimeout(this.textTimer);
    this.textTimer = undefined;
    this.current = text;
    this.losses = [];
    this.drawingUnread = false;
    this.text.value = text;
    this.page = drawing?.trim() ? drawing : null;
    this.generation++;
    this.show(text);
  }

  onChange(listener: (text: string) => void): void {
    this.listeners.push(listener);
  }

  onDrawing(listener: (drawing: string | null) => void): void {
    this.drawingListeners.push(listener);
  }

  onShown(listener: () => void): void {
    this.shownListeners.push(listener);
  }

  focus(): void {
    this.host.querySelector<HTMLElement>(':scope > crux-sketch')?.focus(); // the box holds the loading cover too
  }

  drawingExport(): DrawingExport {
    // the drawing exported only while it holds the text (not one too large, or refused: the text is exported); nor
    // while a Text typed has not reached it yet, or a kept drawing is being read (else the right pane, drawn 400 ms
    // after a keystroke, before the Text's 500 ms, drew the previous drawing's texts over the new text's arrows)
    const sketch = this.tooLarge || this.refused !== null || this.textTimer !== undefined || this.opening !== 0 ?
      null : this.sketch;
    if (!sketch)
      return {};
    return {cdxml: () => drawingCdxml(sketch),
      textsLeftOut: warningsOf(sketch).some((w) => w.kind === 'format' && w.code === 'texts'),
      stepTexts: stepTextsOfSketch(sketch), asDrawn: sketcherPicture(sketch)};
  }

  dispose(): void {
    this.disposed = true;
    clearTimeout(this.textTimer);
    this.listeners.length = 0;
    this.drawingListeners.length = 0;
    this.shownListeners.length = 0;
    this.text.dispose();
    this.sketch?.destroy();
    this.sketch = this.area.sketch = null;
    this.area.onRendered.complete();
    this.area.detach();
  }
}

/** The shell's panels beside the view that the app shuts as it opens, as other Datagrok apps do (SequenceTranslator's
 * apps, Peptides, BiostructureViewer): the context panel and the help panel. At 1440x900 with both open
 * and Browse, the view was 667 px wide, the sketcher 359 px and the opened route's bonds 8 px; with both shut,
 * 1137 px, 635 px and 15 px. A person opens them again as ever; the app does not put them back. */
export function shutSidePanels(): void {
  grok.shell.windows.showContextPanel = false;
  grok.shell.windows.showHelp = false;
}

/** The app's view: the editing area on the left, the reaction drawn and told of on the right. */
export async function reactionsAppView(): Promise<DG.View> {
  shutSidePanels();
  const view = DG.View.create();
  view.name = APP_NAME;
  view.box = true;
  view.helpUrl = reactionsHelpUrl('the-reaction-sketcher');
  try {
    // Crux draws every picture, in either engine mode
    new ReactionsApp(view, await loadCrux('reactions app'), loadDocument() ?? defaultExample());
  } catch (e) {
    console.error('Chem | the Reaction Sketcher engine did not load', e);
    view.append(ui.divText(ENGINE_FAILED));
  }
  return view;
}

/** The app on its view: the document, its picture, and the elements that show them. */
class ReactionsApp {
  private picture: ReactionPicture;
  private readonly area: EditingArea = new SketcherEditingArea((a, b) => sameReactionText(this.crux, a, b),
    (text) => lastProduct(this.crux, text));
  private readonly drawing = ui.divV([], 'chem-reactions-app-drawing');
  // True while the inputs are set from the document, so that their change events edit nothing.
  private filling = false;
  private drawingTimer: ReturnType<typeof setTimeout> | undefined;
  // The Reaction pane shown on the right, detached (its pictures' subscriptions ended) when replaced or closed.
  private panel: DG.Widget | null = null;
  // The right pane's page shown (its title, then the Reaction pane or a line), or null before the first.
  private page: HTMLElement | null = null;
  // The page drawn to replace it, hidden over it until its pictures are drawn (else the pane goes blank as each
  // picture is made again), and its wait; null when none is.
  private next: {page: HTMLElement, panel: DG.Widget | null, timer: ReturnType<typeof setTimeout>} | null = null;
  // A button held down in the editing area (a drag in the drawing): the right pane is drawn again once it is let go,
  // the drawing settled.
  private pressed = false;
  private drawWaits = false;
  // The example's title the app dropped (the Text emptied, or a reaction with none of its molecules), put back when
  // its molecules come back (an undo); null once a person sets a title or a document is loaded.
  private dropped: string | null = null;
  // The step texts the right pane was last drawn with (JSON): a drawing whose texts alone changed is drawn again.
  private drawnTexts = 'null';

  private readonly titleInput = ui.input.string('Title', {onValueChanged: (v) => this.edit({title: v ?? ''})});
  private readonly mapInput = ui.input.bool('Map numbers', {onValueChanged: (v) => this.edit({mapNumbers: !!v})});
  private readonly status = ui.divText('', 'chem-reactions-app-status');
  // Under the status: what the text the drawing wrote does not keep of it, one line a loss.
  private readonly losses = ui.divText('', LOSSES_CLASS);

  constructor(private readonly view: DG.View, private readonly crux: ChemWasm, private doc: DocumentReaction) {
    this.picture = new ReactionPicture(crux, doc);
    this.area.onChange((text) => this.edit({text}));
    this.area.onDrawing((drawing) => this.keep(drawing ?? undefined));
    this.area.onShown(() => this.drawingShown());
    Object.assign(this.losses.style, {color: 'var(--orange-2)', whiteSpace: 'pre-line', display: 'none'});
    (this.titleInput.input as HTMLInputElement).placeholder = UNTITLED;
    const editor = ui.divV([ui.form([this.titleInput]), this.area.root, this.status, this.losses],
      'chem-reactions-app-editor');
    const left = ui.div([editor], {style: {overflow: 'auto', height: '100%', padding: '8px 12px',
      boxSizing: 'border-box'}});
    Object.assign(this.drawing.style, {background: 'white', padding: `${PAGE_PADDING}px`, boxSizing: 'border-box',
      minHeight: '100%', position: 'relative'});
    // the preview's own controls first, over the pages swapped in under them (swapIn): the Map numbers switch, right
    // of its title's line as the Reaction pane's icons are under it
    const controls = ui.divH([this.mapInput.root], PREVIEW_CONTROLS_CLASS);
    Object.assign(controls.style, {justifyContent: 'flex-end', margin: '0 0 4px 0'});
    this.drawing.append(controls);
    const right = ui.div([this.drawing], {style: {overflow: 'auto', height: '100%', background: 'var(--grey-1)'}});
    const split = ui.splitH([left, right], {style: {width: '100%', height: '100%'}}, true);
    view.append(split);
    const shares = keepShares(split, left.parentElement!, right.parentElement!); // the split's boxes around them
    const menu = view.ribbonMenu.group(EXAMPLES_MENU);
    for (const example of examples())
      menu.item(example.title, () => this.load(example));
    view.setRibbonPanels([[
      ui.button('Reset to example', () => this.confirmReset(), 'The default example again; your changes are lost'),
    ], [
      // the Export panel (in place of two one-click buttons, with no choice of width, dpi or style)
      ui.button(EXPORT_BUTTON, () => this.openExport(), 'The reaction as SVG or PNG, previewed at print size'),
      ui.button(ADD_TO_TABLE_BUTTON, () => this.putInTable(), ADD_TO_TABLE_TOOLTIP),
    ]]);
    // closed: nothing left behind (the old app kept its schemes, nodes and listeners). The platform ends no plain
    // view's `subs` when it closes (View.close() does not detach it; one scheme was left per open), so the app lets go
    // when the platform says this view was removed, or on a detach, whichever comes first.
    // a drag in the drawing: the right pane keeps its picture until the button is let go
    const press = () => this.pressed = true;
    const release = () => {
      if (!this.pressed)
        return;
      this.pressed = false;
      if (this.drawWaits) {
        this.drawWaits = false;
        this.drawSoon();
      }
    };
    this.area.root.addEventListener('pointerdown', press, true);
    window.addEventListener('pointerup', release, true);
    window.addEventListener('pointercancel', release, true);
    const letGo = new Subscription(() => {
      clearTimeout(this.drawingTimer);
      window.removeEventListener('pointerup', release, true);
      window.removeEventListener('pointercancel', release, true);
      this.dropNext();
      this.picture.free();
      this.panel?.detach();
      this.panel = null;
      this.area.dispose();
    });
    view.subs.push(letGo);
    letGo.add(shares);
    letGo.add(grok.events.onViewRemoved.subscribe((removed) => {
      if (removed?.root === view.root)
        letGo.unsubscribe();
    }));
    this.fill();
  }

  // The inputs, the status and the drawing set from the document.
  private fill(): void {
    this.filling = true;
    try {
      this.titleInput.value = this.doc.title;
      this.area.open(this.doc.text, this.doc.drawing);
      this.mapInput.value = this.doc.mapNumbers;
    } finally {
      this.filling = false;
    }
    this.showStatus();
    this.draw();
  }

  // What Crux reads from the text, under the editing area.
  private showStatus(): void {
    const steps = this.picture.steps().length;
    this.status.textContent = (this.picture.error ?? readsAs(steps)) +
      (this.area.drawingUnread ? ` ${DRAWING_UNREAD}` : '');
    // a lone molecule is a reaction not yet made, no error: said in grey
    this.status.style.color = this.picture.error === NO_ARROW ? 'var(--grey-5)' :
      this.picture.error ? 'var(--red-3)' : 'var(--green-2)';
    const lost = this.doc.text.trim() ? this.area.losses : [];
    this.losses.textContent = flagsText(lost, '\n');
    this.losses.style.display = lost.length ? '' : 'none';
  }

  // The right pane: the title, then the Reaction pane of the text (the reaction drawn large, full screen, the steps,
  // the template's notes, Map atoms), the map numbers as the switch says.
  private draw(): void {
    clearTimeout(this.drawingTimer);
    this.dropNext();
    const title = ui.divText(this.doc.title, 'chem-reactions-app-title');
    Object.assign(title.style, {fontWeight: '500', fontSize: '16px', margin: '0 0 8px 0'});
    const page = ui.div([title], 'chem-reactions-app-page');
    let panel: DG.Widget | null = null;
    if (!this.doc.text.trim())
      page.append(ui.divText('Type or paste a reaction on the left, or choose one of the Examples.'));
    else if (tooLargeToDraw(this.doc.text)) // seconds of drawing, the tab frozen: the cells' rule
      page.append(ui.divText(TOO_LARGE, 'chem-reactions-app-too-large'));
    else if (noReactionArrow(this.doc.text)) // a lone molecule: one short line, no error dump
      page.append(ui.divText(NO_ARROW, NO_ARROW_CLASS));
    else {
      // the drawing's texts over its arrows drawn too, while it holds the document's reaction
      const stepTexts = this.stepTexts();
      this.drawnTexts = JSON.stringify(stepTexts);
      // its Export icon offers the drawing as drawn, as the ribbon's Export... does
      panel = reactionPanel(this.doc.text, {mapNumbers: this.doc.mapNumbers ? 'show' : 'hide', stepTexts,
        asDrawn: () => this.area.drawingExport().asDrawn ?? null});
      page.append(panel.root);
    }
    this.swapIn(page, panel);
  }

  // The right pane drawn again a moment after the last change ({@link DRAWING_DELAY_MS}), or once a button held in
  // the editing area is let go (a drag in the drawing).
  private drawSoon(): void {
    clearTimeout(this.drawingTimer);
    this.drawingTimer = setTimeout(() => {
      if (this.pressed)
        this.drawWaits = true;
      else
        this.draw();
    }, DRAWING_DELAY_MS);
  }

  // `page` in place of the page shown: at once when it has no pictures (or none is shown yet); else drawn hidden over
  // the page shown, which stays until its pictures are drawn (at most {@link SWAP_MS}), so that the pane never goes
  // blank between the two.
  private swapIn(page: HTMLElement, panel: DG.Widget | null): void {
    const show = () => {
      if (this.next?.page === page)
        this.next = null;
      this.page?.remove();
      this.panel?.detach();
      Object.assign(page.style, {position: '', visibility: '', left: '', right: '', top: '', padding: '',
        pointerEvents: ''});
      [this.page, this.panel] = [page, panel];
    };
    const canvases = [...page.querySelectorAll('canvas')] as HTMLCanvasElement[];
    if (this.page === null || canvases.length === 0) {
      this.drawing.append(page);
      show();
      return;
    }
    Object.assign(page.style, {position: 'absolute', visibility: 'hidden', left: '0', right: '0', top: '0',
      padding: `${PAGE_PADDING}px`, boxSizing: 'border-box', pointerEvents: 'none'});
    this.drawing.append(page);
    const started = Date.now();
    const wait = () => {
      if (canvases.every((c) => c.style.height !== '') || Date.now() - started > SWAP_MS)
        show();
      else if (this.next?.page === page)
        this.next.timer = setTimeout(wait, 20);
    };
    this.next = {page, panel, timer: setTimeout(wait, 0)};
  }

  // The page being drawn to replace the one shown, let go (a newer one is drawn, or the app closed).
  private dropNext(): void {
    const next = this.next;
    if (!next)
      return;
    this.next = null;
    clearTimeout(next.timer);
    next.panel?.detach();
    next.page.remove();
  }

  // The document changed by a person: kept, read again, its drawing made again (a moment after the last keystroke).
  private edit(change: Partial<DocumentReaction>): void {
    if (this.filling)
      return;
    if (change.title !== undefined) // the person's title: theirs, never put back
      this.dropped = null;
    const title = change.text !== undefined && change.title === undefined ? this.titleWith(change.text) : null;
    if (title !== null) {
      change = {...change, title};
      this.filling = true;
      try {
        this.titleInput.value = title;
      } finally {
        this.filling = false;
      }
    }
    if (Object.entries(change).every(([key, value]) => (this.doc as any)[key] === value)) {
      // the drawing alone changed (a text typed over an arrow, an arrow's style): the line of what the text does not
      // keep says so (not waiting for the next change of the text), and the right pane is drawn again with the
      // drawing's texts
      this.showStatus();
      if (change.text !== undefined && JSON.stringify(this.stepTexts()) !== this.drawnTexts)
        this.drawSoon();
      return;
    }
    // a new text drops the kept drawing; the drawing area keeps it again when the drawing holds more than the text
    if (change.text !== undefined && change.text !== this.doc.text && this.doc.drawing !== undefined)
      change = {...change, drawing: undefined};
    this.doc = {...this.doc, ...change};
    this.picture.free();
    this.picture = new ReactionPicture(this.crux, this.doc);
    this.save();
    this.showStatus();
    if (change.text !== undefined)
      this.drawSoon();
    else if (change.mapNumbers !== undefined)
      this.draw();
    else // the title alone: the drawing as it is (and the page being drawn)
      this.drawing.querySelectorAll('.chem-reactions-app-title').forEach((t) => t.textContent = this.doc.title);
  }

  // The document's kept drawing opened in the drawing area (after its CDXML module came in): the status and the line
  // of what the text does not keep of it said again.
  private drawingShown(): void {
    this.showStatus();
    // its texts over the arrows drawn on the right too (a reopened drawing's were not, until an edit)
    if (!this.picture.error && JSON.stringify(this.stepTexts()) !== this.drawnTexts)
      this.drawSoon();
  }

  // The drawing kept with the document (`drawing`, a CDXML page: DocumentReaction.drawing), or none: saved, nothing
  // drawn again (the text is unchanged).
  private keep(drawing: string | undefined): void {
    if (this.filling || this.doc.drawing === drawing)
      return;
    this.doc = {...this.doc, drawing};
    if (drawing === undefined)
      delete this.doc.drawing;
    this.save();
  }

  // The drawing's texts over and under its steps' arrows, while it holds the document's reaction; else null.
  private stepTexts(): StepTexts | null {
    return this.picture.error ? null : this.area.drawingExport().stepTexts ?? null;
  }

  // The title the document takes with `text` in it, or null to keep its own. An example's title stays through every
  // edit of its reaction (a step added by Next step, the drawing pasted back as CDXML, had wiped it to
  // "Untitled"); it is dropped once the Text is emptied or holds a reaction with none of the example's reactants and
  // products (it had stayed over Clear and a new drawing), and put back when they come back (Clear undone). A
  // title of the person's own always stays; so does an example's when Crux cannot tell (a template's SMARTS).
  private titleWith(text: string): string | null {
    const exampleTitled = (title: string) => [defaultExample(), ...examples()].find((e) => e.title === title);
    const example = this.doc.title ? exampleTitled(this.doc.title) : undefined;
    if (example) {
      if (text.trim() !== '' && shareAMolecule(this.crux, example.text, text) !== false)
        return null;
      this.dropped = example.title;
      return '';
    }
    const back = !this.doc.title && this.dropped !== null && text.trim() !== '' ? exampleTitled(this.dropped) :
      undefined;
    if (back && shareAMolecule(this.crux, back.text, text) === true) {
      this.dropped = null;
      return back.title;
    }
    return null;
  }

  // An example in place of the document: kept as any change is.
  private load(example: DocumentReaction): void {
    this.replace({...example});
    this.save();
  }

  // `doc` in place of the document, the inputs and the drawing with it.
  private replace(doc: DocumentReaction): void {
    this.doc = doc;
    this.dropped = null;
    this.picture.free();
    this.picture = new ReactionPicture(this.crux, this.doc);
    this.fill();
  }

  private save(): void {
    if (!saveDocument(this.doc))
      grok.shell.warning('The browser keeps no copy of the reaction: its storage is blocked or full');
  }

  private confirmReset(): void {
    ui.dialog({title: 'Reset to example', helpUrl: reactionsHelpUrl('the-reaction-sketcher')})
      .add(ui.divText('The default example replaces your reaction, and your changes are lost.'))
      .onOK(() => this.reset())
      .show();
  }

  // The default example again, and nothing kept in the browser: the app opens on it next time too.
  private reset(): void {
    forgetDocument();
    this.replace(defaultExample());
  }

  // The reaction into a table ({@link addToTable}), unless there is none to add.
  private putInTable(): void {
    if (!this.doc.text.trim() || (this.picture.error && this.picture.error !== TOO_LARGE)) {
      grok.shell.warning(NOTHING_TO_ADD);
      return;
    }
    addToTable(this.doc.text, this.doc.title).catch((e) => {
      console.error('Chem | the reaction was not added to a table', e);
      grok.shell.warning('The reaction was not added to a table');
    });
  }

  // The Export panel on the reaction, its map numbers as the app shows them.
  private openExport(): void {
    // a ChemDraw page as the Text: exported as a CDXML cell is, its reaction pictured, its CDXML the page, and what its
    // steps leave out said
    if (isCdxmlText(this.doc.text)) {
      exportCellValue(this.doc.text).catch((e) => console.error('Chem | the Export panel did not open', e));
      return;
    }
    // CDXML as drawn, the texts kept (Export had dropped the text over the arrow): the
    // sketcher's document, while it holds the document's reaction (not while the Text holds what does not read). The
    // SVG and PNG are drawn from the reaction text, which holds no texts: the panel says so.
    const drawn = this.picture.error ? {} : this.area.drawingExport();
    openReactionExport(this.doc.text, {mapNumbers: this.doc.mapNumbers ? 'show' : 'hide', ...drawn});
  }
}
