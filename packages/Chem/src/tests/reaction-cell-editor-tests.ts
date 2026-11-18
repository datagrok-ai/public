/**
 * Category `Reactions: editor`: the reaction cell editor (utils/reactions/reaction-cell-editor.ts). In Crux mode it
 * opens the Crux sketcher's reaction document on the cell's reaction, drawn; OK unchanged writes nothing (the
 * sketcher's value is the cell's string byte for byte, for the corpus and Chem's templates); an edit is written;
 * Cancel writes nothing; closing lets go of the sketcher. RDKit mode's is the same, its footer's flags too.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {EDITOR_FULL_HEIGHT, NOT_READ_NOW, ReactionCellEditor, arrowRows, cellReadFormat, editReactionCell,
  editorBoxHeight, openReactionCellEditor, plainRefusal}
  from '../utils/reactions/reaction-cell-editor';
import {EDITOR_FULL_SCREEN_CLASS, REACTION_EDITOR_CLASS, SKETCH_LOADING_CLASS}
  from '../utils/reactions/reaction-cell-editor';
import {EDIT_ITEM, addReactionExportItem} from '../widgets/reaction-export';
import {CruxSketchLoaderForTests, setCruxSketchLoaderForTests, sketcherOwnsKey} from '../crux/crux-sketcher';
import {loadCruxSketch} from '../crux/crux-sketcher';
import {drawingLossFlags, drawingLosses} from '../utils/reactions/reaction-losses';
import {ChemicalEngine} from '../engine/chem-engine';
import {DEFAULT_REACTIONS_FILE} from '../utils/reactions/consts';
import {inEngine} from './engine-utils';
import {loadFileAsText} from './utils';
import {REACTIONS_HELP_FILE} from '../utils/reactions/reactions-help';
import {SET_VALUE_ARRANGES} from '../utils/reactions/reaction-arrange';
import {defaultExample} from '../utils/reactions/reaction-document';
import {watchReactionEditor} from '../utils/footguns/editors';
import {REFUSED_ROUTE} from './reactions-app-tests';

const ESTER = 'CC(=O)O.OCC>>CC(=O)OCC';
const ROUTE = 'CC(=O)O.Cl>>CC(=O)Cl>>CC(=O)Cl.N>>CC(=O)N';
const TEMPLATE = '[c:1]-B(O)O.[Cl,Br,I][c:2]>>[c:1]-[c:2]';
/** Ethanol and an atom labelled HATU (a MOL alias) to acetaldehyde, as an RXN block. */
const LABELLED_RXN = ['$RXN', '', '      Crux', '', '  2  1',
  '$MOL', '', '  Crux', '', '  3  2  0  0  0  0  0  0  0  0999 V2000',
  '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    2.5981    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  1  0', '  2  3  1  0', 'M  END',
  '$MOL', '', '  Crux', '', '  1  0  0  0  0  0  0  0  0  0999 V2000',
  '    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  'A    1', 'HATU', 'M  END',
  '$MOL', '', '  Crux', '', '  3  2  0  0  0  0  0  0  0  0999 V2000',
  '    8.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    9.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   10.5981    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  1  0', '  2  3  2  0', 'M  END', ''].join('\n');
const inCrux = <T>(run: () => Promise<T>): Promise<T> => inEngine(ChemicalEngine.Crux, run);
/** Crux Sketch, for the tests' loader, as a module whose sketchers do not start (createSketcher rejects). */
const NOT_STARTING = {
  createSketcher: async () => {
    throw new Error('test: the sketcher does not start');
  },
} as unknown as Awaited<ReturnType<CruxSketchLoaderForTests>>;

/** A keydown of `key` ('Ctrl+a', 'Alt+a', 'Delete'), its modifiers as named. */
function keyOf(key: string): KeyboardEvent {
  const parts = key.split('+');
  const name = parts.pop()!;
  return new KeyboardEvent('keydown', {key: name, ctrlKey: parts.includes('Ctrl'), altKey: parts.includes('Alt')});
}

/** A table view of one reaction column holding `values`. */
function gridOf(values: string[]): {view: DG.TableView, grid: DG.Grid, col: DG.Column} {
  const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('reaction', values)]);
  const col = df.col('reaction')!;
  col.semType = 'ChemicalReaction';
  const view = grok.shell.addTableView(df);
  return {view, grid: view.grid, col};
}

/** The editor opened on row `row`'s cell, once its sketcher has read it. */
async function open(grid: DG.Grid, row: number): Promise<ReactionCellEditor> {
  const editor = editReactionCell(grid.cell('reaction', row));
  await editor.ready;
  return editor;
}

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

/** The reaction arrows Crux draws in the editor. */
function arrowsOf(editor: ReactionCellEditor): number {
  return editor.dialog.root.querySelector('crux-sketch')?.shadowRoot?.querySelectorAll('g[data-arrow]').length ?? 0;
}

/** The sketcher's box, its crux-sketch element and its canvas, on the screen (CSS px). */
function editorSizes(editor: ReactionCellEditor): {box: DOMRect, element: DOMRect | null, canvas: number[]} {
  const box = editor.dialog.root.querySelector(`.${REACTION_EDITOR_CLASS}`) as HTMLElement;
  const element = box.querySelector('crux-sketch') as HTMLElement | null;
  const canvas = editor.sketch?.layout.canvas;
  return {box: box.getBoundingClientRect(), element: element?.getBoundingClientRect() ?? null,
    canvas: canvas ? [canvas.width, canvas.height] : [0, 0]};
}

/** files/tests/reactions-crux.csv row 37 (a convergent route) opened in the editor once said the same refusal twice,
 * above the drawing and in the footer's slot, over an empty sketcher inviting a drawing ("Draw reactants, a Reaction
 * arrow, then products"). As the platform opens it (the editor and its foot-gun slot, package.ts editReactionCell):
 * said once, in the box, no sketcher left; OK off; Cancel leaves the cell as it was. Crux Sketch on crux-core 15b1bbd65
 * (non-chain-routes) and later opens a convergent route: row 37 now opens, drawn, OK on; the refusal is shown on a
 * route it still does not open, a step with no product (the app's {@link REFUSED_ROUTE}). */
async function refusedOnce(): Promise<void> {
  const corpus = DG.DataFrame.fromCsv(await loadFileAsText('tests/reactions-crux.csv'));
  const convergent = corpus.col('reaction')!.get(37) as string;
  expect(corpus.col('id')!.get(37), 'drug-sildenafil-commercial');
  const route = REFUSED_ROUTE;
  const {view, grid, col} = gridOf([route, convergent]);
  col.meta.units = 'route';
  try {
    const opened = await open(grid, 1);
    try {
      expect(opened.refused, null, 'a convergent route opens');
      await awaitCheck(() => !!opened.sketch && !opened.sketch.isPending && !opened.sketch.isEmpty,
        'the convergent route was not drawn', 10000);
      expect(opened.dialog.getButton('OK').disabled, false, 'a convergent route: OK on');
    } finally {
      opened.dialog.close();
    }
    expect(col.get(1), convergent, 'closed unedited: the cell as it was');
    const editor = editReactionCell(grid.cell('reaction', 0));
    const slot = watchReactionEditor(editor);
    editor.slot = slot;
    await editor.ready;
    await awaitCheck(() => editor.dialog.getButton('OK').disabled, 'OK off', 3000);
    try {
      expect(editor.refused !== null, true, 'refused');
      const text = editor.dialog.root.textContent ?? '';
      const said = text.match(/cannot open this reaction|does not open this reaction/g) ?? [];
      expect(said.length, 1, `the refusal said once: ${text.slice(0, 400)}`);
      const box = editor.dialog.root.querySelector(`.${REACTION_EDITOR_CLASS}`)!;
      const line = box.querySelector('.chem-reaction-cell-editor-refused') as HTMLElement | null;
      expect(line !== null && line.getBoundingClientRect().height > 0, true, 'said in the box, where the sketcher was');
      expect(box.querySelector('crux-sketch'), null, 'no sketcher left to draw in');
      expect(editor.sketch, null);
      // settled, for what waits on the box (the BDD's OK and Cancel steps: a demo scenario timed out on it)
      const status = DG.Widget.find(box as HTMLElement)?.getWidgetStatus()?.values ?? {};
      expect(status.ready === true && status.pending === false && status.refused === true, true,
        `the box's status settled: ${JSON.stringify(status)}`);
      expect(slot?.shown.length ?? 0, 0, 'the footer\'s slot says nothing more');
      const cancel = editor.dialog.getButton('CANCEL');
      if (cancel)
        cancel.click();
      else
        editor.dialog.close();
      expect(col.get(0), route, 'the cell as it was');
    } finally {
      editor.dialog.close();
    }
  } finally {
    view.close();
  }
}

/** The editor's dialog made larger, the sketcher kept its size, empty space below and
 * right of it. Its box fills the dialog and follows it, larger and smaller, the sketcher's canvas with it. */
async function followsDialog(): Promise<void> {
  const {view, grid} = gridOf([ESTER]);
  try {
    const editor = await open(grid, 0);
    try {
      const opened = editorSizes(editor);
      const contents = editor.dialog.root.querySelector('.d4-dialog-contents')!.getBoundingClientRect();
      expect(Math.abs(contents.right - opened.box.right) < 16 && Math.abs(contents.bottom - opened.box.bottom) < 16,
        true, `the box fills the dialog's contents: box ${JSON.stringify(opened.box)}, contents ` +
        JSON.stringify(contents));
      const r = editor.dialog.root.getBoundingClientRect();
      editor.dialog.root.style.width = `${Math.round(r.width + 300)}px`;
      editor.dialog.root.style.height = `${Math.round(r.height + 150)}px`;
      await awaitCheck(() => {
        const now = editorSizes(editor);
        return now.box.width > opened.box.width + 250 && now.box.height > opened.box.height + 100 &&
          now.canvas[0] > opened.canvas[0] + 200;
      }, `larger: ${JSON.stringify(opened)} -> ${JSON.stringify(editorSizes(editor))}`, 5000);
      const larger = editorSizes(editor);
      expect(larger.element !== null && Math.abs(larger.element.width - larger.box.width) < 4 &&
        Math.abs(larger.element.height - larger.box.height) < 4, true, 'the sketcher fills its box');
      editor.dialog.root.style.width = `${Math.round(r.width - 200)}px`;
      editor.dialog.root.style.height = `${Math.round(r.height - 100)}px`;
      await awaitCheck(() => {
        const now = editorSizes(editor);
        return now.box.width < opened.box.width - 150 && now.canvas[0] < opened.canvas[0] - 100;
      }, `smaller: ${JSON.stringify(editorSizes(editor))}`, 5000);
    } finally {
      editor.dialog.close();
    }
  } finally {
    view.close();
  }
}

category('Reactions: editor', () => {
  // a template that is valid SMILES too (15 of the enumerator's templates: 3, 18, 29 ...) opened read as SMILES, its
  // bracket atoms without an H count drawn as radicals. The editor reads a SMARTS column's value (units `smarts`) as
  // SMARTS: a query, no radical; OK unchanged writes nothing. Red on a crux-sketch that does not read a reaction in the
  // format a write names (its readReaction).
  test('a template in a SMARTS column opens as a query, its bracket atoms no radicals; OK unchanged writes nothing',
    () => inCrux(async () => {
      const template = '[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4]';
      const {view, grid, col} = gridOf([template]);
      col.meta.units = 'smarts';
      try {
        expect(cellReadFormat(grid.cell('reaction', 0)), 'smarts', 'read as SMARTS: the column\'s units');
        const editor = await open(grid, 0);
        try {
          await awaitCheck(() => !editor.sketch!.isPending && !editor.sketch!.isEmpty, 'not drawn', 10000);
          const radicals = (editor.sketch!.svg.match(/data-radical='[1-9]'/g) ?? []).length;
          expect(radicals, 0, `bracket atoms drawn as radicals: ${radicals}`);
          expect(editor.value, template, 'OK unchanged writes the cell\'s string');
        } finally {
          editor.dialog.close();
        }
      } finally {
        view.close();
      }
      // the value's own notation: SMARTS syntax outside a CX part; else as Crux Sketch detects it
      const {view: v2, grid: g2} = gridOf(['[#6:1]>>[#6:1]O', ESTER]);
      try {
        expect(`${cellReadFormat(g2.cell('reaction', 0))} ${cellReadFormat(g2.cell('reaction', 1))}`, 'smarts null');
      } finally {
        v2.close();
      }
    }));

  // at 1024 x 700 the dialog's OK, Cancel and its x were off-screen, Escape the only way out. The
  // dialog is kept in the window: no larger than it less its margins, moved back in, the sketcher's box giving way
  test('in a 1024 x 700 window the dialog fits: its x, OK and Cancel in the window, the sketcher taking the rest',
    () => inCrux(async () => {
      const {view, grid} = gridOf([ROUTE]);
      try {
        const editor = await open(grid, 0);
        try {
          const [w, h] = [1024, 700];
          editor.fitToWindow(w, h);
          const inside = (e: Element | null, what: string) => {
            const r = e?.getBoundingClientRect();
            expect(!!r && r.width > 0 && r.left >= 0 && r.top >= 0 && r.right <= w && r.bottom <= h, true,
              `${what} inside ${w} x ${h}: ${JSON.stringify(r)}`);
          };
          inside(editor.dialog.root, 'the dialog');
          inside(editor.dialog.getButton('OK'), 'OK');
          inside(editor.dialog.getButton('CANCEL'), 'Cancel');
          // Datagrok's close icon (its header's grok-font-icon-close)
          inside(editor.dialog.root.querySelector('.d4-dialog-header .grok-font-icon-close'), 'its x');
          // the room the window leaves (its margins, 24 px each) goes to the dialog, the sketcher's box the rest
          const dialog = editor.dialog.root.getBoundingClientRect();
          expect(dialog.height >= h - 2 * 24 - 2, true, `the dialog takes the window's height: ${dialog.height} px`);
          const box = editor.dialog.root.querySelector(`.${REACTION_EDITOR_CLASS}`)!.getBoundingClientRect();
          expect(box.height >= 320, true, `the sketcher's box: ${box.height} px`);
        } finally {
          editor.dialog.close();
        }
      } finally {
        view.close();
      }
    }));

  // a double-click on a reaction cell threw FormatException: 622.4 in Datagrok's Modal.clipToScreen: on a scaled screen
  // (125%) the window and the dialog measure fractional pixels, the dialog moved back into the window was placed at
  // left 622.4px, and Datagrok reads a dialog's left and top back as integers. Its place and size are whole pixels
  test('the dialog moved back into the window, and full screen, are placed at whole pixels: Datagrok reads integers',
    () => inCrux(async () => {
      const {view, grid} = gridOf([ESTER]);
      try {
        const editor = await open(grid, 0);
        try {
          const root = editor.dialog.root;
          const whole = (when: string) => {
            for (const side of ['left', 'top', 'width', 'height'] as const) {
              const v = root.style[side];
              expect(v === '' || /^-?\d+px$/.test(v), true, `${when}: its ${side} ${v}, not whole pixels`);
            }
          };
          // hanging out at the right and the bottom of a window measured fractional, as a scaled screen measures it
          Object.assign(root.style, {left: `${window.innerWidth - 100}px`, top: `${window.innerHeight - 100}px`});
          editor.fitToWindow(window.innerWidth - 0.4, window.innerHeight - 0.6);
          whole('moved back in');
          editor.toggleFullScreen();
          whole('full screen');
          editor.toggleFullScreen();
          whole('back from full screen');
        } finally {
          editor.dialog.close();
        }
      } finally {
        view.close();
      }
    }));

  // there was no full screen for the dialog (not in its title bar, its menus or the cell's), and double-click the only
  // way in: the header's switch, beside its help and x, and a double-click on its title fill the window and give the
  // dialog back; the cell's menu has "Edit reaction...", which opens the editor
  test('the dialog goes full screen from its header and back; the cell\'s menu opens it: Edit reaction...', () =>
    inCrux(async () => {
      const {view, grid} = gridOf([ESTER]);
      try {
        const editor = await open(grid, 0);
        try {
          const root = editor.dialog.root;
          const icon = root.querySelector(`.d4-dialog-header .${EDITOR_FULL_SCREEN_CLASS}`) as HTMLElement | null;
          expect(!!icon && icon.getBoundingClientRect().width > 0, true, 'the header\'s full-screen switch shows');
          const rect = () => root.getBoundingClientRect();
          const boxHeight = () => root.querySelector(`.${REACTION_EDITOR_CLASS}`)!.getBoundingClientRect().height;
          const [before, boxBefore] = [rect(), boxHeight()];
          icon!.click();
          const [w, h] = [window.innerWidth, window.innerHeight];
          await awaitCheck(() => rect().width >= w - 2 * 24 - 2 && rect().height >= h - 2 * 24 - 2,
            `full screen: the window less its margins: ${JSON.stringify(rect())} in ${w} x ${h}`, 3000);
          const full = rect();
          expect(full.left >= 0 && full.top >= 0 && full.right <= w && full.bottom <= h, true,
            `inside the window: ${JSON.stringify(full)}`);
          const ok = editor.dialog.getButton('OK')!.getBoundingClientRect();
          expect(ok.bottom <= h && ok.width > 0, true, `OK in the window: ${JSON.stringify(ok)}`);
          expect(boxHeight() >= boxBefore - 1, true, `the sketcher takes the room: ${boxBefore} -> ${boxHeight()} px`);
          // the switch kept "Full screen" as its tip while full screen
          expect(`${icon!.getAttribute('aria-label')} | ${editor.fullScreenTip()}`,
            'Exit full screen | Exit full screen', 'full screen: the switch says it exits');
          icon!.click();
          await awaitCheck(() => Math.abs(rect().width - before.width) <= 2 &&
            Math.abs(rect().height - before.height) <= 2,
          `back to its size: ${JSON.stringify(before)} -> ${JSON.stringify(rect())}`, 3000);
          expect(Math.abs(rect().left - before.left) <= 2 && Math.abs(rect().top - before.top) <= 2, true,
            `back in its place: ${JSON.stringify(before)} -> ${JSON.stringify(rect())}`);
          expect(icon!.getAttribute('aria-label'), 'Full screen', 'back: the switch says full screen again');
          root.querySelector('.d4-dialog-title')!.dispatchEvent(new MouseEvent('dblclick', {bubbles: true}));
          expect(editor.fullScreen, true, 'a double-click on its title: full screen');
        } finally {
          editor.dialog.close();
        }
        const menu = DG.Menu.popup();
        expect(addReactionExportItem({args: {menu, item: grid.cell('reaction', 0)}}), true, 'the cell\'s menu');
        const editors = () => document.querySelectorAll(`.${REACTION_EDITOR_CLASS}`);
        const opened = editors().length;
        menu.find(EDIT_ITEM).click();
        await awaitCheck(() => editors().length > opened, `${EDIT_ITEM} opens the editor`, 5000);
        const dialog = editors()[editors().length - 1].closest('.d4-dialog');
        expect(dialog?.querySelector('.d4-dialog-title')?.textContent, 'Reaction', 'the cell editor\'s dialog');
        (dialog?.querySelector('.d4-dialog-header .grok-font-icon-close') as HTMLElement | null)?.click();
      } finally {
        view.close();
      }
    }));

  // crux-sketch's host-look tokens: the sketcher's own dialogs and tooltips in Datagrok's look, its
  // tokens Datagrok's colour variables, so that a theme follows, and its dialogs' titles Datagrok's
  test('the sketcher takes Datagrok\'s look: its host-look tokens are Datagrok\'s colours', () => inCrux(async () => {
    const {view, grid} = gridOf([ESTER]);
    try {
      const editor = await open(grid, 0);
      try {
        const el = editor.dialog.root.querySelector(`.${REACTION_EDITOR_CLASS} crux-sketch`) as HTMLElement;
        const token = (name: string) => getComputedStyle(el).getPropertyValue(name).trim().toLowerCase();
        const root = (name: string) => getComputedStyle(document.documentElement).getPropertyValue(name).trim()
          .toLowerCase();
        for (const [token_, datagrok] of [['--crux-sketch-accent', '--blue-1'], ['--crux-sketch-line', '--grey-2'],
          ['--crux-sketch-tip-face', '--beige-1'], ['--crux-sketch-tip-text', '--warm-grey-4'],
          ['--crux-sketch-button-text', '--blue-1']])
          expect(token(token_), root(datagrok) || token(token_), `${token_} is Datagrok's ${datagrok}`);
        expect(token('--crux-sketch-accent'), '#2083d5', 'Datagrok\'s blue');
        expect(token('--crux-sketch-button-case'), 'uppercase');
        expect(token('--crux-sketch-font').startsWith('roboto'), true, token('--crux-sketch-font'));
        // its dialogs' titles as Datagrok's own (the dialog round it): 18px, weight 400
        const title = getComputedStyle(editor.dialog.root.querySelector('.d4-dialog-title') as Element);
        expect(token('--crux-sketch-title-size'), title.fontSize, 'the size of Datagrok\'s dialog title');
        expect(token('--crux-sketch-title-weight'), title.fontWeight, 'the weight of Datagrok\'s dialog title');
        expect(token('--crux-sketch-title-size'), '18px');
        // its dialogs' form rows as Datagrok's (the owner, 2026-11-17: "inputs' labels are right aligned"): a label
        // right-aligned in Datagrok's light grey, 8px from its field; a focused field's underline Datagrok's blue, no
        // box; Settings' check boxes with their names in the label column (no switch set)
        expect(token('--crux-sketch-label-align'), 'right');
        expect(token('--crux-sketch-label-text'), root('--text-color-light') || '#9497a0', 'Datagrok\'s label grey');
        expect(token('--crux-sketch-label-gap'), '8px');
        expect(token('--crux-sketch-field-focus-ring'), 'none');
        expect(token('--crux-sketch-field-focus-line'), root('--blue-1') || '#2083d5', 'Datagrok\'s focused underline');
        expect(token('--crux-sketch-settings-checks'), '', 'Settings\' check boxes: names in the label column');
      } finally {
        editor.dialog.close();
      }
    } finally {
      view.close();
    }
  }));

  test('opens the cell\'s reaction drawn, in a reaction document with the reaction tools', () => inCrux(async () => {
    const {view, grid} = gridOf([ESTER]);
    try {
      const editor = await open(grid, 0);
      try {
        expect(!!editor.sketch, true, 'the sketcher');
        expect(editor.sketch!.document, 'reaction');
        expect(editor.sketch!.toolset, 'reaction');
        expect(editor.value, ESTER, 'what OK would write unchanged: the cell\'s string, byte for byte');
        if (SET_VALUE_ARRANGES) // opened arranged by the write itself: the value as read
          expect(editor.sketch!.value, ESTER, 'the sketcher\'s value, byte for byte');
        await awaitCheck(() => arrowsOf(editor) === 1, 'its arrow is not drawn', 5000);
        // the dialog's help is Chem's Reactions page, at its section
        expect(editor.dialog.helpUrl?.endsWith(`${REACTIONS_HELP_FILE}#editing-a-reaction-cell`), true,
          `the help: ${editor.dialog.helpUrl}`);
      } finally {
        editor.dialog.close();
      }
    } finally {
      view.close();
    }
  }));

  // the dialog made larger, the sketcher kept its size, empty space below and right
  test('the sketcher fills the dialog and follows it as the dialog is resized, larger and smaller', () =>
    inCrux(followsDialog));

  // a click on the sketcher's tools before it had read the reaction was lost (seen in the Reaction Sketcher)
  test('the sketcher\'s box is covered until it has read the cell: a click as it opens reaches no tool',
    () => inCrux(async () => {
      const {view, grid} = gridOf([ROUTE]);
      try {
        const editor = editReactionCell(grid.cell('reaction', 0)); // not awaited: as the dialog opens
        try {
          const host = editor.dialog.root.querySelector(`.${REACTION_EDITOR_CLASS}`) as HTMLElement;
          const cover = host.querySelector(`:scope > .${SKETCH_LOADING_CLASS}`) as HTMLElement | null;
          expect(!!cover && getComputedStyle(cover).display !== 'none', true, 'covered as it opens');
          const r = host.getBoundingClientRect();
          const hit = () => document.elementFromPoint(r.left + r.width / 2, r.top + 20); // the sketcher's top bar
          expect(hit() === cover, true, `a click as it opens reaches ${hit()?.className}`);
          await editor.ready;
          expect(getComputedStyle(cover!).display, 'none', 'uncovered once the cell is read');
          expect(hit() === host.querySelector('crux-sketch'), true, `a click then reaches ${hit()?.tagName}`);
        } finally {
          editor.dialog.close();
        }
      } finally {
        view.close();
      }
    }));

  // a sketcher that does not start must not leave its box under the cover, taking every click
  test('a sketcher that does not start: the cover goes, the dialog says so, and OK writes nothing',
    () => inCrux(async () => {
      const {view, grid, col} = gridOf([ESTER]);
      setCruxSketchLoaderForTests(async () => NOT_STARTING);
      try {
        const editor = editReactionCell(grid.cell('reaction', 0));
        await editor.ready;
        expect(editor.sketch, null, 'no sketcher');
        const host = editor.dialog.root.querySelector(`.${REACTION_EDITOR_CLASS}`) as HTMLElement;
        expect((host.textContent ?? '').includes('could not start'), true, `the dialog says: ${host.textContent}`);
        const cover = host.querySelector(`:scope > .${SKETCH_LOADING_CLASS}`) as HTMLElement;
        expect(getComputedStyle(cover).display, 'none', 'the cover gone');
        const r = host.getBoundingClientRect();
        const hit = document.elementFromPoint(r.left + r.width / 2, r.top + r.height / 2);
        expect(hit !== cover && host.contains(hit), true, `a click on the box reaches ${hit?.className}`);
        editor.dialog.getButton('OK').click();
        expect(col.get(0), ESTER, 'OK wrote nothing');
      } finally {
        setCruxSketchLoaderForTests(null);
        view.close();
      }
    }));

  // at 520 px the left bar scrolled and hid the arrow tool
  // eslint-disable-next-line max-len
  test('the sketcher is tall enough for its left bar: in a 1080-px window every group shows, the arrow tool too, no scroll button',
    () => inCrux(async () => {
      const {view, grid} = gridOf([ESTER]);
      try {
        const editor = await open(grid, 0);
        try {
          const host = editor.dialog.root.querySelector('.chem-reaction-cell-editor') as HTMLElement;
          expect(Math.round(host.getBoundingClientRect().height), editorBoxHeight(),
            'the box\'s height for this window');
          // a shorter window than the package tests' 1920x1080: the box as tall as it allows
          if (window.innerHeight < 1080)
            return;
          expect(editorBoxHeight(), EDITOR_FULL_HEIGHT, 'the full height in a 1080-px window');
          const el = host.querySelector('crux-sketch') as HTMLElement;
          for (const toolset of ['reaction', 'publication'] as const) {
            editor.sketch!.toolset = toolset; // the gear's choice; no edit
            await awaitCheck(() => leftBarWhole(el),
              `${toolset} tools: the left bar scrolls or hides the arrow tool`, 3000);
          }
        } finally {
          editor.dialog.close();
        }
      } finally {
        view.close();
      }
    }));

  test('OK unchanged writes nothing: the value is the cell\'s string for the corpus\'s reactions and Chem\'s templates',
    () => inCrux(async () => {
      const corpus = DG.DataFrame.fromCsv(await loadFileAsText('tests/reactions-crux.csv')).col('reaction')!
        .toList() as string[];
      const templates = (JSON.parse(await loadFileAsText(DEFAULT_REACTIONS_FILE)) as {reactionSmarts: string}[])
        .map((r) => r.reactionSmarts);
      const texts = [...corpus, ...templates, ESTER, ROUTE, TEMPLATE].filter((t) => t?.trim());
      // one sketcher for all: what OK would write unchanged is the sketcher's value as read
      const host = document.createElement('div');
      Object.assign(host.style, {width: '600px', height: '400px', position: 'absolute', left: '-9999px'});
      document.body.append(host);
      const sketch = await (await loadCruxSketch()).createSketcher(host, {document: 'reaction', toolset: 'reaction'});
      // a text Crux Sketch refuses (a route neither a chain nor a tree of chains, a template's `${...}`) is the
      // editor's refused case: it says so and OK writes nothing (the RXN-06 test below); every text it reads must come
      // back as it was, byte for byte
      const changed: string[] = [];
      const refused: string[] = [];
      try {
        for (const text of texts) {
          sketch.value = text;
          if (sketch.lastError !== null)
            refused.push(text);
          else if (sketch.value !== text)
            changed.push(text);
        }
      } finally {
        sketch.destroy();
        host.remove();
      }
      expect(changed.length, 0, `${changed.length} of ${texts.length - refused.length} read not kept: ` +
        changed.slice(0, 3).join(' | '));
      expect(refused.length < texts.length / 10, true, `${refused.length} of ${texts.length} refused: ` +
        refused.slice(0, 3).map((r) => r.slice(0, 80)).join(' | '));
      // and through the dialog: OK without an edit leaves the cell as it was
      const {view, grid, col} = gridOf([TEMPLATE]);
      try {
        const editor = await open(grid, 0);
        editor.dialog.getButton('OK').click();
        expect(col.get(0), TEMPLATE, 'the cell after OK');
        expect(editor.edited, false);
      } finally {
        view.close();
      }
    }));

  // a reagent name drawn as an atom label came back as "*" in a cell, without a word. What the
  // editor's footer says of it comes from Crux Sketch's warnings (the slot's own test is in "Foot-guns: Crux"): on the
  // real engine, a reaction holding an atom labelled HATU loses the label in its reaction SMILES, not in its RXN block.
  test('a label such as HATU: Crux Sketch says reaction SMILES writes it as a plain atom, an RXN block keeps it',
    () => inCrux(async () => {
      const host = document.createElement('div');
      Object.assign(host.style, {width: '600px', height: '400px', position: 'absolute', left: '-9999px'});
      document.body.append(host);
      const sketch = await (await loadCruxSketch()).createSketcher(host, {document: 'reaction', toolset: 'reaction'});
      try {
        sketch.value = LABELLED_RXN;
        expect(sketch.lastError, null, 'the RXN read');
        const smiles = sketch.smiles;
        expect(smiles.includes('>') && !smiles.includes('HATU'), true, `its reaction SMILES: ${smiles}`);
        const asSmiles = drawingLosses(sketch.warnings, smiles);
        expect(asSmiles.labels, 1, `the label lost in reaction SMILES: ${JSON.stringify(sketch.warnings)}`);
        expect(drawingLossFlags(asSmiles).map((f) => f.rule).join(), 'labelAsStar', 'its flag');
        expect(drawingLosses(sketch.warnings, sketch.value).labels, 0, 'the RXN block read keeps it');
      } finally {
        sketch.destroy();
        host.remove();
      }
    }));

  test('a cell Crux Sketch does not read: the editor says why, and OK writes nothing', () => inCrux(async () => {
    const placeholder = '[c:1][Br]>>[c:1]${smiles}'; // a template with a placeholder, as Chem's reaction files hold
    const {view, grid, col} = gridOf([placeholder]);
    try {
      const editor = await open(grid, 0);
      expect(editor.refused !== null, true, 'refused');
      const notice = editor.dialog.root.querySelector('.chem-reaction-cell-editor-refused') as HTMLElement;
      expect(notice.style.display !== 'none' && notice.textContent!.includes('cannot open this reaction'), true,
        'said');
      editor.edited = true; // whatever is drawn meanwhile
      editor.dialog.getButton('OK').click();
      expect(col.get(0), placeholder, 'OK wrote nothing');
      expect(/crux|rdkit|reload/i.test(notice.textContent!), false, `the engine named: ${notice.textContent}`);
    } finally {
      view.close();
    }
  }));

  // a convergent route (the corpus's drug-sildenafil-commercial, row 37) was refused twice, above the drawing and in
  // the footer, over an empty sketcher that invited a drawing; the vendored Crux Sketch opens it
  test('a convergent route opens; a route the editor does not open: said once, in place of the sketcher, ' +
    'no sketcher left; OK off', () => inCrux(refusedOnce));

  // the vendored sketcher's engine stop says "Crux Sketch: the engine failed (...); reload the page ..."
  test('a reason that names the engine or asks for a reload is said in plain words; others as they are', async () => {
    const plain = 'a route that is not a chain of steps opens in a later version';
    expect(plainRefusal(plain), plain);
    const stop = 'Crux Sketch: the engine failed (RuntimeError: unreachable); reload the page to use the sketcher ' +
      'again';
    expect(plainRefusal(stop), NOT_READ_NOW);
    expect(plainRefusal('RDKit could not parse it'), NOT_READ_NOW);
  });

  test('an edit is written by OK; Cancel writes nothing', () => inCrux(async () => {
    const {view, grid, col} = gridOf([ESTER, ESTER]);
    const edit = 'CC(=O)Cl.OCC>>CC(=O)OCC';
    try {
      let editor = await open(grid, 0);
      editor.sketch!.value = edit; // the drawing changed; a person's edit is the BDD scenario's
      editor.edited = true;
      editor.dialog.getButton('OK').click();
      expect(col.get(0), edit, 'written');
      editor = await open(grid, 1);
      editor.sketch!.value = edit;
      editor.edited = true;
      editor.dialog.getButton('CANCEL').click();
      expect(col.get(1), ESTER, 'Cancel wrote nothing');
    } finally {
      view.close();
    }
  }));

  // a route cell with no coordinates of its own opens arranged, as no person's edit
  test('a route cell opens arranged: in rows, its bonds long enough to read; OK without an edit writes nothing',
    () => inCrux(async () => {
      const route = defaultExample().text;
      const {view, grid, col} = gridOf([route]);
      try {
        const editor = await open(grid, 0);
        const el = editor.dialog.root.querySelector('crux-sketch') as HTMLElement;
        const bond = () => editor.sketch!.positions.bondLength;
        await awaitCheck(() => arrowRows(el, bond()) >= 2, `the route in ${arrowRows(el, bond())} row(s)`, 5000);
        expect(bond() >= 16, true, `bonds drawn ${bond().toFixed(1)} px long`);
        expect(editor.edited, false, 'opening arranged is no person\'s edit');
        expect(editor.value, route, 'what OK would write: the cell\'s string');
        if (SET_VALUE_ARRANGES)
          expect(editor.sketch!.canUndo, false, 'opening arranged is no edit to undo');
        editor.dialog.getButton('OK').click();
        expect(col.get(0), route, 'the cell after OK');
      } finally {
        view.close();
      }
    }));

  test('a person\'s edit undone back to the drawing as opened: OK writes nothing', () => inCrux(async () => {
    const {view, grid, col} = gridOf([ROUTE]);
    try {
      const editor = await open(grid, 0);
      editor.sketch!.addMolecule('CCO'); // the drawing changed; a person's edit is the BDD scenario's
      editor.edited = true;
      expect(editor.sketch!.undo(), true, 'the addition undone');
      editor.dialog.getButton('OK').click();
      expect(col.get(0), ROUTE, 'the cell after OK');
    } finally {
      view.close();
    }
  }));

  test('closing lets go: onChange completes and the sketcher is destroyed', () => inCrux(async () => {
    const {view, grid} = gridOf([ROUTE]);
    try {
      const editor = await open(grid, 0);
      let completed = false;
      let passed = 0;
      editor.onChange.subscribe({next: () => passed++, complete: () => completed = true});
      editor.dialog.close();
      await awaitCheck(() => completed, 'onChange did not complete', 2000);
      editor.onChange.next('CC>>CO'); // a completed subject passes nothing on
      expect(passed, 0, 'nothing reaches a subscriber after the dialog closed');
      expect(editor.dialog.root.querySelector('[data-u2-name="cruxSketch"]') === null ||
        !editor.dialog.root.isConnected, true, 'the dialog is gone');
    } finally {
      view.close();
    }
  }));

  // the reaction editor works in both modes: RDKit mode's is Crux Sketch too (it was the reaction's text in a plain
  // dialog), and behaves the same, its footer's checks too (erasing the arrow blocks OK, a molecule written into a
  // reaction column is flagged)
  test('RDKit mode: the same Crux Sketch editor, its help and its flags; an edit is written, Cancel writes nothing',
    () => inEngine(ChemicalEngine.RDKit, async () => {
      const {view, grid, col} = gridOf([ESTER, ESTER]);
      const edit = 'CC(=O)Cl.OCC>>CC(=O)OCC';
      try {
        let editor = await open(grid, 0);
        expect(!!editor.sketch, true, 'Crux Sketch in RDKit mode');
        expect(editor.sketch!.document, 'reaction');
        expect(editor.dialog.root.querySelector('textarea'), null, 'no plain text box');
        await awaitCheck(() => arrowsOf(editor) === 1, 'its arrow is not drawn', 5000);
        expect(editor.dialog.helpUrl?.endsWith(`${REACTIONS_HELP_FILE}#editing-a-reaction-cell`), true,
          `the help: ${editor.dialog.helpUrl}`);
        // the arrow gone: a molecule would be written into the reaction column, a block (RXN-16), OK off
        const slot = watchReactionEditor(editor);
        editor.slot = slot;
        editor.sketch!.value = 'CC(=O)O.OCC.CC(=O)OCC';
        editor.edited = true;
        editor.onChange.next(editor.value);
        await awaitCheck(() => slot.blocked, 'RDKit mode: no arrow left blocks OK', 3000);
        const said = editor.dialog.root.querySelector('.chem-footgun-slot')?.getAttribute('aria-label') ?? '';
        expect(/arrow/i.test(said), true, `the slot names the arrow: ${said}`);
        editor.sketch!.value = edit;
        editor.onChange.next(editor.value);
        await awaitCheck(() => !slot.blocked, 'a reaction again: OK on', 3000);
        editor.edited = true;
        editor.dialog.getButton('OK').click();
        expect(col.get(0), edit, 'written');
        editor = await open(grid, 1);
        editor.sketch!.value = edit;
        editor.edited = true;
        editor.dialog.getButton('CANCEL').click();
        expect(col.get(1), ESTER, 'Cancel wrote nothing');
      } finally {
        view.close();
      }
    }));

  // a double-click on any reaction cell opens the editor, in both modes (it had been kept to the columns of CDXML
  // pages, units=cdxml, RDKit mode's double-click the grid's inline edit)
  test('a double-click opens the editor on every reaction cell, in RDKit mode and in Crux mode', async () => {
    const editor = DG.Func.find({package: 'Chem', name: 'editReactionCell'})
      .find((f) => f.options['role'] === 'cellEditor');
    expect(!!editor, true, 'Chem\'s reaction cell editor');
    expect(editor!.options['columnTags'], 'quality=ChemicalReaction',
      'every reaction column (no units asked: SMILES, SMARTS, RXN, routes, CDXML pages)');
    const kinds: [string, string][] =
      [['smiles', ESTER], ['smarts', TEMPLATE], ['route', ROUTE], ['rxn', LABELLED_RXN]];
    for (const engine of [ChemicalEngine.RDKit, ChemicalEngine.Crux]) {
      await inEngine(engine, async () => {
        for (const [units, value] of kinds) {
          const {view, grid, col} = gridOf([value]);
          col.meta.units = units;
          try {
            // what the double-click calls: the function's body (package.ts editReactionCell)
            const opened = await openReactionCellEditor(grid.cell('reaction', 0));
            try {
              await opened.ready;
              expect(!!opened.sketch, true, `${engine} mode, units ${units}: Crux Sketch opened on the cell`);
              expect(opened.refused, null, `${engine} mode, units ${units}: read (${opened.refused})`);
              expect(opened.value, value, `${engine} mode, units ${units}: OK unchanged writes the cell's string`);
            } finally {
              opened.dialog.close();
            }
          } finally {
            view.close();
          }
        }
      });
    }
  });

  // Ctrl+A over the editor's canvas (on reactions-cdxml.csv) also selected every row of the grid behind the dialog: the
  // dialog kept the focus as it opened, so the key went to the window, where Datagrok's shortcuts listen (core
  // command.dart: Select all on the table); the canvas clicked, Crux kept it. The editor gives its drawing the keyboard
  // as it opens, and its dialog keeps the sketcher's keys pressed anywhere in it (but a text field) from the page;
  // Escape, Enter and Tab go on (Cancel, OK, focus)
  test('the sketcher\'s keys stay in the editor: the drawing focused as it opens; Ctrl+A, Ctrl+Z, Delete, the arrows ' +
    'never reach the platform\'s shortcuts, Tab and Escape do', () => inCrux(async () => {
    expect(['Ctrl+a', 'Ctrl+A', 'Ctrl+z', 'Ctrl+y', 'Delete', 'Backspace', 'ArrowLeft'].every((k) =>
      sketcherOwnsKey(keyOf(k))), true, 'the sketcher\'s keys');
    expect(['Escape', 'Enter', 'Tab', 'a', 'Alt+a'].some((k) => sketcherOwnsKey(keyOf(k))), false, 'the page\'s keys');
    const {view, grid} = gridOf(['CCO>>CC=O', 'CC>>C=C', 'C>>CO']);
    const reached: string[] = [];
    const listen = (e: KeyboardEvent) => reached.push(`${e.ctrlKey ? 'Ctrl+' : ''}${e.key}`);
    window.addEventListener('keydown', listen);
    const editor = await open(grid, 0);
    try {
      const element = editor.dialog.root.querySelector('crux-sketch') as HTMLElement;
      await awaitCheck(() => element.shadowRoot?.activeElement?.localName === 'svg',
        `the drawing focused as the editor opens: ${document.activeElement?.className}`, 5000);
      const keys: [string, boolean][] = [['A', true], ['a', true], ['z', true], ['y', true], ['Delete', false],
        ['ArrowLeft', false]];
      for (const target of [editor.dialog.root, element.shadowRoot!.activeElement!]) {
        for (const [key, ctrl] of keys)
          target.dispatchEvent(new KeyboardEvent('keydown', {key, ctrlKey: ctrl, bubbles: true, composed: true}));
      }
      expect(reached.join(), '', 'the sketcher\'s keys kept from the page, from the dialog and from the drawing');
      expect(grid.dataFrame.selection.trueCount, 0, 'no row of the table selected');
      editor.dialog.root.dispatchEvent(new KeyboardEvent('keydown', {key: 'Tab', bubbles: true, composed: true}));
      expect(reached.join(), 'Tab', 'Tab goes on: focus moves');
      editor.dialog.root.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true, composed: true}));
      expect(reached.join(), 'Tab,Escape', 'Escape goes on: the dialog\'s Cancel');
    } finally {
      window.removeEventListener('keydown', listen);
      editor.dialog.close();
      view.close();
    }
  }));
});
