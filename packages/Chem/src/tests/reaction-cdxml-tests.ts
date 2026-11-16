/**
 * Category `Reactions: CDXML`: reaction cells that hold CDXML (rendering/cdxml-reaction-cells.ts), on Chem's demo
 * table files/reactions-cdxml.csv (its notice: rendering/NOTICE.reactions-cdxml): its column is a reaction column of
 * units `cdxml` in Crux mode and stays text in RDKit mode; each page is read by Crux Sketch and pictured (cut to its
 * reaction when it holds more), its reaction SMILES read from it; a cell edited is written back as CDXML, its texts
 * kept.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {CDXML_TOO_SMALL, CDXML_UNITS, cdxmlFit, cdxmlPictureForTests, cdxmlReactionSmiles, detectCdxmlColumn,
  cdxmlReadCounts, drawCdxmlCell, isCdxmlText, legibleLines, moleculesOnNoStep, noStepWords, queueCdxmlReadForTests,
  reactionPage, resetCdxmlReadsForTests, svgBondLength} from '../rendering/cdxml-reaction-cells';
import {NO_STEP_CLASS, reactionPanel} from '../widgets/reaction-panel';
import {NOTE_CLASS, exportCellValue} from '../widgets/reaction-export';
import {editReactionCell} from '../utils/reactions/reaction-cell-editor';
import {REACTION_CELL_TEXT_COLOUR} from '../rendering/crux-reaction-drawing';
import {ChemicalEngine} from '../engine/chem-engine';
import {inEngine} from './engine-utils';
import {loadFileAsText} from './utils';

/** Chem's demo table of reaction schemes as CDXML. */
export const CDXML_DEMO_FILE = 'reactions-cdxml.csv';

/** The demo table's CDXML column, as text. */
async function demoColumn(): Promise<DG.Column> {
  const df = DG.DataFrame.fromCsv(await loadFileAsText(CDXML_DEMO_FILE));
  return df.col('cdxml')!;
}

category('Reactions: CDXML', () => {
  test('the demo table\'s pages are CDXML; a SMILES, an RXN block and other XML are not', async () => {
    const pages = (await demoColumn()).toList() as string[];
    expect(pages.length, 8, 'the demo table\'s rows');
    expect(pages.every(isCdxmlText), true, 'every page');
    for (const s of ['CC>>CO', '$RXN\n\n  Crux\n\n  1  1\n', '<?xml version="1.0"?><svg/>', '<html><CDXML/>', ''])
      expect(isCdxmlText(s), false, s.slice(0, 20));
  });

  test('a column of CDXML pages is a reaction column of units cdxml in Crux mode, text in RDKit mode (12A)',
    async () => {
      await inEngine(ChemicalEngine.RDKit, async () => {
        const col = await demoColumn();
        detectCdxmlColumn(col);
        expect(col.semType !== DG.SEMTYPE.CHEMICAL_REACTION, true, 'RDKit mode: not a reaction column');
      });
      await inEngine(ChemicalEngine.Crux, async () => {
        const col = await demoColumn();
        detectCdxmlColumn(col);
        expect(col.semType, DG.SEMTYPE.CHEMICAL_REACTION);
        expect(col.meta.units, CDXML_UNITS);
      });
    });

  test('an empty cell of a CDXML column draws nothing, not "Drawing…"', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('cdxml', [''])]);
    df.col('cdxml')!.meta.units = CDXML_UNITS;
    const grid = DG.Viewer.grid(df);
    const canvas = document.createElement('canvas');
    [canvas.width, canvas.height] = [200, 100];
    const g = canvas.getContext('2d')!;
    drawCdxmlCell(g, 0, 0, 200, 100, grid.cell('cdxml', 0));
    const inked = g.getImageData(0, 0, 200, 100).data.some((v, i) => i % 4 === 3 && v > 0);
    expect(inked, false, 'nothing drawn');
  });

  test('every page of the demo table is read: a reaction with its steps', async () => inEngine(ChemicalEngine.Crux,
    async () => {
      await resetCdxmlReadsForTests();
      const pages = (await demoColumn()).toList() as string[];
      const unread: string[] = [];
      for (const [i, page] of pages.entries()) {
        const smiles = await cdxmlReactionSmiles(page);
        if (!smiles?.includes('>'))
          unread.push(`row ${i + 1}: ${smiles}`);
      }
      expect(unread.length, 0, unread.join(' | '));
    }));

  // each cell of the demo table is its page's picture, not a line
  test('every cell of the demo table is pictured once read: no "Drawing…" left and no unread line', async () =>
    inEngine(ChemicalEngine.Crux, async () => {
      await resetCdxmlReadsForTests();
      const col = await demoColumn();
      detectCdxmlColumn(col);
      const grid = DG.Viewer.grid(col.dataFrame);
      const pages = col.toList() as string[];
      const lines: string[] = [];
      for (const [i, page] of pages.entries()) {
        const cell = grid.cell('cdxml', i);
        drawCdxmlCell(cellCanvas().getContext('2d')!, 0, 0, CELL.width, CELL.height, cell); // starts its read
        await cdxmlReactionSmiles(page); // the read done, whatever it reads as
        const canvas = cellCanvas();
        drawCdxmlCell(canvas.getContext('2d')!, 0, 0, CELL.width, CELL.height, cell);
        const height = inkedHeight(canvas);
        if (height <= LINE_HEIGHT)
          lines.push(`row ${i + 1}: ${height} px of ink`);
      }
      expect(lines.length, 0, `a line, not a picture: ${lines.join(' | ')}`);
    }));

  // at a browser zoom over 100% the cells' pictures slid against a sideways scroll, and
  // row 5's route was cut off right and below
  test('at a device pixel ratio of 1.5 a page is drawn in its cell, whole, and moves with it', async () =>
    inEngine(ChemicalEngine.Crux, scaledGrid));

  test('a cell too small for bonds of 4 px says the page is too large to show there, never squeezed; a larger one ' +
    'pictures it', async () => inEngine(ChemicalEngine.Crux, tooSmall));

  // the cell says "double-click to open": the cell editor (what a double-click opens) opens the page
  test('a page too large for its cell opens in the cell editor, drawn, and Cancel leaves it', async () =>
    inEngine(ChemicalEngine.Crux, tooLargeOpens));

  // an edited CDXML cell is written back as CDXML (the sketcher's `cdxml`), not as reaction
  // SMILES, so the page's texts are kept with the edit
  test('a CDXML cell edited: OK writes a CDXML page, the molecule added and the page\'s texts kept', async () =>
    inEngine(ChemicalEngine.Crux, async () => {
      const page = ((await demoColumn()).toList() as string[])[EDITED_ROW];
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('reaction', [page])]);
      const col = df.col('reaction')!;
      col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
      col.meta.units = CDXML_UNITS;
      const view = grok.shell.addTableView(df);
      try {
        const editor = editReactionCell(view.grid.cell('reaction', 0));
        await editor.ready;
        expect(editor.holdsCdxml, true, 'the editor knows the cell holds CDXML');
        expect(editor.refused, null, 'the page read');
        await awaitCheck(() => !editor.sketch!.isPending && !editor.sketch!.isEmpty, 'the page was not drawn', 15000);
        editor.sketch!.addMolecule('CCO'); // the drawing changed; a person's edit is the BDD scenario's
        editor.edited = true;
        editor.dialog.getButton('OK').click();
        const written = col.get(0) as string;
        expect(isCdxmlText(written), true, `a CDXML page: ${written.slice(0, 60)}`);
        expect(written !== page, true, 'the edit written');
        expect(fragments(written), fragments(page) + 1, 'the page\'s molecules and the one added');
        expect(captions(page).length > 0, true, 'the demo page has texts');
        const lost = captions(page).filter((t) => !captions(written).includes(t));
        expect(lost.length, 0, `texts not kept: ${lost.join(' | ')}`);
      } finally {
        view.close();
      }
    }));

  // a reaction text made from a ChemDraw page holds its steps' molecules only; the
  // molecules on no step (a key of structures, a starting material drawn apart) were left out with nothing said
  test('a page with molecules on no step says how many were left out: the Reaction panel and the Export panel', () =>
    inEngine(ChemicalEngine.Crux, async () => {
      const pages = (await demoColumn()).toList() as string[];
      // counted from each page's own <fragment> and <step> ids (the Suzuki page's key: 15 structures)
      expect(pages.map(moleculesOnNoStep).join(' '), '0 0 15 5 2 2 0 0', 'the demo pages\' molecules on no step');
      expect(noStepWords(1), '1 molecule on no step was left out');
      const page = pages[EDITED_ROW];
      const widget = reactionPanel(page);
      try {
        await awaitCheck(() => widget.root.querySelector(`.${NO_STEP_CLASS}`) !== null,
          'the Reaction panel says nothing of the molecules on no step', 15000);
        expect(widget.root.querySelector(`.${NO_STEP_CLASS}`)!.textContent,
          '2 molecules on no step were left out of the reaction');
      } finally {
        widget.detach();
      }
      const dialog = await exportCellValue(page);
      try {
        const note = dialog.root.querySelector(`.${NOTE_CLASS}`)?.textContent ?? '';
        expect(note.includes('2 molecules on no step were left out of the SVG and PNG; CDXML keeps them'), true, note);
      } finally {
        dialog.close();
      }
      // a page whose steps hold every molecule says nothing
      const whole = reactionPanel(pages[0]);
      try {
        await awaitCheck(() => whole.root.querySelector('canvas') !== null, 'the first page was not read', 15000);
        expect(whole.root.querySelector(`.${NO_STEP_CLASS}`), null, 'nothing on no step: nothing said');
      } finally {
        whole.detach();
      }
    }));

  // at the opening row height the Suzuki and N-sulfonylation cells drew their whole page, the key box too
  test('a page with a key box beside its reaction is pictured cut to its reaction, readable at the opening row ' +
    'height; the cell\'s value stays the whole page', async () => inEngine(ChemicalEngine.Crux, keyBoxLeftOut));

  // a fast scroll through a 1,500-page table queued a read per row passed, with no cancel, and the
  // queue took ~40 s of CPU to drain after the scroll stopped. A queued read whose cells have all left the screen when
  // its turn comes is skipped (a cell that comes back on screen asks again); a read someone needs is never skipped
  test('a queued read whose cells have left the screen is skipped when its turn comes; those on screen are read',
    async () => inEngine(ChemicalEngine.Crux, async () => {
      const pages = (await demoColumn()).toList() as string[];
      const before = {...cdxmlReadCounts};
      const off = Array.from({length: 500}, (_, i) => `${pages[0]}<!-- off screen ${i} -->`);
      const on = pages.slice(0, 3).map((p, i) => `${p}<!-- on screen ${i} -->`);
      const answers = await Promise.all([...off.map((v) => queueCdxmlReadForTests(v, () => false)),
        ...on.map((v) => queueCdxmlReadForTests(v, () => true))]);
      expect(answers.slice(0, off.length).filter((a) => a).length, 0, 'no page off screen read');
      expect(answers.slice(off.length).every((a) => a), true, 'the pages on screen read');
      expect(cdxmlReadCounts.read - before.read, on.length, 'the reads made: those on screen');
      expect(cdxmlReadCounts.skipped - before.skipped, off.length, 'the reads skipped: those off screen');
      expect(await cdxmlPictureForTests(off[0]) !== null, true, 'a page skipped is read when someone needs it');
    }));
});

/** The demo table's pages whose reaction sits over a key box (the Suzuki and N-sulfonylation pages, rows 3 and 4), a
 * cell's height as the table opens (about 150 px), and the least bond a cell then draws them with (CSS px:
 * the whole pages' were about 4). */
const KEYED_ROWS = [2, 3];
const OPENING_ROW_HEIGHT = 150;
const KEYED_MIN_BOND = 6;

/** The ids a CDXML page's steps list (their molecules, arrows, texts over and under the arrows). */
function stepIds(page: string): string[] {
  const doc = new DOMParser().parseFromString(page, 'text/xml');
  return Array.from(doc.getElementsByTagName('step')).flatMap((step) => ['ReactionStepReactants',
    'ReactionStepProducts', 'ReactionStepArrows', 'ReactionStepObjectsAboveArrow', 'ReactionStepObjectsBelowArrow']
    .flatMap((name) => (step.getAttribute(name) ?? '').split(/\s+/).filter((id) => id)));
}

/** The Suzuki and N-sulfonylation pages drew their key box (a framed list of reagents and roles) under the
 * reaction, so that at the opening row height the reaction was a few pixels high and the key's labels overlapped. A
 * cell pictures such a page cut to its reaction: the objects its steps list and the texts beside them, the key left
 * out; its bonds readable in a 400 x 150 cell. A page with nothing away from its reaction is pictured whole, and the
 * cell's reaction is the whole page's (the editor opens the page, Export writes it). */
async function keyBoxLeftOut(): Promise<void> {
  await resetCdxmlReadsForTests();
  const pages = (await demoColumn()).toList() as string[];
  for (const row of KEYED_ROWS) {
    const page = pages[row];
    const cut = reactionPage(page);
    expect(cut !== null, true, `row ${row + 1}: cut to its reaction`);
    const doc = new DOMParser().parseFromString(cut!, 'text/xml');
    const texts = Array.from(doc.getElementsByTagName('t')).map((t) => (t.textContent ?? '').trim());
    expect(texts.includes('Key') || Array.from(doc.getElementsByTagName('graphic'))
      .some((g) => g.getAttribute('GraphicType') === 'Rectangle'), false, `row ${row + 1}: the key box left out`);
    const ids = new Set(Array.from(doc.querySelectorAll('[id]')).map((e) => e.getAttribute('id')));
    const lost = stepIds(page).filter((id) => !ids.has(id));
    expect(lost.length, 0, `row ${row + 1}: what its steps list kept; lost ${lost.join(' ')}`);
    expect(moleculesOnNoStep(cut!), 0, `row ${row + 1}: no molecule of the key left`);
    const picture = await cdxmlPictureForTests(page);
    expect(picture !== null && picture.height / picture.width < 0.45, true,
      `row ${row + 1}: the picture a reaction's row, not the page: ${JSON.stringify(picture)}`);
    const fit = cdxmlFit(picture!, CELL.width, OPENING_ROW_HEIGHT);
    expect(!fit.tooSmall && picture!.bond * fit.scale >= KEYED_MIN_BOND, true,
      `row ${row + 1}: bonds of ${(picture!.bond * fit.scale).toFixed(1)} px in a ${CELL.width} x ` +
      `${OPENING_ROW_HEIGHT} cell`);
    expect((await cdxmlReactionSmiles(page))?.includes('>') ?? false, true, `row ${row + 1}: its reaction read`);
  }
  const whole = pages.map((p, i) => KEYED_ROWS.includes(i) || reactionPage(p) === null ? '' : `row ${i + 1}`)
    .filter((r) => r);
  expect(whole.length, 0, `nothing away from their reaction, pictured whole: ${whole.join(', ')}`);
}

/** The cell the picture tests paint (CSS px), and the most ink a one-line message makes in it (px high). */
const CELL = {width: 400, height: 200};
const LINE_HEIGHT = 24;
/** The demo row the edit and export tests use: a three-step route over two rows with texts over its arrows (DCM,
 * (COCl)2 DMF (cat.), Fragment 1, ...: 9 captions; the sketcher's CDXML keeps them all, an addition too). It
 * was row 1, the Boc page, whose structures over the arrow are no captions: "texts kept" held of none there, and the
 * app's export test could not start (no texts said). Rows 0, 2 and 3 have texts that crux-sketch d0b8bb9's CDXML drops
 * in part ("Rxn" labels, a "+" text, a long catalyst name): crux-sketch's to fix. */
export const EDITED_ROW = 4;

/** The box of `canvas`'s ink (device px), or null when it has none. */
function inkBox(canvas: HTMLCanvasElement): {left: number, top: number, right: number, bottom: number} | null {
  const {data, width, height} = canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height);
  let [left, top, right, bottom] = [width, height, -1, -1];
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      if (data[(y * width + x) * 4 + 3] > 128)
        [left, top, right, bottom] = [Math.min(left, x), Math.min(top, y), Math.max(right, x), Math.max(bottom, y)];
    }
  }
  return right < 0 ? null : {left, top, right, bottom};
}

/** The demo table's row 5 (index {@link EDITED_ROW}) as a grid cell, its page read. */
async function readCell(): Promise<DG.GridCell> {
  const col = await demoColumn();
  detectCdxmlColumn(col);
  const grid = DG.Viewer.grid(col.dataFrame);
  const cell = grid.cell('cdxml', EDITED_ROW);
  drawCdxmlCell(cellCanvas().getContext('2d')!, 0, 0, CELL.width, CELL.height, cell); // starts its read
  await cdxmlReactionSmiles(col.get(EDITED_ROW));
  return cell;
}

/** Zoomed in, the cells' pictures slid against a sideways scroll and row 5's route was
 * cut off. The grid hands a cell renderer a context scaled by the device pixel ratio and the cell's box in CSS px,
 * which the reaction renderer turns into canvas px: the picture is drawn in canvas px (the context's scale undone),
 * whole inside its cell, and moves with the cell as the grid scrolls. */
async function scaledGrid(): Promise<void> {
  await resetCdxmlReadsForTests();
  const cell = await readCell();
  const ratio = 1.5;
  const canvas = document.createElement('canvas');
  [canvas.width, canvas.height] = [Math.round(900 * ratio), Math.round(400 * ratio)];
  const at = (cssX: number) => {
    const g = canvas.getContext('2d')!;
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.clearRect(0, 0, canvas.width, canvas.height);
    g.setTransform(ratio, 0, 0, ratio, 0, 0); // as the grid's context comes at a device pixel ratio of 1.5
    drawCdxmlCell(g, ratio * cssX, ratio * 50, ratio * CELL.width, ratio * CELL.height, cell, ratio);
    return inkBox(canvas);
  };
  const box = at(300);
  const cellBox = {left: 300 * ratio, top: 50 * ratio, right: (300 + CELL.width) * ratio,
    bottom: (50 + CELL.height) * ratio};
  expect(box !== null, true, 'drawn');
  expect(box!.left >= cellBox.left && box!.top >= cellBox.top && box!.right < cellBox.right &&
    box!.bottom < cellBox.bottom, true, `inside its cell ${JSON.stringify(cellBox)}: ${JSON.stringify(box)}`);
  // whole: the ink clear of the cell's edges (a picture cut off by them would touch them)
  expect(box!.left > cellBox.left + 2 && box!.right < cellBox.right - 2 && box!.top > cellBox.top + 2 &&
    box!.bottom < cellBox.bottom - 2, true, `whole, not cut by the cell's edges: ${JSON.stringify(box)}`);
  // the grid scrolled 120 CSS px sideways: the cell 120 px to the left, its picture with it
  const scrolled = at(180);
  expect(Math.abs(scrolled!.left - (box!.left - 120 * ratio)) <= 2, true,
    `moved with its cell: ${box!.left} -> ${scrolled!.left} (${-120 * ratio} expected)`);
}

/** A page fitted to its cell, shrunk until it fits; below bonds of 4 CSS px the molecules
 * are not drawn, the cell says the page is too large to show at that size. */
async function tooSmall(): Promise<void> {
  await resetCdxmlReadsForTests();
  const cell = await readCell();
  const read = {width: 0, height: 0, bond: 0};
  // what a cell draws, read from its context's calls (a colour test misread the grey line's faint edge pixels)
  const painted = (w: number, h: number) => {
    const canvas = document.createElement('canvas');
    [canvas.width, canvas.height] = [w, h];
    const g = canvas.getContext('2d')!;
    const calls = {images: 0, texts: [] as string[], looks: [] as {colour: string, align: string, x: number}[],
      squeezed: [] as string[], wider: [] as string[]};
    const [image, text] = [g.drawImage.bind(g), g.fillText.bind(g)];
    g.drawImage = ((...args: Parameters<typeof g.drawImage>) => {
      calls.images++;
      return (image as (...a: unknown[]) => void)(...args);
    }) as typeof g.drawImage;
    g.fillText = (t: string, x: number, y: number, max?: number) => {
      calls.texts.push(t);
      calls.looks.push({colour: String(g.fillStyle), align: g.textAlign, x});
      const width = g.measureText(t).width;
      if (max !== undefined && width > max)
        calls.squeezed.push(t);
      if (width > w - 8) // the cell's room, its 4 px padding kept
        calls.wider.push(t);
      return text(t, x, y, max);
    };
    drawCdxmlCell(g, 0, 0, w, h, cell);
    return calls;
  };
  const big = painted(CELL.width, CELL.height);
  expect(`${big.images} ${big.texts.join()}`, '1 ', 'a 400 x 200 cell: the page pictured');
  // the message squeezed to the box's width read as a smear in Datagrok's view header: a box too narrow
  // for it shows its first words and an ellipsis, or nothing, never squeezed
  const small = painted(90, 40);
  expect(small.images, 0, 'a 90 x 40 cell: not drawn');
  expect(small.texts.length === 1 && small.texts[0].endsWith('…') &&
    CDXML_TOO_SMALL.startsWith(small.texts[0].slice(0, -1)), true, `a 90 x 40 cell: shortened: ${small.texts}`);
  const header = painted(150, 20); // Datagrok's view header: one line high
  expect(header.texts.every((t) => t.endsWith('…')), true, `a 150 x 20 box: shortened: ${header.texts}`);
  expect(painted(150, 10).texts.join(), '', 'a box lower than a line: nothing said');
  for (const said of [small, header])
    expect(`${said.squeezed.join()}|${said.wider.join()}`, '|', 'no text squeezed, none wider than its box');
  const width = (t: string) => 7 * t.length;
  expect(legibleLines(width, CDXML_TOO_SMALL, 400, 40, 12, 15.6).join('|'), CDXML_TOO_SMALL, 'one line');
  expect(legibleLines(width, CDXML_TOO_SMALL, 200, 40, 12, 15.6).join('|'),
    'Too large for this cell:|double-click to open', 'broken after its colon');
  expect(legibleLines(width, CDXML_TOO_SMALL, 100, 40, 12, 15.6).join('|'), 'Too large for…', 'its first words');
  expect(legibleLines(width, CDXML_TOO_SMALL, 5, 40, 12, 15.6).join('|'), '', 'no room even for an ellipsis');
  // reactions-cdxml.csv row 6, the mechanisms page: the cell said only "Too large to show at this size";
  // it says how to see the page, in the reaction cells' muted grey, centred in the cell
  expect(CDXML_TOO_SMALL, 'Too large for this cell: double-click to open');
  const wide = painted(300, 40);
  expect(`${wide.images} ${wide.texts.join(' ')}`, `0 ${CDXML_TOO_SMALL}`, 'a 300 x 40 cell: one line');
  for (const [said, w] of [[wide, 300], [small, 90]] as const) {
    expect(said.looks.every((l) => l.colour === REACTION_CELL_TEXT_COLOUR), true,
      `the cells' muted grey: ${JSON.stringify(said.looks)}`);
    expect(said.looks.every((l) => l.align === 'center' && Math.abs(l.x - w / 2) < 1), true,
      `centred in a ${w} px cell: ${JSON.stringify(said.looks)}`);
  }
  expect(cdxmlFit({...read, width: 800, height: 200, bond: 20}, 90, 40).tooSmall, true, 'bonds of 2 px: too small');
  expect(cdxmlFit({...read, width: 800, height: 200, bond: 20}, 400, 200).tooSmall, false, 'bonds of 9 px: drawn');
  expect(cdxmlFit({...read, width: 100, height: 50, bond: 20}, 400, 200).scale, 1, 'never larger than natural');
  expect(svgBondLength('<svg xmlns="http://www.w3.org/2000/svg"><g id="p1-bond-0"><line x1="0" y1="0" x2="30" ' +
    'y2="40"/></g></svg>'), 50, 'a bond\'s length read from the picture');
  // a line stops short of an atom's label: its atoms' centres give the length (the pages keep ChemDraw's labels
  // contracted, and the lines' median put the mechanisms page, 4.20 px bonds, under 4 px in a 400 x 200 cell)
  expect(svgBondLength('<svg xmlns="http://www.w3.org/2000/svg"><g id="p1-atom-0" data-atom="0" data-x="0" ' +
    'data-y="0"/><g id="p1-atom-1" data-atom="1" data-x="30" data-y="40"/><g id="p1-bond-0" data-atoms="0 1">' +
    '<line x1="6" y1="8" x2="30" y2="40"/></g></svg>'), 50, 'a bond to a label: atom centre to atom centre');
}

/** The demo table's mechanisms page (row 6: 18 arrows on 7 rows), which a grid's default reaction cell is too small
 * for. */
const MECHANISMS_ROW = 5;

/** The cell says "double-click to open" of a page too large for it: the mechanisms page, said so in a 300 x 80 cell,
 * opens in the cell editor (what the grid's double-click opens), read and drawn;
 * Cancel writes nothing. */
async function tooLargeOpens(): Promise<void> {
  await resetCdxmlReadsForTests();
  const page = ((await demoColumn()).toList() as string[])[MECHANISMS_ROW];
  const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('reaction', [page])]);
  const col = df.col('reaction')!;
  col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
  col.meta.units = CDXML_UNITS;
  const view = grok.shell.addTableView(df);
  try {
    const cell = view.grid.cell('reaction', 0);
    const canvas = document.createElement('canvas');
    [canvas.width, canvas.height] = [300, 80];
    drawCdxmlCell(canvas.getContext('2d')!, 0, 0, 300, 80, cell); // starts its read
    await cdxmlReactionSmiles(page);
    const said: string[] = [];
    const g = canvas.getContext('2d')!;
    const text = g.fillText.bind(g);
    g.fillText = (t: string, x: number, y: number, max?: number) => {
      said.push(t);
      return text(t, x, y, max);
    };
    drawCdxmlCell(g, 0, 0, 300, 80, cell);
    expect(said.join(' '), CDXML_TOO_SMALL, 'a 300 x 80 cell says the page is too large for it');
    const editor = editReactionCell(cell);
    await editor.ready;
    expect(editor.refused, null, 'the page read');
    await awaitCheck(() => !editor.sketch!.isPending && !editor.sketch!.isEmpty, 'the page was not drawn', 15000);
    editor.dialog.getButton('CANCEL').click();
    await awaitCheck(() => !editor.dialog.root.isConnected, 'the dialog did not close', 3000);
    expect(col.get(0), page, 'Cancel wrote nothing');
  } finally {
    view.close();
  }
}

/** A blank canvas the size of {@link CELL}. */
function cellCanvas(): HTMLCanvasElement {
  const canvas = document.createElement('canvas');
  [canvas.width, canvas.height] = [CELL.width, CELL.height];
  return canvas;
}

/** How many rows of `canvas` hold any ink, from the first to the last (px). */
function inkedHeight(canvas: HTMLCanvasElement): number {
  const {data, width, height} = canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height);
  let first = -1;
  let last = -1;
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      if (data[(y * width + x) * 4 + 3] > 0) {
        first = first < 0 ? y : first;
        last = y;
        break;
      }
    }
  }
  return first < 0 ? 0 : last - first + 1;
}

/** A CDXML page's molecules (its `fragment` elements, outside any other fragment). */
function fragments(cdxml: string): number {
  const doc = new DOMParser().parseFromString(cdxml, 'application/xml');
  return [...doc.getElementsByTagName('fragment')].filter((f) => f.parentElement?.closest('fragment') == null).length;
}

/** A CDXML page's texts (its `t` elements outside the atoms' labels), each its words with spaces collapsed. */
export function captions(cdxml: string): string[] {
  const doc = new DOMParser().parseFromString(cdxml, 'application/xml');
  return [...doc.getElementsByTagName('t')].filter((t) => t.closest('n') == null)
    .map((t) => (t.textContent ?? '').replace(/\s+/g, ' ').trim()).filter((t) => t !== '');
}
