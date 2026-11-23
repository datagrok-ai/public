/**
 * Category `Reactions: Crux`: reaction cells, the Reaction panel and the reaction previews (the reaction dialogs', full
 * screen included, the reaction browser's cards, the reaction editor's: utils/reactions/reaction-drawing.ts)
 * drawn by Crux in Crux mode, and by RDKit as before in RDKit mode (rendering/crux-reaction-drawing.ts). The rows are
 * 41 of a reaction corpus (files/tests/reactions-crux.csv, the corpus's ids kept: single steps of every
 * category RDKit reads, routes of every kind) and Chem's own test-reactions.csv. The renderer drawing them is this
 * bundle's, so that the engine's test override reaches it (engine-utils.ts), on the cells of a real grid.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {RDKitReactionRenderer, parseMultiStepReaction} from '../rendering/rdkit-reaction-renderer';
import {CHEM_CENTRE, CruxPictureOptions, chemDrawOptions, cruxDrawsSteps, cruxReactionCounts, cruxReactionNotes,
  cruxReactionSteps, cruxReadsReaction, drawReactionByCrux, productAtomsMissing}
  from '../rendering/crux-reaction-drawing';
import {EXPORT_ICON_CLASS, FULL_SCREEN_CENTRE, FULL_SCREEN_CLASS, FULL_SCREEN_FILL, FULL_SCREEN_LEGEND_GAP,
  FULL_SCREEN_LINE_SHARE, FULL_SCREEN_MAX_BOND, MAPPER_FAILED_WORDS, MAPPING_CLASS, MAP_LOAD_TIMEOUT_MS, MAP_TIMEOUT_MS,
  NOT_READ_WORDS, noSideWords, openReactionFullScreen, reactionPanel, setMapperLoaderForTests,
  setReactionMapperForTests} from '../widgets/reaction-panel';
import {LEGEND_CLASS, PICTURE_LEGEND_GAP} from '../widgets/reaction-centre-legend';
import {EXPORT_TITLE, PREVIEW_CLASS, exportInputsOf} from '../widgets/reaction-export';
import {defaultExample, examples} from '../utils/reactions/reaction-document';
import {agentsBackOverTheArrow, mapReaction, mappedInTheValue, mapperExported}
  from '../utils/reactions/reaction-mapping';
import {ChemicalEngine} from '../engine/chem-engine';
import {ChemWasm, cruxModuleSync} from '../engine/crux-module';
import {setSketchMapperForTests, sketchReactionMapper} from '../utils/reactions/reaction-sketch-mapper';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {RDKIT_COMMON_RENDER_OPTS, createRenderingOpts, getRdKitModule} from '../utils/chem-common-rdkit';
import {UNSPECIFIED_STEREO_IS_UNKNOWN, renderingOptions, renderingValue, setRenderingValue}
  from '../rendering/rendering-settings';
import {PackageFunctions} from '../package';
import {_package} from '../package-test';
import {expectFallback, expectNoFallback, inEngine} from './engine-utils';
import {loadFileAsText} from './utils';
import {PREVIEW_CALL, drawReactionPreview, setCruxUnreadablePreviewsForTests}
  from '../utils/reactions/reaction-drawing';
import {ReactionBrowser} from '../utils/reactions/reaction-browser';
import {openReactionEditor} from '../utils/reactions/reaction-editor';
import {transformationReactionsUI, twoComponentReactionUI} from '../utils/reactions/ui';
import {NamedReaction} from '../utils/reactions/types';
import {setDefaultReactionsForTests} from '../utils/reactions/reaction-storage';

/** The cell the grid gives a reaction column by default (CSS px). */
const [CELL_W, CELL_H] = [600, 150];

/** Pixels with any ink in the canvas, or in the box x0..x1, y0..y1 (device px) of it. */
function ink(canvas: HTMLCanvasElement, box?: [number, number, number, number]): number {
  const [x0, y0, x1, y1] = box ?? [0, 0, canvas.width, canvas.height];
  const data = canvas.getContext('2d')!.getImageData(x0, y0, x1 - x0, y1 - y0).data;
  let n = 0;
  for (let i = 3; i < data.length; i += 4)
    n += data[i] > 0 ? 1 : 0;
  return n;
}

/** The width and height (device px) of the box the canvas's ink spans; [0, 0] for none. */
function inkExtent(canvas: HTMLCanvasElement): [number, number] {
  const data = canvas.getContext('2d')!.getImageData(0, 0, canvas.width, canvas.height).data;
  let [x0, y0, x1, y1] = [Infinity, Infinity, -1, -1];
  for (let i = 3; i < data.length; i += 4) {
    if (data[i] === 0)
      continue;
    const p = (i - 3) / 4;
    const [x, y] = [p % canvas.width, Math.floor(p / canvas.width)];
    [x0, y0, x1, y1] = [Math.min(x0, x), Math.min(y0, y), Math.max(x1, x), Math.max(y1, y)];
  }
  return x1 < 0 ? [0, 0] : [x1 - x0 + 1, y1 - y0 + 1];
}

/** The cell `row` of `column` drawn by `renderer` on a canvas of its own, at the page's device pixel ratio. */
function drawCell(renderer: RDKitReactionRenderer, grid: DG.Grid, column: string, row: number): HTMLCanvasElement {
  const r = window.devicePixelRatio;
  const canvas = ui.canvas(CELL_W * r, CELL_H * r);
  renderer.render(canvas.getContext('2d')!, 0, 0, CELL_W, CELL_H, grid.cell(column, row));
  return canvas;
}

/** `run()` in Crux mode, then the server's engine setting again. */
const inCrux = <T>(run: () => Promise<T>): Promise<T> => inEngine(ChemicalEngine.Crux, run);

/** The Reaction panel's Mapping section once it has answered (its `data-state` no longer waiting, loading or mapping),
 * waited for at most `ms`: the state and the section's words. */
async function mappingAnswer(panel: DG.Widget, ms: number = 20000): Promise<{state: string, text: string}> {
  const status = () => panel.root.querySelector(`.${MAPPING_CLASS}`) as HTMLElement | null;
  const going = () => ['waiting', 'loading', 'mapping'].includes(status()?.dataset.state ?? 'waiting');
  for (const start = Date.now(); going() && Date.now() - start < ms;)
    await DG.delay(50);
  return {state: status()?.dataset.state ?? 'no Mapping section', text: status()?.textContent ?? ''};
}

/** The options Crux's reaction drawing is sent while `run()` runs (Crux mode), each draw's parsed, in order. */
async function sentOptions(run: () => Promise<void>): Promise<any[]> {
  const scheme: any = (cruxModuleSync() as unknown as ChemWasm).ChemReactionScheme.prototype;
  const draw = scheme.draw_to_canvas;
  const sent: any[] = [];
  scheme.draw_to_canvas = function(this: any, canvas: any, options: string) {
    sent.push(JSON.parse(options));
    return draw.call(this, canvas, options);
  };
  try {
    await run();
  } finally {
    scheme.draw_to_canvas = draw;
  }
  return sent;
}

/** A reaction drawing by Crux: the canvas it drew on and the options it was sent. */
interface CruxDraw {
  canvas: HTMLCanvasElement;
  options: any;
}

/** What drew reactions while `run()` ran: Crux's drawings (none when Chem's Crux file is not loaded), and how many
 * times RDKit's reaction objects were asked to draw (a picture RDKit's renderer kept is copied, not counted); with
 * `run()`'s result. */
async function reactionDraws<T>(run: () => Promise<T>): Promise<{result: T, crux: CruxDraw[], rdkit: number}> {
  const crux: CruxDraw[] = [];
  const scheme: any = (cruxModuleSync() as unknown as ChemWasm | null)?.ChemReactionScheme.prototype ?? null;
  const cruxDraw = scheme?.draw_to_canvas;
  if (scheme) {
    scheme.draw_to_canvas = function(this: any, canvas: HTMLCanvasElement, options: string) {
      crux.push({canvas, options: JSON.parse(options)});
      return cruxDraw.call(this, canvas, options);
    };
  }
  const rxn = getRdKitModule().get_rxn('C>>C');
  const rdkitProto: any = Object.getPrototypeOf(rxn);
  rxn.delete();
  const rdkitCalls = ['draw_to_canvas', 'draw_to_canvas_with_highlights', 'get_svg', 'get_svg_with_highlights']
    .filter((name) => typeof rdkitProto[name] === 'function');
  const kept = rdkitCalls.map((name) => rdkitProto[name]);
  let rdkit = 0;
  rdkitCalls.forEach((name, i) => rdkitProto[name] = function(this: any, ...args: any[]) {
    rdkit++;
    return kept[i].apply(this, args);
  });
  try {
    const result = await run();
    return {result, crux, rdkit};
  } finally {
    if (scheme)
      scheme.draw_to_canvas = cruxDraw;
    rdkitCalls.forEach((name, i) => rdkitProto[name] = kept[i]);
  }
}

/** The open dialog titled `title`, waited for at most `ms`. */
async function dialogTitled(title: string, ms: number = 10000): Promise<DG.Dialog> {
  await awaitCheck(() => DG.Dialog.getOpenDialogs().some((d) => d.title === title), `no "${title}" dialog`, ms);
  return DG.Dialog.getOpenDialogs().find((d) => d.title === title)!;
}

/** Clicks the reaction browser's card of the reaction `name` in `root` (its name label; the click reaches the card),
 * once the card is there (waited for at most `ms`). */
async function clickCard(root: HTMLElement, name: string, ms: number = 15000): Promise<void> {
  const label = () => [...root.querySelectorAll('label.ui-label')].find((l) => l.textContent === name) as
    HTMLElement | undefined;
  await awaitCheck(() => !!label(), `no "${name}" card`, ms);
  label()!.click();
}

/** The canvas the dialog's previews are drawn on: a click on it opens the reaction full screen. */
function previewCanvasOf(dialog: DG.Dialog): HTMLCanvasElement {
  const canvas = [...dialog.root.querySelectorAll('canvas')].find((c) => c.style.cursor === 'pointer');
  expect(!!canvas, true, `the preview canvas of "${dialog.title}"`);
  return canvas!;
}

category('Reactions: Crux', () => {
  let corpus: DG.DataFrame;
  let samples: DG.DataFrame;
  let view: DG.TableView;

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) { // this bundle's RDKit loaded: the category may run alone
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    corpus = DG.DataFrame.fromCsv(await loadFileAsText('tests/reactions-crux.csv'));
    samples = DG.DataFrame.fromCsv(await loadFileAsText('test-reactions.csv'));
    // (test-reactions.csv ends in blank lines: its empty cells are left out)
    const reactions = [...corpus.col('reaction')!.toList(), ...samples.col('reaction')!.toList()].filter((s) => s);
    const table = DG.DataFrame.fromColumns([DG.Column.fromStrings('reaction', reactions)]);
    await grok.data.detectSemanticTypes(table);
    view = grok.shell.addTableView(table);
  });

  after(async () => view?.close());

  test('every cell of the corpus and the samples draws through Crux', async () => {
    const grid = view.grid;
    expect(grid.col('reaction')!.cellType, 'ChemicalReaction');
    await inCrux(async () => {
      const renderer = new RDKitReactionRenderer(getRdKitModule());
      const before = {...cruxReactionCounts};
      for (let row = 0; row < view.dataFrame.rowCount; row++) {
        const canvas = drawCell(renderer, grid, 'reaction', row);
        const value = grid.cell('reaction', row).cell.value;
        expect(ink(canvas) > 200, true, `row ${row} drew ${ink(canvas)} pixels: ${value}`);
      }
      expect(cruxReactionCounts.asText, before.asText, 'values drawn as text');
      // (a value drawn twice at one size is drawn once and copied)
      const values = new Set(view.dataFrame.col('reaction')!.toList()).size;
      expect(cruxReactionCounts.drawn >= before.drawn + values, true, 'reactions Crux drew');
      expect(renderer.canvasCounter, 0, 'pictures RDKit drew');
    });
  });

  test('a cell painted again is a copy of its first picture; another size is drawn anew', async () => {
    const route = view.dataFrame.col('reaction')!.toList().find((s: string) => s.split('>>').length > 2)!;
    const row = view.dataFrame.col('reaction')!.toList().indexOf(route);
    await inCrux(async () => {
      const renderer = new RDKitReactionRenderer(getRdKitModule());
      for (const r of [0, row]) {
        const first = drawCell(renderer, view.grid, 'reaction', r);
        const drawn = cruxReactionCounts.drawn;
        const again = drawCell(renderer, view.grid, 'reaction', r);
        expect(cruxReactionCounts.drawn, drawn, `drawings for the second paint of row ${r}`);
        const [a, b] = [first, again].map((c) => c.getContext('2d')!.getImageData(0, 0, c.width, c.height).data);
        expect(a.every((v, i) => v === b[i]), true, `row ${r}: the second paint's pixels are the first's`);
        const smaller = ui.canvas(420 * window.devicePixelRatio, 120 * window.devicePixelRatio);
        renderer.render(smaller.getContext('2d')!, 0, 0, 420, 120, view.grid.cell('reaction', r));
        expect(cruxReactionCounts.drawn > drawn, true, `row ${r} drawn anew at 420 x 120`);
      }
      expect(renderer.canvasCounter, 0, 'pictures RDKit drew');
    });
  });

  test('a cell\'s molecules get Chem\'s draw options, its column\'s Rendering settings among them', async () => {
    const table = DG.DataFrame.fromColumns([DG.Column.fromStrings('reaction', ['ClCC(=O)O.OCC>>ClCC(=O)OCC'])]);
    const col = table.col('reaction')!;
    col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
    const grid = table.plot.grid();
    await inCrux(async () => {
      const renderer = new RDKitReactionRenderer(getRdKitModule());
      const stereo = renderingValue(UNSPECIFIED_STEREO_IS_UNKNOWN, col) === 'true';
      const sent = await sentOptions(async () => {
        drawCell(renderer, grid, 'reaction', 0);
        drawCell(renderer, grid, 'reaction', 0); // a copy
        setRenderingValue(UNSPECIFIED_STEREO_IS_UNKNOWN, col, String(!stereo)); // drawn anew
        drawCell(renderer, grid, 'reaction', 0);
      });
      expect(sent.length, 2, 'drawings');
      expect(sent[0].drawOptions.bondLineWidth, RDKIT_COMMON_RENDER_OPTS.bondLineWidth, 'the line width');
      expect(JSON.stringify(sent[0].drawOptions.atomColourPalette),
        JSON.stringify(RDKIT_COMMON_RENDER_OPTS.atomColourPalette), 'the palette');
      expect(sent[0].drawOptions.unspecifiedStereoIsUnknown, stereo, 'the first drawing\'s setting');
      expect(sent[1].drawOptions.unspecifiedStereoIsUnknown, !stereo, 'the column\'s new setting');
      expect(JSON.stringify(sent[1].drawOptions), JSON.stringify(createRenderingOpts(renderingOptions(col))),
        'the options of the column\'s molecule cells');
    });
  });

  /** The line widths Crux strokes `reaction` with in a cell's box, one per stroke. */
  const strokeWidths = (reaction: string, options: CruxPictureOptions): number[] => {
    const canvas = ui.canvas(CELL_W, CELL_H);
    const g = canvas.getContext('2d')!;
    const widths: number[] = [];
    const stroke = g.stroke;
    g.stroke = function(this: CanvasRenderingContext2D, ...args: any[]) {
      widths.push(this.lineWidth);
      return (stroke as any).apply(this, args);
    };
    drawReactionByCrux(g, reaction, 0, 0, CELL_W, CELL_H, options);
    return widths;
  };

  test('Crux strokes a reaction at Chem\'s molecule line width (choice 1A)', async () => {
    // (without Chem's draw options Crux strokes bonds, arrow and plus signs at 2 px)
    await inCrux(async () => {
      const widths = strokeWidths('ClCC(=O)O.OCC>>ClCC(=O)OCC', {drawOptions: chemDrawOptions()});
      const width = RDKIT_COMMON_RENDER_OPTS.bondLineWidth;
      expect(widths.length > 10, true, `strokes: ${widths.length}`);
      expect(widths.every((w) => Math.abs(w - width) < 0.05), true, `stroke widths ${[...new Set(widths)].join(', ')}`);
    });
  });

  // the reaction centre was drawn at Crux's 2x and 1x of Chem's 0.7 px line, hairlines in Datagrok (the owner: "still
  // too thin"): Chem asks for 4x and 2x (CHEM_CENTRE: a bond of one line, each line of a double one), full screen for
  // 3x and 1.5x of its thicker line. Red on an engine before chem-wasm's centreWidth and centreLinesWidth
  test('Chem strokes a reaction centre at 4x its line, each line of a double centre bond at 2x; full screen at 3x ' +
    'and 1.5x', async () => {
    await inCrux(async () => {
      const dielsAlder = examples()[1].text; // mapped: bonds made (single) and bonds whose order changed (double)
      const has = (widths: number[], w: number) => widths.some((x) => Math.abs(x - w) < 0.05);
      const line = RDKIT_COMMON_RENDER_OPTS.bondLineWidth;
      const cell = strokeWidths(dielsAlder, {drawOptions: chemDrawOptions()});
      expect(`${has(cell, line)} ${has(cell, CHEM_CENTRE.centreWidth * line)} ` +
        `${has(cell, CHEM_CENTRE.centreLinesWidth * line)}`, 'true true true',
      `a cell: plain, centre and centre double-bond lines; stroke widths ${[...new Set(cell)].join(', ')}`);
      const fullLine = 2.8; // full screen's line at a 100 px bond (FULL_SCREEN_LINE_SHARE)
      const full = strokeWidths(dielsAlder, {drawOptions: {...chemDrawOptions(), bondLineWidth: fullLine},
        centre: FULL_SCREEN_CENTRE});
      expect(`${has(full, FULL_SCREEN_CENTRE.centreWidth * fullLine)} ` +
        `${has(full, FULL_SCREEN_CENTRE.centreLinesWidth * fullLine)}`, 'true true',
      `full screen: stroke widths ${[...new Set(full)].join(', ')}`);
    });
  });

  test('the Reaction panel\'s molecules get Chem\'s draw options', async () => {
    await inCrux(async () => {
      const sent = await sentOptions(async () => {
        const panel = reactionPanel(corpus.col('reaction')!.get(0));
        const host = ui.dialog('Reaction panel test').add(panel.root).show();
        try {
          await DG.delay(300);
        } finally {
          host.close();
        }
      });
      expect(sent.length > 0, true, 'drawings');
      for (const options of sent)
        expect(JSON.stringify(options.drawOptions), JSON.stringify(createRenderingOpts({})), 'the package\'s options');
    });
  });

  test('a route Crux reads is drawn whole, as one chain (decision 9b); another in step panels', async () => {
    // (Crux's binding reads one step per value until the route reader reaches it: then every route is drawn whole)
    const reactions: string[] = view.dataFrame.col('reaction')!.toList();
    const routes = [...new Set(reactions.filter((s) => s.split('>>').length > 2))];
    expect(routes.length > 10, true, 'routes in the table');
    await inCrux(async () => {
      const renderer = new RDKitReactionRenderer(getRdKitModule());
      for (const route of routes) {
        const steps = new Set(parseMultiStepReaction(route).flat()).size;
        const before = cruxReactionCounts.drawn;
        drawCell(renderer, view.grid, 'reaction', reactions.indexOf(route));
        const drawings = cruxReactionCounts.drawn - before;
        if (cruxReadsReaction(route))
          expect(drawings, 1, `one drawing for the chain: ${route}`);
        else
          expect(drawings > 1 || steps === 1, true, `${drawings} drawings for ${steps} panels: ${route}`);
      }
    });
  });

  test('a value Crux reads no reaction from is written as text, not a cross', async () => {
    await inCrux(async () => {
      const renderer = new RDKitReactionRenderer(getRdKitModule());
      const table = DG.DataFrame.fromColumns([DG.Column.fromStrings('reaction', ['CC(=O)O.OCC>>CC(=O)OCC('])]);
      table.col('reaction')!.semType = DG.SEMTYPE.CHEMICAL_REACTION;
      const before = cruxReactionCounts.asText;
      const canvas = drawCell(renderer, table.plot.grid(), 'reaction', 0);
      expect(cruxReactionCounts.asText, before + 1);
      // the text sits on the cell's middle line, the cross would reach its corners
      const r = window.devicePixelRatio;
      expect(ink(canvas, [0, 0, 10 * r, 10 * r]), 0, 'ink in the top left corner');
      expect(ink(canvas) > 50, true, 'the text');
      expect(renderer.canvasCounter, 0, 'pictures RDKit drew');
    });
  });

  test('RDKit mode draws reactions with RDKit, as before', async () => {
    await inEngine(ChemicalEngine.RDKit, async () => {
      const renderer = new RDKitReactionRenderer(getRdKitModule());
      const before = cruxReactionCounts.drawn;
      const canvas = drawCell(renderer, view.grid, 'reaction', 0);
      expect(ink(canvas) > 200, true);
      expect(renderer.canvasCounter, 1, 'pictures RDKit drew');
      expect(cruxReactionCounts.drawn, before, 'reactions Crux drew');
    });
  });

  test('the Reaction panel shows in Crux mode only', async () => {
    const reaction = corpus.col('reaction')!.get(0);
    expect(await inCrux(async () => PackageFunctions.reactionPanelShown(reaction)), true, 'Crux mode');
    expect(await inEngine(ChemicalEngine.RDKit, async () => PackageFunctions.reactionPanelShown(reaction)), false,
      'RDKit mode');
  });

  test('the Reaction panel draws the reaction, opens full screen and closes', async () => {
    const reaction = corpus.col('reaction')!.get(0);
    await inCrux(async () => {
      const panel = reactionPanel(reaction);
      const host = ui.dialog('Reaction panel test').add(panel.root).show();
      try {
        await DG.delay(300);
        const canvas = panel.root.querySelector('canvas')!;
        expect(canvas.width > 0 && ink(canvas) > 200, true, 'the panel\'s picture');
        const full = openReactionFullScreen(reaction);
        try {
          await DG.delay(300);
          const big = full.root.querySelector('canvas')!;
          expect(big.clientWidth > canvas.clientWidth, true, 'larger than the panel\'s');
          expect(ink(big) > 2 * ink(canvas), true, 'the full-screen picture, drawn at twice the bond length');
          // not a right-shifted dialog, CANCEL its only button: it fills the window, a view with no buttons (closed
          // by its x or Escape, as Datagrok's picture views)
          const r = full.root.getBoundingClientRect();
          const [ww, wh] = [window.innerWidth, window.innerHeight];
          expect(r.left >= -1 && r.top >= -1 && r.right <= ww + 1 && r.bottom <= wh + 1 && r.width >= 0.9 * ww &&
            r.height >= 0.9 * wh, true,
          `the window filled: ${JSON.stringify(r)} in ${window.innerWidth} x ${window.innerHeight}`);
          const buttons = [...full.root.querySelectorAll('.d4-dialog-footer button')]
            .filter((b) => (b as HTMLElement).offsetParent !== null).map((b) => b.textContent);
          expect(buttons.join(), '', 'no OK or Cancel: a view');
          const pictureBox = big.getBoundingClientRect();
          expect(pictureBox.right <= r.right + 1, true, `the picture inside the dialog: ${JSON.stringify(pictureBox)}`);
        } finally {
          full.close();
        }
        expect(DG.Dialog.getOpenDialogs().some((d) => d.title === 'Reaction'), false, 'the full-screen view closed');
      } finally {
        host.close();
      }
    });
  });

  // in a 2549 x 1304 window the full-screen view drew the reaction small in a large empty area, its body scrolled and
  // its legend was under the fold, over a band of nothing and oversized buttons; then the legend sat far from a centred
  // drawing still small for the window, its lines hairlines. Now the drawing and its legend are one centred group, the
  // drawing filling 80% of its room's limiting dimension, its lines scaled
  test('the full-screen view fits its window: no scroll, the picture and its legend one centred group, the picture ' +
    'filling its room, its lines scaled; a route\'s steps in a list of their own', async () => inCrux(async () => {
    const route = defaultExample().text;
    const cases: [string, string, string[]][] = [['a mapped reaction', examples()[1].text, []],
      ['a route', route, cruxReactionSteps(route) ?? []]];
    for (const [name, text, steps] of cases) {
      const full = openReactionFullScreen(text, steps);
      try {
        const host = full.root.querySelector(`.${FULL_SCREEN_CLASS}`) as HTMLElement;
        const picture = host.querySelector('canvas') as HTMLCanvasElement;
        const stage = picture.parentElement!;
        await awaitCheck(() => ink(picture) > 200, `${name}: the picture drawn`, 5000);
        await DG.delay(200); // a resize as the dialog settles draws it again
        const contents = full.root.querySelector('.d4-dialog-contents') as HTMLElement;
        for (const box of [contents, host]) {
          expect(box.scrollHeight <= box.clientHeight + 1, true,
            `${name}: no scroll in ${box.className}: ${box.scrollHeight} > ${box.clientHeight}`);
        }
        const r = full.root.getBoundingClientRect();
        const c = contents.getBoundingClientRect();
        expect(r.bottom - c.bottom <= 2, true, `${name}: the contents fill the dialog, no band under them: ` +
          `${c.bottom} in ${r.bottom}`);
        const legend = picture.nextElementSibling as HTMLElement;
        expect(legend.classList.contains(LEGEND_CLASS) && legend.style.display !== 'none', true,
          `${name}: the legend shown, right after the picture`);
        const [p, l, s] = [picture, legend, stage].map((e) => e.getBoundingClientRect());
        // a route too long for its room even at Crux's own bond is drawn taller, its room scrolling (from its top)
        const overflows = stage.scrollHeight > stage.clientHeight + 1;
        expect(l.height > 0 && l.top >= p.bottom - 1 && l.top - p.bottom <= FULL_SCREEN_LEGEND_GAP + 1 &&
          (overflows || (l.bottom <= c.bottom + 1 && l.bottom <= window.innerHeight)), true,
        `${name}: the legend right under the picture, in view: ${JSON.stringify(l)} under ${JSON.stringify(p)}`);
        const mid = (a: number, b: number) => (a + b) / 2;
        expect(Math.abs(mid(p.left, p.right) - mid(s.left, s.right)) <= 2 &&
          Math.abs(mid(l.left, l.right) - mid(p.left, p.right)) <= 2 &&
          (overflows || Math.abs(mid(p.top, l.bottom) - mid(s.top, s.bottom)) <= 2), true,
        `${name}: the picture and its legend one group centred in its room: ${JSON.stringify([p, l, s])}`);
        const bond = Number(picture.dataset.bondLength);
        expect(bond > 0 && bond <= FULL_SCREEN_MAX_BOND, true, `${name}: a fitted bond, ${bond}`);
        const [roomWidth, roomHeight] = [s.width, s.height - l.height - FULL_SCREEN_LEGEND_GAP];
        expect(bond === FULL_SCREEN_MAX_BOND || overflows || p.width >= 0.7 * roomWidth ||
          p.height >= 0.7 * roomHeight, true,
        `${name}: the picture fills its room (${p.width} x ${p.height} in ${roomWidth} x ${roomHeight}, ` +
        `bond ${bond})`);
        const [iw, ih] = inkExtent(picture);
        expect(iw >= 0.7 * picture.width || ih >= 0.7 * picture.height, true,
          `${name}: the drawing fills its canvas (${iw} x ${ih} in ${picture.width} x ${picture.height})`);
        const line = Number(picture.dataset.lineWidth);
        expect(line >= bond * FULL_SCREEN_LINE_SHARE - 0.1, true,
          `${name}: its lines scaled with its bond: ${line} px at a ${bond} px bond`);
        if (steps.length > 1) {
          // the default route was fitted as one row, at a 4-11 px bond. A route is broken into
          // rows, at a bond of 44 px (the view's old one) or more where its rows at 44 px fit the room; at Crux's own
          // 22 px or more where only those fit; else at 44 px, its room scrolling
          const column = Number(picture.dataset.column);
          expect(column > 0, true, `${name}: drawn as a column, its rows wrapping (${picture.dataset.column})`);
          const scheme = new (cruxModuleSync() as unknown as ChemWasm).ChemReactionScheme(text, '');
          try {
            const at = (o: object) => Array.from(scheme.get_size(JSON.stringify(
              {drawOptions: chemDrawOptions(), footnotes: false, ...o}))) as number[];
            // the room the view fits a route's rows to (FULL_SCREEN_FILL of it), less a margin: where this says a bond
            // fits, the view's own search (to within 3 px) finds it fits too
            const roomFor = (b: number) => {
              const [, height, percent] = at({bondLength: b, column: true, width: column});
              return (percent === undefined || percent >= 99.5) && height <= 0.95 * FULL_SCREEN_FILL * roomHeight;
            };
            const least = roomFor(44) ? 44 : roomFor(22) ? 22 : 0;
            const said = least ? `at least ${least} px, the room allowing it` : '44 px, its room scrolling';
            expect(least === 0 ? !overflows || bond >= 43 : bond >= least - 3, true,
              `${name}: a bond of ${said}: ${bond} px in ${roomWidth} x ${roomHeight}` +
              `${overflows ? ', scrolling' : ''}`);
            const [oneRowWidth, oneRowHeight] = at({bondLength: bond});
            if (oneRowWidth > column) {
              expect(picture.clientHeight > 1.5 * oneRowHeight, true, `${name}: its rows wrap: ` +
                `${picture.clientHeight} px tall, one row ${Math.round(oneRowHeight)} x ${Math.round(oneRowWidth)}`);
            }
          } finally {
            scheme.free();
          }
          const list = host.querySelector(':scope > div:last-child') as HTMLElement;
          expect(list.querySelectorAll('canvas').length, steps.length, `${name}: the steps in their list`);
          expect(list.getBoundingClientRect().bottom <= c.bottom + 1, true, `${name}: the list inside the dialog`);
        }
      } finally {
        full.close();
      }
    }
  }));

  test('the Reaction panel offers atom mapping exactly when Chem\'s Crux file has the mapper', async () => {
    await inCrux(async () => {
      const panel = reactionPanel(corpus.col('reaction')!.get(0));
      const headers = [...panel.root.querySelectorAll('h3')].map((h) => h.textContent);
      expect(headers.includes('Mapping'), mapperExported());
    });
  });

  test('the Reaction panel maps the reaction by itself, on a worker, and draws it mapped: no button', async () => {
    const reaction = 'CC(=O)Cl.OCC>>CCOC(C)=O';
    const maps = (s: string) => [...s.matchAll(/:(\d+)\]/g)].map((m) => m[1]);
    await inCrux(async () => {
      expect(mapperExported(), true, 'the mapper in Chem\'s Crux file');
      const {mapped, confidence} = await mapReaction(reaction);
      const [reactants, products] = mapped.split('>>');
      expect(maps(products).length, 6, `every heavy atom of the ester mapped: ${mapped}`);
      expect(maps(products).every((m) => maps(reactants).includes(m)), true, `each from a reactant: ${mapped}`);
      expect(confidence > 0 && confidence <= 1, true, `confidence ${confidence}`);
      const panel = reactionPanel(reaction);
      const host = ui.dialog('Reaction panel test').add(panel.root).show();
      try {
        // the panel runs the mapper by itself, with no button
        expect([...panel.root.querySelectorAll('button')].some((b) => /map atoms/i.test(b.textContent ?? '')), false,
          'a Map atoms button');
        const {state, text} = await mappingAnswer(panel);
        expect(state, 'done', text);
        expect(text.includes('Mapped automatically (the mapper is'), true, text);
        const picture = panel.root.querySelector('canvas')!;
        expect(/:\d+\]/.test(picture.dataset.reaction ?? ''), true,
          `the picture drawn mapped: ${picture.dataset.reaction}`);
        // drawn on its first size change, a frame after the swap
        await awaitCheck(() => ink(picture) > 200, 'the mapped picture', 3000);
      } finally {
        host.close();
      }
    });
  });

  test('the Reaction panel maps a route step by step: RXNMapper takes one reaction at a time', async () => {
    // Map atoms on a route said "not a valid reaction string" (the whole route was sent to the mapper); each step is
    // mapped on its own and drawn mapped under Steps (by itself, with no button)
    const route = 'CC(=O)O.Cl>>CC(=O)Cl>>CC(=O)Cl.N>>CC(=O)N';
    await inCrux(async () => {
      expect(mapperExported(), true, 'the mapper in Chem\'s Crux file');
      const panel = reactionPanel(route);
      const host = ui.dialog('Reaction panel test').add(panel.root).show();
      try {
        // until both steps are answered, or the whole route was refused at once (the bug: "Not mapped: …")
        const {text} = await mappingAnswer(panel);
        expect(text.includes('Step 1: mapped automatically (the mapper is'), true, text);
        expect(text.includes('Step 2: mapped automatically (the mapper is'), true, text);
        const canvases = [...panel.root.querySelectorAll('canvas')];
        expect(canvases.length, 3, 'the route\'s picture and the two mapped steps');
        // a mapped step's canvas draws on its first size change, a frame after the swap (read sooner, it is blank)
        await awaitCheck(() => canvases.every((c) => ink(c) > 100), 'every picture drawn', 3000);
      } finally {
        host.close();
      }
    });
  });

  // the panel maps by itself (5 s a reaction, no button; errors caught; never a template); each reaction's 5 s counts
  // from the mapper's one-time load, which has a cap of its own (30 s)

  /** `run(calls)` with the panels' mapper replaced by `answer` (each call's reaction kept in `calls`), and no load. */
  const withMapper = async (answer: (reaction: string) => Promise<{mapped: string, confidence: number}>,
    run: (calls: string[]) => Promise<void>) => {
    const calls: string[] = [];
    setReactionMapperForTests((reaction) => {
      calls.push(reaction);
      return answer(reaction);
    });
    try {
      await inCrux(() => run(calls));
    } finally {
      setReactionMapperForTests(null);
    }
  };
  /** A panel of `reaction` shown in a dialog while `run(panel)` runs; the dialog closed and the widget let go after. */
  const shownPanel = async (reaction: string, run: (panel: DG.Widget) => Promise<void>) => {
    const panel = reactionPanel(reaction);
    const host = ui.dialog('Reaction panel test').add(panel.root).show();
    try {
      await run(panel);
    } finally {
      host.close();
      panel.detach();
    }
  };
  const PLAIN = 'CC(=O)Cl.OCC>>CCOC(C)=O';

  // a route the panel maps step by step was passed on unmapped: its picture stayed uncoloured, its legend sat at the
  // panel's end under the steps, full screen's legend at the end of its list of steps, and the Export panel's Map
  // numbers and Reaction centre colours did nothing. Now the route's picture is drawn again from its mapped steps, its
  // one legend right under it, centred on it, and full screen and Export are given the mapped route
  test('a route mapped step by step: its picture drawn mapped, one legend right under it; full screen and Export get ' +
    'it mapped', async () => {
    let mappedSteps: string[] = [];
    let asked = 0;
    await withMapper(async () => ({mapped: mappedSteps[asked++], confidence: 0.9}), async () => {
      mappedSteps = cruxReactionSteps(defaultExample().text) ?? [];
      expect(mappedSteps.length > 1, true, 'the default route\'s steps, mapped');
      const route = mappedSteps.map((s) => s.replace(/:\d+\]/g, ']')).join('--**--');
      const isMapped = (text?: string) => /:\d+\]/.test(text ?? '');
      expect(isMapped(route), false, 'the route asked for, unmapped');
      // the legend right under the first picture, centred on it, and no other; the picture drawn mapped
      const underPicture = (root: HTMLElement, where: string) => {
        const legends = (Array.from(root.querySelectorAll(`.${LEGEND_CLASS}`)) as HTMLElement[])
          .filter((l) => l.style.display !== 'none');
        expect(legends.length, 1, `${where}: one legend shown`);
        const picture = legends[0].previousElementSibling as HTMLCanvasElement | null;
        expect(picture !== null && picture === root.querySelector('canvas'), true,
          `${where}: right after the route's picture, the first`);
        expect(isMapped(picture!.dataset.reaction), true, `${where}: the route's picture drawn mapped`);
        const [p, l] = [picture!, legends[0]].map((e) => e.getBoundingClientRect());
        expect(l.top >= p.bottom - 1 && l.top - p.bottom <= PICTURE_LEGEND_GAP + 1 &&
          Math.abs((l.left + l.right) / 2 - (p.left + p.right) / 2) <= 2, true,
        `${where}: right under the picture, centred on it: ${JSON.stringify(l)} under ${JSON.stringify(p)}`);
      };
      await shownPanel(route, async (panel) => {
        const {state, text} = await mappingAnswer(panel);
        expect(`${state} ${asked}`, `done ${mappedSteps.length}`, text);
        await awaitCheck(() => isMapped(panel.root.querySelector('canvas')?.dataset.reaction), 'drawn mapped', 3000);
        underPicture(panel.root, 'the panel');
        (panel.root.querySelector('.fa-expand-alt') as HTMLElement).click();
        const full = DG.Dialog.getOpenDialogs().find((d) => d.title === 'Reaction')!;
        try {
          await awaitCheck(() => ((full.root.querySelector(`.${LEGEND_CLASS}`) as HTMLElement | null)?.offsetHeight ??
            0) > 0, 'full screen\'s legend', 5000);
          underPicture(full.root, 'full screen');
        } finally {
          full.close();
        }
        (panel.root.querySelector(`.${EXPORT_ICON_CLASS}`) as HTMLElement).click();
        const exported = DG.Dialog.getOpenDialogs().find((d) => d.title === EXPORT_TITLE)!;
        try {
          const inputs = exportInputsOf(exported.root)!;
          expect(`${inputs.colours.enabled} ${inputs.mapNumbers.enabled}`, 'true true',
            'Export: the mapped route\'s Reaction centre colours and Map numbers apply');
          inputs.colours.value = true;
          const preview = exported.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
          await awaitCheck(() => preview.dataset.legend === 'true', 'Export: the colours on, the legend drawn', 5000);
        } finally {
          exported.close();
        }
      });
    });
  });

  // C1CCNCC1>>CC(=O)N1CCCCC1, an N-acetylation with no acetyl source drawn, was mapped with a ring
  // carbon on the acetyl CH3, 63% sure, over the 50% warning: nothing said. Products holding heavy atoms that no
  // reactant or reagent drawn holds are said, however sure the mapper is; with the reagent over the arrow, nothing more
  test('a mapping whose products hold atoms no reactant or reagent drawn holds says check it, however sure',
    async () => {
      const unbalanced = 'C1CCNCC1>>CC(=O)N1CCCCC1';
      await withMapper(async (r) => ({mapped: r, confidence: 0.63}), async () => {
        expect(productAtomsMissing(unbalanced), 'C 2, O 1', 'the acetyl\'s atoms');
        expect(productAtomsMissing('C1CCNCC1>CC(=O)OC(C)=O>CC(=O)N1CCCCC1'), null, 'acetic anhydride over the arrow');
        expect(productAtomsMissing(PLAIN), null, 'balanced');
        for (const [reaction, warned] of [[unbalanced, true], [PLAIN, false]] as const) {
          await shownPanel(reaction, async (panel) => {
            const {state, text} = await mappingAnswer(panel);
            expect(state, 'done', text);
            expect(text.includes('Check this mapping: the reactants drawn don\'t hold every product atom (C 2, O 1 ' +
              'more in the products)'), warned, `${reaction}: ${text}`);
            expect(panel.root.querySelector(`.${MAPPING_CLASS} [data-level="warn"]`) !== null, warned,
              `${reaction}: an orange line`);
          });
        }
      });
    });

  // reaction SMILES writes no stoichiometry, so a reactant used twice (a
  // self-condensation, a homocoupling, a double acylation) leaves the products with its atoms once more: the line asks
  // whether a reactant is used more than once, not for a reagent that is drawn already. The N-acetylation above (the
  // atoms more are no molecule's drawn) still asks for the reagent.
  test('products holding one reactant\'s atoms once more ask whether it is used more than once, not for a reagent',
    async () => {
      const aldol = 'CC=O>>CC(O)CC=O';
      expect(productAtomsMissing(aldol), 'C 2, O 1|repeated', 'a self-condensation: the aldehyde\'s atoms again');
      expect(productAtomsMissing('Brc1ccccc1>>c1ccc(-c2ccccc2)cc1'), 'C 6|repeated',
        'a homocoupling: the aryl again, its bromide left');
      expect(productAtomsMissing('NCCN.CC(=O)Cl>>CC(=O)NCCNC(C)=O'), 'C 2, O 1|repeated',
        'a double acylation: the acyl chloride again');
      expect(productAtomsMissing('C1CCNCC1>>CC(=O)N1CCCCC1'), 'C 2, O 1', 'EVAL5: the acetyl is no molecule drawn');
      expect(productAtomsMissing('CO>>CCO'), 'C 1', 'one carbon more is no methanol again: its O is not');
      await withMapper(async (r) => ({mapped: r, confidence: 0.9}), async () => {
        for (const [reaction, words] of [[aldol, 'is a reactant used more than once?'],
          ['C1CCNCC1>>CC(=O)N1CCCCC1', 'add the reagent that brings them']] as const) {
          await shownPanel(reaction, async (panel) => {
            const {state, text} = await mappingAnswer(panel);
            expect(state, 'done', text);
            expect(text.includes('Check this mapping: the reactants drawn don\'t hold every product atom (C 2, O 1 ' +
              `more in the products): ${words}`), true, `${reaction}: ${text}`);
            expect(text.includes('|'), false, `${reaction}: no detail mark in the words: ${text}`);
          });
        }
      });
    });

  test('the Reaction panel does not map a template: the mapper is never asked', async () => {
    const template = corpus.col('reaction')!.get(0); // dgdef-suzuki-coupling, a SMIRKS with map numbers
    await withMapper(async (r) => ({mapped: r, confidence: 1}), async (calls) => {
      await shownPanel(template, async (panel) => {
        const {state, text} = await mappingAnswer(panel, 2000);
        await DG.delay(1000);
        expect(state, 'refused', text);
        expect(text.includes('Not mapped: a template keeps the map numbers written in it'), true, text);
        expect(calls.length, 0, `the mapper asked for: ${calls.join(' | ')}`);
      });
    });
  });

  test('the Reaction panel keeps a mapping the value has: drawn as written, the mapper never asked', async () => {
    // automatic mapping must not replace a mapping the value already has (any mapped atom)
    const own = '[CH3:1][C:2](=[O:3])Cl.[OH:4][CH2:5][CH3:6]>>[CH3:1][C:2](=[O:3])[O:4][CH2:5][CH3:6]';
    const partly = 'CC(=O)Cl.[OH:4]CC>>CC(=O)[O:4]CC';
    expect(mappedInTheValue(own), true, 'a mapped SMILES');
    expect(mappedInTheValue(partly), true, 'one mapped atom is enough');
    expect(mappedInTheValue(PLAIN), false, 'an unmapped SMILES');
    expect(mappedInTheValue('[CH3:0]C>>CC |(0,0,;1,0,;0,0,;1,0,)|'), false, 'map number 0, a CX extension');
    await withMapper(async (r) => ({mapped: `${r} |mapper|`, confidence: 1}), async (calls) => {
      for (const value of [own, partly]) {
        await shownPanel(value, async (panel) => {
          const {state, text} = await mappingAnswer(panel, 2000);
          await DG.delay(1000); // past the panel's delay: still not asked
          expect(state, 'refused', text);
          expect(text, 'Mapped in the value itself');
          expect(calls.length, 0, `the mapper asked for: ${calls.join(' | ')}`);
          const picture = panel.root.querySelector('canvas')!;
          expect(picture.dataset.reaction, value, 'the value\'s own mapping drawn');
          expect(ink(picture) > 200, true, 'the picture drawn');
        });
      }
    });
  });

  test('the Reaction panel shows a failing mapper\'s words; nothing thrown, the picture kept', async () => {
    const thrown: string[] = [];
    const onError = (e: ErrorEvent | PromiseRejectionEvent) =>
      thrown.push(String((e as any).reason ?? (e as any).message));
    window.addEventListener('error', onError);
    window.addEventListener('unhandledrejection', onError);
    try {
      await withMapper(() => Promise.reject(new Error('simulated: the mapper failed')), async (calls) => {
        await shownPanel(PLAIN, async (panel) => {
          const {state, text} = await mappingAnswer(panel);
          expect(state, 'failed', text);
          // in words, the mapper's own text to the console (no raw engine text in the UI)
          expect(text, MAPPER_FAILED_WORDS);
          expect(text.includes('simulated'), false, `the mapper's text in the panel: ${text}`);
          expect(calls.length, 1, 'mappings asked');
          expect(panel.root.querySelector('canvas')!.dataset.reaction, PLAIN, 'the picture kept unmapped');
        });
      });
    } finally {
      window.removeEventListener('error', onError);
      window.removeEventListener('unhandledrejection', onError);
    }
    expect(thrown.length, 0, `errors that escaped: ${thrown.join(' | ')}`);
  });

  test('the Reaction panel says in words that a step with no product or reactant cannot be mapped, never asking the ' +
    'mapper for it', async () => {
    // after a step whose product side was empty, the panel said `Step 1: not mapped: ""
    // is not a valid SMILES string`, the mapper's own text; such a step is skipped, said in words, and no error
    const empty = '"" is not a valid SMILES string'; // what the engine's mapper says of an empty side
    const engineLike = (r: string) => r.split(' |')[0].split('>').some((side, i) => i !== 1 && side.trim() === '') ?
      Promise.reject(new Error(empty)) : Promise.resolve({mapped: r, confidence: 1});
    await withMapper(engineLike, async (calls) => {
      // a route whose last step has no product (a chain ending in ">>": Crux reads `CCO>>CC=O` and `CC=O>>`)
      await shownPanel('CCO>>CC=O>>', async (panel) => {
        const {state, text} = await mappingAnswer(panel);
        expect(text.includes('Step 1: mapped automatically'), true, text);
        expect(text.includes(noSideWords('product', 2)), true, text);
        expect(text.includes(empty), false, `the engine's text in the panel: ${text}`);
        expect(state, 'done', `a step skipped is no failure: ${text}`);
        expect(calls.join(' | '), 'CCO>>CC=O', 'the mapper asked for step 1 alone');
      });
      calls.length = 0;
      for (const [reaction, side] of [['CCO>>', 'product'], ['CCO>O>', 'product'], ['>>CC=O', 'reactant']] as const) {
        await shownPanel(reaction, async (panel) => {
          const {state, text} = await mappingAnswer(panel, 2000);
          expect(text, noSideWords(side), reaction);
          expect(state, 'refused', `${reaction}: ${text}`);
        });
      }
      // a value Crux reads no reaction from (its |f:| groups name components out of range): the mapper was once
      // asked and the panel said "Not mapped: Cannot read properties of undefined (reading 'split')"
      await shownPanel('CC>>CC |f:0.1.2|', async (panel) => {
        const {state, text} = await mappingAnswer(panel, 2000);
        expect(text.startsWith(`${NOT_READ_WORDS}: its salt and mixture groups`) && !text.includes('f:['), true,
          `unread, and why in words: ${text}`);
        expect(state, 'refused', text);
      });
      // a text typed wrong was called a template ("This is a template: ... Map atoms maps
      // reactions, not templates"); it is unread, and the panel says so and why (Crux's words, as the app's status)
      await shownPanel('C1CC(>>N#Q[Zz]', async (panel) => {
        const {state, text} = await mappingAnswer(panel, 2000);
        expect(text.startsWith(`${NOT_READ_WORDS}: `) && text.length > NOT_READ_WORDS.length + 3, true,
          `unread, and why: ${text}`);
        expect(/template/i.test(panel.root.textContent ?? ''), false, `no template words: ${text}`);
        expect(state, 'refused', text);
      });
      expect(calls.length, 0, `the mapper never asked for a reaction with an empty side: ${calls.join(' | ')}`);
    });
  });

  test('the Reaction panel gives up a mapping after 5 s and ignores its late answer', async () => {
    const late = MAP_TIMEOUT_MS + 1500;
    await withMapper((r) => new Promise((resolve) => setTimeout(() => resolve({mapped: `${r} |late|`, confidence: 1}),
      late)), async () => {
      await shownPanel(PLAIN, async (panel) => {
        const start = Date.now();
        const {state, text} = await mappingAnswer(panel, MAP_TIMEOUT_MS + 3000);
        const ms = Date.now() - start;
        expect(state, 'timeout', text);
        expect(text, `Mapping took longer than ${MAP_TIMEOUT_MS / 1000} s`);
        expect(ms < late, true, `given up after ${ms} ms (the answer comes at ${late} ms)`);
        await DG.delay(late - ms + 500);
        expect(panel.root.querySelector('canvas')!.dataset.reaction, PLAIN, 'the late answer ignored');
        expect(panel.root.querySelector(`.${MAPPING_CLASS}`)!.textContent, text, 'the line kept');
      });
    });
  });

  test('the Reaction panel says the mapper could not be loaded when its load takes over 30 s', async () => {
    let loads = 0;
    setMapperLoaderForTests(() => {
      loads++;
      return new Promise(() => {}); // never loads
    });
    try {
      await withMapper(async (r) => ({mapped: r, confidence: 1}), async (calls) => {
        await shownPanel(PLAIN, async (panel) => {
          await DG.delay(1000);
          const status = panel.root.querySelector(`.${MAPPING_CLASS}`) as HTMLElement;
          expect(status.dataset.state, 'loading', status.textContent ?? '');
          expect(status.textContent, 'Loading the mapper…');
          const {state, text} = await mappingAnswer(panel, MAP_LOAD_TIMEOUT_MS + 5000);
          expect(state, 'failed', text);
          expect(text, 'The mapper could not be loaded');
          expect(calls.length, 0, 'mappings asked of a mapper not loaded');
          expect(loads, 1, 'loads asked');
        });
      });
    } finally {
      setMapperLoaderForTests(null);
    }
  }, {timeout: MAP_LOAD_TIMEOUT_MS + 30000});

  test('rows passed by map nothing, and a panel gone mid-mapping is not updated', async () => {
    let release: () => void = () => {};
    await withMapper((r) => new Promise((resolve) => {
      release = () => resolve({mapped: `${r} |mapped|`, confidence: 1});
    }), async (calls) => {
      // twenty rows passed by, each panel let go before its turn (the context panel lets go of the last one's)
      for (let i = 0; i < 20; i++)
        reactionPanel(`C${'C'.repeat(i)}O>>C${'C'.repeat(i)}=O`).detach();
      await DG.delay(1000);
      expect(calls.length, 0, `mappings asked for rows passed by: ${calls.join(' | ')}`);
      // a panel let go while its mapping runs: the answer updates nothing
      const panel = reactionPanel(PLAIN);
      for (let i = 0; i < 40 && calls.length === 0; i++)
        await DG.delay(50);
      expect(calls.length, 1, 'the panel\'s mapping asked');
      const status = panel.root.querySelector(`.${MAPPING_CLASS}`) as HTMLElement;
      const before = `${status.dataset.state}: ${status.textContent}`;
      panel.detach();
      release();
      await DG.delay(300);
      expect(`${status.dataset.state}: ${status.textContent}`, before, 'the Mapping section after its panel went');
      expect(panel.root.querySelector('canvas')!.dataset.reaction, PLAIN, 'the picture after its panel went');
    });
  });

  test('full screen shows the reaction as the panel shows it: mapped once mapped', async () => {
    // after Map atoms, the full-screen icon opened the reaction unmapped
    const mapped = (s: string | undefined) => /:\d+\]/.test(s ?? '');
    const fullScreenShows = async (panel: DG.Widget): Promise<string[]> => {
      (panel.root.querySelector('.fa-expand-alt') as HTMLElement).click();
      await DG.delay(300);
      const full = DG.Dialog.getOpenDialogs().find((d) => d.title === 'Reaction')!;
      try {
        return [...full.root.querySelectorAll('canvas')].map((c) => (c as HTMLCanvasElement).dataset.reaction ?? '');
      } finally {
        full.close();
      }
    };
    const mapAll = async (panel: DG.Widget, done: RegExp) => {
      const {text} = await mappingAnswer(panel); // mapped by itself
      expect(done.test(text), true, text);
    };
    await inCrux(async () => {
      expect(mapperExported(), true, 'the mapper in Chem\'s Crux file');
      // one reaction: full screen draws the mapped reaction
      const single = reactionPanel('CC(=O)Cl.OCC>>CCOC(C)=O');
      let host = ui.dialog('Reaction panel test').add(single.root).show();
      try {
        await DG.delay(300);
        await mapAll(single, /Mapped automatically|Not mapped:/); // the panel's words: no "confidence"
        const drawn = await fullScreenShows(single);
        expect(mapped(drawn[0]), true, `full screen draws: ${drawn[0]}`);
      } finally {
        host.close();
      }
      // a route: full screen draws the route and its steps, each step mapped as in the panel
      const route = reactionPanel('CC(=O)O.Cl>>CC(=O)Cl>>CC(=O)Cl.N>>CC(=O)N');
      host = ui.dialog('Reaction panel test').add(route.root).show();
      try {
        await DG.delay(300);
        await mapAll(route, /Step 2: (mapped|not mapped)|Not mapped:/);
        const drawn = await fullScreenShows(route);
        expect(drawn.length, 3, `the route and its two steps: ${drawn.join(' | ')}`);
        expect(drawn.slice(1).every(mapped), true, `the steps mapped: ${drawn.join(' | ')}`);
      } finally {
        host.close();
      }
    });
  });

  test('the Reaction panel lists a route\'s steps', async () => {
    const route = samples.col('reaction')!.toList().find((s: string) => s.includes('>>CC(=O)Cl>>'))!;
    await inCrux(async () => {
      const panel = reactionPanel(route);
      const host = ui.dialog('Reaction panel test').add(panel.root).show();
      try {
        await DG.delay(300);
        const labels = [...panel.root.querySelectorAll('label')].map((l) => l.textContent);
        expect(labels.filter((l) => l?.startsWith('Step ')).length > 1, true, labels.join(', '));
        const canvases = [...panel.root.querySelectorAll('canvas')];
        expect(canvases.every((c) => ink(c) > 100), true, 'every step drawn');
      } finally {
        host.close();
      }
    });
  });

  // the same intermediate turned over from step to step, the cards at scales of their own
  test('a route\'s step cards are drawn from the route, its step each, in one box, where Crux draws a step alone; ' +
    'else each step by itself', async () => {
    const route = 'CC(=O)O.Cl>>CC(=O)Cl>>CC(=O)Cl.N>>CC(=O)N';
    await inCrux(async () => {
      setReactionMapperForTests(() => new Promise(() => {})); // the cards as first drawn: no mapping swaps them
      try {
        const own = cruxDrawsSteps();
        let cards: HTMLCanvasElement[] = [];
        const sent = await sentOptions(async () => {
          const panel = reactionPanel(route);
          const host = ui.dialog('Reaction panel test').add(panel.root).show();
          try {
            await awaitCheck(() => [...panel.root.querySelectorAll('canvas')].every((c) => ink(c) > 100),
              'every picture drawn', 5000);
            cards = [...panel.root.querySelectorAll('canvas')].slice(1) as HTMLCanvasElement[];
          } finally {
            host.close();
            panel.detach();
          }
        });
        expect(cards.map((c) => c.dataset.reaction).join(' | '), 'CC(=O)O.Cl>>CC(=O)Cl | CC(=O)Cl.N>>CC(=O)N',
          'each card says its step');
        const stepped = sent.filter((o) => o.step !== undefined);
        if (own) {
          expect(stepped.map((o) => o.step).join(), '0,1', 'the cards: the route\'s steps 0 and 1');
          expect(new Set(stepped.map((o) => `${o.width}x${o.height}`)).size, 1, 'one box for every card');
          expect(stepped.every((o) => o.oneScale === undefined), true, 'one scale (the engine\'s default)');
        } else
          expect(stepped.length, 0, 'an engine before the step option: each step drawn by itself');
      } finally {
        setReactionMapperForTests(null);
      }
    });
  });

  // where every card is drawn from the route mapped so far, drawing all the cards again after each
  // mapped step drew a 5-step route's cards 5 times (~85 ms); they are drawn again once. An engine before it draws
  // only the step mapped, by itself: one card a step, as before. The route's own picture is drawn again once with them,
  // from its mapped steps
  test('a mapping pass draws each step card once: the cards once at its end where Crux draws a step alone, else ' +
    'each mapped step\'s own', async () => {
    const route = 'CCO>>CC=O>>CC(=O)O>>CC(=O)OC';
    let release = () => {};
    const gate = new Promise<void>((resolve) => release = resolve);
    await withMapper(async (r) => {
      await gate;
      return {mapped: r, confidence: 0.9};
    }, async (calls) => {
      const steps = cruxReactionSteps(route) ?? [];
      expect(steps.length, 3, `the route's steps: ${steps.join(' | ')}`);
      const scheme: any = (cruxModuleSync() as unknown as ChemWasm).ChemReactionScheme.prototype;
      const draw = scheme.draw_to_canvas;
      const drawn: HTMLCanvasElement[] = []; // the cards drawn after the first pictures, a draw each
      let counting = false;
      scheme.draw_to_canvas = function(this: any, canvas: any, options: string) {
        if (counting && canvas instanceof HTMLCanvasElement && canvas.dataset.reaction !== route)
          drawn.push(canvas);
        return draw.call(this, canvas, options);
      };
      try {
        await shownPanel(route, async (panel) => {
          const canvases = () => [...panel.root.querySelectorAll('canvas')] as HTMLCanvasElement[];
          await awaitCheck(() => canvases().length === 4 && canvases().every((c) => ink(c) > 100),
            'the route\'s picture and its 3 cards drawn', 5000);
          await awaitCheck(() => calls.length === 1, 'the mapper asked for step 1', 5000);
          counting = true;
          release();
          const {state, text} = await mappingAnswer(panel);
          expect(state, 'done', text);
          await awaitCheck(() => canvases().slice(1).every((c) => ink(c) > 100), 'the mapped cards drawn', 3000);
          await DG.delay(500); // nothing more to come: a late redraw would be counted
          const [picture, ...cards] = canvases();
          // the route's picture is drawn from its steps joined (`--**--`), a card from its step
          const routeDraws = drawn.filter((c) => (c.dataset.reaction ?? '').includes('--**--'));
          const cardDraws = drawn.filter((c) => !routeDraws.includes(c));
          expect(new Set(cardDraws).size, 3,
            `cards drawn in the pass: ${cardDraws.length} draws of ${new Set(cardDraws).size}`);
          expect(cardDraws.every((c) => cards.includes(c)), true, 'only the cards shown drawn');
          expect(`${routeDraws.length} ${routeDraws.includes(picture)}`, `${cruxDrawsSteps() ? 1 : 3} true`,
            'the route\'s picture drawn again from its mapped steps: once, or once a step mapped');
        });
      } finally {
        scheme.draw_to_canvas = draw;
      }
      expect(calls.length, 3, 'each step mapped');
    });
  });

  test('the Reaction panel lists the steps Crux reads: a restated intermediate is no step', async () => {
    // test-reactions.csv's two-step chain: cut pairwise (parseMultiStepReaction) the panel listed three steps, the
    // middle one "acyl chloride → acyl chloride + NH3", under a picture of two
    const route = 'CC(=O)O.Cl>>CC(=O)Cl>>CC(=O)Cl.N>>CC(=O)N';
    expect(samples.col('reaction')!.toList().includes(route), true, 'the chain in test-reactions.csv');
    await inCrux(async () => {
      expect(JSON.stringify(cruxReactionSteps(route)), JSON.stringify(['CC(=O)O.Cl>>CC(=O)Cl', 'CC(=O)Cl.N>>CC(=O)N']),
        'Crux\'s steps as the value writes them (ChemReactionScheme.get_steps; verbatim since crux-core a9fc04786)');
      const panel = reactionPanel(route);
      const host = ui.dialog('Reaction panel test').add(panel.root).show();
      try {
        await DG.delay(300);
        const labels = [...panel.root.querySelectorAll('label')].map((l) => l.textContent ?? '');
        expect(labels.filter((l) => l.startsWith('Step ')).join(', '), 'Step 1, Step 2');
        const canvases = [...panel.root.querySelectorAll('canvas')];
        expect(canvases.length, 3, 'the picture and the two steps');
        expect(canvases.every((c) => ink(c) > 100), true, 'every step drawn');
      } finally {
        host.close();
      }
    });
  });

  test('a mapped reaction keeps its agents over the arrow', async () => {
    // once mapped, a step drew H2O, H2, Pd and EtOH as reactants with plus signs (RXNMapper
    // writes every left-hand component as a reactant; one with no mapped atom is a reagent)
    expect(agentsBackOverTheArrow('[CH3:1][C:2](=[O:3])[OH:4].O.[Pd]>>[CH3:1][C:2](=[O:3])[O-:4]'),
      '[CH3:1][C:2](=[O:3])[OH:4]>O.[Pd]>[CH3:1][C:2](=[O:3])[O-:4]', 'unmapped components to the agents');
    expect(agentsBackOverTheArrow('[CH3:1][OH:2].O>[H+]>[CH3:1][OH:2]'), '[CH3:1][OH:2]>[H+].O>[CH3:1][OH:2]',
      'after the agents the text has');
    expect(agentsBackOverTheArrow('CC>>CC'), 'CC>>CC', 'nothing mapped: as it is');
    expect(agentsBackOverTheArrow('[CH3:1][CH3:2]>>[CH3:1][CH3:2]'), '[CH3:1][CH3:2]>>[CH3:1][CH3:2]', 'no agents');
    // the app's Isotopes example, mapped: its salt's fragment group kept f:0.1.2 once the salt went over
    // the arrow, a group over the phenol and the iodomethane, which Crux refused, and the app's pane wrote the
    // reaction as text; the groups are numbered with their components
    const salt = 'O=C([O-])[O-].[K+].[K+].[OH:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1.[I:1][C:3]([2H:2])([2H:4])' +
      '[2H:5]>>[IH:1].[2H:2][C:3]([2H:4])([2H:5])[O:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1 |f:0.1.2|';
    const back = agentsBackOverTheArrow(salt);
    expect(back, '[OH:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1.[I:1][C:3]([2H:2])([2H:4])[2H:5]>O=C([O-])[O-].' +
      '[K+].[K+]>[IH:1].[2H:2][C:3]([2H:4])([2H:5])[O:6][c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1 |f:2.3.4|',
    'the salt over the arrow, its group with it');
    expect(cruxReadsReaction(back), true, 'Crux reads it');
    expect(agentsBackOverTheArrow('[Na+].[Cl-].[CH3:1][OH:2]>>[CH3:1][OH:2].[Na+].[Cl-] |f:0.1,3.4|'),
      '[CH3:1][OH:2]>[Na+].[Cl-]>[CH3:1][OH:2].[Na+].[Cl-] |f:1.2,3.4|', 'two groups, the products\' as they were');
    expect(agentsBackOverTheArrow('O.[CH3:1][OH:2]>>[CH3:1][OH:2] |c:1|'), 'O.[CH3:1][OH:2]>>[CH3:1][OH:2] |c:1|',
      'a CX section numbering atoms: as it is');
    if (!mapperExported())
      return;
    // a step of the demo's paracetamol route: the nitro group reduced over H2 and Pd
    const {mapped} = await mapReaction('O=[N+]([O-])c1ccc(O)cc1.[H][H].[Pd]>>Nc1ccc(O)cc1');
    const [left, middle] = mapped.split('>');
    expect(left.split('.').every((c) => /:[1-9]\d*\]/.test(c)), true, `every reactant mapped: ${mapped}`);
    expect(middle.length > 0, true, `the agents over the arrow: ${mapped}`);
  });

  // the sketchers' Auto-map is Crux Sketch's mapReaction hook answered by Chem's mapper, one step at a time, null
  // for a step it does not map (a failure, a late answer), as the hook asks; a mapping under 50 % sure is said
  test('Auto-map in the sketchers: Chem\'s mapper maps a step; a failing or late mapper maps nothing', async () => {
    const ester = 'CC(=O)O.OCC>>CC(=O)OCC';
    if (mapperExported()) {
      const hook = sketchReactionMapper();
      expect(hook !== null, true, 'a hook while Chem\'s Crux file has the mapper');
      const mapped = await hook!(ester);
      expect(/:[1-9]\d*\]/.test(mapped ?? '') && (mapped ?? '').includes('>'), true, `mapped: ${mapped}`);
    }
    const warned: string[] = [];
    const balloon = grok.shell.warning;
    grok.shell.warning = (s: string | HTMLElement) => warned.push(String(s));
    try {
      setSketchMapperForTests(async () => {
        throw new Error('the mapper failed (test)');
      });
      expect(await sketchReactionMapper()!(ester), null, 'a failing mapper');
      setSketchMapperForTests(() => new Promise((resolve) =>
        setTimeout(() => resolve({mapped: '[CH3:1][OH:2]>>[CH3:1][OH:2]', confidence: 0.9}), MAP_TIMEOUT_MS + 500)));
      expect(await sketchReactionMapper()!(ester), null, 'a late mapper');
      setSketchMapperForTests(async () => ({mapped: '[CH3:1][OH:2]>>[CH2:1]=[O:2]', confidence: 0.42}));
      expect(await sketchReactionMapper()!(ester), '[CH3:1][OH:2]>>[CH2:1]=[O:2]', 'an unsure mapping is kept');
      expect(warned.length === 1 && warned[0].includes('42% sure'), true, `said: ${warned.join(' | ')}`);
    } finally {
      grok.shell.warning = balloon;
      setSketchMapperForTests(null);
    }
  });

  test('the Reaction panel draws its picture larger than a cell', async () => {
    // a small reaction's "large picture" was no larger than its cell (Crux draws at 22 px
    // a bond unless asked, and never enlarges)
    await inCrux(async () => {
      let host: DG.Dialog | null = null;
      const sent = await sentOptions(async () => {
        const panel = reactionPanel('CC(=O)O.OCC>>CC(=O)OCC');
        host = ui.dialog('Reaction panel test').add(panel.root).show();
        await DG.delay(300);
      });
      (host as DG.Dialog | null)?.close();
      expect(sent.length > 0, true, 'the picture drawn');
      expect(sent[0].bondLength > 22, true, `the picture's bond length: ${sent[0].bondLength}`);
    });
  });

  test('a template\'s note words are in the Reaction panel; its cell keeps only the letters', async () => {
    // a wide, low template cell wrapped its footnotes to the scheme's width; that information belongs in the context
    // panel when the reaction is clicked, not in the cell: the small letters beside the atoms stay in the cell, only
    // the words go to the panel (crux-core 3f4d88a89).
    const template = '[C:1](=[O:2])[OH].[N;H1,H2:3]>>[C:1](=[O:2])[N:3]';
    const table = DG.DataFrame.fromColumns([DG.Column.fromStrings('reaction', [template])]);
    table.col('reaction')!.semType = DG.SEMTYPE.CHEMICAL_REACTION;
    const grid = table.plot.grid();
    await inCrux(async () => {
      const notes = cruxReactionNotes(template);
      expect(notes.length > 0, true, 'Crux has notes for the template (ChemReactionScheme.get_notes)');
      expect(cruxReactionNotes('CC(=O)O.OCC>>CC(=O)OCC').length, 0, 'notes of a plain reaction');
      const sent = await sentOptions(async () => {
        drawCell(new RDKitReactionRenderer(getRdKitModule()), grid, 'reaction', 0);
      });
      expect(sent.length > 0 && sent.every((o) => o.footnotes === false), true,
        `the cell's options: ${JSON.stringify(sent.map((o) => o.footnotes))}`);
      // what Crux draws with the cell's options: the letters beside the atoms (atom-note), no footnote words
      // (reaction-note)
      const crux = cruxModuleSync() as unknown as ChemWasm;
      const scheme = new crux.ChemReactionScheme(template, '');
      try {
        const svg: string = scheme.get_svg(JSON.stringify({...sent[0], offsetx: 0, offsety: 0}));
        expect(svg.includes('reaction-note'), false, 'a footnote\'s words in the cell');
        expect(svg.includes('atom-note'), true, 'the footnote letter beside its atom in the cell');
      } finally {
        scheme.free();
      }
      const panel = reactionPanel(template);
      const host = ui.dialog('Reaction panel test').add(panel.root).show();
      try {
        const list = panel.root.querySelector('.chem-reaction-notes');
        expect(!!list, true, 'the panel\'s notes');
        const shown = [...list!.children].map((d) => d.textContent);
        expect(JSON.stringify(shown), JSON.stringify(notes), 'the notes the panel lists');
        expect(!!reactionPanel('CC(=O)O.OCC>>CC(=O)OCC').root.querySelector('.chem-reaction-notes'), false,
          'notes for a plain reaction');
      } finally {
        host.close();
      }
    });
  });

  // the reaction previews were RDKit's in both modes, from when Crux drew no reactions; they follow the engine in
  // effect now, as the cells do: Crux in Crux mode (RDKit behind it for a template Crux reads none from), RDKit in
  // RDKit mode

  /** A table of acids and amines (Molecule columns) the reaction dialogs open on. */
  const reactantsTable = (): DG.DataFrame => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('acids', ['OC(=O)c1ccccc1', 'CC(=O)O', 'OC(=O)CNC(=O)OCC1c2ccccc2-c2ccccc21']),
      DG.Column.fromStrings('amines', ['Nc1ccccc1', 'NCc1ccccc1', 'CCN'])]);
    for (const col of df.columns.toList()) {
      col.semType = DG.SEMTYPE.MOLECULE;
      col.meta.units = DG.UNITS.Molecule.SMILES;
    }
    return df;
  };

  /** `run(library)` with the reaction library (files/reactions/default-reactions.json) as the dialogs' and the
   * browser's built-in reactions: this bundle's package object reads no package file (reaction-storage.ts). */
  async function withLibrary<T>(run: (library: NamedReaction[]) => Promise<T>): Promise<T> {
    const library = JSON.parse(await loadFileAsText('reactions/default-reactions.json')) as NamedReaction[];
    setDefaultReactionsForTests(library);
    try {
      return await run(library);
    } finally {
      setDefaultReactionsForTests(null);
    }
  }

  test('the reaction browser\'s cards are drawn by Crux in Crux mode, RDKit drawing none; by RDKit in RDKit mode',
    async () => withLibrary(async (library) => {
      // every card of the library (and the user's own reactions), each template with its variables' defaults
      const cards = async (): Promise<HTMLCanvasElement[]> => {
        const browser = new ReactionBrowser(getRdKitModule(), {allowAdd: false, allowEdit: false});
        await awaitCheck(() => browser.root.querySelectorAll('canvas').length >= library.length,
          `fewer cards than the library's ${library.length} reactions`, 15000);
        return [...browser.root.querySelectorAll('canvas')] as HTMLCanvasElement[];
      };
      await inCrux(async () => {
        const draws = await reactionDraws(() => expectNoFallback([PREVIEW_CALL], cards));
        const drawn = new Set(draws.crux.map((d) => d.canvas));
        expect(draws.result.filter((c) => !drawn.has(c)).length, 0,
          `cards Crux did not draw, of ${draws.result.length}`);
        expect(draws.rdkit, 0, 'RDKit\'s reaction drawings');
        expect(draws.result.filter((c) => ink(c) < 100).length, 0, 'cards with no picture');
      });
      await inEngine(ChemicalEngine.RDKit, async () => {
        const draws = await reactionDraws(cards);
        expect(draws.crux.length, 0, 'reactions Crux drew');
        expect(draws.result.filter((c) => ink(c) < 100).length, 0, 'cards with no picture');
      });
    }));

  test('the reaction editor\'s preview is drawn by Crux in Crux mode, RDKit drawing none; by RDKit in RDKit mode',
    async () => {
      const preset: NamedReaction = {id: 'test-preview', name: 'Preview test', category: 'Test', mode: 'two-component',
        reactionSmarts: '[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]'};
      const preview = async (): Promise<HTMLCanvasElement> => {
        await openReactionEditor(getRdKitModule(), {preset});
        const dialog = await dialogTitled('Edit Reaction');
        try {
          const canvas = dialog.root.querySelector('canvas') as HTMLCanvasElement | null;
          expect(!!canvas, true, 'the editor\'s preview');
          return canvas!;
        } finally {
          dialog.close();
        }
      };
      await inCrux(async () => {
        const draws = await reactionDraws(() => expectNoFallback([PREVIEW_CALL], preview));
        expect(draws.crux.some((d) => d.canvas === draws.result), true, 'Crux drew the preview');
        expect(draws.rdkit, 0, 'RDKit\'s reaction drawings');
        expect(ink(draws.result) > 100, true, 'the preview\'s picture');
      });
      await inEngine(ChemicalEngine.RDKit, async () => {
        const draws = await reactionDraws(preview);
        expect(draws.crux.length, 0, 'reactions Crux drew');
        expect(ink(draws.result) > 100, true, 'the preview\'s picture');
      });
    });

  test('the reaction dialogs\' previews and the full-screen preview are drawn by Crux in Crux mode, RDKit drawing ' +
    'none; by RDKit in RDKit mode', async () => withLibrary(async () => {
    const table = reactantsTable();
    // a template with a variable, drawn with its default (`${smiles}`: CCCC)
    const template = 'Fmoc Substitution';
    // Transformation...: the preview of the picked reaction, then full screen (a click on the preview)
    const transformation = async (): Promise<{preview: HTMLCanvasElement, full: HTMLCanvasElement}> => {
      await transformationReactionsUI(table, table.col('acids'));
      const dialog = await dialogTitled('Run Reaction');
      try {
        const preview = previewCanvasOf(dialog);
        await clickCard(dialog.root, template);
        await awaitCheck(() => ink(preview) > 100, 'the transformation preview drawn', 5000);
        preview.click();
        const full = await dialogTitled(template);
        try {
          const canvas = full.root.querySelector('canvas') as HTMLCanvasElement | null;
          expect(!!canvas && ink(canvas) > 100, true, 'the full-screen picture');
          return {preview, full: canvas!};
        } finally {
          full.close();
        }
      } finally {
        dialog.close();
      }
    };
    // Two-Component Reaction...: the preview of the picked reaction
    const twoComponent = async (): Promise<HTMLCanvasElement> => {
      await twoComponentReactionUI(table);
      const dialog = await dialogTitled('Two-Component Reaction');
      try {
        const preview = previewCanvasOf(dialog);
        await clickCard(dialog.root, 'Amide Coupling');
        await awaitCheck(() => ink(preview) > 100, 'the two-component preview drawn', 5000);
        return preview;
      } finally {
        dialog.close();
      }
    };
    await inCrux(async () => {
      const one = await reactionDraws(() => expectNoFallback([PREVIEW_CALL], transformation));
      const on = (canvas: HTMLCanvasElement) => one.crux.filter((d) => d.canvas === canvas);
      expect(on(one.result.preview).length > 0, true, 'Crux drew the transformation preview');
      const full = on(one.result.full);
      expect(full.length > 0, true, 'Crux drew the full-screen preview');
      // fitted to its window, as the Reaction panel's full screen (Crux's own bond is 22 px)
      expect(full[0].options.bondLength > 22, true, `the full-screen bond length: ${full[0].options.bondLength}`);
      expect(one.rdkit, 0, 'RDKit\'s reaction drawings, Transformation');
      const two = await reactionDraws(() => expectNoFallback([PREVIEW_CALL], twoComponent));
      expect(two.crux.some((d) => d.canvas === two.result), true, 'Crux drew the two-component preview');
      expect(two.rdkit, 0, 'RDKit\'s reaction drawings, Two-Component Reaction');
    });
    await inEngine(ChemicalEngine.RDKit, async () => {
      const one = await reactionDraws(transformation);
      expect(one.crux.length, 0, 'reactions Crux drew, Transformation');
      const two = await reactionDraws(twoComponent);
      expect(two.crux.length, 0, 'reactions Crux drew, Two-Component Reaction');
    });
  }));

  test('a template Crux reads no reaction from is previewed by RDKit in Crux mode, as the fallback "reaction preview"',
    async () => {
      const template = '[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]';
      const r = window.devicePixelRatio;
      await inCrux(async () => {
        const canvas = ui.canvas(400 * r, 100 * r);
        setCruxUnreadablePreviewsForTests([template]);
        try {
          // one preview RDKit drew for Crux (the fallback's reason names the template; the platform's log gets a real
          // one as "Crux fallback: reaction preview", not this test's own)
          const {result: draws, reason} = await expectFallback(PREVIEW_CALL,
            () => reactionDraws(async () => drawReactionPreview(canvas, template)));
          expect(draws.result, true, 'drawn');
          expect(draws.crux.length, 0, 'reactions Crux drew');
          expect(reason, `Crux reads no reaction from ${template}`, 'why');
        } finally {
          setCruxUnreadablePreviewsForTests(null);
        }
        expect(ink(canvas) > 100, true, 'RDKit\'s picture');
      });
    });
});
