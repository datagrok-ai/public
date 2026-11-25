/**
 * Category `Reactions: export`: the Export panel (widgets/reaction-export.ts), in place of two one-click buttons with
 * the colours baked in and no choice of width, dpi or style. Here:
 * - the panel offers the presets, the widths, the formats and Copy as picture, and no PDF;
 * - ACS 1996 has the colours off;
 * - a single column at 300 dpi is a PNG 980 px wide, at print size, with no colour and no legend;
 * - with the colours on, the picture carries the legend;
 * - the SVG is sized in points;
 * - the Reaction panel and a reaction cell's menu open it, the menu in both modes.
 */
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {ChemicalEngine} from '../engine/chem-engine';
import {CRUX_OFF, defaultExample, withResolution} from '../utils/reactions/reaction-document';
import {ACS_1996, AS_DRAWN, CDXML, CDXML_FILE, CDXML_NOT_YET, CDXML_WITH_MAPPING, CDXML_WITHOUT_MAPPING, CDXML_WRITTEN,
  CHOICES_CLASS,
  COPY_AS_PICTURE, DOUBLE_COLUMN, EDIT_ITEM, EXPORT_ITEM, EXPORT_TITLE, ExportChoice, ExportInputs, ExportOptions,
  DRAWN_OWN_LOOK, FIT_STEPS, PNG_FILE, SVG_FILE,
  FREE_WIDTH, MAP_NUMBERS, NOTE_CLASS, NOTE_WARNING_CLASS, PNG_300, PNG_600, PREVIEW_CLASS, ReactionExport, SCREEN,
  SINGLE_COLUMN, STEP_NUMBER_CLASS, SVG, TEXTS_CDXML_ONLY, addReactionExportItem, cdxmlMapped, cdxmlMapping, exportNote,
  MAP_NUMBERS_SMALL, cruxDrawsSmallMapNumbers, exportCellValue, exportInputsOf, exportPanel, exportWarning,
  openReactionExport} from '../widgets/reaction-export';
import {DrawnPicture, DrawnPictureOptions, drawnWithCoordinates, isCdxmlText, reactionCdxml, sketcherPicture}
  from '../rendering/cdxml-reaction-cells';
import {DrawnStepTexts, cruxDrawsStepTexts, stepTextsOf} from '../rendering/crux-reaction-drawing';
import {EXPORT_ICON_CLASS, reactionPanel} from '../widgets/reaction-panel';
import {CENTRE_COLOURS, LEGEND_TITLE, LEGEND_WORDS} from '../widgets/reaction-centre-legend';
import {inEngine} from './engine-utils';
import {loadFileAsText} from './utils';
import {CDXML_DEMO_FILE, EDITED_ROW} from './reaction-cdxml-tests';

/** Two steps of small molecules (files/test-reactions.csv, row 5): 272 pt wide at ACS 1996's size on one row, so two
 * rows in 8.3 cm, at print size. */
const TWO_STEPS = 'CC(=O)O.Cl>>CC(=O)Cl>>CC(=O)Cl.N>>CC(=O)N';
/** Ten steps of small molecules, no agents: 857 pt wide on one row at ACS 1996's size (chem-wasm 3330cc1e6 in node),
 * so at print size five rows in a single column (160 pt tall), two in a double one (78 pt), one when free (41 pt). */
const TEN_STEPS = 'C=C>>CCBr>>CCO>>CC=O>>CC(=O)O>>CC(=O)Cl>>CC(=O)N>>CC#N>>CCN>>CCNC>>CCN(C)C';
/** An amide coupling, mapped: one bond made (C-N), one broken (C-O). */
const MAPPED = '[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH2:6][CH3:7]>>[CH3:1][C:2](=[O:3])[NH:5][CH2:6][CH3:7]';

const choice = (preset: string, width: string, centre: boolean, transparent = false): ExportChoice =>
  ({preset, width, centre, transparent});
/** {@link choice} with each step fitted on a line of its own. */
const fitted = (width: string): ExportChoice => ({...choice(ACS_1996, width, false), fitSteps: true});

/** Whether each band of `png` (top and bottom in units, the picture `unitsHigh` tall) has ink: dark pixels. */
async function bandsInked(png: Blob, unitsHigh: number, bands: [number, number][]): Promise<boolean[]> {
  const bitmap = await createImageBitmap(png);
  const canvas = new OffscreenCanvas(bitmap.width, bitmap.height);
  const g = canvas.getContext('2d')!;
  g.drawImage(bitmap, 0, 0);
  const d = g.getImageData(0, 0, canvas.width, canvas.height).data;
  const k = canvas.height / unitsHigh;
  return bands.map(([top, bottom]) => {
    for (let y = Math.floor(top * k); y < Math.min(canvas.height, Math.ceil(bottom * k)); y++) {
      for (let x = 0; x < canvas.width; x++) {
        const i = (y * canvas.width + x) * 4;
        if (d[i + 3] > 128 && d[i] < 100 && d[i + 1] < 100 && d[i + 2] < 100)
          return true;
      }
    }
    return false;
  });
}

/** The text Save writes from the Export panel of `reaction` (opened with `options`), once `set` made its choices (an
 * SVG; a CDXML page once the note says it is written), and the note it was saved under. */
async function savedSvg(reaction: string, options: ExportOptions,
  set: (inputs: ExportInputs) => void): Promise<string> {
  return (await savedText(reaction, options, set)).text;
}

async function savedText(reaction: string, options: ExportOptions,
  set: (inputs: ExportInputs) => void): Promise<{text: string, note: string}> {
  const download = DG.Utils.download;
  const got: [string, BlobPart][] = [];
  DG.Utils.download = (name: string, content: BlobPart) => {
    got.push([name, content]);
  };
  const {dialog, inputs} = exportPanel(reaction, options);
  const note = () => dialog.root.querySelector(`.${NOTE_CLASS}`)!.textContent ?? '';
  try {
    set(inputs);
    if (inputs.format.value === CDXML)
      await awaitCheck(() => note().includes(CDXML_WRITTEN), `the CDXML page written: ${note()}`, 15000);
    dialog.getButton('Save').click();
    await awaitCheck(() => got.length === 1, 'one download', 5000);
    return {text: got[0][1] as string, note: note()};
  } finally {
    DG.Utils.download = download;
    dialog.close();
  }
}

/** The Export panel's "Map numbers" box. Its default is the preset's, as the centre colours' are: off for ACS 1996; for
 * Screen as the caller shows them, else off; a preset chosen sets it again. It applies to what Save writes (the map
 * numbers drawn as Crux's `map-number` paths), CDXML too (a mapped page's ReactionStepAtomMap kept with the box on,
 * left out with it off), and the note says which. */
async function mapNumbersBox(): Promise<void> {
  const drawn = (svg: string) => svg.includes('class=\'map-number');
  expect(drawn(await savedSvg(MAPPED, {mapNumbers: 'show'}, () => {})), false,
    'ACS 1996 of a mapped reaction: no map number, whatever the caller shows');
  expect(drawn(await savedSvg(MAPPED, {}, (inputs) => {
    inputs.mapNumbers.value = true;
  })), true, 'the box on: the map numbers drawn');
  const {dialog, inputs} = exportPanel(MAPPED, {mapNumbers: 'show'});
  try {
    expect((dialog.root.textContent ?? '').includes(MAP_NUMBERS), true, 'offered');
    expect(inputs.mapNumbers.value, false, 'ACS 1996: off');
    inputs.style.value = SCREEN;
    await awaitCheck(() => inputs.mapNumbers.value === true, 'Screen: on, as the caller shows them', 2000);
    inputs.style.value = ACS_1996;
    await awaitCheck(() => inputs.mapNumbers.value === false, 'back to ACS 1996: off again', 2000);
  } finally {
    dialog.close();
  }
  // CDXML: the box governs the page's atom mapping, and the note says which
  const mapped = (page: string) => page.includes('ReactionStepAtomMap');
  const off = await savedText(MAPPED, {}, (inputs) => {
    inputs.format.value = CDXML;
  });
  expect(mapped(off.text), false, `ACS 1996, the box off: no atom mapping in the page: ${off.text.slice(0, 300)}`);
  expect(off.note.includes(CDXML_WITHOUT_MAPPING), true, `the note says so: ${off.note}`);
  expect(off.text.includes('ReactionStepReactants'), true, 'still a reaction step');
  const on = await savedText(MAPPED, {}, (inputs) => {
    inputs.mapNumbers.value = true;
    inputs.format.value = CDXML;
  });
  expect(mapped(on.text), true, `the box on: the page keeps its atom mapping: ${on.text.slice(0, 300)}`);
  expect(on.note.includes(CDXML_WITH_MAPPING), true, `the note says so: ${on.note}`);
  expect(cdxmlMapped(on.text) && !cdxmlMapped(cdxmlMapping(on.text, false)) && cdxmlMapping(on.text, true) === on.text,
    true, 'cdxmlMapping leaves the mapping out, or keeps the page as written');
  const plainPage = await savedText('CC(=O)Cl.NCc1ccccc1>>CC(=O)NCc1ccccc1', {}, (inputs) => {
    inputs.mapNumbers.value = true;
    inputs.format.value = CDXML;
  });
  expect(plainPage.note.includes(CDXML_WITH_MAPPING) || plainPage.note.includes(CDXML_WITHOUT_MAPPING), false,
    `an unmapped reaction: no word on the mapping: ${plainPage.note}`);
  const plain = exportPanel(MAPPED);
  try {
    plain.inputs.style.value = SCREEN;
    expect(plain.inputs.mapNumbers.value, false, 'Screen, a caller that does not say: off');
  } finally {
    plain.dialog.close();
  }
}

/** The preview's canvas on the screen, and whether it fits its box (padding 12 px) and fills it one way. */
function previewFit(preview: HTMLElement): {width: number, height: number, px: number, fills: boolean} {
  const canvas = preview.querySelector('canvas') as HTMLCanvasElement | null;
  const r = canvas?.getBoundingClientRect();
  if (!canvas || !r)
    return {width: 0, height: 0, px: 0, fills: false};
  const [w, h] = [preview.clientWidth - 24, preview.clientHeight - 24];
  const fills = r.width <= w + 1 && r.height <= h + 1 && (Math.abs(r.width - w) < 2 || Math.abs(r.height - h) < 2);
  return {width: r.width, height: r.height, px: canvas.width, fills};
}

/** The Export dialog made larger, its preview stayed 720 x 420. The box fills the dialog and
 * the picture is drawn anew to fit it as the dialog is resized, larger and smaller; Save still writes print size. */
async function previewFollowsDialog(): Promise<void> {
  const download = DG.Utils.download;
  const got: [string, BlobPart][] = [];
  DG.Utils.download = (name: string, content: BlobPart) => {
    got.push([name, content]);
  };
  const {dialog, inputs} = exportPanel(TWO_STEPS);
  try {
    const preview = dialog.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
    await awaitCheck(() => previewFit(preview).fills, `fitted to the box as it opens: ${JSON.stringify(
      previewFit(preview))} in ${preview.clientWidth} x ${preview.clientHeight}`, 3000);
    const opened = previewFit(preview);
    const box = preview.getBoundingClientRect();
    const r = dialog.root.getBoundingClientRect();
    // larger, as the dialog's resize leaves it
    dialog.root.style.width = `${Math.round(r.width + 300)}px`;
    dialog.root.style.height = `${Math.round(r.height + 150)}px`;
    await awaitCheck(() => previewFit(preview).width > opened.width + 100 && previewFit(preview).fills,
      `drawn larger to fit the larger box: ${opened.width} -> ${JSON.stringify(previewFit(preview))}`, 3000);
    const larger = previewFit(preview);
    expect(preview.getBoundingClientRect().width > box.width + 250, true,
      `the box fills the dialog: ${box.width} -> ${preview.getBoundingClientRect().width}`);
    expect(Math.abs(larger.px / larger.width - window.devicePixelRatio) < 0.05, true,
      `drawn anew at the screen's pixels, not stretched: ${larger.px} px for ${larger.width} CSS px`);
    // smaller
    dialog.root.style.width = `${Math.round(r.width - 250)}px`;
    await awaitCheck(() => previewFit(preview).width < opened.width - 50 && previewFit(preview).fills,
      `drawn smaller to fit the smaller box: ${JSON.stringify(previewFit(preview))}`, 3000);
    // what Save writes is the print size, whatever the preview's zoom
    inputs.format.value = PNG_300;
    dialog.getButton('Save').click();
    await awaitCheck(() => got.length === 1, 'one download', 5000);
    const facts = await pngFacts(got[0][1] as Blob);
    const l = await withExport(TWO_STEPS, async (e) => e.layout(choice(ACS_1996, SINGLE_COLUMN, false))!);
    expect(`${facts.width} ${facts.height} ${facts.dpi}`, `980 ${Math.round(l.height * 980 / 235)} 300`,
      'Save: 8.3 cm at 300 dpi, as before');
  } finally {
    DG.Utils.download = download;
    dialog.close();
  }
}

/** The preview followed the dialog only as it shrank. The default route in a single column is a tall
 * picture, held by the box's height: a dialog made only wider widened the grey around it. Made wider, the choices take
 * two columns, and the rows they free go to the preview: the picture grows. */
async function previewFollowsWidening(): Promise<void> {
  const {dialog} = exportPanel(defaultExample().text);
  try {
    const preview = dialog.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
    const choices = dialog.root.querySelector(`.${CHOICES_CLASS}`) as HTMLElement;
    const columns = () => getComputedStyle(choices).gridTemplateColumns.trim().split(/\s+/).length;
    await awaitCheck(() => previewFit(preview).fills, `fitted to the box as it opens: ${JSON.stringify(
      previewFit(preview))} in ${preview.clientWidth} x ${preview.clientHeight}`, 3000);
    const opened = previewFit(preview);
    expect(columns(), 1, 'the choices in one column as it opens');
    expect(opened.height >= preview.clientHeight - 24 - 2, true,
      `held by the box's height as it opens: ${JSON.stringify(opened)} in ${preview.clientHeight} px`);
    const r = dialog.root.getBoundingClientRect();
    dialog.root.style.width = `${Math.round(r.width + 400)}px`; // only wider, as a resize leaves it
    dialog.root.style.height = `${Math.round(r.height)}px`;
    await awaitCheck(() => columns() === 2, `the choices in two columns: ${columns()}`, 3000);
    await awaitCheck(() => previewFit(preview).height > opened.height + 40 && previewFit(preview).fills,
      `the picture grown with the dialog's width: ${opened.height} -> ${JSON.stringify(previewFit(preview))}`, 3000);
    expect(dialog.root.getBoundingClientRect().bottom <= window.innerHeight + 1, true, 'the dialog in the window');
  } finally {
    dialog.close();
  }
}

/** Neither column fits a step of the demo route on one line (a column breaks the route at its arrows); so "Fit each
 * step on one line", each step on its row, its number kept (left of the row), every step at one scale, the largest at
 * which the widest fits its line (one bond length across the figure: rows at scales of their own looked wrong), under
 * half size the note saying the text's size. */
async function stepsFitted(): Promise<void> {
  /** The bond length a step's SVG draws (its commonest full-length bond path: a C-C single bond runs centre to
   * centre), in its own units. */
  const drawnBond = (svg: string): number => {
    const counts = new Map<string, number>();
    for (const m of svg.matchAll(/<path class='bond-\d+[^']*' d='M ([\d.]+),([\d.]+) L ([\d.]+),([\d.]+)'/g)) {
      const l = Math.hypot(Number(m[3]) - Number(m[1]), Number(m[4]) - Number(m[2])).toFixed(1);
      counts.set(l, (counts.get(l) ?? 0) + 1);
    }
    return Number([...counts.entries()].sort((a, b) => b[1] - a[1] || Number(b[0]) - Number(a[0]))[0]?.[0] ?? 0);
  };
  await withExport(defaultExample().text, async (e) => {
    const single = e.layout(fitted(SINGLE_COLUMN))!;
    expect(single !== null && single.rows.length, 5, `five steps, five rows: ${e.error}`);
    expect(single.rows.map((r) => r.label).join(), '1,2,3,4,5', 'each row numbered as its arrow is in the route');
    expect(single.width, 235, 'the single column\'s width');
    const percent = single.rows[0].percent;
    expect(single.rows.every((r) => r.percent === percent && r.bond === single.rows[0].bond) && percent < 100, true,
      `every step of the demo route at one scale, made smaller to fit 8.3 cm: ${single.rows.map((r) => r.percent)}`);
    for (const [i, r] of single.rows.entries()) {
      expect(r.x + r.width <= single.width, true, `step ${r.label} within the line`);
      expect(i === 0 || r.y >= single.rows[i - 1].y + single.rows[i - 1].height, true, `step ${r.label} below`);
    }
    const widest = Math.max(...single.rows.map((r) => r.x + r.width));
    expect(widest >= single.width - 0.03 * single.width, true,
      `the largest scale that fits: the widest step reaches the line (${widest} of ${single.width})`);
    const note = exportNote(single);
    expect(percent < 50 && note.includes(`each step on one line, drawn at ${percent}% of print size`), true,
      `the note says the scale: ${note}`);
    // "its text 3.7 pt", grey in the note, read as no warning: a warning of its own says it
    const warned = exportWarning(single) ?? '';
    expect(warned.startsWith('Too small to read in print') && warned.includes(`its text is ` +
      `${(10 * percent / 100).toFixed(1)} pt`) && warned.includes('double column'), true,
    `the warning says the text's size and a wider width: ${warned}`);
    // the drawing: one bond length on every row (Crux's bondLength; a smaller box shrank each step's bonds unevenly)
    const bonds = single.rows.map((r) => drawnBond(e.svg(fitted(SINGLE_COLUMN))!.split(`<svg x='${r.x}' y='${r.y}'`)[1]
      ?.split('</svg>')[0] ?? ''));
    expect(bonds.every((b) => Math.abs(b - single.rows[0].bond) <= 0.15), true,
      `each row's bonds ${single.rows[0].bond.toFixed(2)} long: ${bonds.join(' ')}`);
    // the double column: one scale too, under 100 (steps 4 and 5 are wider than 17.8 cm), over half size
    const double = e.layout(fitted(DOUBLE_COLUMN))!;
    const percents = double.rows.map((r) => r.percent);
    expect(new Set(percents).size === 1 && percents[0] < 100 && percents[0] >= 50, true,
      `one scale, under print size: ${percents.join(' ')}`);
    expect(exportWarning(double), null, `over half size: no warning (${exportNote(double)})`);
    // what Save writes is what is laid out: the SVG's rows and numbers, the PNG at the column's width, each row inked
    const svg = e.svg(fitted(SINGLE_COLUMN))!;
    const numbers = Array.from(new DOMParser().parseFromString(svg, 'image/svg+xml')
      .querySelectorAll(`.${STEP_NUMBER_CLASS}`)).map((t) => t.textContent);
    expect(numbers.join(), '1,2,3,4,5', 'the SVG numbers its rows');
    expect(svg.includes(`width='235pt' height='${single.height}pt'`), true, svg.slice(0, 300));
    const png = (await e.png(fitted(SINGLE_COLUMN), 300))!;
    const facts = await pngFacts(png.png);
    expect(`${facts.width} ${facts.height}`, `980 ${Math.round(single.height * 980 / 235)}`, 'the PNG as laid out');
    const inked = await bandsInked(png.png, single.height, single.rows.map((r) => [r.y, r.y + r.height]));
    expect(inked.join(), 'true,true,true,true,true', 'each step drawn on its row');
    expect(e.layout(choice(ACS_1996, SINGLE_COLUMN, false))!.rows.length, 0, 'off: the column as before');
  });
  await withExport(TWO_STEPS, async (e) => {
    const l = e.layout(fitted(SINGLE_COLUMN))!;
    expect(l.rows.map((r) => r.percent).join(' '), '100 100', 'two small steps fit a single column at print size');
    expect(exportNote(l).endsWith('each step on one line, at print size'), true, exportNote(l));
  });
}

/** A stand-in drawing as drawn ({@link DrawnPicture}): a line 160 px long in a picture 200 x 100 px, white behind it,
 * whatever the options; what each picture was asked with, and each PNG's scale. */
function drawnStandIn(takesOptions = true): DrawnPicture & {asked: DrawnPictureOptions[], scales: number[]} {
  const asked: DrawnPictureOptions[] = [];
  const scales: number[] = [];
  const svg = `<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='200px' height='100px' ` +
    `viewBox='0 0 200 100'>\n<g data-layer='background'>\n<rect width='200' height='100' fill='#FFFFFF'/>\n</g>\n` +
    `<g data-layer='bonds'><path class='drawn-line' d='M 20,50 L 180,50' style='stroke:#000000;stroke-width:6px'/>` +
    `</g>\n</svg>`;
  return {takesOptions, asked, scales,
    svg: async (options) => {
      asked.push({...options});
      return svg;
    },
    png: async (options) => {
      const {scale, ...picture} = options;
      asked.push(picture);
      scales.push(scale);
      const canvas = new OffscreenCanvas(Math.round(200 * scale), Math.round(100 * scale));
      const g = canvas.getContext('2d')!;
      g.fillStyle = '#FFFFFF';
      g.fillRect(0, 0, canvas.width, canvas.height);
      g.fillStyle = '#000000';
      g.fillRect(20 * scale, 47 * scale, 160 * scale, 6 * scale);
      return withResolution(await canvas.convertToBlob({type: 'image/png'}), 96 * scale);
    }};
}

/** What Save writes from an Export panel on TWO_STEPS with `drawing` as drawn, {@link AS_DRAWN} on and `set`'s
 * choices made, once the preview shows the drawing; and the note. */
async function savedDrawn(drawing: DrawnPicture, set: (inputs: ExportInputs) => void = () => {}):
  Promise<{name: string, content: BlobPart, note: string}> {
  const download = DG.Utils.download;
  const got: [string, BlobPart][] = [];
  DG.Utils.download = (name: string, content: BlobPart) => {
    got.push([name, content]);
  };
  const {dialog, inputs} = exportPanel(TWO_STEPS, {asDrawn: drawing});
  const preview = dialog.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
  const note = () => dialog.root.querySelector(`.${NOTE_CLASS}`)!.textContent ?? '';
  try {
    inputs.asDrawn!.value = true;
    set(inputs);
    await awaitCheck(() => preview.dataset.drawn === 'true' && note().includes(', as drawn, ') &&
      !!preview.querySelector('canvas'), `previewed as drawn: ${note()}`, 5000);
    dialog.getButton('Save').click();
    await awaitCheck(() => got.length === 1, 'one download', 5000);
    return {name: got[0][0], content: got[0][1], note: note()};
  } finally {
    DG.Utils.download = download;
    dialog.close();
  }
}

/** "As drawn" beside the Export panel's arrangements exports the reaction as it is drawn on the canvas, not arranged,
 * with the sketcher's picture API (crux-sketch's svgWith and png with options): the
 * style's print size (ACS 1996's 14.4 pt bond), its colours, the Map numbers box. Offered only with a drawing; on, the
 * width, the arrangements and the centre's colours do not apply; the preview, Save (SVG, PNG) and Copy are the
 * drawing's. */
async function asDrawnExported(): Promise<void> {
  const plain = exportPanel(TWO_STEPS);
  try {
    expect(plain.inputs.asDrawn, null, 'no drawing: not offered');
    expect((plain.dialog.root.textContent ?? '').includes(AS_DRAWN), false, 'not in the panel');
  } finally {
    plain.dialog.close();
  }
  // offered, off; on: the drawing previewed, the arrangements off, its note
  const drawing = drawnStandIn();
  const {dialog, inputs} = exportPanel(TWO_STEPS, {asDrawn: drawing});
  const preview = dialog.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
  const note = () => dialog.root.querySelector(`.${NOTE_CLASS}`)!.textContent ?? '';
  const clipboard = navigator.clipboard;
  const copied: Blob[] = [];
  Object.defineProperty(clipboard, 'write', {configurable: true, value: async (items: ClipboardItem[]) => {
    copied.push(await items[0].getType('image/png'));
  }});
  try {
    expect(!!inputs.asDrawn && (dialog.root.textContent ?? '').includes(AS_DRAWN), true, 'offered');
    expect(inputs.asDrawn!.value, false, 'off by default');
    expect(preview.dataset.drawn, 'false', 'the reaction arranged');
    inputs.asDrawn!.value = true;
    await awaitCheck(() => preview.dataset.drawn === 'true' && preview.dataset.width === '200' &&
      !!preview.querySelector('canvas'), `previewed as drawn: ${note()}`, 5000);
    expect(JSON.stringify(drawing.asked[drawing.asked.length - 1]),
      JSON.stringify({bondLength: 14.4, colour: false, mapNumbers: false}),
      'ACS 1996: its 14.4 pt bond, black, the map numbers off as its box');
    expect(note().startsWith('5.3 × 2.6 cm, as drawn, at print size'), true, note());
    expect([inputs.width, inputs.fitSteps, inputs.colours].map((i) => i.enabled).join(), 'false,false,false',
      'the width, the arrangements and the centre\'s colours do not apply');
    // the preview draws the picture: ink along the middle of its canvas
    const canvas = preview.querySelector('canvas') as HTMLCanvasElement;
    const row = canvas.getContext('2d')!.getImageData(0, Math.round(canvas.height / 2), canvas.width, 1).data;
    let ink = 0;
    for (let i = 0; i < row.length; i += 4)
      ink += row[i] < 100 ? 1 : 0;
    expect(ink > canvas.width * 0.6, true, `the drawing in the preview: ${ink} of ${canvas.width} px inked`);
    // Copy: the drawing's PNG, 300 dpi for SVG
    dialog.getButton(COPY_AS_PICTURE).click();
    await awaitCheck(() => copied.length === 1, 'copied', 5000);
    expect(drawing.scales.join(), '3.125', 'Copy: the drawing\'s PNG at 300 dpi');
    // off again: the arrangements back
    inputs.asDrawn!.value = false;
    await awaitCheck(() => preview.dataset.drawn === 'false', 'arranged again', 3000);
    expect([inputs.width, inputs.fitSteps].every((i) => i.enabled), true, 'the arrangements back');
    expect(inputs.colours.enabled, false, 'the centre\'s colours still off: the route has no map numbers');
  } finally {
    delete (clipboard as unknown as {write?: unknown}).write;
    dialog.close();
  }
  // Save: the SVG sized in points at print size, its view box kept; transparent: no background
  const svg = await savedDrawn(drawnStandIn());
  expect(svg.name, SVG_FILE);
  const text = String(svg.content);
  expect(text.includes(`width='150pt'`) && text.includes(`height='75pt'`) && text.includes(`viewBox='0 0 200 100'`),
    true, `sized in points: ${text.slice(0, 240)}`);
  expect(text.includes('drawn-line') && text.includes('<rect'), true, 'the drawing on its white');
  const clear = String((await savedDrawn(drawnStandIn(), (i) => i.transparent.value = true)).content);
  expect(clear.includes('drawn-line') && !clear.includes('<rect'), true, 'transparent: the background left out');
  // Screen: its 22 px bonds (16.5 pt), the element colours; the box's map numbers
  const screen = drawnStandIn();
  await savedDrawn(screen, (i) => {
    i.style.value = SCREEN;
    i.mapNumbers.value = true;
  });
  expect(JSON.stringify(screen.asked[screen.asked.length - 1]),
    JSON.stringify({bondLength: 16.5, colour: true, mapNumbers: true}), 'Screen, Map numbers on');
  // PNG at 300 and 600 dpi: the drawing's PNG at 3.125 and 6.25 times its picture, saying so
  for (const [format, scale, dpi] of [[PNG_300, 3.125, 300], [PNG_600, 6.25, 600]] as const) {
    const png = drawnStandIn();
    const saved = await savedDrawn(png, (i) => i.format.value = format);
    expect(`${saved.name} ${png.scales.join()}`, `${PNG_FILE} ${scale}`, format);
    const facts = await pngFacts(saved.content as Blob);
    expect(`${facts.width} ${facts.height} ${facts.dpi}`,
      `${Math.round(200 * scale)} ${Math.round(100 * scale)} ${dpi}`, `${format}: the picture at ${dpi} dpi`);
  }
  // a transparent PNG: the SVG drawn on no background
  const clearPng = drawnStandIn();
  const saved = await savedDrawn(clearPng, (i) => {
    i.format.value = PNG_300;
    i.transparent.value = true;
  });
  expect(clearPng.scales.length, 0, 'not the sketcher\'s PNG, which is white');
  const bitmap = await createImageBitmap(saved.content as Blob);
  const corner = new OffscreenCanvas(bitmap.width, bitmap.height).getContext('2d')!;
  corner.drawImage(bitmap, 0, 0);
  expect(corner.getImageData(1, 1, 1, 1).data[3], 0, 'transparent round the drawing');
  expect((await pngFacts(saved.content as Blob)).dpi, 300, 'at 300 dpi');
}

/** A crux-sketch without its picture API: the canvas's own picture, the options left out, the note saying so; the
 * options passed on to one with it (svgWith, png with options). And which cells are offered "As drawn": those whose
 * value has coordinates of its own (a CDXML page, an RXN block drawn), not a reaction SMILES. */
async function asDrawnFallback(): Promise<void> {
  const calls: string[] = [];
  const before = {svg: '<svg width=\'10\' height=\'5\'></svg>',
    png: async (o?: {scale?: number}) => {
      calls.push(`png ${JSON.stringify(o)}`);
      return null;
    }};
  const old = sketcherPicture(before);
  expect(old.takesOptions, false, 'no svgWith: takes no options');
  expect(await old.svg({bondLength: 14.4, mapNumbers: false}), before.svg, 'its svg');
  await old.png({bondLength: 14.4, scale: 2});
  expect(calls.pop(), 'png {"scale":2}', 'png with its scale alone');
  const now = sketcherPicture({...before, svgWith: (o?: DrawnPictureOptions) => {
    calls.push(`svgWith ${JSON.stringify(o)}`);
    return '<svg/>';
  }});
  expect(now.takesOptions, true, 'svgWith: takes the options');
  await now.svg({bondLength: 14.4, colour: false});
  expect(calls.pop(), 'svgWith {"bondLength":14.4,"colour":false}');
  await now.png({bondLength: 14.4, scale: 2});
  expect(calls.pop(), 'png {"bondLength":14.4,"scale":2}');
  // the note of a picture at the canvas's own size
  const saved = await savedDrawn(drawnStandIn(false));
  expect(saved.note.includes(DRAWN_OWN_LOOK), true, saved.note);
  // the cells offered it
  const rxn = (x: string) => `$RXN\n\n  Chem\n\n  1  1\n$MOL\n\n  Chem\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n` +
    `    ${x}    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$MOL\n\n  Chem\n\n` +
    `  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n` +
    `M  END\n`;
  expect(drawnWithCoordinates('CCO>>CC=O'), false, 'a reaction SMILES: laid out, not drawn');
  expect(drawnWithCoordinates(rxn('0.0000')), false, 'an RXN block with no coordinates');
  expect(drawnWithCoordinates(rxn('1.2990')), true, 'an RXN block drawn');
  expect(drawnWithCoordinates('<?xml version="1.0"?><CDXML><page/></CDXML>'), true, 'a CDXML page');
}

/** A CDXML cell's Export ("Export reaction..."): "As drawn" offered, the page pictured by the cells' reader, Save
 * writing that picture; a reaction SMILES cell's: not offered (it has no drawing). */
async function cellAsDrawn(): Promise<void> {
  const page = (DG.DataFrame.fromCsv(await loadFileAsText(CDXML_DEMO_FILE)).col('cdxml')!
    .toList() as string[])[EDITED_ROW];
  const panel = () => DG.Dialog.getOpenDialogs().find((d) => d.title === EXPORT_TITLE);
  const download = DG.Utils.download;
  const got: [string, BlobPart][] = [];
  DG.Utils.download = (name: string, content: BlobPart) => {
    got.push([name, content]);
  };
  try {
    await exportCellValue('CCO>>CC=O');
    await awaitCheck(() => !!panel(), 'the SMILES cell\'s panel', 15000);
    expect(exportInputsOf(panel()!.root)?.asDrawn ?? null, null, 'a SMILES cell: not offered');
    panel()!.close();
    await exportCellValue(page);
    await awaitCheck(() => !!panel(), 'the CDXML cell\'s panel', 15000);
    const inputs = exportInputsOf(panel()!.root)!;
    expect(!!inputs.asDrawn, true, 'a CDXML cell: offered');
    inputs.asDrawn!.value = true;
    const preview = panel()!.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
    const note = () => panel()?.root.querySelector(`.${NOTE_CLASS}`)?.textContent ?? '';
    await awaitCheck(() => note().includes(', as drawn, ') && !!preview.querySelector('canvas'),
      `the page previewed as drawn: ${note()}`, 15000);
    panel()!.getButton('Save').click();
    await awaitCheck(() => got.length === 1, 'one download', 5000);
    const svg = String(got[0][1]);
    expect(svg.includes('data-layer') && /width='[\d.]+pt'/.test(svg), true,
      `the page's own picture, in points: ${svg.slice(0, 300)}`);
  } finally {
    DG.Utils.download = download;
    panel()?.close();
  }
}

/** A PNG's size (its IHDR) and the dpi its pHYs chunk (right after IHDR, as withResolution writes it) says. */
async function pngFacts(png: Blob): Promise<{width: number, height: number, dpi: number}> {
  const b = new DataView(await png.arrayBuffer());
  return {width: b.getUint32(16), height: b.getUint32(20), dpi: Math.round(b.getUint32(41) * 0.0254)};
}

/** Whether a pixel of `png` has colour (the centre's green, red or blue): its channels far apart. */
async function hasColour(png: Blob): Promise<boolean> {
  const bitmap = await createImageBitmap(png);
  const canvas = new OffscreenCanvas(bitmap.width, bitmap.height);
  const g = canvas.getContext('2d')!;
  g.drawImage(bitmap, 0, 0);
  const d = g.getImageData(0, 0, canvas.width, canvas.height).data;
  for (let i = 0; i < d.length; i += 4) {
    if (d[i + 3] > 128 && Math.max(Math.abs(d[i] - d[i + 1]), Math.abs(d[i + 1] - d[i + 2]),
      Math.abs(d[i] - d[i + 2])) > 80)
      return true;
  }
  return false;
}

/** "HATU, DIPEA" typed over the arrow was in no picture. A drawing's step texts (crux-sketch's reactionSteps) become
 * Crux's `stepTexts` by step, and the export draws them: the scheme with all of them,
 * each fitted step with its own alone (red on chem-wasm 3330cc1e6, which draws no host texts). */
async function stepTextsDrawn(): Promise<void> {
  const drawn: DrawnStepTexts[] = [{step: 1, texts: [{place: 'above', text: 'HATU, DIPEA'}]},
    {step: 3, texts: [{place: 'below', text: 'rt, 12 h'}, {place: 'above', text: 'DIPEA', drawing: {kind: 'text'}}]}];
  const texts = stepTextsOf(drawn)!;
  expect(JSON.stringify(texts), JSON.stringify([{above: ['HATU, DIPEA'], below: []}, null,
    {above: [{kind: 'text'}], below: ['rt, 12 h']}]), 'by step, each text\'s drawing object as it is');
  expect(stepTextsOf([]), null, 'none');
  expect(stepTextsOf([{step: 2, texts: []}]), null, 'no step with a text');
  expect(cruxDrawsStepTexts(), true, 'Chem\'s Crux file draws a host\'s step texts');
  const route = defaultExample().text;
  // each fitted row's own SVG, by the layout's rows, its place left out (the SVG's first `<svg x='0' y='0'` is the
  // scheme's clipping frame round the rows, none of them: a run once counted it as a sixth row)
  const inner = (e: ReactionExport) => {
    const svg = e.svg(fitted(SINGLE_COLUMN))!;
    return e.layout(fitted(SINGLE_COLUMN))!.rows.map((r) =>
      svg.split(`<svg x='${r.x}' y='${r.y}'`)[1]?.split('</svg>')[0] ?? '');
  };
  await withExport(route, async (plain) => {
    const withTexts = new ReactionExport(route, texts);
    try {
      const double = choice(ACS_1996, DOUBLE_COLUMN, false);
      expect(withTexts.svg(double) !== plain.svg(double), true, 'the scheme drawn with its step texts');
      const [a, b] = [inner(plain), inner(withTexts)];
      expect(a.length === 5 && b.length === 5 && [...a, ...b].every((s) => s.includes('<path')), true,
        `five rows, each drawn: ${a.length}, ${b.length}`);
      // steps 2, 4 and 5 drawn as without texts, at the steps' one scale (which the texts may change)
      const sameScale = plain.layout(fitted(SINGLE_COLUMN))!.rows[0].bond ===
        withTexts.layout(fitted(SINGLE_COLUMN))!.rows[0].bond;
      expect(a.map((s, k) => s === b[k] ? 'same' : 'texts').join(), sameScale ? 'texts,same,texts,same,same' :
        'texts,texts,texts,texts,texts', 'each fitted step drawn with its own texts alone');
    } finally {
      withTexts.free();
    }
  });
}

/** `run` with an export of `text`, freed after. */
async function withExport<T>(text: string, run: (e: ReactionExport) => Promise<T>): Promise<T> {
  const e = new ReactionExport(text);
  try {
    return await run(e);
  } finally {
    e.free();
  }
}

category('Reactions: export', () => {
  test('the panel offers the presets, widths and formats, no PDF, and the colours off for ACS 1996', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const {dialog, inputs} = exportPanel(defaultExample().text);
      try {
        const text = dialog.root.textContent ?? '';
        for (const offered of [EXPORT_TITLE, ACS_1996, SCREEN, SINGLE_COLUMN, DOUBLE_COLUMN, FREE_WIDTH, SVG, PNG_300,
          PNG_600, CDXML, 'Transparent background', COPY_AS_PICTURE, 'Save'])
          expect(text.includes(offered), true, `the panel offers "${offered}"`);
        expect(/PDF/i.test(text), false, 'no PDF (17A)');
        expect(inputs.colours.value, false, 'ACS 1996: the centre colours off (P9)');
        inputs.style.value = SCREEN;
        await awaitCheck(() => inputs.colours.value === true, 'Screen: the colours on', 2000);
        inputs.style.value = ACS_1996;
        await awaitCheck(() => inputs.colours.value === false, 'back to ACS 1996: off again', 2000);
        const preview = dialog.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
        expect(preview.querySelector('canvas') !== null, true, 'a preview is drawn');
        expect(preview.dataset.width, '235', 'the preview is a single column: 8.3 cm in points');
      } finally {
        dialog.close();
      }
    });
  });

  // Map numbers and Reaction centre colours looked as if they did nothing: a reaction without map numbers has neither
  // to show. Such a reaction's boxes are off (their tooltip says why, NOT_MAPPED_TOOLTIP); a mapped one's apply
  test('a reaction without map numbers: Map numbers and Reaction centre colours off; a mapped one\'s on', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      for (const [reaction, mapped] of [[TWO_STEPS, false], [MAPPED, true]] as const) {
        const {dialog, inputs} = exportPanel(reaction);
        try {
          expect(`${inputs.colours.enabled} ${inputs.mapNumbers.enabled}`, `${mapped} ${mapped}`, reaction);
        } finally {
          dialog.close();
        }
      }
    });
  });

  test('ACS 1996 in a single column at 300 dpi: 980 px wide, at print size, no colour, no legend', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      await withExport(TWO_STEPS, async (e) => {
        const l = e.layout(choice(ACS_1996, SINGLE_COLUMN, false))!;
        expect([l.width, l.column, l.percent].join(' '), '235 true 100',
          'two steps keep ACS 1996 size in 8.3 cm, on the rows they need (the column option of the drawing)');
        expect(exportNote(l), `8.3 × ${(l.height / 72 * 2.54).toFixed(1)} cm, at print size`, 'the note');
      });
      // the default route, mapped: black with the colours off
      await withExport(defaultExample().text, async (e) => {
        const acs = choice(ACS_1996, SINGLE_COLUMN, false);
        const l = e.layout(acs);
        expect(l !== null, true, `drawn: ${e.error}`);
        const png = await e.png(acs, 300);
        const facts = await pngFacts(png!.png);
        expect(`${facts.width} ${facts.dpi}`, '980 300', '8.3 cm at 300 dpi');
        expect(facts.height, Math.round(l!.height * 980 / 235), 'as tall as the preview, at the same scale');
        expect(await hasColour(png!.png), false, 'ACS 1996: black, the centre colours off');
        expect(l!.legend.length, 0, 'no legend without colours');
        expect(e.svg(acs)!.includes(LEGEND_TITLE), false, 'no legend in the SVG either');
      });
    });
  });

  test('the colours on: the mapped reaction is coloured and carries the legend', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      await withExport(MAPPED, async (e) => {
        const on = choice(ACS_1996, SINGLE_COLUMN, true);
        const png = await e.png(on, 300);
        expect(await hasColour(png!.png), true, 'the bond made green, the bond broken red');
        const svg = e.svg(on)!;
        for (const words of [LEGEND_TITLE, LEGEND_WORDS[0], LEGEND_WORDS[1]])
          expect(svg.includes(words), true, `the legend says "${words}"`);
        expect(svg.includes(LEGEND_WORDS[2]), false, 'no bond changes its order');
        expect(svg.includes(`fill='${CENTRE_COLOURS[0]}'`), true, 'the green swatch');
        const off = e.svg(choice(ACS_1996, SINGLE_COLUMN, false))!;
        expect(off.includes(LEGEND_TITLE), false, 'colours off: no legend');
        expect(e.layout(on)!.height > e.layout(choice(ACS_1996, SINGLE_COLUMN, false))!.height, true,
          'the legend under the scheme');
      });
    });
  });

  // the legend was oversized: Screen's 11 px against its 8.7 px atom labels, and at full size under
  // a scheme drawn smaller. Now 0.8 of the look's label size (ACS 1996: 8 pt under 10 pt labels; Screen: 7 px),
  // smaller with the scheme (to 60% of that)
  test('the legend is sized to the figure: under the labels\' size, smaller with a scheme drawn smaller', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      await withExport(MAPPED, async (e) => {
        const acs = e.layout(choice(ACS_1996, FREE_WIDTH, true))!;
        const screen = e.layout(choice(SCREEN, FREE_WIDTH, true))!;
        expect(acs.legend.length > 0 && screen.legend.length > 0, true, 'legends');
        expect(acs.legendFont.size, 8, 'ACS 1996: 8 pt at print size');
        expect(Math.abs(screen.legendFont.size - 0.8 * 8.7) < 0.01, true, `Screen: ${screen.legendFont.size} px`);
      });
      await withExport(defaultExample().text, async (e) => {
        const small = e.layout({...fitted(SINGLE_COLUMN), centre: true})!;
        const percent = Math.min(...small.rows.map((r) => r.percent));
        expect(percent < 50 && small.legendFont.size < 8 && small.legendFont.size >= 0.6 * 8 - 0.01, true,
          `at ${percent}%: the legend ${small.legendFont.size} pt`);
      });
    });
  });

  test('the SVG is sized in points at print size; transparent leaves the white out', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      await withExport(MAPPED, async (e) => {
        const white = e.svg(choice(ACS_1996, SINGLE_COLUMN, false))!;
        const h = e.layout(choice(ACS_1996, SINGLE_COLUMN, false))!.height;
        expect(white.includes(`width='235pt' height='${h}pt' viewBox='0 0 235 ${h}'`), true, white.slice(0, 300));
        expect(white.includes(`<rect x='0' y='0' width='235' height='${h}' fill='#FFFFFF'/>`), true, 'white');
        const clear = e.svg(choice(ACS_1996, SINGLE_COLUMN, false, true))!;
        expect(clear.includes(`fill='#FFFFFF'/>`), false, 'transparent: no white');
        const screen = e.svg(choice(SCREEN, SINGLE_COLUMN, true))!;
        expect(screen.includes(`width='314px'`), true, 'Screen: 8.3 cm in CSS px');
      });
    });
  });

  test('a double column takes fewer rows than a single one; free is one row; PNG 600 dpi', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      // the rows, on a route whose every piece fits a single column: all three at print size, so their heights are
      // their rows (the default route takes five rows in both columns, a step a row, the single drawn at 66%:
      // heights at print size 412 and 414 said nothing of rows)
      await withExport(TEN_STEPS, async (e) => {
        const [single, double, free] = [SINGLE_COLUMN, DOUBLE_COLUMN, FREE_WIDTH]
          .map((w) => e.layout(choice(ACS_1996, w, false))!);
        expect(double.width, 505, '17.8 cm in points (504.6, rounded)');
        expect([single, double, free].map((l) => l.percent).join(' '), '100 100 100', 'all at print size');
        expect(single.height > double.height && double.height > free.height, true,
          `heights ${single.height}, ${double.height}, ${free.height}: fewer rows in the double column, one when free`);
      });
      await withExport(defaultExample().text, async (e) => {
        const single = e.layout(choice(ACS_1996, SINGLE_COLUMN, false))!;
        const double = e.layout(choice(ACS_1996, DOUBLE_COLUMN, false))!;
        const free = e.layout(choice(ACS_1996, FREE_WIDTH, false))!;
        expect(free.percent, 100, 'free: the natural size');
        // its agents are drawn as structures over the arrows, so a piece of it is wider than a single column: there it
        // is drawn smaller, and the note says how much
        expect(single.percent !== null && single.percent < 100, true, `the single column at ${single.percent}%`);
        for (const l of [single, double]) {
          expect(l.percent !== null && l.percent <= 100, true, `${l.width} pt: ${l.percent}%`);
          expect(exportNote(l).includes('% of print size'), l.percent! < 100, exportNote(l));
        }
        const png = await e.png(choice(ACS_1996, DOUBLE_COLUMN, false), 600);
        const facts = await pngFacts(png!.png);
        expect(`${facts.width} ${facts.dpi}`, `${Math.round(17.8 / 2.54 * 600)} 600`, '17.8 cm at 600 dpi');
      });
    });
  });

  test('the Reaction panel and a reaction cell\'s menu open it, the menu in both modes', async () => {
    const cellArgs = (menu: DG.Menu) => ({args: {menu, item: {isTableCell: true,
      tableColumn: {semType: DG.SEMTYPE.CHEMICAL_REACTION}, cell: {valueString: MAPPED}}}});
    await inEngine(ChemicalEngine.Crux, async () => {
      const widget = reactionPanel(MAPPED);
      try {
        expect(widget.root.querySelector(`.${EXPORT_ICON_CLASS}`) !== null, true, 'the panel\'s Export icon');
      } finally {
        widget.detach();
      }
      const menu = DG.Menu.popup();
      expect(addReactionExportItem(cellArgs(menu)), true, `Crux mode: ${EXPORT_ITEM} in the menu`);
      expect(!!menu.find(EXPORT_ITEM), true, `the menu has ${EXPORT_ITEM}`);
      // double-click was the only way into the editor
      expect(!!menu.find(EDIT_ITEM), true, `the menu has ${EDIT_ITEM}`);
      const empty = DG.Menu.popup();
      const emptyArgs = cellArgs(empty);
      emptyArgs.args.item.cell.valueString = '';
      expect(addReactionExportItem(emptyArgs), true, 'an empty reaction cell');
      expect(`${!!empty.find(EDIT_ITEM)} ${!!empty.find(EXPORT_ITEM)}`, 'true false',
        `an empty cell: ${EDIT_ITEM}, nothing to export`);
      const molecule = {args: {menu: DG.Menu.popup(), item: {isTableCell: true,
        tableColumn: {semType: DG.SEMTYPE.MOLECULE}, cell: {valueString: 'CCO'}}}};
      expect(addReactionExportItem(molecule), false, 'not on a molecule cell');
    });
    // the reaction features work in RDKit mode as in Crux mode: the menu's two items in RDKit mode too, and its Export
    // draws (Chem's Crux engine, loaded for it if RDKit mode had not)
    await inEngine(ChemicalEngine.RDKit, async () => {
      const menu = DG.Menu.popup();
      expect(addReactionExportItem(cellArgs(menu)), true, 'RDKit mode: the reaction cell\'s menu');
      expect(`${!!menu.find(EDIT_ITEM)} ${!!menu.find(EXPORT_ITEM)}`, 'true true', `RDKit mode: ${EDIT_ITEM} and ` +
        `${EXPORT_ITEM}`);
      const dialog = await exportCellValue(MAPPED);
      try {
        const note = dialog.root.querySelector(`.${NOTE_CLASS}`)?.textContent ?? '';
        expect(note.includes(CRUX_OFF), false, `RDKit mode: the panel draws (${note})`);
        expect(dialog.getButton('Save').classList.contains('disabled'), false, 'RDKit mode: Save on');
      } finally {
        dialog.close();
      }
    });
  });

  // a CDXML cell's "Export reaction..." read the page as reaction text, which Crux does
  // not read (SVG and PNG off), and its Reaction panel's Export wrote the CDXML anew from the page's SMILES (its texts
  // and arrow styles lost)
  test('a CDXML cell\'s Export, from its menu and its Reaction panel: the page\'s reaction drawn, its CDXML the page',
    async () => {
      await inEngine(ChemicalEngine.Crux, async () => {
        const page = (DG.DataFrame.fromCsv(await loadFileAsText(CDXML_DEMO_FILE)).col('cdxml')!
          .toList() as string[])[EDITED_ROW];
        expect(isCdxmlText(page), true, 'the demo page');
        const download = DG.Utils.download;
        const got: [string, BlobPart][] = [];
        DG.Utils.download = (name: string, content: BlobPart) => {
          got.push([name, content]);
        };
        const panel = () => DG.Dialog.getOpenDialogs().find((d) => d.title === EXPORT_TITLE);
        const note = () => panel()?.root.querySelector(`.${NOTE_CLASS}`)?.textContent ?? '';
        const savedCdxml = async (how: string): Promise<string> => {
          await awaitCheck(() => !!panel(), `${how}: the Export panel opens`, 15000);
          expect(panel()!.root.querySelector(`.${PREVIEW_CLASS} canvas`) !== null, true,
            `${how}: the page's reaction previewed: ${note()}`);
          expect(panel()!.getButton('Save').classList.contains('disabled'), false, `${how}: an SVG to save`);
          panel()!.input('Format').value = CDXML;
          await awaitCheck(() => note().includes(CDXML_WRITTEN), `${how}: the CDXML page: ${note()}`, 15000);
          got.length = 0;
          panel()!.getButton('Save').click();
          await awaitCheck(() => got.length === 1, `${how}: one download`, 5000);
          return String(got[0][1]);
        };
        try {
          const menu = DG.Menu.popup();
          expect(addReactionExportItem({args: {menu, item: {isTableCell: true,
            tableColumn: {semType: DG.SEMTYPE.CHEMICAL_REACTION}, cell: {valueString: page}}}}), true, 'the menu');
          menu.find(EXPORT_ITEM).click();
          expect(await savedCdxml('the menu') === cdxmlMapping(page, false), true, 'the menu: the page itself');
          panel()?.close();
          const widget = reactionPanel(page);
          try {
            await awaitCheck(() => widget.root.querySelector(`.${EXPORT_ICON_CLASS}`) !== null,
              'the Reaction panel did not read the page', 15000);
            (widget.root.querySelector(`.${EXPORT_ICON_CLASS}`) as HTMLElement).click();
            expect(await savedCdxml('the Reaction panel') === cdxmlMapping(page, false), true,
              'the Reaction panel: the page itself');
          } finally {
            widget.detach();
          }
        } finally {
          DG.Utils.download = download;
          panel()?.close();
        }
      });
    });

  test('a text that reads as no reaction is said so, and Save is off; CDXML is written when the sketcher writes it',
    async () => {
      await inEngine(ChemicalEngine.Crux, async () => {
        const dialog = openReactionExport('this is not a reaction');
        try {
          const note = dialog.root.querySelector(`.${NOTE_CLASS}`)!.textContent ?? '';
          expect(note.includes('does not read as a reaction'), true, note);
          expect(dialog.getButton('Save').classList.contains('disabled'), true, 'nothing to save');
          expect(/crux|rdkit/i.test(note), false, 'P10: no engine named');
        } finally {
          dialog.close();
        }
        const mapped = exportPanel(MAPPED);
        try {
          mapped.inputs.format.value = CDXML;
          const note = () => mapped.dialog.root.querySelector(`.${NOTE_CLASS}`)!.textContent ?? '';
          // the page is written by Crux Sketch (its CDXML module loaded on the way): said once it is, or that this
          // version of the sketcher writes no reaction as CDXML (Save off)
          await awaitCheck(() => note().includes(CDXML_WRITTEN) || note().includes(CDXML_NOT_YET), note(), 15000);
          const saveOff = mapped.dialog.getButton('Save').classList.contains('disabled');
          if (note().includes(CDXML_NOT_YET))
            expect(saveOff, true, 'CDXML not written: not saved');
          else {
            expect(saveOff, false, 'CDXML written: Save on');
            expect(isCdxmlText(await reactionCdxml(MAPPED) ?? ''), true, 'the page written is CDXML');
          }
        } finally {
          mapped.dialog.close();
        }
      });
    });

  // the app's Export wrote CDXML from the reaction text, the text over the arrow lost
  test('a drawing\'s own CDXML page is what Save writes; texts the pictures leave out are said to be in CDXML only',
    async () => {
      await inEngine(ChemicalEngine.Crux, async () => {
        const page = '<?xml version="1.0"?><CDXML><page><t><s>DMF, rt, 12 h</s></t></page></CDXML>';
        const download = DG.Utils.download;
        const got: [string, BlobPart][] = [];
        DG.Utils.download = (name: string, content: BlobPart) => {
          got.push([name, content]);
        };
        let asked = 0;
        const own = exportPanel(MAPPED, {textsLeftOut: true, cdxml: async () => {
          asked++;
          return page;
        }});
        const note = () => own.dialog.root.querySelector(`.${NOTE_CLASS}`)!.textContent ?? '';
        try {
          expect(note().includes(TEXTS_CDXML_ONLY), true, note());
          // the line says which formats leave the texts out, not only which keeps them
          expect(note().includes('SVG and PNG leave them out'), true, `the formats named: ${note()}`);
          expect(asked, 0, 'the page is asked for once CDXML is chosen, not before');
          own.inputs.format.value = CDXML;
          await awaitCheck(() => note().includes(CDXML_WRITTEN), note(), 5000);
          expect(note().includes(TEXTS_CDXML_ONLY), true, `still said: ${note()}`);
          own.dialog.getButton('Save').click();
          await awaitCheck(() => got.length === 1, 'one download', 5000);
          expect(got[0][0], CDXML_FILE);
          expect(got[0][1], page, 'the drawing\'s own page, as given');
          expect(asked, 1, 'asked once');
        } finally {
          DG.Utils.download = download;
          own.dialog.close();
        }
        const plain = exportPanel(MAPPED);
        try {
          const said = plain.dialog.root.querySelector(`.${NOTE_CLASS}`)!.textContent ?? '';
          expect(said.includes(TEXTS_CDXML_ONLY), false, `no texts, no line: ${said}`);
        } finally {
          plain.dialog.close();
        }
      });
    });

  test('the preview fills the dialog and is drawn again to fit it as the dialog is resized; Save keeps print size',
    async () => inEngine(ChemicalEngine.Crux, previewFollowsDialog));

  test('the preview grows as the dialog is made only wider: the choices in two columns, their rows the preview\'s ' +
    '(EVAL6 #15)', async () => inEngine(ChemicalEngine.Crux, previewFollowsWidening));

  test('fit each step on one line: the demo route\'s steps on rows 1-5, at one scale, the largest that fits the line',
    async () => inEngine(ChemicalEngine.Crux, stepsFitted));

  // a "Map numbers" box, its default the preset's as the centre colours'
  test('map numbers: off for ACS 1996, as the caller shows them for Screen; the box puts them in what Save writes',
    async () => inEngine(ChemicalEngine.Crux, mapNumbersBox));

  // Map numbers on, the default route in a single column (each step on one line, at 37% of print size)
  // showed no map number in the preview or the SVG, and nothing said why: Crux leaves map numbers out of a scheme
  // drawn under 10 units a bond (crux-core scheme-draw), 70% of ACS 1996's print size. So the box on asks
  // for them at any size (`mapNumbersMinBond: 0`); where the engine draws them so, the note says they are small, else
  // the warning says they are left out
  test('Map numbers on, a figure too small for them: drawn small and said so, or left out and warned of, by the ' +
    'engine; the free width draws them (EVAL6 #3)', async () => inEngine(ChemicalEngine.Crux, async () => {
    const route = defaultExample().text;
    const any = cruxDrawsSmallMapNumbers();
    const {dialog, inputs} = exportPanel(route, {});
    const warned = () => dialog.root.querySelector(`.${NOTE_WARNING_CLASS}`)?.textContent ?? '';
    const noted = () => dialog.root.querySelector(`.${NOTE_CLASS}`)?.textContent ?? '';
    const said = () => any ? noted() : warned();
    try {
      inputs.fitSteps.value = true;
      inputs.mapNumbers.value = true;
      expect(ReactionExport.options({...choice(ACS_1996, SINGLE_COLUMN, false), mapNumbers: 'show'})
        .mapNumbersMinBond, 0, 'the box on asks for map numbers at any size');
      if (any) {
        await awaitCheck(() => noted().includes(MAP_NUMBERS_SMALL), `said: ${noted()}`, 3000);
        expect(warned().includes('Map numbers'), false, `no warning: ${warned()}`);
      } else {
        await awaitCheck(() => warned().includes('Map numbers are left out'), `said: ${warned()}`, 3000);
        expect(warned().includes('drawn from 70% of print size') && warned().includes('free width'), true, warned());
      }
      inputs.width.value = FREE_WIDTH;
      await awaitCheck(() => !said().includes('Map numbers'), `the free width: nothing said (${said()})`, 3000);
      inputs.width.value = SINGLE_COLUMN;
      inputs.mapNumbers.value = false;
      await awaitCheck(() => !said().includes('Map numbers'), `the box off: nothing said (${said()})`, 3000);
    } finally {
      dialog.close();
    }
    const unmapped = exportPanel(TEN_STEPS, {});
    try {
      unmapped.inputs.mapNumbers.value = true;
      const said = unmapped.dialog.root.querySelector(`.${NOTE_WARNING_CLASS}`)?.textContent ?? '';
      expect(said.includes('Map numbers'), false, `a reaction with no map numbers: nothing said (${said})`);
    } finally {
      unmapped.dialog.close();
    }
    const drawn = (svg: string) => svg.includes('class=\'map-number');
    expect(drawn(await savedSvg(route, {}, (i) => {
      i.fitSteps.value = true;
      i.mapNumbers.value = true;
      i.width.value = FREE_WIDTH;
    })), true, 'the free width: the map numbers drawn');
    expect(drawn(await savedSvg(route, {}, (i) => {
      i.fitSteps.value = true;
      i.mapNumbers.value = true;
    })), any, any ? 'a single column: drawn small, as asked' : 'a single column: left out by this engine');
  }));

  // unticked, a single column dropped every reagent and condition: the engine's levels of detail, a cell's (agents
  // left out under 9 px bonds), applied to the page; the export asks for none (crux-core's `levels` option)
  test('a single column keeps the agents over the arrows: the export asks for no levels of detail', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      expect(ReactionExport.options(choice(ACS_1996, SINGLE_COLUMN, false)).levels, false, 'levels off');
      await withExport(defaultExample().text, async (e) => {
        for (const width of [SINGLE_COLUMN, DOUBLE_COLUMN]) {
          const svg = e.svg(choice(ACS_1996, width, false)) ?? '';
          const names = (svg.match(/class='reaction-agent-name'/g) ?? []).length;
          expect(names > 0, true, `${width}: the agents' names over the arrows (${names} paths; ${e.error})`);
        }
      });
    });
  });

  test('a drawing\'s step texts are drawn: the scheme with all, each fitted step with its own', async () =>
    inEngine(ChemicalEngine.Crux, stepTextsDrawn));

  test('the panel offers "Fit each step on one line", off; checked, the preview has a row a step', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const {dialog, inputs} = exportPanel(defaultExample().text);
      try {
        const preview = dialog.root.querySelector(`.${PREVIEW_CLASS}`) as HTMLElement;
        expect((dialog.root.textContent ?? '').includes(FIT_STEPS), true, 'offered');
        expect(inputs.fitSteps.value, false, 'off by default');
        expect(preview.dataset.rows, '0', 'off: the column\'s rows');
        inputs.fitSteps.value = true;
        await awaitCheck(() => preview.dataset.rows === '5', `a row a step: ${preview.dataset.rows}`, 3000);
        const note = dialog.root.querySelector(`.${NOTE_CLASS}`)!.textContent ?? '';
        expect(note.includes('each step on one line'), true, note);
        // under half size, a warning seen as one: shown, its sign, the warnings' orange, its words
        const warning = dialog.root.querySelector(`.${NOTE_WARNING_CLASS}`) as HTMLElement | null;
        expect(!!warning && warning.style.display !== 'none' && warning.getBoundingClientRect().height > 0, true,
          'the warning shown');
        expect((warning!.textContent ?? '').includes('Too small to read in print'), true, warning!.textContent ?? '');
        expect(warning!.querySelector('.fa-exclamation-triangle') !== null, true, 'its sign');
        const colourOf = (token: string) => {
          const probe = document.createElement('div');
          probe.style.color = `var(${token})`;
          document.body.append(probe);
          const c = getComputedStyle(probe).color;
          probe.remove();
          return c;
        };
        const sign = warning!.querySelector('.fa-exclamation-triangle') as HTMLElement;
        expect(getComputedStyle(sign).color, colourOf('--orange-2'), 'its sign in the warnings\' orange');
        const [shown, expected] = [getComputedStyle(warning!).color, colourOf('--orange-3')];
        expect(shown, expected, 'its words in the dark orange');
        inputs.fitSteps.value = false;
        inputs.width.value = DOUBLE_COLUMN;
        await awaitCheck(() => warning!.style.display === 'none', 'no warning at the double column\'s scale', 3000);
      } finally {
        dialog.close();
      }
    });
  });

  // in "Fit each step on one line" the same intermediate turned over between its two steps
  test('fit each step on one line draws the steps from the route\'s own layout, each fitted to its row (needs the ' +
    'vendored chem-wasm, ENGINE-7)', async () => inEngine(ChemicalEngine.Crux, async () => {
    await withExport(defaultExample().text, async (e) => {
      const l = e.layout(fitted(SINGLE_COLUMN))!;
      expect(l.rows.length, 5, 'a row a step');
      expect(new Set(l.rows.map((r) => r.percent)).size, 1, 'one scale');
      expect((e.svg(fitted(SINGLE_COLUMN)) ?? '').length > 0, true, 'drawn');
      expect(e.drawsRouteSteps, true, 'the steps drawn from the route: needs the vendored chem-wasm (ENGINE-7)');
    });
  }));

  // "As drawn" beside the arrangements
  test('As drawn: offered with a drawing only; on, the preview, Save and Copy are the drawing\'s, at the style\'s ' +
    'print size, its map numbers and colours as the panel says', async () =>
    inEngine(ChemicalEngine.Crux, asDrawnExported));

  test('As drawn before crux-sketch\'s picture API: the canvas\'s own picture, said; cells with coordinates only',
    async () => inEngine(ChemicalEngine.Crux, asDrawnFallback));

  test('a CDXML cell\'s Export offers As drawn, the page\'s own picture; a SMILES cell\'s does not', async () =>
    inEngine(ChemicalEngine.Crux, cellAsDrawn));
});
