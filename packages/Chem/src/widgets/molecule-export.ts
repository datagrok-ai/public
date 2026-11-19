/**
 * The Export dialog for a drawn molecule (Apps | Chem | Sketcher's ribbon, Export...), laid out as the reactions'
 * Export panel (widgets/reaction-export.ts): its choices over a preview fitted to the resizable dialog, a note under
 * it, Save and Copy.
 * - Two styles, for the pictures and CDXML:
 *   - "ACS Document 1996": the molecule laid out afresh, as Crux lays out a molecule read from its text, and drawn in
 *     ACS 1996's look: black, bonds 0.508 cm (14.4 pt), lines 0.6 pt, labels 10 pt; the SVG sized in points, at print
 *     size. The reactions' panel lays its scheme out the same way in this style;
 *   - "As drawn": the drawing as it is on the canvas, the sketcher's own coordinates and look (its atom colours, its
 *     20 pt bond); the SVG sized in CSS px.
 *   Both are Crux Sketch's pictures (`svgWith`, `png`) and its CDXML (ChemDraw's ACS Document 1996 document style,
 *   whichever the coordinates): of the drawing itself, or of a copy laid out afresh on the page's one offscreen
 *   sketcher.
 * - Formats: SVG; PNG at 300 dpi (its pHYs chunk saying so); CDXML; MOL V2000 and V3000, SDF, SMILES, CXSMILES, SMARTS,
 *   InChI and InChIKey, which are the drawing's whatever the style (its own coordinates in the MOL and SDF). A format
 *   the drawing has no text in (a query has no SMILES, InChI holds no query or R label) says why, and Save and Copy are
 *   off.
 * - Save downloads the format's file; Copy puts it on the clipboard (the PNG as a picture, every other format as text).
 * - The background white or transparent (SVG and PNG).
 */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {SketcherHandle, SketcherFormat} from '../../vendor/crux-sketch/crux-sketch';
import {loadCruxSketchForReactions} from '../crux/crux-sketcher';
import {ACS_1996, AS_DRAWN, previewZoom, svgPixels, withoutBackground} from './reaction-export';
import {withResolution} from '../utils/reactions/reaction-document';
import {sketcherHelpUrl} from '../utils/sketcher-app/sketcher-help';

export {ACS_1996, AS_DRAWN};
/** The styles, in the order the dialog offers them; the first is chosen as it opens, as in the reactions' panel. */
export const STYLES = [ACS_1996, AS_DRAWN];

/** The formats, in the order the dialog offers them. */
export const SVG = 'SVG';
export const PNG = 'PNG (300 dpi)';
export const CDXML = 'CDXML';
export const MOL_V2000 = 'MOL V2000';
export const MOL_V3000 = 'MOL V3000';
export const SDF = 'SDF';
export const SMILES = 'SMILES';
export const CXSMILES = 'CXSMILES';
export const SMARTS = 'SMARTS';
export const INCHI = 'InChI';
export const INCHIKEY = 'InChIKey';
export const FORMATS = [SVG, PNG, CDXML, MOL_V2000, MOL_V3000, SDF, SMILES, CXSMILES, SMARTS, INCHI, INCHIKEY];
/** The formats the style applies to; the others are the drawing's text. */
export const STYLED_FORMATS = [SVG, PNG, CDXML];
/** The PNG's resolution. */
export const PNG_DPI = 300;

/** Each format's file and media type. */
export const FORMAT_FILES: {[format: string]: {file: string, type: string}} = {
  [SVG]: {file: 'molecule.svg', type: 'image/svg+xml'},
  [PNG]: {file: 'molecule.png', type: 'image/png'},
  [CDXML]: {file: 'molecule.cdxml', type: 'chemical/x-cdxml'},
  [MOL_V2000]: {file: 'molecule.mol', type: 'chemical/x-mdl-molfile'},
  [MOL_V3000]: {file: 'molecule.mol', type: 'chemical/x-mdl-molfile'},
  [SDF]: {file: 'molecule.sdf', type: 'chemical/x-mdl-sdfile'},
  [SMILES]: {file: 'molecule.smi', type: 'chemical/x-daylight-smiles'},
  [CXSMILES]: {file: 'molecule.cxsmiles', type: 'text/plain'},
  [SMARTS]: {file: 'molecule.sma', type: 'text/plain'},
  [INCHI]: {file: 'molecule.inchi', type: 'chemical/x-inchi'},
  [INCHIKEY]: {file: 'molecule.inchikey', type: 'text/plain'},
};

/** ACS 1996's bond (points) and the PNG's scale against the pictures' natural size (one px a CSS px, 96 dpi). */
const ACS_BOND_PT = 14.4;
const PNG_SCALE = PNG_DPI / 96;

/** The dialog's title, its classes (the tests find the preview, the text and the note by them), and its buttons. */
export const EXPORT_TITLE = 'Export molecule';
export const PREVIEW_CLASS = 'chem-molecule-export-preview';
export const TEXT_CLASS = 'chem-molecule-export-text';
export const NOTE_CLASS = 'chem-molecule-export-note';
export const CHOICES_CLASS = 'chem-molecule-export-choices';
export const SAVE = 'Save';
export const COPY = 'Copy';

/** What the note says of a format the drawing has no text in. */
export const QUERY_HAS_NO_SMILES = 'A drawing with query features has no SMILES: choose SMARTS';
export const SMILES_LOSES = 'SMILES leaves out what CXSMILES keeps of this molecule';
export const V2000_AS_V3000 = 'Written as V3000: MOL V2000 cannot hold this molecule';
export const NOTHING_DRAWN = 'There is nothing drawn to export';
export const NOT_LAID_OUT = 'The molecule could not be laid out afresh';

/** The preview's box as the dialog opens (CSS px), its least height, its padding, the most it enlarges a picture. */
const PREVIEW_OPEN_WIDTH = 640;
const PREVIEW_OPEN_HEIGHT = 380;
const PREVIEW_MIN_HEIGHT = 120;
const PREVIEW_DIALOG_CHROME = 300;
const PREVIEW_PAD = 12;

/** A picture: its SVG (CSS px at its natural size) and that size. */
export interface MoleculePicture {
  svg: string;
  width: number;
  height: number;
}

/** What one choice exports: the content and its file, or why it cannot. `note`: what the dialog says of it. */
export type ExportResult = {content: string | Blob, file: string, type: string, note: string, text: string | null} |
  {error: string};

/** The page's one offscreen sketcher, which lays a molecule out afresh for the ACS style (made on first use, again
 * after it failed to start). */
let arranger: Promise<SketcherHandle> | null = null;
function theArranger(): Promise<SketcherHandle> {
  arranger ??= (async () => {
    const host = document.createElement('div');
    Object.assign(host.style, {position: 'absolute', left: '-10000px', top: '0', width: '400px', height: '300px'});
    document.body.append(host);
    try {
      return await (await loadCruxSketchForReactions()).createSketcher(host, {mode: 'query', toolbar: {}});
    } catch (e) {
      host.remove();
      throw e;
    }
  })().catch((e) => {
    arranger = null;
    throw e;
  });
  return arranger;
}

/** A Crux Sketch error's or an exception's words. */
function message(e: unknown): string {
  return String((e as {message?: string})?.message ?? e);
}

/** The text of a drawing without coordinates for it to be laid out afresh: its CXSMILES, or (a query, which has no
 * SMILES) its SMARTS; null for an empty drawing. */
export function textToLayOut(sketch: SketcherHandle): {text: string, format: SketcherFormat} | null {
  if (sketch.isEmpty)
    return null;
  const cxsmiles = sketch.cxsmiles;
  return cxsmiles !== '' ? {text: cxsmiles, format: 'cxsmiles'} : {text: sketch.smarts, format: 'smarts'};
}

/** A drawn molecule exported: its pictures and texts in each style, made once each until the drawing changes. */
export class MoleculeExport {
  private laidOut: Promise<SketcherHandle | null> | null = null;
  private readonly made = new Map<string, Promise<MoleculePicture | null>>();

  constructor(readonly sketch: SketcherHandle) {}

  /** The sketcher whose drawing style `style` pictures: the drawing itself, or (ACS 1996) the page's offscreen
   * sketcher holding it laid out afresh; null when that could not be made. */
  sketcherOf(style: string): Promise<SketcherHandle | null> {
    if (style === AS_DRAWN)
      return Promise.resolve(this.sketch);
    this.laidOut ??= (async () => {
      const text = textToLayOut(this.sketch);
      if (text === null)
        return null;
      const copy = await theArranger();
      copy.setValue(text.text, {format: text.format});
      return copy.lastError === null && !copy.isEmpty ? copy : null;
    })().catch((e) => {
      console.warn('Chem | the molecule was not laid out afresh', e);
      return null;
    });
    return this.laidOut;
  }

  /** The picture options of `style`: ACS 1996's bond, black; the sketcher's own look as drawn. */
  static pictureOptions(style: string): {bondLength?: number, colour?: boolean} {
    return style === AS_DRAWN ? {} : {bondLength: ACS_BOND_PT, colour: false};
  }

  /** The picture of the drawing in `style` (CSS px; its background left out when `transparent`), or null. The
   * offscreen sketcher is read at once, before another export can set it. */
  picture(style: string, transparent: boolean): Promise<MoleculePicture | null> {
    const key = JSON.stringify([style, transparent]);
    let job = this.made.get(key);
    if (!job) {
      job = this.sketcherOf(style).then((s) => {
        if (s === null || s.isEmpty)
          return null;
        const drawn = s.svgWith(MoleculeExport.pictureOptions(style));
        const svg = transparent ? withoutBackground(drawn) : drawn;
        const size = svgPixels(svg);
        return size ? {svg, ...size} : null;
      }).catch((e) => {
        console.warn('Chem | the molecule was not pictured', e);
        return null;
      });
      this.made.set(key, job);
    }
    return job;
  }

  /** The SVG file of `style`: sized in points for ACS 1996 (placed at print size), CSS px as drawn; an XML prolog. */
  async svgFile(style: string, transparent: boolean): Promise<string | null> {
    const p = await this.picture(style, transparent);
    if (!p)
      return null;
    const unit = style === AS_DRAWN ? {k: 1, unit: 'px'} : {k: 0.75, unit: 'pt'};
    const round = (v: number) => Math.round(v * unit.k * 100) / 100;
    const root = /<svg\b[^>]*>/.exec(p.svg)![0];
    let sized = root.replace(/\swidth\s*=\s*['"][^'"]*['"]/, ` width='${round(p.width)}${unit.unit}'`)
      .replace(/\sheight\s*=\s*['"][^'"]*['"]/, ` height='${round(p.height)}${unit.unit}'`);
    if (!/\sviewBox\s*=/.test(sized))
      sized = sized.replace(/^<svg\b/, `<svg viewBox='0 0 ${p.width} ${p.height}'`);
    return `<?xml version='1.0' encoding='UTF-8'?>\n${p.svg.replace(root, sized)}`;
  }

  /** The PNG of `style` at {@link PNG_DPI} (Crux Sketch's `png`, white; a transparent one, the SVG drawn with no
   * background), or null. */
  async png(style: string, transparent: boolean): Promise<Blob | null> {
    if (!transparent) {
      const s = await this.sketcherOf(style);
      return s === null || s.isEmpty ? null : s.png({...MoleculeExport.pictureOptions(style), scale: PNG_SCALE});
    }
    const p = await this.picture(style, true);
    if (!p)
      return null;
    const image = await svgImage(p.svg);
    const canvas = new OffscreenCanvas(Math.round(p.width * PNG_SCALE), Math.round(p.height * PNG_SCALE));
    canvas.getContext('2d')!.drawImage(image, 0, 0, canvas.width, canvas.height);
    return withResolution(await canvas.convertToBlob({type: 'image/png'}), PNG_DPI);
  }

  /** The CDXML of `style` (Crux Sketch's writer, its module loaded first), or null. */
  async cdxml(style: string): Promise<string | null> {
    const s = await this.sketcherOf(style);
    if (s === null || s.isEmpty)
      return null;
    await s.loadFormat('cdxml');
    return s.cdxml || null;
  }

  /** The drawing as text in `format` (one of the unstyled formats) and what the note adds, or why it has none. */
  async text(format: string): Promise<{text: string, note: string} | {error: string}> {
    const s = this.sketch;
    if (s.isEmpty)
      return {error: NOTHING_DRAWN};
    try {
      switch (format) {
      case MOL_V2000: {
        const molfile = s.molfile;
        return {text: molfile, note: molfile.includes('V3000') ? V2000_AS_V3000 : ''};
      }
      case MOL_V3000:
        return {text: s.molV3000, note: ''};
      case SDF: {
        const molfile = s.molfile;
        return {text: `${molfile}${molfile.endsWith('\n') ? '' : '\n'}$$$$\n`,
          note: molfile.includes('V3000') ? V2000_AS_V3000 : ''};
      }
      case SMILES: {
        const smiles = s.plainSmiles;
        if (smiles === '')
          return {error: QUERY_HAS_NO_SMILES};
        return {text: smiles, note: s.smiles !== smiles ? `${SMILES_LOSES}: choose CXSMILES` : ''};
      }
      case CXSMILES: {
        const cxsmiles = s.cxsmiles;
        return cxsmiles === '' ? {error: QUERY_HAS_NO_SMILES} : {text: cxsmiles, note: ''};
      }
      case SMARTS:
        return {text: s.smarts, note: ''};
      case INCHI:
      case INCHIKEY: {
        await s.loadFormat('inchi');
        const text = format === INCHI ? s.inchi : s.inchiKey;
        if (text === '') {
          const why = s.warnings.find((w) => w.kind === 'format' && w.format === 'inchi');
          return {error: `${format}: InChI cannot hold this drawing${why ? ` (${why.message})` : ''}`};
        }
        return {text, note: ''};
      }
      }
    } catch (e) {
      return {error: `${format}: ${message(e)}`};
    }
    return {error: `${format} is not a format of this dialog`};
  }

  /** What `format` in `style` exports, and what the note says. */
  async result(format: string, style: string, transparent: boolean): Promise<ExportResult> {
    const {file, type} = FORMAT_FILES[format];
    if (this.sketch.isEmpty)
      return {error: NOTHING_DRAWN};
    if (!STYLED_FORMATS.includes(format)) {
      const t = await this.text(format);
      return 'error' in t ? t : {content: t.text, file, type, note: t.note, text: t.text};
    }
    const p = await this.picture(style, transparent);
    if (p === null)
      return {error: style === AS_DRAWN ? 'The drawing could not be pictured' : NOT_LAID_OUT};
    const size = sizeNote(p, style);
    try {
      if (format === SVG) {
        const svg = await this.svgFile(style, transparent);
        return svg === null ? {error: 'The SVG was not written'} : {content: svg, file, type, note: size, text: svg};
      }
      if (format === PNG) {
        const png = await this.png(style, transparent);
        return png === null ? {error: 'The PNG was not written'} :
          {content: png, file, type, note: `${size}; ${Math.round(p.width * PNG_SCALE)} × ` +
            `${Math.round(p.height * PNG_SCALE)} px at ${PNG_DPI} dpi`, text: null};
      }
      const cdxml = await this.cdxml(style);
      return cdxml === null ? {error: 'The CDXML was not written'} : {content: cdxml, file, type, text: cdxml,
        note: `CDXML for ChemDraw (its ACS Document 1996 document style), ${style === AS_DRAWN ?
          'the molecule as drawn' : 'the molecule laid out afresh'}`};
    } catch (e) {
      return {error: `${format}: ${message(e)}`};
    }
  }
}

/** What the note says of picture `p` in `style`: its size on paper. */
export function sizeNote(p: MoleculePicture, style: string): string {
  const cm = (px: number) => (px / 96 * 2.54).toFixed(1);
  return `${cm(p.width)} × ${cm(p.height)} cm, ${style === AS_DRAWN ? 'as drawn on the canvas' :
    'laid out afresh, ACS Document 1996 (black, bonds 0.508 cm), at print size'}`;
}

/** `svg` as an image, once it is decoded. */
function svgImage(svg: string): Promise<HTMLImageElement> {
  return new Promise((resolve, reject) => {
    const image = new Image();
    image.onload = () => resolve(image);
    image.onerror = () => reject(new Error('the picture could not be decoded'));
    image.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
  });
}

/** The dialog's inputs (the tests set them as a person would). */
export interface MoleculeExportInputs {
  style: DG.InputBase<string | null>;
  format: DG.InputBase<string | null>;
  transparent: DG.InputBase<boolean | null>;
}

/** An open dialog's inputs and export, hung on its root by a registered symbol (the tests reach them by it; Chem's
 * package tests are a bundle of their own, so a module's WeakMap there is not the dialog's). */
type OpenExport = {inputs: MoleculeExportInputs, exported: MoleculeExport, shown: () => Promise<ExportResult>};
const EXPORT_KEY = Symbol.for('datagrok.chem.molecule-export');
type ExportRoot = HTMLElement & {[EXPORT_KEY]?: OpenExport};

/** Tests: the inputs of the Export dialog whose root is `root`, its export, and what it shows now (once made). */
export function moleculeExportOf(root: HTMLElement): OpenExport | undefined {
  return (root as ExportRoot)[EXPORT_KEY];
}

/** Opens the Export dialog on the drawing `sketch` holds. */
export function openMoleculeExport(sketch: SketcherHandle): DG.Dialog {
  const exported = new MoleculeExport(sketch);
  let style = ACS_1996;
  let format = SVG;
  let transparent = false;
  const openHeight = Math.max(PREVIEW_MIN_HEIGHT,
    Math.min(PREVIEW_OPEN_HEIGHT, window.innerHeight - PREVIEW_DIALOG_CHROME));
  const preview = ui.div([], {classes: PREVIEW_CLASS, style: {flex: '1 1 auto', height: `${openHeight}px`,
    minHeight: `${PREVIEW_MIN_HEIGHT}px`, position: 'relative', overflow: 'hidden', padding: `${PREVIEW_PAD}px`,
    background: 'var(--grey-1)', boxSizing: 'border-box'}});
  const note = ui.divText('', NOTE_CLASS);
  note.style.color = 'var(--grey-5)';
  const styleInput = ui.input.choice('Style', {items: STYLES, value: style, onValueChanged: (v) => {
    style = v ?? ACS_1996;
    refresh();
  }});
  styleInput.setTooltip('ACS Document 1996: laid out afresh, black, bonds 0.508 cm. As drawn: the canvas\'s ' +
    'drawing, its coordinates and colours. For SVG, PNG and CDXML');
  const formatInput = ui.input.choice('Format', {items: FORMATS, value: format, onValueChanged: (v) => {
    format = v ?? SVG;
    refresh();
  }});
  const transparentInput = ui.input.bool('Transparent background', {value: false, onValueChanged: (v) => {
    transparent = !!v;
    refresh();
  }});
  const dialog = ui.dialog({title: EXPORT_TITLE, helpUrl: sketcherHelpUrl('export')});
  /** The result shown now; the turn of the last one asked for (a later choice overtakes it). */
  let shown: Promise<ExportResult> = Promise.resolve({error: NOTHING_DRAWN});
  let turn = 0;
  /** The picture shown, fitted to the box again as the dialog is resized. */
  let picture: MoleculePicture | null = null;
  let fittedTo = '';
  const enable = (on: boolean) => {
    for (const name of [SAVE, COPY]) {
      const button = dialog.getButton(name);
      if (button)
        ui.setClass(button, 'disabled', !on);
    }
  };
  const fit = () => {
    ui.empty(preview);
    fittedTo = `${preview.clientWidth}x${preview.clientHeight}`;
    if (!picture)
      return;
    const zoom = previewZoom(picture.width, picture.height, preview.clientWidth - 2 * PREVIEW_PAD,
      preview.clientHeight - 2 * PREVIEW_PAD);
    const image = document.createElement('img');
    image.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(picture.svg)}`;
    // a transparent picture over a chequer
    const chequer = transparent && format !== CDXML;
    Object.assign(image.style, {position: 'absolute', inset: '0', margin: 'auto', display: 'block',
      width: `${picture.width * zoom}px`, height: `${picture.height * zoom}px`, backgroundSize: '16px 16px',
      backgroundImage: chequer ? 'repeating-conic-gradient(#e6e6e6 0% 25%, #ffffff 0% 50%)' : 'none',
      backgroundColor: chequer ? 'transparent' : '#FFFFFF'});
    preview.append(image);
    preview.dataset.zoom = zoom.toFixed(3);
  };
  const showText = (text: string) => {
    ui.empty(preview);
    picture = null;
    const box = ui.element('textarea') as HTMLTextAreaElement;
    box.readOnly = true;
    box.value = text;
    box.classList.add(TEXT_CLASS);
    const inner = `calc(100% - ${2 * PREVIEW_PAD}px)`;
    Object.assign(box.style, {position: 'absolute', inset: `${PREVIEW_PAD}px`, width: inner, height: inner,
      resize: 'none', fontFamily: 'monospace', fontSize: '12px', boxSizing: 'border-box', whiteSpace: 'pre',
      overflow: 'auto'});
    preview.append(box);
  };
  const refresh = () => {
    const styled = STYLED_FORMATS.includes(format);
    styleInput.enabled = styled;
    transparentInput.enabled = format === SVG || format === PNG;
    const mine = ++turn;
    Object.assign(preview.dataset, {format, style: styled ? style : '', ready: 'false'});
    note.textContent = 'Writing…';
    enable(false);
    shown = exported.result(format, style, transparent);
    const pictureJob = styled ? exported.picture(style, format === CDXML ? false : transparent) : null;
    void Promise.all([shown, pictureJob]).then(([r, p]) => {
      if (mine !== turn || !dialog.root.isConnected)
        return;
      if ('error' in r) {
        ui.empty(preview);
        picture = null;
        note.textContent = r.error;
      } else {
        note.textContent = [r.note].filter((t) => t).join('. ');
        if (styled) {
          picture = p;
          fit();
        } else
          showText(r.text ?? '');
      }
      enable(!('error' in r));
      preview.dataset.ready = 'true';
    });
  };
  const save = async () => {
    const r = await shown;
    if ('error' in r) {
      grok.shell.warning(`Not saved: ${r.error}`);
      return;
    }
    DG.Utils.download(r.file, r.content, r.type);
    dialog.close(); // as the reactions' panel: Save is done with the dialog; Copy leaves it open
  };
  const copy = async () => {
    const r = await shown;
    try {
      if ('error' in r)
        throw new Error(r.error);
      if (r.content instanceof Blob)
        await navigator.clipboard.write([new ClipboardItem({[r.type]: r.content})]);
      else
        await navigator.clipboard.writeText(r.content);
      grok.shell.info(`Copied as ${format}`);
    } catch (e) {
      grok.shell.warning(`Not copied: ${message(e)}`);
    }
  };
  const resized = new ResizeObserver(() => {
    if (preview.isConnected && picture && fittedTo !== `${preview.clientWidth}x${preview.clientHeight}`)
      fit();
  });
  const choices = ui.div([styleInput.root, formatInput.root, transparentInput.root], CHOICES_CLASS);
  Object.assign(choices.style, {display: 'grid', columnGap: '16px',
    gridTemplateColumns: 'repeat(auto-fill, minmax(300px, 1fr))'});
  dialog.add(choices).add(preview).add(note)
    .addButton(SAVE, () => save(), 0, 'Download the file (and close the dialog)')
    .addButton(COPY, () => copy(), 1, 'On the clipboard: the PNG as a picture, every other format as text');
  dialog.onClose.subscribe(() => resized.disconnect());
  dialog.show({resizable: true, width: PREVIEW_OPEN_WIDTH + 4 * PREVIEW_PAD});
  const contents = dialog.root.querySelector('.d4-dialog-contents') as HTMLElement | null;
  if (contents) {
    contents.style.maxHeight = 'none';
    Object.assign(contents.style, {display: 'flex', flexDirection: 'column'});
  }
  refresh();
  resized.observe(preview);
  (dialog.root as ExportRoot)[EXPORT_KEY] = {inputs: {style: styleInput, format: formatInput,
    transparent: transparentInput}, exported, shown: () => shown};
  return dialog;
}
