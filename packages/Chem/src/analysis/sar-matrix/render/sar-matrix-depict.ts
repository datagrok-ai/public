import * as DG from 'datagrok-api/dg';

import {withEngineMol} from '../sar-matrix-engine';
import {engineGeneration} from '../../../engine/chem-engine';
import {getMoleculeRenderer} from '../../../package';
import {RDKitCellRenderer} from '../../../rendering/rdkit-cell-renderer';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {SarMatrix} from '../sar-matrix-types';
import {attachmentNumbers} from '../build/sar-matrix-link';
import {MatrixGridState} from '../sar-matrix-ui-common';

/** Drawing molecules onto a grid canvas, aligned so a shared core points the same way in every cell.
 *
 *  Separate from the DOM-canvas helpers in `sar-matrix-ui-common`: those re-straighten a depiction,
 *  which destroys the alignment this path exists to preserve. The alignment layouts and the scratch
 *  canvas are module state, so everything that draws through them lives here. */

// RDKit's own depictor, not CoordGen: the free layout alignOnly fits to the template is ~20x faster per cell
// (CoordGen's clash minimization: 9 ms per cell, a 1,000-cell matrix ~9 s on RDKit.js 2024.09) and the templates are
// RDKit's too
const ALIGN_OPTS = JSON.stringify({useCoordGen: false, allowRGroups: true, acceptFailure: true, alignOnly: true});

/** Aligned-molblock cache keyed by template + molecule; cleared wholesale past the cap to avoid leaks. */
const alignCache = new Map<string, string>();
const ALIGN_CACHE_MAX = 4000;
/** Per-core alignment templates by core SMILES. `null` is cached too, so failures aren't retried. */
const templateCache = new Map<string, string | null>();
const TEMPLATE_CACHE_MAX = 4000;

/** Canonical molblock template every cell/core aligns to, so the shared core points the same way.
 *  `[*:n]` dummies are swapped for H so the template is a plain substructure. */
export function buildAlignmentTemplate(coreSmiles: string): string | null {
  const cached = templateCache.get(coreSmiles);
  if (cached !== undefined)
    return cached;
  // "Chemical engine": laid out by the main thread's module (Crux's in Crux mode, which draws the cells)
  const template = withEngineMol('SAR template', coreSmiles.replace(/\[\*:\d+\]/g, '[H]'), (mol) => {
    mol.set_new_coords();
    return mol.get_molblock();
  });
  if (templateCache.size >= TEMPLATE_CACHE_MAX)
    templateCache.clear();
  templateCache.set(coreSmiles, template);
  return template;
}

/** Aligned molblock for `molStr` oriented to the template, or the input unchanged if it can't align. */
export function alignToTemplate(molStr: string, templateMolblock: string): string {
  const key = `${templateMolblock}${molStr}`;
  const cached = alignCache.get(key);
  if (cached !== undefined)
    return cached;
  // "Chemical engine": aligned by the main thread's module (Crux's in Crux mode, which draws the cells), RDKit's where
  // Crux reads neither or throws
  const result = withEngineMol('SAR alignment', molStr, (mol, module) => {
    const template = module.get_mol(templateMolblock);
    try {
      if (!template)
        return null;
      mol.generate_aligned_coords(template, ALIGN_OPTS);
      return mol.get_molblock() || molStr;
    } finally {
      template?.delete();
    }
  }) ?? molStr;
  if (alignCache.size >= ALIGN_CACHE_MAX)
    alignCache.clear();
  alignCache.set(key, result);
  return result;
}

/** Drop the per-dataset core-alignment layouts on rebuild; the shared renderer's raster cache is LRU. */
export function clearDepictionCaches(): void {
  alignCache.clear();
  images.clear();
  imageBytes = 0;
}

const coreBlockCache = new DG.LruCache<string, string>(TEMPLATE_CACHE_MAX);

/** Writes every attachment point except the column's as `*` instead of `R`. Matrices keying on one
 *  scaffold differ only in which position they vary, so that position must be legible at the size a
 *  core is drawn — two distinct symbols carry it where a superscript on one letter does not.
 *  `columnOrdinal` counts attachment points in atom order. */
function starNonColumnSites(molblock: string, columnOrdinal: number): string {
  const lines = molblock.split('\n');
  const atomCount = Number.parseInt((lines[3] ?? '').slice(0, 3), 10);
  if (!Number.isFinite(atomCount))
    return molblock;
  let seen = 0;
  for (let i = 0; i < atomCount; i++) {
    const line = lines[4 + i] ?? '';
    // Columns 31-33 of a V2000 atom line hold the element symbol; RDKit writes a dummy there as `R`.
    if (line.slice(31, 34).trim() !== 'R')
      continue;
    if (seen !== columnOrdinal)
      lines[4 + i] = `${line.slice(0, 31)}*  ${line.slice(34)}`;
    seen++;
  }
  return lines.join('\n');
}

/** A core as a laid-out molblock, so its attachment points draw as labels at all: RDKit names a dummy
 *  only when it reads one from a molblock, and from SMILES the same atom draws as an unlabelled stub.
 *  Unparseable input is returned as-is, so a bad key still draws something. */
export function coreDepictionBlock(smiles: string, columnOrdinal: number): string {
  if (!smiles)
    return smiles;
  // "Chemical engine": laid out by the main thread's module (Crux's in Crux mode, which draws it)
  return coreBlockCache.getOrCreate(`${columnOrdinal} ${smiles}`, () => withEngineMol('SAR core depiction', smiles,
    (mol) => {
      mol.set_new_coords();
      return starNonColumnSites(mol.get_molblock() || smiles, columnOrdinal);
    }) ?? smiles);
}

/**
 * The scaffold a matrix varies, with the R at the position its columns enumerate and a `*` at
 * every position its rows differ at. The axis is read off the attachment points the column
 * substituents carry, since position names need not be numbered.
 */
export function matrixCore(matrix: SarMatrix): string {
  const core = matrix.rows[0]?.coreSmiles ?? '';
  // Rows of a single-position matrix are each their own core, so none of them is the matrix's; its
  // key is, and the one attachment it carries is where the columns hang.
  if (!core || Object.keys(matrix.refValues).length === 0)
    return coreDepictionBlock((matrix.siteKey || core).replace(/\[\d+\*\]|\[\*:\d+\]/g, '[*]'), 0);
  const coreNumbers = attachmentNumbers(core);
  const columnNumber = matrix.columns.flatMap((c) => [...attachmentNumbers(c.substSmiles)])
    .find((n) => coreNumbers.has(n));
  let seen = 0;
  let ordinal = -1;
  const smiles = core.replace(/\[\*:(\d+)\]/g, (_match, number) => {
    if (Number(number) === columnNumber)
      ordinal = seen;
    seen++;
    return '[*]';
  });
  return coreDepictionBlock(smiles, ordinal);
}

/** The molblock a row's cells align to: the row key aligned to the pane's shared core (so cells
 *  match the core printed beside them), else laid out on its own. */
export function rowTemplate(state: MatrixGridState, rowIdx: number): string | null {
  const cached = state.rowTemplates[rowIdx];
  if (cached !== null)
    return cached;
  const paneRow = state.rows[rowIdx];
  const key = paneRow.matrix.rows[paneRow.rowIndex].keySmiles;
  const aligned = state.paneTemplate !== null ? alignToTemplate(key, state.paneTemplate) : '';
  // A key that can't align to the shared core comes back untouched, so lay it out on its own.
  const template = aligned.includes('V2000') ? aligned : buildAlignmentTemplate(key);
  state.rowTemplates[rowIdx] = template;
  return template;
}

/** Depiction work one task may do before the rest waits for a repaint: the first paint of a matrix of
 *  large molecules aligns and rasterizes every visible cell, seconds in which the page answers nothing. */
const PAINT_BUDGET_MS = 50;
let spent = 0;
let resetQueued = false;
/** The cells' rasters (`null`: one the renderer cannot draw), by engine, size, template and molecule: a paint
 *  blits these without spending its budget. Kept here because the shared renderer keeps 100 images and a zoomed-out
 *  matrix shows hundreds of cells: with only its cache behind them, every paint redrew the evicted cells, spent its
 *  budget on them and deferred the cells painted last for good. Bounded by bytes, not by count (so a Map, not a
 *  DG.LruCache): a raster is its cell in device pixels, 126 KB for a 104x76 cell on a Retina screen and a few KB
 *  zoomed far out, while the visible cells' rasters add up to at most the screen's device pixels at any zoom. */
const images = new Map<string, ImageData | null>();
/** A 5K screen covered in cells (5120 x 2880 x 4 bytes, 59 MB), and some to scroll back to. */
const IMAGES_MAX_BYTES = 64 * 2 ** 20;
let imageBytes = 0;

function cachedImage(key: string): ImageData | null | undefined {
  const image = images.get(key);
  if (image !== undefined) { // most recently used last, as the eviction below reads them oldest first
    images.delete(key);
    images.set(key, image);
  }
  return image;
}

function cacheImage(key: string, image: ImageData | null): void {
  images.set(key, image);
  imageBytes += image?.data.byteLength ?? 0;
  for (const [oldKey, old] of images) {
    if (imageBytes <= IMAGES_MAX_BYTES)
      break;
    images.delete(oldKey);
    imageBytes -= old?.data.byteLength ?? 0;
  }
}

function charge(ms: number): void {
  spent += ms;
  if (!resetQueued) {
    resetQueued = true;
    setTimeout(() => {
      spent = 0;
      resetQueued = false;
    });
  }
}

/** Whether this task has used its depiction budget; a painter defers what is not drawn yet past it. */
export function depictionBudgetSpent(): boolean {
  return spent >= PAINT_BUDGET_MS;
}

/** Runs depiction work, charging its time to this task's budget. */
export function charged<T>(work: () => T): T {
  const t0 = performance.now();
  try {
    return work();
  } finally {
    charge(performance.now() - t0);
  }
}

/** The package's molecule renderer, or one of this module's own where the package's init has not run (the test
 *  bundle), so the grid draws either way. */
let ownRenderer: RDKitCellRenderer | null = null;
function moleculeRenderer(): RDKitCellRenderer | null {
  return getMoleculeRenderer() ?? (ownRenderer ??= new RDKitCellRenderer(getRdKitModule()));
}

/** Scratch canvas the cached `ImageData` is blitted through so the grid's clip is respected
 *  (`putImageData` ignores the clip and would paint over the pinned core column or past the edge). */
let blitCanvas: OffscreenCanvas | null = null;
function ensureBlitCanvas(w: number, h: number): OffscreenCanvas {
  if (!blitCanvas || blitCanvas.width < w || blitCanvas.height < h)
    blitCanvas = new OffscreenCanvas(Math.max(w, blitCanvas?.width ?? 0), Math.max(h, blitCanvas?.height ?? 0));
  return blitCanvas;
}

/** Draw a depiction onto the grid canvas in device pixels, from this module's image cache or the shared renderer.
 *  The bitmap is produced at device size and blitted 1:1 at a whole-pixel offset so bond lines stay
 *  crisp (drawing through the grid's scaled transform bilinear-filters every bond). The device rect
 *  comes from the context transform, not `devicePixelRatio`, so a grid translate can't displace cells.
 *
 *  The molecule raster is rendered on a TRANSPARENT background (so it caches by molecule + size, not
 *  by the per-cell tint); `argb` is laid down here as the background and the molecule composited over
 *  it, which blends the anti-aliased bonds into the tint with no pale fringe. Both the fill and the
 *  blit run at device pixels through the scratch canvas, so the grid's clip is respected.
 *
 *  Past this task's budget, a depiction not in the cache is left out and `deferred` called instead, to
 *  paint it on a later pass. The map-stripped string is drawn, aligned to `template` when given. */
export function drawDepiction(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
  molStr: string, template: string | null, argb: number, deferred: () => void): void {
  const renderer = moleculeRenderer();
  if (!renderer || !molStr)
    return;
  const m = g.getTransform();
  const dw = Math.round(m.a * w);
  const dh = Math.round(m.d * h);
  if (dw < 1 || dh < 1)
    return;
  const plain = molStr.replace(/\[\*:\d+\]/g, '[*]');
  const key = `${engineGeneration()}|${dw}x${dh}|${template ?? ''}|${plain}`;
  let image = cachedImage(key);
  if (image === undefined) {
    if (depictionBudgetSpent()) {
      deferred();
      return;
    }
    const t0 = performance.now();
    try {
      image = renderer.getCachedMolImageData(template ? alignToTemplate(plain, template) : plain, dw, dh);
    } catch {
      image = null; // malformed: cached too, or the grid would be repainted for it as long as it shows
    } finally {
      charge(performance.now() - t0);
    }
    cacheImage(key, image);
  }
  if (image === null)
    return; // left to its background
  const scratch = ensureBlitCanvas(dw, dh);
  scratch.getContext('2d')!.putImageData(image, 0, 0);
  g.save();
  g.setTransform(1, 0, 0, 1, 0, 0);
  const destX = Math.round(m.a * x + m.e);
  const destY = Math.round(m.d * y + m.f);
  g.fillStyle = DG.Color.toHtml(argb);
  g.fillRect(destX, destY, dw, dh);
  g.drawImage(scratch, 0, 0, dw, dh, destX, destY, dw, dh);
  g.restore();
}
