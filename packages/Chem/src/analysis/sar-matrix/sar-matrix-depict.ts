import * as DG from 'datagrok-api/dg';

import {getRdKitModule} from '../../utils/chem-common-rdkit';
import {getMoleculeRenderer} from '../../package';
import {SarMatrix} from './sar-matrix-types';
import {MatrixGridState} from './sar-matrix-ui-common';

/** Drawing molecules onto a grid canvas, aligned so a shared core points the same way in every cell.
 *
 *  Separate from the DOM-canvas helpers in `sar-matrix-ui-common`: those re-straighten a depiction,
 *  which destroys the alignment this path exists to preserve. The alignment layouts and the scratch
 *  canvas are module state, so everything that draws through them lives here. */

const ALIGN_OPTS = JSON.stringify({useCoordGen: true, allowRGroups: true, acceptFailure: true, alignOnly: true});

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
  const rdkit = getRdKitModule();
  let mol = null;
  let template: string | null = null;
  try {
    mol = rdkit.get_mol(coreSmiles.replace(/\[\*:\d+\]/g, '[H]'));
    if (mol && mol.is_valid()) {
      mol.set_new_coords();
      template = mol.get_molblock();
    }
  } catch {
    template = null;
  } finally {
    mol?.delete();
  }
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
  const rdkit = getRdKitModule();
  let template = null;
  let mol = null;
  let result = molStr;
  try {
    template = rdkit.get_mol(templateMolblock);
    mol = rdkit.get_mol(molStr);
    if (template && mol && mol.is_valid()) {
      mol.generate_aligned_coords(template, ALIGN_OPTS);
      result = mol.get_molblock() || molStr;
    }
  } catch {
    result = molStr;
  } finally {
    template?.delete();
    mol?.delete();
  }
  if (alignCache.size >= ALIGN_CACHE_MAX)
    alignCache.clear();
  alignCache.set(key, result);
  return result;
}

/** Drop the per-dataset core-alignment layouts on rebuild; the shared renderer's raster cache is LRU. */
export function clearDepictionCaches(): void {
  alignCache.clear();
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
  return coreBlockCache.getOrCreate(`${columnOrdinal} ${smiles}`, () => {
    let mol = null;
    let block = smiles;
    try {
      mol = getRdKitModule().get_mol(smiles);
      if (mol?.is_valid()) {
        mol.set_new_coords();
        block = starNonColumnSites(mol.get_molblock() || smiles, columnOrdinal);
      }
    } catch {
      block = smiles;
    } finally {
      mol?.delete();
    }
    return block;
  });
}

/**
 * The scaffold a matrix varies, with the R at the position its columns enumerate and a `*` at
 * every position its rows differ at.
 *
 * Every row of a decomposed matrix shares one core, numbered to match {@link SarMatrix.positions},
 * and the first of those is the column axis. Read off the core rather than guessed from the site
 * key's isotopes: which mark is the axis is not fixed, and guessing puts the R where the matrix
 * holds constant. Rows keep a mark, not a hydrogen, which would claim they are unsubstituted.
 */
export function matrixCore(matrix: SarMatrix): string {
  const core = matrix.rows[0]?.coreSmiles ?? '';
  const columnNumber = Number.parseInt((matrix.positions[0] ?? '').replace(/^\D+/, ''), 10);
  // Rows of a single-position matrix are each their own core, so none of them is the matrix's; its
  // key is, and the one attachment it carries is where the columns hang.
  if (!core || !Number.isFinite(columnNumber) || Object.keys(matrix.refValues).length === 0)
    return coreDepictionBlock((matrix.siteKey || core).replace(/\[\d+\*\]|\[\*:\d+\]/g, '[*]'), 0);
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

/** The molecule string a cell/core/header is drawn from: aligned to `template` when given, else the
 *  map-stripped string for the shared renderer. A molblock arrives ready to draw and passes through.
 *  `null` if empty. */
function preparedDepiction(molStr: string, template: string | null): string | null {
  if (!molStr)
    return null;
  const plain = molStr.replace(/\[\*:\d+\]/g, '[*]');
  return template ? alignToTemplate(plain, template) : plain;
}

/** Scratch canvas the cached `ImageData` is blitted through so the grid's clip is respected
 *  (`putImageData` ignores the clip and would paint over the pinned core column or past the edge). */
let blitCanvas: OffscreenCanvas | null = null;
function ensureBlitCanvas(w: number, h: number): OffscreenCanvas {
  if (!blitCanvas || blitCanvas.width < w || blitCanvas.height < h)
    blitCanvas = new OffscreenCanvas(Math.max(w, blitCanvas?.width ?? 0), Math.max(h, blitCanvas?.height ?? 0));
  return blitCanvas;
}

/** Draw a depiction onto the grid canvas in device pixels, reusing the shared renderer's LRU caches.
 *  The bitmap is produced at device size and blitted 1:1 at a whole-pixel offset so bond lines stay
 *  crisp (drawing through the grid's scaled transform bilinear-filters every bond). The device rect
 *  comes from the context transform, not `devicePixelRatio`, so a grid translate can't displace cells.
 *
 *  The molecule raster is rendered on a TRANSPARENT background (so it caches by molecule + size, not
 *  by the per-cell tint); `argb` is laid down here as the background and the molecule composited over
 *  it, which blends the anti-aliased bonds into the tint with no pale fringe. Both the fill and the
 *  blit run at device pixels through the scratch canvas, so the grid's clip is respected. */
export function drawDepiction(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
  molStr: string, template: string | null, argb: number): void {
  const renderer = getMoleculeRenderer();
  if (!renderer)
    return;
  const molblock = preparedDepiction(molStr, template);
  if (molblock === null)
    return;
  const m = g.getTransform();
  const dw = Math.round(m.a * w);
  const dh = Math.round(m.d * h);
  if (dw < 1 || dh < 1)
    return;
  let image: ImageData;
  try {
    image = renderer.getCachedMolImageData(molblock, dw, dh);
  } catch {
    return; // malformed structure — leave the cell to its background
  }
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
