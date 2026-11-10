/**
 * OligoNucleotide grid cell renderer.
 *
 * - Cell value is HELM (under the hood). The renderer parses once per
 *   distinct value (cached by reference) and draws a duplex view.
 * - No image cache: drawing is cheap; reparses only when the value changes.
 * - Tooltip on hover shows monomer details + an async-loaded RDKit
 *   structure for the hovered monomer.
 *
 * Self-contained: does not require `_package.initLibData()` to draw.
 * The library is consulted lazily inside the tooltip for structure rendering.
 */

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {drawDuplex, hitTest, DuplexLayout, RenderOpts, DEFAULT_OPTS} from './canvas-renderer';
import {looksLikeHelm, parseHelmDuplex} from './helm-parser';
import {ParsedDuplex} from './types';
import {showMonomerTooltip} from './tooltip';

const CELL_TYPE = 'OligoNucleotide';
/** The grid column's temp key the renderer keeps itself under, where hover tests find {@link getMonomerPoint}. */
export const OLIGO_RENDERER_TEMP = 'oligoNucleotideRenderer';

export class OligoNucleotideCellRenderer extends DG.GridCellRenderer {
  /** WeakMap-by-value cache of parsed HELM. Avoids reparsing on redraw. */
  private modelCache = new DG.LruCache<string, ParsedDuplex>();
  /** Last-rendered layout per cell key, for hit-testing on subsequent moves. */
  private layoutCache = new DG.LruCache<string, DuplexLayout>();

  get name(): string { return CELL_TYPE; }
  get cellType(): string { return CELL_TYPE; }
  get defaultWidth(): number | null { return 320; }
  get defaultHeight(): number | null { return 64; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle,
  ): void {
    const value = gridCell.cell.value as string | null;
    if (!value || !looksLikeHelm(value)) {
      // Render raw value as plain text — same fallback as default DG cell.
      g.save();
      g.beginPath();
      g.rect(x, y, w, h);
      g.clip();
      g.fillStyle = '#8b949e';
      g.font = '11px ui-monospace, Menlo, monospace';
      g.textBaseline = 'middle';
      g.textAlign = 'left';
      g.fillText(value ?? '', x + 4, y + h / 2);
      g.restore();
      return;
    }

    const model = this.getOrParse(value);
    if (!gridCell.cell?.dart || !gridCell.cell.column) {
      w = g.canvas.width;
      h = g.canvas.height;
      x = 0;
      y = 0;
      // this happens in forms....
    }
    else if (gridCell.gridColumn && gridCell.gridColumn.temp[OLIGO_RENDERER_TEMP] !== this)
      gridCell.gridColumn.temp[OLIGO_RENDERER_TEMP] = this;
    const layout = drawDuplex(g, x, y, w, h, model);
    // Cache for hit-test on subsequent mouse moves over the same cell.
    this.layoutCache.set(this.cellKey(value, w, h, DEFAULT_OPTS), layout);
  }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const value = gridCell.cell.value as string | null;
    if (!value || !looksLikeHelm(value)) {
      ui.tooltip.hide();
      return;
    }
    const model = this.getOrParse(value);
    const layout = this.layoutOf(gridCell, value, model);
    if (!layout) return;

    const bounds = gridCell.bounds;
    const localX = e.offsetX - bounds.x;
    const localY = e.offsetY - bounds.y;
    const hit = hitTest(localX, localY, model, layout);
    if (!hit) {
      ui.tooltip.hide();
      return;
    }

    showMonomerTooltip(hit, e.clientX, e.clientY);
  }

  override onMouseLeave(_gridCell: DG.GridCell, _e: MouseEvent): void {
    ui.tooltip.hide();
  }

  /** Where the monomer at {@link pos} (0-based over the sense strand's monomers, then the antisense's, as the
   * duplex parses) sits in the cell as laid out, in px from the cell's top left: the middle of a nucleotide's or a
   * conjugate's chip, of a standalone linker's arc. The point the pointer rests on for its tooltip. */
  getMonomerPoint(gridCell: DG.GridCell, pos: number): {x: number, y: number} | null {
    return this.pointOf(gridCell, pos, false);
  }

  /** The middle of the arc of the 3' linkage after the monomer at {@link pos} (counted as in
   * {@link getMonomerPoint}): where the pointer rests for the linkage's tooltip. */
  getLinkagePoint(gridCell: DG.GridCell, pos: number): {x: number, y: number} | null {
    return this.pointOf(gridCell, pos, true);
  }

  private pointOf(gridCell: DG.GridCell, pos: number, linkage: boolean): {x: number, y: number} | null {
    const value = gridCell.cell.value as string | null;
    if (!value || !looksLikeHelm(value))
      return null;
    const model = this.getOrParse(value);
    const layout = this.layoutOf(gridCell, value, model);
    const senseCount = model.sense.monomers.length;
    const sense = pos < senseCount;
    const monomer = (sense ? model.sense : model.antisense)?.monomers[sense ? pos : pos - senseCount];
    if (!layout || layout.textOnlyFallback || !monomer)
      return null;
    const chip = linkage ? undefined : (sense ? layout.senseChips : layout.antiChips).find((c) => c.monomer === monomer);
    if (chip)
      return {x: chip.x + chip.w / 2, y: (sense ? layout.senseY : layout.antiY) + layout.chipH / 2};
    const arc = (sense ? layout.senseLinks : layout.antiLinks).find((l) => l.ownerOrigIdx === monomer.position);
    if (!arc)
      return null;
    return {x: arc.x + arc.w / 2, y: sense ? arc.y - layout.apexH / 2 : arc.y + arc.h + layout.apexH / 2};
  }

  private layoutOf(gridCell: DG.GridCell, value: string, model: ParsedDuplex): DuplexLayout | null {
    const {width, height} = gridCell.bounds;
    return this.layoutCache.getOrCreate(this.cellKey(value, width, height, DEFAULT_OPTS),
      () => drawDuplex(null as unknown as CanvasRenderingContext2D, 0, 0, width, height, model, DEFAULT_OPTS, true));
  }

  private getOrParse(helm: string): ParsedDuplex {
    return this.modelCache.getOrCreate(helm, (h) => parseHelmDuplex(h));
  }

  /** Cache key for a cell's layout. Includes the column's `version` so any
   * edit to the column (which bumps version) orphans previous cache entries
   * — preventing onMouseMove from hit-testing a stale layout that was cached
   * before the edit and not yet replaced by a fresh render(). */
  private cellKey(value: string, width: number, height: number, opts: RenderOpts): string {
    return `${value}::${Math.floor(width)}x${Math.floor(height)}::${JSON.stringify(opts)}`;
  }
}
