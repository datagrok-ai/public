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
    const layout = this.layoutCache.getOrCreate(
      this.cellKey(value, gridCell.bounds.width, gridCell.bounds.height, DEFAULT_OPTS), (v) => drawDuplex(
        null as unknown as CanvasRenderingContext2D, 0, 0, gridCell.bounds.width, gridCell.bounds.height,
        model, DEFAULT_OPTS, true,
      ));
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
