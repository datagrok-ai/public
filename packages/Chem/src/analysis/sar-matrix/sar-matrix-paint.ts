import * as DG from 'datagrok-api/dg';

import * as TextUtils from '@datagrok-libraries/gridext/src/utils/TextUtils';
import {SarMatrix, SarMatrixCell} from './sar-matrix-types';
import {drawDepiction, rowTemplate} from './sar-matrix-depict';
import {ACTIVITY_SCHEME, CHIP_H, CHIP_MARGIN, CHIP_PAD, cssColor, GRID_FONT, HEADER_ARGB,
  HEADER_W, MatrixGridState, PaneRow, REAL_ALPHA, VIRTUAL_ALPHA, VIRTUAL_ALPHA_MIN, WHITE_ARGB,
  FULL_SUPPORT} from './sar-matrix-ui-common';

/** How a cell corner chip is drawn. The two corners carry different facts and sit on opposite
 *  diagonals so they never collide however wide either gets. */
export interface ChipStyle {
  corner: 'top-left' | 'bottom-right';
  color: string;
  italic?: boolean;
  faint?: boolean;
}

/** What the painters read back from the viewer. Everything else they need arrives as an argument,
 *  because a paint runs per visible cell and must not go looking for state. */
export interface PaintHost {
  readonly root: HTMLElement;
  readonly higherIsBetter: boolean;
  readonly dataFrame: DG.DataFrame;
  formatActivity(value: number): string;
  cellVisible(matrix: SarMatrix, ri: number, ci: number): boolean;
  cellIdText(cell: SarMatrixCell): string | null;
  bestRowValue(paneRow: PaneRow): number | null;
}

/** Canvas painting for a pane grid: header, core column, body cells and their corner chips. Split
 *  from the viewer because these run per visible cell on every scroll and are the one part of the
 *  viewer with no DOM or lifecycle of its own. */
export class MatrixPainter {
  constructor(private readonly host: PaintHost) {}

  /** Potency tint composited over white and returned opaque, for a depiction background: RDKit
   *  anti-aliases bonds against the draw's background, so a translucent tint leaves a pale fringe. */
  private flatTint(matrix: SarMatrix, value: number, alpha: number): number {
    const c = this.tint(matrix, value, alpha);
    const a = alpha / 255;
    const overWhite = (channel: number): number => Math.round(channel * a + 255 * (1 - a));
    return DG.Color.argb(255, overWhite(DG.Color.r(c)), overWhite(DG.Color.g(c)), overWhite(DG.Color.b(c)));
  }
  /** Potency tint at explicit alpha (0-255), green = more potent. `scaleColor`'s own alpha arg is
   *  unusable (transparent mid-range), so RGB is taken and alpha set here; mirrored for lower-is-better. */
  private tint(matrix: SarMatrix, value: number, alpha: number): number {
    // Zero-width range would divide by zero in scaleColor; paint the green end.
    if (matrix.maxActivity === matrix.minActivity)
      return DG.Color.argb(alpha, DG.Color.r(DG.Color.green), DG.Color.g(DG.Color.green), DG.Color.b(DG.Color.green));
    const potency = this.host.higherIsBetter ? value : matrix.minActivity + matrix.maxActivity - value;
    const base = DG.Color.scaleColor(potency, matrix.minActivity, matrix.maxActivity, undefined, ACTIVITY_SCHEME);
    return DG.Color.argb(alpha, DG.Color.r(base), DG.Color.g(base), DG.Color.b(base));
  }
  /** Paint an R-group column header (position band + substituent depiction + sort caption) or the
   *  'Core' column's "Aligned core" label. */
  paintHeader(g: CanvasRenderingContext2D, b: DG.Rect, name: string, state: MatrixGridState): void {
    const grey6 = cssColor(this.host.root, '--grey-6', '#4a4a4a');
    const grey5 = cssColor(this.host.root, '--grey-5', '#7d7d7d');
    if (name === 'Core') {
      // Each header paints its own background (default paint suppressed), else it keeps stale pixels.
      g.fillStyle = DG.Color.toHtml(HEADER_ARGB);
      g.fillRect(b.x, b.y, b.width, b.height);
      const captionH = 14;
      if (state.headerDepiction !== null) {
        drawDepiction(g, b.x, b.y, b.width, Math.max(1, b.height - captionH), state.headerDepiction,
          state.paneTemplate, HEADER_ARGB);
      }
      g.fillStyle = grey5;
      g.font = `italic 11px ${GRID_FONT}`;
      g.textAlign = 'center';
      g.textBaseline = state.headerDepiction !== null ? 'top' : 'middle';
      const captionY = state.headerDepiction !== null ? b.y + b.height - captionH : b.y + b.height / 2;
      g.fillText('Aligned core', b.x + b.width / 2, captionY, b.width - 6);
      return;
    }
    const idx = state.colKeyToIdx.get(name);
    if (idx === undefined)
      return;
    const column = state.columns[idx];
    g.fillStyle = DG.Color.toHtml(HEADER_ARGB);
    g.fillRect(b.x, b.y, b.width, b.height);

    // Position label only on the first column of a group (per-cell clipping can't span columns).
    // 'R' matches the unnumbered mark the core draws at the site these columns attach to; the core's
    // numbered marks are the positions its rows differ at, which no column here enumerates.
    const posBandH = 16;
    if (state.firstOfGroup.has(idx)) {
      g.fillStyle = grey6;
      g.font = `600 11px ${GRID_FONT}`;
      g.textAlign = 'left';
      g.textBaseline = 'top';
      g.fillText('R', b.x + 5, b.y + 2, b.width - 10);
    }

    const caption = column.caption;
    const captionBandH = caption ? 14 : 0;
    const depH = Math.max(1, b.height - posBandH - captionBandH);
    const depW = Math.min(b.width, HEADER_W * 2);
    drawDepiction(g, b.x + (b.width - depW) / 2, b.y + posBandH, depW, depH,
      column.substSmiles, null, HEADER_ARGB);

    if (caption) {
      g.fillStyle = grey5;
      g.font = `10px ${GRID_FONT}`;
      g.textAlign = 'center';
      g.textBaseline = 'top';
      g.fillText(caption, b.x + b.width / 2, b.y + posBandH + depH, b.width);
    }
  }
  /** Paint the pinned-left core cell: the core aligned to its own template + the row label beneath. */
  paintCore(g: CanvasRenderingContext2D, b: DG.Rect, paneRow: PaneRow, rowIdx: number,
    state: MatrixGridState): void {
    const row = paneRow.matrix.rows[paneRow.rowIndex];
    const template = rowTemplate(state, rowIdx);
    const labelH = 16;
    const molH = Math.max(1, b.height - labelH);
    // Fill in the depiction's baked-in colour so anti-aliased bonds don't fringe.
    g.fillStyle = DG.Color.toHtml(WHITE_ARGB);
    g.fillRect(b.x, b.y, b.width, b.height);
    drawDepiction(g, b.x, b.y, b.width, molH, row.keySmiles, template, WHITE_ARGB);
    g.fillStyle = cssColor(this.host.root, '--grey-6', '#4a4a4a');
    g.font = `600 11px ${GRID_FONT}`;
    g.textAlign = 'center';
    g.textBaseline = 'top';
    g.fillText(paneRow.label, b.x + b.width / 2, b.y + molH + 1, b.width);
  }
  /** Paint one core×substituent cell: potency tint, aligned molecule, value/id chips, markers, and
   *  selection ring. Tint uses the row's OWN matrix, so transfer sides colour against their own range. */
  paintBodyCell(g: CanvasRenderingContext2D, b: DG.Rect, paneRow: PaneRow, rowIdx: number,
    colIndex: number, state: MatrixGridState): void {
    const matrix = paneRow.matrix;
    const ri = paneRow.rowIndex;
    const cell = matrix.cells[ri][paneRow.colIdxs[colIndex]];
    // A cell below the threshold is blanked, not removed, so the grid stays a readable table.
    if (cell.kind === 'empty' || cell.value === null || !this.host.cellVisible(matrix, ri, paneRow.colIdxs[colIndex])) {
      g.fillStyle = cssColor(this.host.root, '--white', '#ffffff');
      g.fillRect(b.x, b.y, b.width, b.height);
      return;
    }
    const value = cell.value;
    const isVirtual = cell.kind === 'virtual';
    // A compound the set holds but never measured: the value is predicted like a virtual one, but it
    // already exists. The two are told apart by the outline below — tint and chip follow the number's
    // provenance, the frame follows whether the compound is real.
    const isPredicted = isVirtual || cell.kind === 'unmeasured';
    const support = cell.support ?? 0;
    // Thin predictions read fainter than well-supported ones; real cells are solid.
    const alpha = isPredicted ?
      Math.round(VIRTUAL_ALPHA_MIN + (VIRTUAL_ALPHA - VIRTUAL_ALPHA_MIN) * Math.min(1, support / FULL_SUPPORT)) :
      REAL_ALPHA;
    // Flattened opaque and used as the depiction background so bonds blend into the tint; fill first
    // for cells with no structure.
    const bg = this.flatTint(matrix, value, alpha);
    g.fillStyle = DG.Color.toHtml(bg);
    g.fillRect(b.x, b.y, b.width, b.height);
    if (cell.smiles !== null)
      drawDepiction(g, b.x, b.y, b.width, b.height, cell.smiles, rowTemplate(state, rowIdx), bg);

    // '~value' chip, top-left; green when it's the row's most potent observation.
    const isBest = paneRow.markBest === true && !isPredicted && value === this.host.bestRowValue(paneRow);
    this.paintChip(g, b, `${isPredicted ? '~' : ''}${this.host.formatActivity(value)}`, {
      corner: 'top-left',
      italic: isPredicted,
      faint: isPredicted && support <= 1,
      color: isBest ? cssColor(this.host.root, '--green-2', '#1a8a3a') :
        cssColor(this.host.root, isPredicted ? '--grey-5' : '--grey-6', isPredicted ? '#7d7d7d' : '#4a4a4a'),
    });

    // Id chip, diagonally opposite the potency it belongs to.
    const idText = this.host.cellIdText(cell);
    if (idText !== null)
      this.paintChip(g, b, idText, {corner: 'bottom-right', color: cssColor(this.host.root, '--grey-5', '#7d7d7d')});

    // Cell outline: dashed for a virtual analog, a light solid frame otherwise.
    g.lineWidth = 1;
    if (isVirtual) {
      g.setLineDash([3, 2]);
      g.strokeStyle = cssColor(this.host.root, '--grey-4', '#b5b5b5');
    } else {
      g.setLineDash([]);
      g.strokeStyle = cssColor(this.host.root, '--grey-2', '#e0e0e0');
    }
    g.strokeRect(b.x + 0.5, b.y + 0.5, b.width - 1, b.height - 1);
    g.setLineDash([]);

    // Host-grid link: a selected row's cell gets a blue ring.
    if (cell.molIdx !== null && this.host.dataFrame.selection.get(cell.molIdx)) {
      g.lineWidth = 2;
      g.strokeStyle = cssColor(this.host.root, '--blue-1', '#2083d5');
      g.strokeRect(b.x + 1, b.y + 1, b.width - 2, b.height - 2);
    }
  }
  /** Draw one corner chip: a translucent white plate carrying a line of text (kept legible over a
   *  bond). Text is ellipsized to the cell rather than condensed by `fillText`'s max width. */
  paintChip(g: CanvasRenderingContext2D, b: DG.Rect, text: string, style: ChipStyle): void {
    g.save();
    g.font = `${style.italic === true ? 'italic ' : ''}600 10px ${GRID_FONT}`;
    g.textAlign = 'left';
    g.textBaseline = 'top';
    const trimmed = TextUtils.trimText(text, g, b.width - (CHIP_MARGIN + CHIP_PAD) * 2);
    const chipW = g.measureText(trimmed).width + CHIP_PAD * 2;
    // Measure the far corner from x/y + size: the rect has no right/bottom accessors (would be NaN).
    const atEnd = style.corner === 'bottom-right';
    const x = atEnd ? b.x + b.width - CHIP_MARGIN - chipW : b.x + CHIP_MARGIN;
    const y = atEnd ? b.y + b.height - CHIP_MARGIN - CHIP_H : b.y + CHIP_MARGIN;
    g.globalAlpha = style.faint === true ? 0.65 : 1;
    g.fillStyle = 'rgba(255, 255, 255, 0.82)';
    g.fillRect(x, y, chipW, CHIP_H);
    g.fillStyle = style.color;
    g.fillText(trimmed, x + CHIP_PAD, y + 2);
    g.restore();
  }
}
