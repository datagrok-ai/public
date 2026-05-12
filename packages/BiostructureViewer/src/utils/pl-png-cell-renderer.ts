import * as DG from 'datagrok-api/dg';

/**
 * Cell renderer for `PL Diagram` column values — base64 PNGs of the
 * captured LigNetwork canvas, produced by `_captureOne` in
 * `prolif-panel.ts`.
 *
 * Functionally equivalent to PowerGrid's `RawPNGRenderer`, but
 * registered under the semType `PL.LigNetwork` — the same semType the
 * context panel `plDiagramInteractionsWidget` binds to. We bundle this
 * renderer locally rather than rely on PowerGrid's `rawPng` for two
 * reasons:
 *
 *   1. Datagrok prioritises semType-based renderer lookup over the
 *      `cell.renderer` column tag. Setting `semType = 'PL.LigNetwork'`
 *      (needed for the panel) makes the platform look for a renderer
 *      registered for `PL.LigNetwork` — if there isn't one (as with
 *      PowerGrid's `rawPng`, which is keyed on its own semType), the
 *      column falls back to default text rendering and the tag is
 *      effectively ignored. Bundling our own renderer here closes that
 *      loop so the column finds a renderer immediately.
 *   2. Removes the PowerGrid dependency for the PL panel feature — BSV
 *      can render `PL Diagram` cells whether or not PowerGrid is
 *      installed on the Datagrok server.
 *
 * Async load: `img.onload` triggers a grid invalidate so that freshly
 * captured PNGs paint as soon as decoding finishes (without this, the
 * cell would stay blank until the user scrolls or otherwise forces a
 * grid repaint).
 */
export class PlLigNetworkPngRenderer extends DG.GridCellRenderer {
  get name(): string { return 'plLigNetworkPng'; }

  get cellType(): string { return 'PL.LigNetwork'; }

  get defaultHeight(): number { return 130; }

  get defaultWidth(): number { return 220; }

  /** Datagrok's true LRU (proper access-order eviction, not just insertion-
   *  order). Maps a stable `(dfId, colName, rowIndex, colVersion)` key to
   *  the decoded HTMLImageElement so subsequent renders are synchronous.
   *  Capacity is 200 — at the configured row height (130px) the visible
   *  grid is rarely more than ~40 rows × a few diagram columns, so 200
   *  comfortably covers the visible viewport plus scroll-back. */
  private readonly _cache: DG.LruCache<string, HTMLImageElement> =
    new DG.LruCache<string, HTMLImageElement>(200);

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle): void {
    if (gridCell.gridRow < 0 || !gridCell.cell.value) return;
    const key = this._cacheKey(gridCell);
    if (key != null) {
      const cached = this._cache.get(key);
      if (cached != null) { this._drawImg(g, x, y, w, h, cached); return; }
    }
    const pngString = String(gridCell.cell.value);
    if (pngString.length < 100) return; // too short to be a real PNG
    const img = new Image();
    img.onload = () => {
      if (key != null) this._cache.set(key, img);
      // Force a grid redraw so the just-decoded image gets drawn. Walk
      // via `gridColumn.grid` (typed on `DG.GridColumn`) rather than the
      // untyped `gridCell.grid` shortcut, which exists at runtime but
      // isn't on the `DG.GridCell` type.
      try { gridCell.gridColumn?.grid?.invalidate?.(); } catch (_e) { /* */ }
    };
    img.onerror = () => { /* malformed PNG — cell stays blank */ };
    img.src = 'data:image/png;base64,' + pngString;
  }

  private _cacheKey(gridCell: DG.GridCell): string | null {
    const c = gridCell.cell;
    if (!c?.column?.name || (c.rowIndex ?? -1) < 0) return null;
    // `dataFrame.id` is null until the DF has been added to dapi; fall
    // back to `dataFrame.name` (stable across re-runs of the batch).
    const dfId = c.dataFrame?.id ?? c.dataFrame?.name ?? '_';
    const ver = c.column.version ?? 0;
    return `${dfId}|${c.column.name}|${c.rowIndex}|${ver}`;
  }

  private _drawImg(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    img: HTMLImageElement): void {
    g.clearRect(x, y, w, h);
    if (img.width < 1 || img.height < 1 || w < 4 || h < 4) return;
    // Fit while preserving aspect ratio, with a 2px inset so neighbour
    // grid lines stay visible.
    const scale = Math.min((w - 4) / img.width, (h - 4) / img.height);
    const dw = img.width * scale;
    const dh = img.height * scale;
    g.drawImage(img,
      Math.floor(x + (w - dw) / 2),
      Math.floor(y + (h - dh) / 2),
      Math.floor(dw),
      Math.floor(dh),
    );
  }
}
