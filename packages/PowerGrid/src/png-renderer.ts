/* eslint-disable max-len */

/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

@grok.decorators.cellRenderer({name: 'rawPng', cellType: 'rawPng',
  // @ts-ignore
  tags: ['cellRenderer']})
export class RawPNGRenderer extends DG.GridCellRenderer {
  get name(): string { return 'rawPng'; }

  get cellType(): string { return 'rawPng'; }

  get defaultHeight(): number { return 200; }

  get defaultWidth(): number { return 200; }

  /**
     * Cell renderer function.
     *
     * @param {CanvasRenderingContext2D} g Canvas rendering context.
     * @param {number} x x coordinate on the canvas.
     * @param {number} y y coordinate on the canvas.
     * @param {number} w width of the cell.
     * @param {number} h height of the cell.
     * @param {DG.GridCell} gridCell Grid cell.
     * @param {DG.GridCellStyle} _cellStyle Cell style.
     */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    _cellStyle: DG.GridCellStyle
  ): void {
    if (gridCell.gridRow < 0 || !gridCell.cell.value) return;

    const cacheKey = this.getCacheKey(gridCell);
    // console.log('Cache key:', cacheKey);
    if (cacheKey && this.imageCache.has(cacheKey)) {
      const img = this.imageCache.get(cacheKey);
      if (img) // can be null if image failed to load or is loading currently
        this.drawInBounds(g, x, y, w, h, gridCell, img);
      return;
    }
    const pngString: string = gridCell.cell.value;
    const img = new Image();

    img.src = 'data:image/png;base64,' + pngString;
    if (cacheKey)
      this.imageCache.set(cacheKey, null); // mark as loading
    img.onload = () => {
      if (cacheKey) {
        this.imageCache.set(cacheKey, img);
        gridCell.render();
      } else { this.drawInBounds(g, x, y, w, h, gridCell, img); }
    };
  }

  getCacheKey(gridCell: DG.GridCell): string | null {
    if (!gridCell?.cell?.value || !gridCell.cell.dataFrame || !gridCell.cell.dataFrame.id ||
      !gridCell.cell.column?.name || (gridCell.cell.rowIndex ?? -1) < 0 || gridCell.cell.column.version == null)
      return null;
    return `${gridCell.cell.dataFrame.id}-${gridCell.cell.column.name}-${gridCell.cell.rowIndex}-${gridCell.cell.column.version}`;
  }

  imageCache: DG.LruCache<string, HTMLImageElement | null> = new DG.LruCache<string, HTMLImageElement | null>(100);

  drawInBounds(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    img: HTMLImageElement): void {
    g.clearRect(x, y, w, h);
    const minAxis = Math.min(w, h, (img.width ?? 0), (img.height ?? 0));
    if (minAxis < 4)
      return;
    // image should be scaled such that it is in the center of the cell
    const scaleBy = w < h ? (w - 4) / img.width : (h - 4) / img.height;
    const imgWidth = img.width * scaleBy;
    const imgHeight = img.height * scaleBy;
    const xOffset = (w - imgWidth) / 2;
    const yOffset = (h - imgHeight) / 2;
    g.drawImage(img, Math.floor(x + xOffset), Math.floor(y + yOffset), Math.floor(imgWidth), Math.floor(imgHeight));
  }

  onDoubleClick(gridCell: DG.GridCell<any>, e: MouseEvent): void {
    if (gridCell.cell.value) {
      const pngString: string = gridCell.cell.value;
      const img = new Image();
      img.src = 'data:image/png;base64,' + pngString;
      const canvas = ui.canvas(600, 600);
      canvas.width = 600;
      canvas.height = 600;

      const ctx = canvas.getContext('2d')!;
      img.onload = () => this.drawInBounds(ctx, 0, 0, 600, 600, gridCell, img);
      ui.dialog({title: 'PNG Image'})
        .add(canvas)
        .show({resizable: true});
      e.preventDefault();
      e.stopPropagation();
      setTimeout(() => {
        canvas.style.width = '100%';
        canvas.style.height = '100%';
      }, 300);
      // d.sub(ui.onSizeChanged(canvas).subscribe(() => console.log(canvas.offsetHeight, canvas.offsetWidth)));
    }
  }
}
