import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';


@grok.decorators.cellRenderer({
  name: 'svgCellRenderer',
  cellType: 'SVG',
  //@ts-ignore
  tags: ['cellRenderer'],
})
export class SvgCellRenderer extends DG.GridCellRenderer {
  get name() { return 'SVG'; }

  get cellType() { return 'SVG'; }

  get defaultWidth(): number { return 100; }

  get defaultHeight(): number { return 100; }

  imageCache: DG.LruCache<string, HTMLImageElement | null> = new DG.LruCache<string, HTMLImageElement | null>(200);

  private getCacheKey(svgString: string): string {
    const len = svgString.length;
    if (len < 200)
      return svgString;
    return `${len}:${svgString.substring(0, 100)}:${svgString.substring(len - 50)}`;
  }

  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle
  ): void {
    if (gridCell.gridRow < 0 || !gridCell.cell.value || w < 5 || h < 5)
      return;

    const svgString: string = gridCell.cell.value;
    const cacheKey = this.getCacheKey(svgString);

    if (this.imageCache.has(cacheKey)) {
      const img = this.imageCache.get(cacheKey);
      if (img)
        this.drawImage(g, x, y, w, h, img);
      return;
    }

    this.imageCache.set(cacheKey, null);
    const img = new Image();
    img.onload = () => {
      this.imageCache.set(cacheKey, img);
      gridCell.render();
    };
    img.onerror = () => {};
    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgString);
  }

  private drawImage(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    img: HTMLImageElement): void {
    const iw = img.naturalWidth || w;
    const ih = img.naturalHeight || h;
    const bounds = new DG.Rect(x, y, w, h).fit(iw, ih);
    g.drawImage(img, bounds.x, bounds.y, bounds.width, bounds.height);
  }

  onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.cell.value)
      return;
    const svgString: string = gridCell.cell.value;
    const img = new Image();
    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgString);
    img.style.maxWidth = '100%';
    img.style.maxHeight = '100%';
    ui.dialog({title: 'SVG Image'})
      .add(img)
      .show({resizable: true});
    e.preventDefault();
    e.stopPropagation();
  }
}
