import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class ImageCellRenderer extends DG.GridCellRenderer {
  get name() { return 'ImageUrl'; }
  get cellType() { return 'ImageUrl'; }

  get defaultWidth(): number | null { return 200; }
  get defaultHeight(): number | null { return 100; }

  static images: DG.LruCache = new DG.LruCache();

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    g.fillText(`image ${gridCell.cell.value}`, x + 3, y + 3);

    if (w < 5 || h < 5) return;
    const url = gridCell.cell.value;
    if (url == null) return;

    var img = ImageCellRenderer.images.getOrCreate(url, (_: any) => {
      let image = new Image();
      image.src = url;
      image.onload = () => gridCell.grid.invalidate();
      return image;
    });

    if (!img.complete || img.naturalWidth == 0)
      return;

    g.save();
    g.rect(x, y, w, h);
    g.clip();

    // const width: number = (img.naturalWidth * (cellStyle.imageScale ?? 1)).toInt();
    // int height = (img.naturalHeight * (cellStyle.imageScale ?? 1)).toInt();

    g.drawImage(img, x, y);

    // if (img.naturalHeight > 0 && cellStyle.imageScale == null)
    //   g.drawImageToRect(img, bounds.fit(new Size(width, height)));
    // else
    //   g.drawImageToRect(img, new Rectangle(bounds.center.x - width ~/ 2, bounds.center.y - height ~/ 2, width, height));

    // if (cellStyle.opacity != null) {
    //   fillRect(g, bounds, Color.rgba(255, 255, 255, (255 * cellStyle.opacity).toInt()));
    //   if (cellStyle.opacity == 1) {
    //     fillRect(g, bounds, Color.green);
    //     print('opacity');
    //   }
    // }

    g.restore();
  }
}