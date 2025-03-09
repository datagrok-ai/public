import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';


const MAX_IMG_PX_WIDTH = 600;
const MAX_IMG_PX_HEIGHT = 600;
const DLG_WINDOW_PADDING_PX = 24;

@grok.decorators.cellRenderer({
  name: 'imageUrlCellRenderer',
  cellType: 'ImageUrl',
})
export class ImageCellRenderer extends DG.GridCellRenderer {
  get name() { return 'ImageUrl'; }

  get cellType() { return 'ImageUrl'; }

  get defaultWidth(): number | null { return 200; }

  get defaultHeight(): number | null { return 100; }

  static images: DG.LruCache = new DG.LruCache();

  async getImage(gridCell: DG.GridCell): Promise<HTMLImageElement | undefined> {
    const url: string = gridCell.cell.value;
    if (url == null)
      return;
    const image = new Image();
    image.onload = () => gridCell.grid.invalidate();
    image.crossOrigin = 'anonymous';
    if (url.toLowerCase().startsWith('http'))
      image.src = url;
    else if (['.png', '.jpg', '.jpeg'].some((end) => url.endsWith(end)) && await grok.dapi.files.exists(url)) {
      const endSlashIdx = url.lastIndexOf('/');
      const folderPath = url.slice(0, endSlashIdx + 1);
      const fileList = await grok.dapi.files.list(folderPath);
      const file = fileList.find((file) => file.nqName === url);
      if (file)
        image.src = file.url;
    }
    return image;
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    if (w < 5 || h < 5) return;
    const url = gridCell.cell.value;
    if (url == null) return;

    const img: Promise<HTMLImageElement> = ImageCellRenderer.images.getOrCreate(url, () => this.getImage(gridCell));
    img.then((img: HTMLImageElement) => {
      if (!img || !img.complete || img.naturalWidth == 0)
        return;

      // g.save();
      g.rect(x, y, w, h);
      const bounds = new DG.Rect(x, y, w, h).fit(img.width, img.height);

      // const width: number = (img.naturalWidth * (cellStyle.imageScale ?? 1)).toInt();
      // int height = (img.naturalHeight * (cellStyle.imageScale ?? 1)).toInt();

      g.drawImage(img, bounds.x, bounds.y, bounds.width, bounds.height);

      // if (img.naturalHeight > 0 && cellStyle.imageScale == null)
      //   g.drawImageToRect(img, bounds.fit(new Size(width, height)));
      // else
      //   g.drawImageToRect(
      //     img, new Rectangle(bounds.center.x - width ~/ 2, bounds.center.y - height ~/ 2, width, height));

      // if (cellStyle.opacity != null) {
      //   fillRect(g, bounds, Color.rgba(255, 255, 255, (255 * cellStyle.opacity).toInt()));
      //   if (cellStyle.opacity == 1) {
      //     fillRect(g, bounds, Color.green);
      //     print('opacity');
      //   }
      // }

      // g.restore();
    });
  }

  onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.cell.value)
      return;

    this.getImage(gridCell).then((img) => {
      if (!img)
        return;

      const srcWidth = img.width;
      const srcHeight = img.height;

      const dlg = ui.dialog({title: 'Image', showFooter: false})
        .add(img)
        .show({resizable: true});

      const dlgContentsBox = dlg.root.getElementsByClassName('d4-dialog-contents dlg-image')[0] as HTMLElement;
      const nonImgHeight = dlg.root.getElementsByClassName('d4-dialog-header')[0].clientHeight +
        dlg.root.getElementsByClassName('d4-dialog-footer')[0].clientHeight;
      const aspectRatio = Math.min((MAX_IMG_PX_WIDTH - DLG_WINDOW_PADDING_PX) / srcWidth,
        (MAX_IMG_PX_HEIGHT - nonImgHeight - DLG_WINDOW_PADDING_PX) / srcHeight);

      dlgContentsBox.style.width = `${srcWidth * aspectRatio + DLG_WINDOW_PADDING_PX}px`;
      dlgContentsBox.style.height = `${srcHeight * aspectRatio + DLG_WINDOW_PADDING_PX}px`;
    });
  }
}
