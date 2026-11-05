import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';


@grok.decorators.cellRenderer({
  name: 'binaryImageCellRenderer',
  cellType: 'BinaryImage',
  //@ts-ignore
  tags: ['cellRenderer'],
})
export class BinaryImageCellRenderer extends DG.GridCellRenderer {
  get name() { return 'BinaryImage'; }

  get cellType() { return 'BinaryImage'; }

  get defaultWidth(): number | null { return 200; }

  get defaultHeight(): number | null { return 100; }

  /** Decoded images by cell value, least recently drawn first, so that a repaint does not decode them again. */
  static images: Map<Uint8Array, HTMLImageElement> = new Map();
  static maxImages = 256;

  static getImage(bytes: Uint8Array): HTMLImageElement {
    const images = BinaryImageCellRenderer.images;
    const img = images.get(bytes) ?? BinaryImageCellRenderer.loadImage(bytes);
    images.delete(bytes);
    images.set(bytes, img);
    if (images.size > BinaryImageCellRenderer.maxImages)
      images.delete(images.keys().next().value!);
    return img;
  }

  static loadImage(bytes: Uint8Array): HTMLImageElement {
    const img = new Image();
    const url = URL.createObjectURL(new Blob([bytes as BlobPart]));
    img.onload = img.onerror = () => URL.revokeObjectURL(url);
    img.src = url;
    return img;
  }

  static drawImage(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, img: HTMLImageElement) {
    if (img.naturalWidth === 0)
      return;
    const r = new DG.Rect(x, y, w, h).inflate(-2, -2).fit(img.naturalWidth, img.naturalHeight);
    g.drawImage(img, r.x, r.y, r.width, r.height);
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle,
  ) {
    const bytes = gridCell.cell.value;
    if (w < 5 || h < 5 || !(bytes instanceof Uint8Array))
      return;
    const img = BinaryImageCellRenderer.getImage(bytes);
    if (img.complete) {
      BinaryImageCellRenderer.drawImage(g, x, y, w, h, img);
      return;
    }
    // hosts that a grid invalidation does not repaint, such as the Form viewer, keep this draw
    img.addEventListener('load', () => {
      BinaryImageCellRenderer.drawImage(g, x, y, w, h, img);
      gridCell.grid?.invalidate();
    }, {once: true});
  }

  onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void {
    // keeps the grid from opening its text editor, which shows the byte count
    e.preventDefault();
    const bytes = gridCell.cell.value;
    const replace = () => DG.Utils.openFileBytes({
      accept: 'image/png, image/jpeg, image/gif, image/webp, image/bmp',
      open: (b) => gridCell.cell.value = b,
    });
    if (!(bytes instanceof Uint8Array)) {
      replace();
      return;
    }
    const img = BinaryImageCellRenderer.loadImage(bytes);
    img.style.maxWidth = '80vw';
    img.style.maxHeight = '70vh';
    const dlg = ui.dialog({title: gridCell.tableColumn?.name ?? 'Image'}).add(img);
    dlg.addButton('Replace...', () => {
      dlg.close();
      replace();
    }).show({resizable: true});
  }
}
