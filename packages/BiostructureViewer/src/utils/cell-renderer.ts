import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {TwinProteinView} from '../viewers/twin-p-cells';
import {initViewer, byData} from '../viewers/molstar-viewer';
import {Rect, TableView} from 'datagrok-api/dg';
import {getNglGlService} from '../package';
import {NglGlTask} from '@datagrok-libraries/bio';

const PDB_RENDERER_IMAGE_CACHE_KEY = 'PdbRendererImageCache';

async function base64ToImg(base64: string): Promise<HTMLImageElement> {
  return new Promise<HTMLImageElement>((resolve, reject) => {
    try {
      const img = new Image();
      img.onload = () => {
        resolve(img);
      };
      img.src = base64;
    } catch (err) {
      reject(err);
    }
  });
}

export class PdbRenderer extends DG.GridCellRenderer {
  get name(): string { return 'xray'; }

  get cellType(): string { return 'xray'; }

  get defaultHeight(): number { return 300 / window.devicePixelRatio; }

  get defaultWidth(): number { return 300 / window.devicePixelRatio; }

  // static imgDG = await base64ToImg('data:image/png;base64,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');

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    // let twin = new TwinProteinView();
    // let ligandSelection: {[key: string]: any} = {};
    // twin.init(gridCell.cell.value, grok.shell.v as DG.TableView, ligandSelection);
    // twin.show(grok.shell.v as DG.TableView);
    const indexStart = gridCell.cell.value.indexOf('TITLE');
    const indexEnd = gridCell.cell.value.indexOf('\n', indexStart + 1);
    const name = gridCell.cell.value.substring(indexStart + 10, indexEnd);
    byData(gridCell.cell.value, name).then((_) => {
      for (const v of grok.shell.views) {
        if (v.name === name) {
          setTimeout(() => {grok.shell.v = v;}, 100);
          break;
        }
      }
    });
  }

  /**
   * Cell renderer function.
   *
   * @param {CanvasRenderingContext2D} g Canvas rendering context.
   * @param {number} x x coordinate on the canvas.
   * @param {number} y y coordinate on the canvas.
   * @param {number} w width of the cell.
   * @param {number} h height of the cell.
   * @param {DG.GridCell} gridCell Grid cell.
   * @param {DG.GridCellStyle} cellStyle Cell style.
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle
  ): void {
    const r = window.devicePixelRatio;
    const service = getNglGlService();

    const rowIdx = gridCell.tableRowIndex;
    console.debug('PTV: PdbRenderer.render() start ' + `rowIdx=${rowIdx}`);

    g.save();
    try {
      if (!gridCell.tableColumn.temp.has(PDB_RENDERER_IMAGE_CACHE_KEY)) {
        gridCell.tableColumn.temp.set(PDB_RENDERER_IMAGE_CACHE_KEY, {});
      }
      const imageCache: { [rowIdx: number]: HTMLImageElement } =
        gridCell.tableColumn.temp.get(PDB_RENDERER_IMAGE_CACHE_KEY);

      const image = rowIdx in imageCache ? imageCache[rowIdx] : null;

      if (!image ||
        Math.abs(image.width / r - gridCell.gridColumn.width) > 0.5 ||
        Math.abs(image.height / r - gridCell.grid.props.rowHeight) > 0.5
      ) {
        // draw image

        // mark rowIdx as in drawing state
        //imageCache[rowIdx] = null;
        //gridCell.tableColumn.temp.set(PDB_RENDERER_IMAGE_CACHE_KEY, imageCache);

        const task: NglGlTask = {
          name: gridCell.cell.rowIndex.toString(),
          backColor: gridCell.grid.props.backColor,
          props: {
            pdb: gridCell.cell.value,
            width: gridCell.gridColumn.width * r,
            height: gridCell.grid.props.rowHeight * r,
          },
          onAfterRender: async (canvas: HTMLCanvasElement) => {
            console.debug('PTV: PdbRenderer.render() onAfterRender() ' + `rowIdx = ${rowIdx}`);
            const imageStr: string = canvas.toDataURL();
            const image: HTMLImageElement = await base64ToImg(imageStr);

            const imageCache = gridCell.tableColumn.temp.get(PDB_RENDERER_IMAGE_CACHE_KEY);
            imageCache[rowIdx] = image;
            gridCell.tableColumn.temp.set(PDB_RENDERER_IMAGE_CACHE_KEY, imageCache);

            service.renderOnGridCell(g, new Rect(x, y, w, h), gridCell, image);
          }
        };

        service.render(task, rowIdx);
      } else {
        const image: HTMLImageElement = imageCache[rowIdx];
        if (image) {
          console.debug('PTV: PdbRenderer.render() ' + `imageCache[${rowIdx}] = <image>`);
          service.renderOnGridCell(g, new Rect(x, y, w, h), gridCell, image);
        } else {
          console.debug('PTV: PdbRenderer.render() ' + `imageCache[${rowIdx}] = null`);
        }
      }
    } finally {
      g.restore();
    }

    return;
  }
}
