/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();


//tags: init
export async function initChem(): Promise<void> {
  return new Promise((resolve, reject) => {
    // @ts-ignore
    dojo.ready(function () { resolve(null); });
  });
}


//name: helmCellRenderer
//tags: cellRenderer,cellRenderer-HELM
//meta.cellType: HELM
//meta-cell-renderer-sem-type: HELM
//output: grid_cell_renderer result
export function helmCellRenderer(): DG.GridCellRenderer {
  return new HelmCellRenderer();
}


class HelmCellRenderer extends DG.GridCellRenderer {

  get name() { return 'helm'; }
  get cellType() { return 'helm'; }
  get defaultWidth(): number | null { return 400; }
  get defaultHeight(): number | null { return 100; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    let host = ui.div([], { style: { width: `${w}px`, height: `${h}px`}});
    host.setAttribute('dataformat', 'helm');
    host.setAttribute('data', gridCell.cell.value);

    gridCell.element = host;
    setTimeout(function() {
      // @ts-ignore
      new JSDraw(host, { width: w, height: h, skin: "w8", viewonly: true });
    }, 200);
  }
}