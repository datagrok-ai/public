import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

@grok.decorators.cellRenderer({
  name: 'Scatter Plot',
  cellType: 'scatterplot',
  virtual: true,
})
export class ScatterPlotCellRenderer extends DG.GridCellRenderer {
  get name() : string { return 'scatterplot'; }
  get cellType() : string { return 'scatterplot'; }
  onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void {
    if (gridCell!.tableRowIndex == null)
      return;

    const plot = gridCell!.gridColumn!.column!.get(gridCell!.tableRowIndex);
    const c = gridCell.grid.root.getElementsByClassName('d4-value-editor d4-value-editor-text');
    if (c.length > 0) {
      const e = (c[0] as HTMLElement);
      gridCell.grid.root.appendChild(plot.root);
      plot.root.style.left = e.offsetLeft + 'px';
      plot.root.style.top = e.offsetTop + 'px';
      plot.root.style.width = e.offsetWidth + 'px';
      plot.root.style.height = e.offsetHeight + 'px';

      const cc = plot.root.getElementsByTagName('CANVAS');
      if (cc.length > 0) {
        const canvas = cc[0] as HTMLCanvasElement;
        canvas.style.width = '100%';
        canvas.style.height = '100%';
      }

      const options = {
        rootMargin: '0px',
        threshold: 1.0
      };
      const observer = new IntersectionObserver((entries, observer) => {
        entries.forEach((entry) => {
          if (entry.intersectionRatio <= 0 && plot.root.parentElement != null) {
            plot.root.parentElement.removeChild(plot.root);
            gridCell.grid.invalidate();
          }
          //console.log(entry.intersectionRatio > 0 ? 'visible' : 'invisible');
        });
      }, options);
      observer.observe(e);
    }
  }

  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    if (gridCell!.tableRowIndex == null)
      return;

    const plot = gridCell!.gridColumn!.column!.get(gridCell!.tableRowIndex) as DG.ScatterPlotViewer;
    if (w < 20 || h < 10 || gridCell.grid.dataFrame === void 0) return;

    const c = plot.root.getElementsByTagName('CANVAS');
    let canvas = null;
    if (c.length > 0) {
      canvas = c[0] as HTMLCanvasElement;
      //const r = window.devicePixelRatio;

      document.body.appendChild(plot.root);
      canvas.width = w * window.devicePixelRatio;
      canvas.height = h * window.devicePixelRatio;
      canvas.style.width = w + 'px';
      canvas.style.height = h + 'px';
    }
    g.translate(x, y);
    //version conflict  plot.render(g);
    g.translate(-x, -y);

    if (canvas !== null)
      document.body.removeChild(plot.root);
  }
}
