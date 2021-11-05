/**
 * OCL-based molecule cell renderer.
 * */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

export class OCLCellRenderer extends DG.GridCellRenderer {
  molCahe: DG.LruCache;
  static _canvas: HTMLCanvasElement = ui.canvas();
  get name() { return 'OCL cell renderer'; }
  get cellType() { return DG.SEMTYPE.MOLECULE; }
  get defaultWidth() { return 200; }
  get defaultHeight() { return 100; }

  constructor() {
    super();
    this.molCahe = new DG.LruCache();
  }

  static _createMol(molString: string) {
    return molString.endsWith('M  END', 10) ? OCL.Molecule.fromMolfile(molString) : OCL.Molecule.fromSmiles(molString);
  }

  //TODO: sdf
  render(
    g: CanvasRenderingContext2D,
    x: number,
    y: number,
    w: number,
    h: number,
    gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle,
  ) {
    if (w < 5 || h < 5) {
      return;
    }
    const molString: string = gridCell.cell.value;
    let mol: OCL.Molecule;
    try {
      if (molString === null) {
        return;
      }
      mol = this.molCahe.getOrCreate(molString, () => OCLCellRenderer._createMol(molString));
      OCLCellRenderer._canvas.width = w;
      OCLCellRenderer._canvas.height = h;
      OCL.StructureView.drawMolecule(OCLCellRenderer._canvas, mol);
      g.drawImage(OCLCellRenderer._canvas, x, y);
    } catch (exception) {
      const midX = x + w / 2;
      const midY = y + h / 2;
      g.strokeStyle = 'red';
      g.beginPath()
      g.moveTo(midX - 20, midY - 20);
      g.lineTo(midX + 20, midY + 20);
      g.moveTo(midX - 20, midY + 20);
      g.lineTo(midX + 20, midY - 20);
      g.stroke();
    }
  }
}
