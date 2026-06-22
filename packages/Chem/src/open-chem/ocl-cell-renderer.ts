/**
 * OCL-based molecule cell renderer.
 * */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {MAX_SMILES_LENGTH} from '../utils/chem-constants';

function oclMol(mol: string): OCL.Molecule {
  const isMolFile = mol.includes('M  END');
  if (!isMolFile && mol.length > MAX_SMILES_LENGTH)
    throw new Error('Invalid molecule string');
  return mol.includes('M  END') ? OCL.Molecule.fromMolfile(mol) : OCL.Molecule.fromSmiles(mol);
}

export class OCLCellRenderer extends DG.GridCellRenderer {
  molCache: DG.LruCache<string, OCL.Molecule> = new DG.LruCache<string, OCL.Molecule>();
  static _canvas: HTMLCanvasElement = ui.canvas();
  get name() {return 'OpenChemLib';}
  get cellType() {return DG.SEMTYPE.MOLECULE;}
  get defaultWidth() {return 200;}
  get defaultHeight() {return 100;}

  constructor() {
    super();
  }

  render(
    g: CanvasRenderingContext2D,
    x: number,
    y: number,
    w: number,
    h: number,
    gridCell: DG.GridCell,
  ): void {
    if (w < 5 || h < 5)
      return;

    const molString: string = gridCell.cell.value;
    if (molString === null)
      return;

    try {
      const mol = this.molCache.getOrCreate(molString, () => oclMol(molString));
      const r = window.devicePixelRatio;
      const dw = Math.floor(w * r);
      const dh = Math.floor(h * r);
      OCLCellRenderer._canvas.width = dw;
      OCLCellRenderer._canvas.height = dh;
      const offCtx = OCLCellRenderer._canvas.getContext('2d', {willReadFrequently: true})!;
      // @ts-ignore
      OCL.StructureView.drawMolecule(OCLCellRenderer._canvas, mol, {suppressChiralText: true});
      g.putImageData(offCtx.getImageData(0, 0, dw, dh), Math.floor(x * r), Math.floor(y * r));
    } catch (exception) {
      const midX = x + w / 2;
      const midY = y + h / 2;
      g.strokeStyle = 'red';
      g.beginPath();
      g.moveTo(midX - 20, midY - 20);
      g.lineTo(midX + 20, midY + 20);
      g.moveTo(midX - 20, midY + 20);
      g.lineTo(midX + 20, midY - 20);
      g.stroke();
    }
  }
}
