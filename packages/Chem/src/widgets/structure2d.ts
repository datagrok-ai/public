import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {drawMoleculeToCanvas} from '../chem_common_rdkit';

export function structure2dWidget(smiles: string) {
  const width = 200;
  const height = 100;
  const host = ui.canvas(width, height);
  drawMoleculeToCanvas(0, 0, width, height, host, smiles);
  return new DG.Widget(host);
}
