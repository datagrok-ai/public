import * as DG from 'datagrok-api/dg';
import {renderMolecule} from '../rendering/render-molecule';

export function structure2dWidget(smiles: string): DG.Widget {
  const width = 200;
  const height = 100;
  // const host = ui.canvas(width, height);
  // drawMoleculeToCanvas(0, 0, width, height, host, smiles);
  const host = renderMolecule(smiles, {renderer: 'OpenChemLib', width: width, height: height});
  return new DG.Widget(host);
}
