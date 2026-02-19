// OCL requires same kind of loading as we do for RdKit
// Therefore, this cannot be currently used from WebWorkers
import * as ui from 'datagrok-api/ui';
import * as OCL from 'openchemlib/full';
import * as DG from 'datagrok-api/dg';

export function renderDescription(description: OCL.IParameterizedString[]): HTMLDivElement {
  const host = ui.div([]);
  const molsHost = ui.div([], 'd4-flex-wrap');
  const width = 200;
  const height = 150;

  let lastMolCanvas: null | HTMLCanvasElement = null;
  let scaffoldMolString: null | string = null;
  description.forEach((entry, idx) => {
    if (entry.type === 2 && ![0, 1].includes(idx) || entry.type === 3) {
      const divElement = ui.div(
        [ui.label(entry.value), lastMolCanvas], lastMolCanvas ? {classes: 'd4-flex-col', style: {margin: '5px'}} : null,
      );
      (lastMolCanvas ? molsHost : host).append(divElement);
      lastMolCanvas = null;
      scaffoldMolString = null;
    }
    if (entry.type === 1) {
      scaffoldMolString = entry.value;
      const mol = OCL.Molecule.fromIDCode(scaffoldMolString);
      lastMolCanvas = _molToCanvas(mol, width, height);
    }
  });
  host.append(molsHost);
  return host;
}

function _molToCanvas(mol: OCL.Molecule, width = 200, height = 100): HTMLCanvasElement {
  const canvas = ui.canvas(width, height);
  canvas!.style.width = width.toString() + 'px';
  canvas!.style.height = height.toString() + 'px';

  if (mol !== null)
    OCL.StructureView.drawMolecule(canvas, mol);

  return canvas;
}

export function oclMol(molStr: string): OCL.Molecule {
  const isMolFile = DG.chem.isMolBlock(molStr);
  if (!isMolFile && molStr.length > 5000)
    throw new Error('Invalid molecule string');
  return isMolFile ? OCL.Molecule.fromMolfile(molStr) : OCL.Molecule.fromSmiles(molStr);
}
