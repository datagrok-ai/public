// OCL requires same kind of loading as we do for RdKit
// Therefore, this cannot be currently used from WebWorkers
import {RdKitService} from './rdkit_service';
import * as ui from 'datagrok-api/ui';
import * as OCL from 'openchemlib/full.js';
import {chemLock, chemUnlock} from './chem_common';

export function renderDescription(description: OCL.IParameterizedString[]) {
  const host = ui.div([], 'd4-flex-wrap');
  const width = 200;
  const height = 150;
  let lastMolCanvas = null;
  for (const entry of description) {
    if (entry.type == 2 || entry.type == 3) {
      host.append(ui.div(
        [ui.label(entry.value), lastMolCanvas],
        lastMolCanvas === null ? {} : {classes: 'd4-flex-col', style: {margin: '5px'}},
      ));
    }
    if (entry.type == 1) {
      const mol = OCL.Molecule.fromIDCode(entry.value);
      lastMolCanvas = _molToCanvas(mol, width, height);
    }
  }
  return host;
}

function _molToCanvas(mol: OCL.Molecule, width=200, height=100) {
  const canvas = ui.canvas(width, height);
  if (mol !== null) {
    OCL.StructureView.drawMolecule(canvas, mol);
  }
  return canvas;
}
