import {RdKitParallel} from './rdkit_parallel';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

let commonRdKitModule: any = null;

export async function setCommonRdKitModule(module: any)
{
  commonRdKitModule = module;
}

export function drawMoleculeToCanvas(
  x: number, y: number, w: number, h: number,
  onscreenCanvas: HTMLCanvasElement, molString: string, scaffoldMolString: string | null = null) {
  let mol = commonRdKitModule.get_mol(molString);
  const molBlock = mol.get_new_coords(true);
  mol.delete();
  mol = commonRdKitModule.get_mol(molBlock);
  mol.normalize_2d_molblock();
  mol.straighten_2d_layout();
  let scaffoldMol = scaffoldMolString == null ? null :
    commonRdKitModule.get_qmol(scaffoldMolString);
  let substructJson = "{}";
  if (scaffoldMol) {
    substructJson = mol.get_substruct_match(scaffoldMol);
    if (substructJson === "") {
      substructJson = "{}";
    }
  }
  // TODO: make this an optional parameter, AND a system-wide setting
  let opts = {
    "clearBackground": false,
    "offsetx": Math.floor(x),
    "offsety": Math.floor(y),
    "width": Math.floor(w),
    "height": Math.floor(h),
    "bondLineWidth": 1,
    "fixedScale": 0.07,
    "minFontSize": 9,
    "highlightBondWidthMultiplier": 12,
    "dummyIsotopeLabels": false,
    "atomColourPalette": {
      16: [0.498, 0.247, 0.0],
      9: [0.0, 0.498, 0.0],
      17: [0.0, 0.498, 0.0],
    }
  };
  if (scaffoldMol) {
    let substruct = JSON.parse(substructJson);
    Object.assign(opts, substruct);
  }
  // we need the offscreen canvas first to not let the molecule scaffold skew on a real canvas
  let offscreenCanvas: OffscreenCanvas | null = new OffscreenCanvas(w, h);
  mol.draw_to_canvas_with_highlights(offscreenCanvas, JSON.stringify(opts));
  let image = offscreenCanvas!.getContext('2d')!.getImageData(0, 0, w, h);
  let context = onscreenCanvas.getContext('2d');
  context!.putImageData(image, x, y);
  offscreenCanvas = null; // ? GC definitely
}

export function renderDescription(description: OCL.IParameterizedString[]) {
  const host = ui.divV([]);
  const width = 200;
  const height = 150;
  for (const entry of description) {
    if (entry.type == 2 || entry.type == 3) {
      host.append(ui.label(entry.value));
    }
    if (entry.type == 1) {
      const mol = OCL.Molecule.fromIDCode(entry.value);
      host.append(_molToCanvas(mol, width, height));
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
