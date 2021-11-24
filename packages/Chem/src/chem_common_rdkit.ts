// This file will be used from Web Workers, so there
// should be no imports from Datagrok or OCL
import {RdKitService} from './rdkit_service';
import {chemLock, chemUnlock} from './chem_common';

let commonRdKitModule: any = null;

export function setCommonRdKitModule(module: any) {
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