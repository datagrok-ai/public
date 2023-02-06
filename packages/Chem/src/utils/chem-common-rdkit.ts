// This file will be used from Web Workers
// There should be no imports from Datagrok or OCL
import {RdKitService} from '../rdkit-service/rdkit-service';
import {convertToRDKit} from '../analysis/r-group-analysis';
//@ts-ignore
import rdKitLibVersion from '../rdkit_lib_version';
//@ts-ignore
import initRDKitModule from '../RDKit_minimal.js';
import { isMolBlock } from './chem-common';
import $ from 'cash-dom';
import {RDModule, RDMol, Reaction} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export let _rdKitModule: RDModule;
export let _rdKitService: RdKitService;
export let _webRoot: string | null;
export let moduleInitialized = false;

const RDKIT_COMMON_RENDER_OPTS: {[key: string]: any} = {
  clearBackground: false,
  offsetx: 0,
  offsety: 0,
  bondLineWidth: 0.7,
  multipleBondOffset: 0.25,
  fixedScale: 0.07,
  baseFontSize: 1.0,
  minFontSize: -1,
  maxFontSize: -1,
  annotationFontScale: 0.7,
  highlightBondWidthMultiplier: 12,
  dummyIsotopeLabels: false,
  atomColourPalette: {
    0: [0.1, 0.1, 0.1],
    1: [0.0, 0.0, 0.0],
    6: [0.0, 0.0, 0.0],
    7: [0.0, 0.0, 1.0],
    8: [1.0, 0.0, 0.0],
    9: [0.0, 0.498, 0.0],
    15: [0.498, 0.0, 0.498],
    16: [0.498, 0.247, 0.0],
    17: [0.0, 0.498, 0.0],
    35: [0.0, 0.498, 0.0],
    53: [0.247, 0.0, 0.498],
  },
  backgroundColour: [1, 1, 1, 1],
};

export function setRdKitWebRoot(webRootValue: string): void {
  _webRoot = webRootValue;
}

export async function initRdKitModuleLocal(): Promise<void> {
  _rdKitModule = await initRDKitModule(
    {locateFile: () => `${_webRoot}/dist/${rdKitLibVersion}.wasm`});
  if (!_rdKitModule)
    throw 'RdKit Module is not loaded';
  _rdKitModule.prefer_coordgen(false);
  console.log('RDKit module package instance was initialized');
  moduleInitialized = true;
  _rdKitService = new RdKitService();
}

export async function initRdKitService(): Promise<void> {
  await _rdKitService!.init(_webRoot!);
}

export function getRdKitModule(): RDModule {
  if (!moduleInitialized)
    throw ('RdKit Module is not initialized');
  return _rdKitModule!;
}

export async function getRdKitService(): Promise<RdKitService> {
  await initRdKitService();
  if (!_rdKitService)
    throw 'RdKit Service isn\'t initialized';
  return _rdKitService;
}

export function getRdKitWebRoot() {
  return _webRoot;
}

export function drawErrorCross(ctx: OffscreenCanvasRenderingContext2D, width: number, height: number) {
  ctx.lineWidth = 1;
      ctx.strokeStyle = '#EFEFEF';
      ctx.beginPath();
      ctx.moveTo(0, 0);
      ctx.lineTo(width, height);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(width, 0);
      ctx.lineTo(0, height);
      ctx.stroke();
}

function createRenderingOpts(addSettings: {[key: string]: any}): {[key: string]: any} {
  const opts: {[key: string]: any} = {};
  Object.keys(RDKIT_COMMON_RENDER_OPTS).forEach((key: string) => opts[key] = RDKIT_COMMON_RENDER_OPTS[key]);
  Object.keys(addSettings).forEach((key: string) => opts[key] = addSettings[key]);
  return opts;
}

export function drawRdKitMoleculeToOffscreenCanvas(
  rdKitMol: RDMol, w: number, h: number, offscreenCanvas: OffscreenCanvas, substruct: Object | null) {
  const opts = createRenderingOpts({width: Math.floor(w), height: Math.floor(h)});
  const g = offscreenCanvas.getContext('2d', {willReadFrequently : true});
  g?.clearRect(0,0, w, h);
  if (substruct)
    Object.assign(opts, substruct);
  rdKitMol.draw_to_canvas_with_highlights((offscreenCanvas as unknown) as HTMLCanvasElement, JSON.stringify(opts));
  // we need the offscreen canvas first to not let the molecule scaffold skew on a real canvas
}

export function drawRdKitReactionToOffscreenCanvas(
  rdKitReaction: Reaction, w: number, h: number, offscreenCanvas: OffscreenCanvas) {
  const opts = createRenderingOpts({width: Math.floor(w), height: Math.floor(h)});
  const g = offscreenCanvas.getContext('2d', {willReadFrequently : true});
  g?.clearRect(0,0, w, h);
  rdKitReaction.draw_to_canvas_with_highlights((offscreenCanvas as unknown) as HTMLCanvasElement, JSON.stringify(opts));
}

export function drawMoleculeToCanvas(
  x: number, y: number, w: number, h: number,
  onscreenCanvas: HTMLCanvasElement, molString: string, scaffoldMolString: string | null = null,
  options = {normalizeDepiction: true, straightenDepiction: true}
) {

  $(onscreenCanvas).addClass('chem-canvas');
  const r = window.devicePixelRatio;
  onscreenCanvas.width = w * r;
  onscreenCanvas.height = h * r;
  onscreenCanvas.style.width = (w).toString() + 'px';
  onscreenCanvas.style.height = (h).toString() + 'px';

  let mol = null;
  try {
    const isMol: boolean = isMolBlock(molString);
    mol = isMol ? getRdKitModule().get_mol(molString) : getRdKitModule().get_mol(convertToRDKit(molString)!);

    if (!isMol)
      mol.set_new_coords(true);

    if (options.normalizeDepiction ?? true)
      !isMol ? mol.normalize_depiction(1) : mol.normalize_depiction(0);

    if (options.straightenDepiction ?? true)
      !isMol ? mol.straighten_depiction(false) : mol.straighten_depiction(true); 

    const scaffoldMol = scaffoldMolString == null ? null :
      (isMolBlock(scaffoldMolString) ? getRdKitModule().get_qmol(scaffoldMolString) : getRdKitModule().get_qmol(convertToRDKit(scaffoldMolString)!));
    let substructJson = '{}';
    if (scaffoldMol) {
      substructJson = mol.get_substruct_match(scaffoldMol);
      if (substructJson === '')
        substructJson = '{}';
    }
    const substruct = JSON.parse(substructJson);
    const offscreenCanvas = new OffscreenCanvas(w, h);
    drawRdKitMoleculeToOffscreenCanvas(mol, w, h, offscreenCanvas, substruct);
    const image = offscreenCanvas!.getContext('2d')!.getImageData(0, 0, w, h);
    const context = onscreenCanvas.getContext('2d')!;
    context.putImageData(image, x, y);
  } finally {
    mol?.delete();
  }
}

export function checkMolEqualSmiles(mol1: any, molfile2: string): boolean {
  const mol2 = checkMoleculeValid(molfile2);
  const result = mol2 ? mol1.get_smiles() === mol2.get_smiles() : false;
  mol2?.delete();
  return result;
}

export function checkMoleculeValid(molecule: string): any {
  let mol;
  try {
    mol = getRdKitModule().get_mol(molecule);
  }
  catch (e: any) {
    mol?.delete();
    return null;
  }
  return mol;
}