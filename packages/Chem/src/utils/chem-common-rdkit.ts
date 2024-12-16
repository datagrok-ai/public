// This file will be used from Web Workers
// There should be no imports from Datagrok or OCL
import {RdKitService} from '../rdkit-service/rdkit-service';
import {convertToRDKit} from '../analysis/r-group-analysis';
//@ts-ignore
import rdKitLibVersion from '../rdkit_lib_version';
//@ts-ignore
import initRDKitModule from '../RDKit_minimal.js';
import {hexToPercentRgb, isMolBlock} from './chem-common';
import $ from 'cash-dom';
import {RDModule, RDMol, RDReaction} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';
import {IMolContext, getMolSafe} from './mol-creation_rdkit';

export let _rdKitModule: RDModule;
export let _rdKitService: RdKitService;
export let _webRoot: string | null;
export let moduleInitialized = false;

export const RDKIT_COMMON_RENDER_OPTS: {[key: string]: any} = {
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
  highlightBondWidthMultiplier: 20,
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
  highlightColour: [1.0, 0.7, 0.7, 1.0],
  backgroundColour: [1, 1, 1, 1],
  queryColour: [0, 0, 0, 1],
};

export function setRdKitWebRoot(webRootValue: string): void {
  _webRoot = webRootValue;
}

export async function initRdKitModuleLocal(): Promise<void> {
  _rdKitModule = await initRDKitModule(
    {locateFile: () => `${_webRoot}/dist/${rdKitLibVersion}.wasm`});
  if (!_rdKitModule)
    throw new Error('RdKit Module is not loaded');
  _rdKitModule.prefer_coordgen(false);
  _rdKitModule.use_legacy_stereo_perception(false);
  console.log('RDKit module package instance was initialized');
  moduleInitialized = true;
  _rdKitService = new RdKitService();
}

export async function initRdKitService(): Promise<void> {
  await _rdKitService!.init(_webRoot!);
}

export function getRdKitModule(): RDModule {
  if (!moduleInitialized)
    throw new Error('RdKit Module is not initialized');
  return _rdKitModule!;
}

export async function getRdKitService(): Promise<RdKitService> {
  await initRdKitService();
  if (!_rdKitService)
    throw new Error('RdKit Service isn\'t initialized');
  return _rdKitService;
}

export function getRdKitWebRoot() {
  return _webRoot;
}

export function drawErrorCross(ctx: OffscreenCanvasRenderingContext2D, width: number, height: number) {
  ctx.lineWidth = 1;
  ctx.strokeStyle = '#DBDCDF';
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
  molCtx: IMolContext, w: number, h: number, offscreenCanvas: OffscreenCanvas, substruct: Object | null) {
  const g = offscreenCanvas.getContext('2d', {willReadFrequently: true}) as OffscreenCanvasRenderingContext2D;
  const rdKitMol: RDMol | null = molCtx.mol;
  if (rdKitMol === null) {
    console.error('Molecule to be rendered cannot be null.');
    drawErrorCross(g, w, h);
    return;
  }


  w = Math.floor(w);
  h = Math.floor(h);

  if (!w || !h) {
    console.warn('Molecule not rendered due to width or height less than 1');
    return;
  }

  const opts = createRenderingOpts({width: w, height: h});

  g?.clearRect(0, 0, w, h);
  if (substruct)
    Object.assign(opts, substruct);
  const kekulize = molCtx.kekulize;
  if (!kekulize)
    Object.assign(opts, {kekulize});

  const useMolBlockWedging = molCtx.useMolBlockWedging;
  const wedgeBonds = false;
  const addChiralHs = false;
  if (useMolBlockWedging)
    Object.assign(opts, {useMolBlockWedging, wedgeBonds, addChiralHs});

  try {
    rdKitMol.draw_to_canvas_with_highlights((offscreenCanvas as unknown) as HTMLCanvasElement, JSON.stringify(opts));
  } catch (e) {
    console.error('Molecule failed to render ' + rdKitMol.get_molblock());
    drawErrorCross(g, w, h);
    return;
  }
  // we need the offscreen canvas first to not let the molecule scaffold skew on a real canvas
}

export function drawRdKitReactionToOffscreenCanvas(
  rdKitReaction: RDReaction, w: number, h: number, offscreenCanvas: OffscreenCanvas) {
  const opts = createRenderingOpts({width: Math.floor(w), height: Math.floor(h)});
  const g = offscreenCanvas.getContext('2d', {willReadFrequently: true});
  g?.clearRect(0, 0, w, h);
  rdKitReaction.draw_to_canvas_with_highlights((offscreenCanvas as unknown) as HTMLCanvasElement, JSON.stringify(opts));
}

export function drawMoleculeToCanvas(x: number, y: number, w: number, h: number,
  onscreenCanvas: HTMLCanvasElement, molString: string, scaffoldMolString: string | null = null,
  options = {normalizeDepiction: true, straightenDepiction: true}, scaffoldStruct: ISubstruct | null = null) {
  if (!w || !h) {
    console.error('Width and height cannot be zero.');
    return;
  }

  $(onscreenCanvas).addClass('chem-canvas');
  const r = window.devicePixelRatio;

  const nW = w * r;
  const nH = h * r;

  onscreenCanvas.width = nW;// w * r;
  onscreenCanvas.height = nH;// h * r;
  onscreenCanvas.style.width = (w).toString() + 'px';
  onscreenCanvas.style.height = (h).toString() + 'px';

  const isMol: boolean = isMolBlock(molString);
  if (!isMol)
    molString = convertToRDKit(molString);

  const offscreenCanvas = new OffscreenCanvas(nW, nH);
  const molCtx = getMolSafe(molString, {}, getRdKitModule());
  const mol : RDMol | null = molCtx.mol;
  if (mol === null) {
    console.error('Molecule is null for ' + molString);
    drawErrorCross(offscreenCanvas.getContext('2d') as OffscreenCanvasRenderingContext2D, w, h);
    return;
  }

  if (!isMol)
    mol.set_new_coords(true);

  const molHasOwnCoords = (mol.has_coords() > 0);
  if (molHasOwnCoords) {
    mol.normalize_depiction(0);
    mol.straighten_depiction(true);
  }

  let scaffoldMol: RDMol | null = null;
  try {
    let substruct: ISubstruct | null = null;
    if (scaffoldStruct)
      substruct = scaffoldStruct;
    else {
      let substructJson = scaffoldStruct ?? '{}';
      scaffoldMol = scaffoldMolString == null ? null :
        (isMolBlock(scaffoldMolString) ? getRdKitModule().get_qmol(scaffoldMolString) :
          getRdKitModule().get_qmol(convertToRDKit(scaffoldMolString)!));
      if (scaffoldMol) {
        try {
          substructJson = mol.get_substruct_match(scaffoldMol);
        } catch (e) {
          console.error(`get_substruct_match failed for ${molString} and ${scaffoldMolString}`);
        }
        if (substructJson === '')
          substructJson = '{}';
      }
      substruct = JSON.parse(substructJson);
    }
    drawRdKitMoleculeToOffscreenCanvas(molCtx, nW, nH, offscreenCanvas, substruct);
    const image = offscreenCanvas!.getContext('2d')!.getImageData(0, 0, nW, nH);
    const context = onscreenCanvas.getContext('2d')!;
    context.putImageData(image, x, y);
  } finally {
    mol?.delete();
    scaffoldMol?.delete();
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
  } catch (e: any) {
    mol?.delete();
    return null;
  }
  return mol;
}