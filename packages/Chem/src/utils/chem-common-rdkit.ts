// This file will be used from Web Workers
// There should be no imports from Datagrok or OCL
// RDKit-only (it goes with RDKit.js, engine/README.md "Removing RDKit"): RDKit's main-thread module and pool
// (`initRdKitModuleLocal`, `initRdKitService`, `getRdKitModule`, `getRdKitService`, `cancelChemOp`),
// `checkMoleculeValid`, `validationMessageByRdKit` and reaction drawing (`drawRdKitReactionToOffscreenCanvas`). The
// rest takes the module or branches by engine, and stays.
import {RdKitService} from '../rdkit-service/rdkit-service';
import {convertToRDKit} from '../analysis/r-group-analysis';
//@ts-ignore
import rdKitLibVersion from '../rdkit_lib_version';
//@ts-ignore
import initRDKitModule from '../RDKit_minimal.js';
import {hasNewLines, isChemOpRunning, isMolBlock} from './chem-common';
import {MAX_SMILES_LENGTH} from './chem-constants';
import $ from 'cash-dom';
import {RDLog, RDModule, RDMol, RDReaction} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';
import {IMolContext, getMolSafe} from './mol-creation_rdkit';
import {initRdKitFrom} from '../rdkit-service/rdkit-service-worker-api';
import {byEngineModule, withMol} from '../engine/chem-engine';
import {noteCruxError} from '../engine/crux-module';
import {cruxReadsCxExtension} from '../crux/crux-service';

export let _rdKitModule: RDModule;
export let _rdKitService: RdKitService;
export let _webRoot: string | null;
export let moduleInitialized = false;
/** The RDKit build compiled once here, so the RDKit workers only instantiate it. */
let _rdKitWasm: WebAssembly.Module | undefined;

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
  // addStereoAnnotation: true,
  simplifiedStereoGroupLabel: true,
  highlightColour: [1.0, 0.7, 0.7, 1.0],
  backgroundColour: [1, 1, 1, 1],
  queryColour: [0, 0, 0, 1],
};

export function setRdKitWebRoot(webRootValue: string): void {
  _webRoot = webRootValue;
}

export async function initRdKitModuleLocal(): Promise<void> {
  const url = `${_webRoot}/dist/${rdKitLibVersion}.wasm`;
  // streaming needs application/wasm; like Emscripten, fall back to the bytes from a server that does not send it
  _rdKitWasm = await WebAssembly.compileStreaming(fetch(url)).catch(async () => {
    const response = await fetch(url);
    if (!response.ok)
      throw new Error(`RDKit wasm is not available: ${response.status} ${response.statusText}`);
    return WebAssembly.compile(await response.arrayBuffer());
  });
  _rdKitModule = await initRdKitFrom(initRDKitModule, _rdKitWasm);
  if (!_rdKitModule)
    throw new Error('RdKit Module is not loaded');
  _rdKitModule.prefer_coordgen(false);
  _rdKitModule.use_legacy_stereo_perception(false);
  console.log('RDKit module package instance was initialized');
  moduleInitialized = true;
  // A second init (a test category's `before`) keeps the worker pool: a new service would leave the old one's workers
  // running for the rest of the page's life, each holding an RDKit instance, and a page has room for only so many
  // WebAssembly memories (instantiation fails with "Out of memory" past it)
  _rdKitService ??= new RdKitService();
}

export async function initRdKitService(): Promise<void> {
  await _rdKitService!.init(_webRoot!, _rdKitWasm!);
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

// Cancels operation `opId` if it's the running one, by restarting the worker(s) it uses (rejecting its
// in-flight call). No-op for a queued operation — that bails on its own canceled check after taking the section.
export async function cancelChemOp(opId: string, workers: number[]): Promise<void> {
  if (!isChemOpRunning(opId))
    return;
  const svc = await getRdKitService();
  if (!isChemOpRunning(opId)) // re-check after the await: the op may have finished, freeing its worker
    return;
  await Promise.all(workers.map((w) => svc.restartWorker(w)));
}

/** Chem's web root (the package's files, the wasm files' folder), set at the package's init: not RDKit's, despite the
 * name; it stays when RDKit goes. */
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
  return {...RDKIT_COMMON_RENDER_OPTS, ...addSettings};
}

/**
 * Draws the molecule (an RDKit one, or Crux's in Crux mode). False when the drawing threw: with `failQuietly` (a Crux
 * molecule, which RDKit then draws) nothing is drawn or logged, else an error cross is drawn.
 */
export function drawRdKitMoleculeToOffscreenCanvas(
  molCtx: IMolContext, w: number, h: number, offscreenCanvas: OffscreenCanvas, substruct: Object | null,
  renderingOptions: {[key: string]: any} = {}, failQuietly: boolean = false): boolean {
  const g = offscreenCanvas.getContext('2d', {willReadFrequently: true}) as OffscreenCanvasRenderingContext2D;
  const rdKitMol: RDMol | null = molCtx.mol;
  if (rdKitMol === null) {
    console.error('Molecule to be rendered cannot be null.');
    drawErrorCross(g, w, h);
    return true;
  }


  w = Math.floor(w);
  h = Math.floor(h);

  if (!w || !h) {
    console.warn('Molecule not rendered due to width or height less than 1');
    return true;
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
  if (renderingOptions)
    Object.assign(opts, renderingOptions);

  try {
    rdKitMol.draw_to_canvas_with_highlights((offscreenCanvas as unknown) as HTMLCanvasElement, JSON.stringify(opts));
  } catch (e) {
    if (failQuietly) { // a Crux drawing that stopped partway: the canvas as RDKit expects it, then RDKit draws
      noteCruxError(e);
      g.setTransform(1, 0, 0, 1, 0, 0);
      g.clearRect(0, 0, w, h);
      return false;
    }
    console.error('Molecule failed to render ' + rdKitMol.get_molblock());
    drawErrorCross(g, w, h);
    return false;
  }
  // we need the offscreen canvas first to not let the molecule scaffold skew on a real canvas
  return true;
}

export function drawRdKitReactionToOffscreenCanvas(
  rdKitReaction: RDReaction, w: number, h: number, offscreenCanvas: OffscreenCanvas) {
  const opts = createRenderingOpts({width: Math.floor(w), height: Math.floor(h), fixedScale: undefined});
  const g = offscreenCanvas.getContext('2d', {willReadFrequently: true});
  g?.clearRect(0, 0, w, h);
  rdKitReaction.draw_to_canvas_with_highlights((offscreenCanvas as unknown) as HTMLCanvasElement, JSON.stringify(opts));
}

export function drawMoleculeToCanvas(x: number, y: number, w: number, h: number,
  onscreenCanvas: HTMLCanvasElement, molString: string, scaffoldMolString: string | null = null,
  options = {normalizeDepiction: true, straightenDepiction: true}, scaffoldStruct: ISubstruct | null = null,
  renderingOptions?: {[key: string]: any}): void {
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

  // in Crux mode ("Chemical engine") Crux draws, and RDKit whatever Crux cannot read or refuses, and a CXSMILES whose
  // query features Crux's reader skips (the molecule's or the scaffold's)
  const cruxReads = (s: string | null): boolean => !s || cruxReadsCxExtension(s);
  byEngineModule('canvasMol', getRdKitModule(), (module, crux) => {
    if (crux && !(cruxReads(molString) && cruxReads(scaffoldMolString)))
      return null;
    const drawn = drawMoleculeWith(module, crux, x, y, w, h, nW, nH, onscreenCanvas, molString, isMol,
      scaffoldMolString, scaffoldStruct, renderingOptions);
    return drawn || !crux ? true : null;
  }, `RDKit draws what Crux cannot read or refuses, e.g. ${molString}`);
}

/** {@link drawMoleculeToCanvas} by one engine's module. Under Crux (`crux`) false, with nothing drawn, for a molecule
 * or a scaffold Crux cannot read or a drawing it refuses. */
function drawMoleculeWith(module: RDModule, crux: boolean, x: number, y: number, w: number, h: number, nW: number,
  nH: number, onscreenCanvas: HTMLCanvasElement, molString: string, isMol: boolean, scaffoldMolString: string | null,
  scaffoldStruct: ISubstruct | null, renderingOptions?: {[key: string]: any}): boolean {
  const offscreenCanvas = new OffscreenCanvas(nW, nH);
  const molCtx = getMolSafe(molString, {}, module, true, crux ? noteCruxError : undefined);
  const mol : RDMol | null = molCtx.mol;
  if (mol === null) {
    if (crux)
      return false;
    console.error('Molecule is null for ' + molString);
    drawErrorCross(offscreenCanvas.getContext('2d') as OffscreenCanvasRenderingContext2D, w, h);
    return true;
  }

  let scaffoldMol: RDMol | null = null;
  try {
    if (!isMol)
      mol.set_new_coords(true);

    const molHasOwnCoords = (mol.has_coords() > 0);
    if (molHasOwnCoords) {
      mol.normalize_depiction(0);
      mol.straighten_depiction(true);
    }

    let substruct: ISubstruct | null = null;
    if (scaffoldStruct)
      substruct = scaffoldStruct;
    else {
      let substructJson = '{}';
      scaffoldMol = scaffoldMolString == null ? null :
        (isMolBlock(scaffoldMolString) ? module.get_qmol(scaffoldMolString) :
          module.get_qmol(convertToRDKit(scaffoldMolString)!));
      if (crux && scaffoldMolString && !scaffoldMol) // ('' is no scaffold, as for RDKit's empty query)
        return false;
      if (scaffoldMol) {
        try {
          substructJson = mol.get_substruct_match(scaffoldMol);
        } catch (e) {
          if (crux)
            throw e;
          console.error(`get_substruct_match failed for ${molString} and ${scaffoldMolString}`);
        }
        if (substructJson === '')
          substructJson = '{}';
      }
      substruct = JSON.parse(substructJson);
    }
    if (!drawRdKitMoleculeToOffscreenCanvas(molCtx, nW, nH, offscreenCanvas, substruct, renderingOptions, crux) && crux)
      return false;
    const image = offscreenCanvas!.getContext('2d')!.getImageData(0, 0, nW, nH);
    const context = onscreenCanvas.getContext('2d')!;
    context.putImageData(image, x, y);
    return true;
  } finally {
    mol?.delete();
    scaffoldMol?.delete();
  }
}

/** `molecule`'s canonical SMILES (null where it does not read: {@link checkMoleculeValid}'s `get_mol`) on `module`;
 * under Crux (`crux`) what the module throws is thrown (the fallback notes a trap). */
export function canonicalWith(module: RDModule, crux: boolean, molecule: string): string | null {
  if (!hasNewLines(molecule) && molecule?.length > MAX_SMILES_LENGTH)
    return null;
  return withMol(module, crux, molecule, (mol) => mol.get_smiles());
}

/** Whether `molecule` reads as a molecule ({@link checkMoleculeValid}'s check) by the main thread's module: "Chemical
 * engine", Crux's in Crux mode, RDKit's where Crux reads none. */
export function readsAsMolecule(molecule: string): boolean {
  return byEngineModule('molecule check', getRdKitModule(),
    (module, crux) => canonicalWith(module, crux, molecule) !== null ? true : null) ?? false;
}

/** Whether two molecule strings are one molecule (the same canonical SMILES) by the main thread's module: "Chemical
 * engine", Crux's in Crux mode, RDKit's where Crux reads either none. A string that does not read is no match. */
export function sameMolecule(molecule1: string, molecule2: string): boolean {
  return byEngineModule('molecule compare', getRdKitModule(), (module, crux) => {
    const smiles1 = canonicalWith(module, crux, molecule1);
    const smiles2 = smiles1 === null ? null : canonicalWith(module, crux, molecule2);
    return smiles2 === null ? null : smiles1 === smiles2;
  }) ?? false;
}

export function checkMoleculeValid(molecule: string): any {
  let mol;
  try {
    if (!hasNewLines(molecule) && molecule?.length > MAX_SMILES_LENGTH)
      return null;
    mol = getRdKitModule().get_mol(molecule);
  } catch (e: any) {
    mol?.delete();
    return null;
  }
  return mol;
}

/** RDKit's message for a molecule string it does not read ('' for one it reads with nothing logged): what
 * `validateMolecule` shows. RDKit-only: Crux keeps no log. */
export function validationMessageByRdKit(rdKitModule: RDModule, molecule: string): string | null {
  let logHandle: RDLog | null = null;
  let mol: RDMol | null = null;
  try {
    logHandle = rdKitModule.set_log_capture('rdApp.error');
    mol = getMolSafe(molecule, {}, rdKitModule, true).mol;
    let logBuffer = logHandle?.get_buffer();
    logHandle?.clear_buffer();
    if (!mol && !logBuffer)
      logBuffer = 'unknown error';
    return logBuffer ?? null;
  } finally {
    logHandle?.delete();
    mol?.delete();
  }
}
