/* eslint-disable max-len */
/**
 * The Molecule cell renderer, for both engines (named "RDKit" for the packages and layouts that know it so: the class,
 * `rdKitCellRenderer`, its `name`). It draws with the engine in effect (`renderingModule`): in Crux mode Crux's module,
 * and RDKit's for a molecule Crux cannot read, lay out or draw ({@link RDKitCellRenderer._fetchMolByEngine},
 * `_rendererGetOrCreate`, `_molToSvg`); its caches hold one engine's molecules and pictures at a time (`_syncEngine`).
 * */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMonomerHover, getSubstructProviders, mergeSubstructs as mergeSubstructsLib}
  from '@datagrok-libraries/chem-meta/src/types';
import {ChemTags, ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

import {
  ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG,
  ALIGN_BY_SCAFFOLD_TAG, FILTER_SCAFFOLD_TAG,
  FIXED_SCALE_TAG,
  HIGHLIGHT_BY_SCAFFOLD_COL, HIGHLIGHT_BY_SCAFFOLD_COL_SYNC,
  HIGHLIGHT_BY_SCAFFOLD_TAG, MIN_MOL_IMAGE_SIZE, PARENT_MOL_COL,
  REGENERATE_COORDS, REGENERATE_COORDS_SYNC,
  SCAFFOLD_COL, SCAFFOLD_COL_SYNC, SCAFFOLD_TREE_HIGHLIGHT,
  getSyncTag,
} from '../constants';
import {hexToPercentRgb} from '../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../utils/chem-constants';
import {drawErrorCross, drawRdKitMoleculeToOffscreenCanvas,
  RDKIT_COMMON_RENDER_OPTS} from '../utils/chem-common-rdkit';
import {IMolContext, getMolSafe} from '../utils/mol-creation_rdkit';
import {ChemicalEngine, engineGeneration, renderingModule} from '../engine/chem-engine';
import {cruxModuleSync, noteCruxError} from '../engine/crux-module';
import {cruxOrRdKit, warnCruxFallback} from '../engine/fallback';
import {cruxReadsCxExtension} from '../crux/crux-service';
import {getGridCellColTemp} from '../utils/ui-utils';
import {errInfo} from '../utils/err-info';
import {extractAtomPositionsFromSvg} from '../utils/chem-atom-picker-utils';
import {AtomPickerController, CellInteractiveInfo} from './atom-picker-controller';

// Re-exported for cross-package consumers (BiostructureViewer, tests).
// `AtomPickerProvider` lives in `@datagrok-libraries/chem-meta/src/types` so
// BSV can read it without depending on Chem; consumers should import it
// directly from chem-meta rather than re-importing here.
export {AtomPickerController} from './atom-picker-controller';
export type {ChemSelectionEvent, CellInteractiveInfo}
  from './atom-picker-controller';

import {_package} from '../package';

interface IMolRenderingInfo {
  molCtx: IMolContext;
  substruct: ISubstruct;
  molString: string;
  /** The engine whose molecule `molCtx.mol` is (the "Chemical engine" setting; Crux mode falls back per molecule). */
  engine: ChemicalEngine;
}

export interface IColoredScaffold {
  molecule: string,
  color?: string,
  priority?: number,
  isSuperstructure?: string,
  align?: boolean,
  highlight?: boolean
}

export interface IHighlightTagInfo {
  scaffolds?: IColoredScaffold[],
  alignByFirstSubstructure: boolean,
}

export const NO_SCAFFOLD_COLOR = 'none';

export function _addColorsToBondsAndAtoms(mainSubstr: ISubstruct, color?: string, tempSubstr?: ISubstruct): void {
  const colorArr = color ? hexToPercentRgb(color) : [1.0, 0.7, 0.7, 1.0];
  const substrToTakeAtomsFrom = tempSubstr ?? mainSubstr;
  if (substrToTakeAtomsFrom) {
    if (substrToTakeAtomsFrom.atoms) {
      for (let j = 0; j < substrToTakeAtomsFrom.atoms.length; j++) {
        mainSubstr.highlightAtomColors ??= {};
        mainSubstr.highlightAtomColors[substrToTakeAtomsFrom.atoms[j]] = colorArr;
      }
    }
    if (substrToTakeAtomsFrom.bonds) {
      for (let j = 0; j < substrToTakeAtomsFrom.bonds.length; j++) {
        mainSubstr.highlightBondColors ??= {};
        mainSubstr.highlightBondColors[substrToTakeAtomsFrom.bonds[j]] = colorArr;
      }
    }
  }
}

export class GridCellRendererProxy extends DG.GridCellRenderer {
  renderer: DG.GridCellRenderer;
  _cellType: string;

  constructor(renderer: DG.GridCellRenderer, cellType: string) {
    super();
    this.renderer = renderer;
    this._cellType = cellType;
  }

  get defaultWidth(): number | null {return this.renderer.defaultWidth;}

  get defaultHeight(): number | null {return this.renderer.defaultHeight;}

  get name(): string {return this.renderer.name;}

  get cellType(): string {return this._cellType;}

  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    this.renderer.render(g, x, y, w, h, gridCell, cellStyle);
  }

  renderInternal(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    this.renderer.renderInternal(g, x, y, w, h, gridCell, cellStyle);
  }

  // -- Interaction event forwarding ----------------------------------------
  // Datagrok dispatches events to this proxy instance, not to the wrapped
  // renderer. Forward all overrides so RDKitCellRenderer's handlers fire.

  override onMouseEnter(gridCell: DG.GridCell, e: MouseEvent): void {
    this.renderer.onMouseEnter(gridCell, e);
  }

  override onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void {
    this.renderer.onMouseLeave(gridCell, e);
  }

  override onMouseDown(gridCell: DG.GridCell, e: MouseEvent): void {
    this.renderer.onMouseDown(gridCell, e);
  }

  override onMouseUp(gridCell: DG.GridCell, e: MouseEvent): void {
    this.renderer.onMouseUp(gridCell, e);
  }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    this.renderer.onMouseMove(gridCell, e);
  }

  override onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    this.renderer.onClick(gridCell, e);
  }

  override onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void {
    this.renderer.onDoubleClick(gridCell, e);
  }

  override onKeyDown(gridCell: DG.GridCell, e: KeyboardEvent): void {
    this.renderer.onKeyDown(gridCell, e);
  }

  override onKeyPress(gridCell: DG.GridCell, e: KeyboardEvent): void {
    this.renderer.onKeyPress(gridCell, e);
  }
}

export class RDKitCellRenderer extends DG.GridCellRenderer {
  readonly WHITE_MOLBLOCK_SUFFIX = `
  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
`;
  /** RDKit's module: RDKit mode's, and Crux mode's fallback per molecule. */
  rdKitModule: RDModule;
  molCache: DG.LruCache<String, IMolRenderingInfo> = new DG.LruCache<String, IMolRenderingInfo>();
  rendersCache: DG.LruCache<String, ImageData> = new DG.LruCache<String, ImageData>();
  sortedScaffoldsCache: DG.LruCache<String, IColoredScaffold[]> = new DG.LruCache<String, IColoredScaffold[]>();
  canvasReused: OffscreenCanvas;
  /** The molecules {@link molCache} holds, freed when it evicts them or when the engine changes. */
  private cachedMols = new Set<IMolRenderingInfo>();
  /** The engine generation ({@link engineGeneration}) the caches were filled under. */
  private cacheGeneration = -1;

  /** Interactive 2D↔3D atom-picker state and event handling. */
  picker: AtomPickerController;

  constructor(rdKitModule: RDModule) {
    super();
    this.rdKitModule = rdKitModule;
    this.canvasReused = new OffscreenCanvas(this.defaultWidth, this.defaultHeight);
    this.molCache.onItemEvicted = (info: IMolRenderingInfo) => this._freeMol(info);
    this.picker = new AtomPickerController({
      rdKitModule: this.rdKitModule,
      _getCellAtomPositions: this._getCellAtomPositions.bind(this),
    });
  }

  /** Frees a molecule the cache lets go of. */
  private _freeMol(info: IMolRenderingInfo): void {
    this.cachedMols.delete(info);
    try {
      info.molCtx.mol?.delete();
    } catch {
      // freed already
    }
    info.molCtx.mol = null;
  }

  /**
   * Empties the caches when the engine in effect changed since they were filled ("Chemical engine": a switch, or the
   * Crux module ready or failed): pictures and molecules are kept per engine, never mixed.
   */
  _syncEngine(): void {
    const generation = engineGeneration();
    if (generation === this.cacheGeneration)
      return;
    this.cacheGeneration = generation;
    for (const info of Array.from(this.cachedMols))
      this._freeMol(info);
    this.molCache = new DG.LruCache<String, IMolRenderingInfo>();
    this.molCache.onItemEvicted = (info: IMolRenderingInfo) => this._freeMol(info);
    this.rendersCache = new DG.LruCache<String, ImageData>();
    this.sortedScaffoldsCache = new DG.LruCache<String, IColoredScaffold[]>();
    this.picker.atomPositionsCache = new DG.LruCache<string, CellInteractiveInfo>();
  }

  /** How many molecules the caches hold now (tests). */
  get cachedMoleculeCount(): number {
    return this.cachedMols.size;
  }

  // -- Test/external accessors for picker state -----------------------------
  // Tests cast `RDKitCellRenderer` to a structural alias that names these
  // members directly — so we proxy them to the controller via getter/setter
  // pairs. Same applies to event handler shims further down.

  get atomPositionsCache(): DG.LruCache<string, CellInteractiveInfo> {
    return this.picker.atomPositionsCache;
  }

  private get _lastHoveredAtom(): {col: string, rowIdx: number, atomIdx: number, erase?: boolean} | null {
    return this.picker._lastHoveredAtom;
  }
  private set _lastHoveredAtom(v: {col: string, rowIdx: number, atomIdx: number, erase?: boolean} | null) {
    this.picker._lastHoveredAtom = v;
  }

  private get _previewAtomIdx(): number | null {return this.picker._previewAtomIdx;}
  private set _previewAtomIdx(v: number | null) {this.picker._previewAtomIdx = v;}

  private get _previewRowIdx(): number {return this.picker._previewRowIdx;}
  private set _previewRowIdx(v: number) {this.picker._previewRowIdx = v;}

  private get _previewFrom3D(): boolean {return this.picker._previewFrom3D;}
  private set _previewFrom3D(v: boolean) {this.picker._previewFrom3D = v;}

  private _trackHoveredAtom(colName: string, rowIdx: number, atomIdx: number,
    erase?: boolean): boolean {
    return this.picker._trackHoveredAtom(colName, rowIdx, atomIdx, erase);
  }

  ensureCanvasSize(w: number, h: number): OffscreenCanvas {
    if (this.canvasReused.width < w || this.canvasReused.height < h)
      this.canvasReused = new OffscreenCanvas(Math.max(this.defaultWidth, w), Math.max(this.defaultHeight, h));
    return this.canvasReused;
  }

  get name(): string {return 'RDKit cell renderer';}

  get cellType(): DG.SemType {return DG.SEMTYPE.MOLECULE;}

  get defaultWidth(): number {return 200;}

  get defaultHeight(): number {return 100;}

  getDefaultSize(): { width: number, height: number } {
    return {width: this.defaultWidth, height: this.defaultHeight};
  }

  /**
   * The molecule of `molString` laid out and matched against `scaffolds`, all by one `engine`. Under Crux, anything
   * Crux cannot read or refuses (the molecule, a scaffold, a layout call) leaves the molecule null, and
   * {@link _fetchMolByEngine} then takes RDKit's for it.
   */
  _fetchMolGetOrCreate(molString: string, scaffolds: IColoredScaffold[], molRegenerateCoords: boolean,
    details: object = {}, alignByFirstSubstr: boolean, engine: ChemicalEngine = ChemicalEngine.RDKit,
  ): IMolRenderingInfo {
    const crux = engine === ChemicalEngine.Crux;
    const module: RDModule = crux ? cruxModuleSync()! : this.rdKitModule;
    let molCtx: IMolContext;
    let mol = null;
    let substruct: ISubstruct = {};
    if ((details as any).isSubstructure) {
      const mappedDummiesAreRGroups = (details as any).mappedDummiesAreRGroups || false;
      if (molString.includes(' H ') || molString.includes('V3000')) {
        molCtx = getMolSafe(molString, {mergeQueryHs: true, mappedDummiesAreRGroups}, module, true,
          crux ? noteCruxError : undefined);
        mol = molCtx.mol;
      } else {
        try {
          if (!molString.includes('\n') && molString.length > MAX_SMILES_LENGTH)
            throw new Error('Invalid molecule string'); // do not attempt to parse very long SMILES, will cause MOB.
          mol = module.get_qmol(molString);
          mol.convert_to_aromatic_form();
        } catch (e) {
          if (crux)
            noteCruxError(e);
          if (mol) {
            mol.delete();
            mol = null;
          }
        }
        molCtx = {mol: mol, kekulize: false, isQMol: true, useMolBlockWedging: false};
      }
    } else {
      molCtx = getMolSafe(molString, details, module, true, crux ? noteCruxError : undefined);
      mol = molCtx.mol;
    }
    /** Under Crux, a scaffold Crux cannot read: the whole molecule goes to RDKit (a scaffold RDKit cannot read either
     * is left out there, as before; an empty one is no scaffold). */
    const cruxMissed = (scaffoldMol: unknown, scaffold: string) => {
      if (crux && !scaffoldMol && scaffold)
        throw new Error('Crux cannot read a scaffold');
    };

    if (mol !== null) {
      try {
        let molHasOwnCoords = (mol.has_coords() > 0);
        let molHasRebuiltCoords = false;
        const scaffoldIsMolBlock = scaffolds.length && alignByFirstSubstr ?
          DG.chem.isMolBlock(scaffolds[0].molecule) : null;
        const alignedByFirstSubstr = scaffoldIsMolBlock && alignByFirstSubstr;
        const {haveReferenceSmarts, parentMolScaffoldMolString} = (details as any);
        if (alignedByFirstSubstr) {
          const scaffoldMol = this._fetchMol(scaffolds[0].molecule,
            parentMolScaffoldMolString ? [{molecule: parentMolScaffoldMolString}] : [], molRegenerateCoords,
            false, {mergeQueryHs: true, isSubstructure: !parentMolScaffoldMolString}, false, engine).molCtx.mol;
          cruxMissed(scaffoldMol, scaffolds[0].molecule);
          if (scaffoldMol) {
            scaffoldMol.normalize_depiction(0);
            if (molHasOwnCoords)
              mol.normalize_depiction(0);
            else if (!haveReferenceSmarts) {
              //need the following 4 rows for smiles with highlights to be rendered in adequate coordinates
              mol.set_new_coords();
              mol.normalize_depiction(1);
              mol.straighten_depiction(false);
              molHasRebuiltCoords = true;
            }
            let substructString = '';
            let useCoordGen = (details as any).useCoordGen;
            useCoordGen = (typeof useCoordGen === 'boolean' ? useCoordGen : true);
            const alignOpts = {
              useCoordGen,
              allowRGroups: true,
              acceptFailure: false,
              alignOnly: molHasOwnCoords || molHasRebuiltCoords,
            };
            let referenceSmarts;
            if (haveReferenceSmarts) {
              try {
                referenceSmarts = mol.get_smiles();
              } catch {
                // do nothing
              }
            }
            if (referenceSmarts)
              (alignOpts as any).referenceSmarts = referenceSmarts;

            try {
              substructString = !scaffolds[0].isSuperstructure ?
                mol.generate_aligned_coords(scaffoldMol, JSON.stringify(alignOpts)) :
                mol.get_substruct_match(mol!);
            } catch (e) {
              // exceptions should not be thrown anymore by RDKit, but let's play safe; one from Crux (a refusal) sends
              // the molecule to RDKit
              if (crux)
                throw e;
            }
            if (substructString === '') {
              substruct = {};
              if (molHasOwnCoords)
                mol.straighten_depiction(true);
            } else {
              if (scaffolds[0].color !== NO_SCAFFOLD_COLOR) {
                substruct = JSON.parse(substructString);
                if (!substruct.atoms)
                  substruct.atoms = [];
                if (!substruct.bonds)
                  substruct.bonds = [];
                _addColorsToBondsAndAtoms(substruct, scaffolds[0].color);
              }
            }
          }
        }
        for (let i = alignedByFirstSubstr ? 1 : 0; i < scaffolds.length; i++) {
          if (scaffolds[i].color !== NO_SCAFFOLD_COLOR) {
            const substructMol = this._fetchMol(scaffolds[i].molecule, [], molRegenerateCoords, false,
              {mergeQueryHs: true, isSubstructure: true}, false, engine).molCtx.mol;
            cruxMissed(substructMol, scaffolds[i].molecule);
            if (substructMol) {
              const matchedAtomsAndBonds: ISubstruct[] = JSON.parse(mol!.get_substruct_matches(substructMol!));
              if (matchedAtomsAndBonds.length) {
                for (let j = 0; j < matchedAtomsAndBonds.length; j++) {
                  if (!substruct.atoms)
                    substruct.atoms = [];
                  if (!substruct.bonds)
                    substruct.bonds = [];
                  this._addAtomsOrBonds(matchedAtomsAndBonds[j].atoms!, substruct.atoms!);
                  this._addAtomsOrBonds(matchedAtomsAndBonds[j].bonds!, substruct.bonds!);
                  _addColorsToBondsAndAtoms(substruct, scaffolds[i].color, matchedAtomsAndBonds[j]);
                }
              }
            }
          }
        }
        molCtx.useMolBlockWedging = molHasOwnCoords;
        if (mol.has_coords() === 0 || molRegenerateCoords || hasNonZeroZCoords(molString, mol.get_num_atoms())) {
          mol.set_new_coords(molRegenerateCoords);
          mol.normalize_depiction(1);
          mol.straighten_depiction();
          molHasOwnCoords = false;
        }
        if (!scaffoldIsMolBlock) {
          mol.normalize_depiction(molHasOwnCoords ? 0 : 1);
          mol.straighten_depiction(molHasOwnCoords);
        } else if (!molHasOwnCoords)
          mol.normalize_depiction(0);
      } catch (e) {
        if (crux)
          noteCruxError(e);
        else {
          console.error(
            'In _fetchMolGetOrCreate: RDKit crashed, possibly a malformed molString molecule: `' + molString + '`');
        }
        if (mol !== null) {
          mol.delete();
          molCtx.mol = null;
        }
      }
    }

    return {
      molCtx: molCtx,
      substruct: substruct,
      molString: molString,
      engine: engine,
    };
  }

  _addAtomsOrBonds(fromAtomsOrBonds: number[], toAtomsOrBonds: number[]): void {
    if (!toAtomsOrBonds) return;
    const seen = new Set(toAtomsOrBonds);
    for (let j = 0; j < fromAtomsOrBonds.length; j++) {
      const item = fromAtomsOrBonds[j];
      if (!seen.has(item)) {
        seen.add(item);
        toAtomsOrBonds.push(item);
      }
    }
  }

  /**
   * The cached molecule of `molString` by `engine`; without an engine, by the engine in effect
   * ({@link _fetchMolByEngine}).
   */
  _fetchMol(molString: string, scaffolds: IColoredScaffold[], molRegenerateCoords: boolean,
    scaffoldRegenerateCoords: boolean, details: object = {}, alignByFirstSubstructure: boolean,
    engine?: ChemicalEngine): IMolRenderingInfo {
    if (engine === undefined) {
      return this._fetchMolByEngine(molString, scaffolds, molRegenerateCoords, scaffoldRegenerateCoords, details,
        alignByFirstSubstructure);
    }
    const name = engine + ' || ' + molString + ' || ' + JSON.stringify(scaffolds) + ' || ' +
      molRegenerateCoords + ' || ' + scaffoldRegenerateCoords + ' || ' +
      ((details as any).isSubstructure || false).toString() +
      ((details as any).haveReferenceSmarts || false).toString();
    return this.molCache.getOrCreate(name, (_: any) => {
      const info = this._fetchMolGetOrCreate(molString, scaffolds, molRegenerateCoords, details,
        alignByFirstSubstructure, engine);
      this.cachedMols.add(info);
      return info;
    });
  }

  /** {@link _fetchMol} by the engine in effect: in Crux mode by RDKit when Crux cannot read or lay the molecule out
   * (logged once). */
  private _fetchMolByEngine(molString: string, scaffolds: IColoredScaffold[], molRegenerateCoords: boolean,
    scaffoldRegenerateCoords: boolean, details: object, alignByFirstSubstructure: boolean): IMolRenderingInfo {
    this._syncEngine();
    const fetch = (engine: ChemicalEngine) => this._fetchMol(molString, scaffolds, molRegenerateCoords,
      scaffoldRegenerateCoords, details, alignByFirstSubstructure, engine);
    if (renderingModule(this.rdKitModule).engine !== ChemicalEngine.Crux)
      return fetch(ChemicalEngine.RDKit);
    // a CXSMILES whose query features Crux's reader skips (the molecule's, a scaffold's): RDKit draws it
    if (!(cruxReadsCxExtension(molString) && scaffolds.every((s) => cruxReadsCxExtension(s.molecule ?? '')))) {
      warnCruxFallback('rendering', `RDKit draws a CXSMILES query Crux's reader skips, e.g. ${molString}`);
      return fetch(ChemicalEngine.RDKit);
    }
    const crux = fetch(ChemicalEngine.Crux);
    if (crux.molCtx.mol !== null)
      return crux;
    const rdkit = fetch(ChemicalEngine.RDKit);
    if (rdkit.molCtx.mol !== null) // (a molecule neither engine reads is no fallback)
      warnCruxFallback('rendering', `RDKit draws what Crux cannot read or lay out, e.g. ${molString}`);
    return rdkit;
  }

  _rendererGetOrCreate(
    width: number, height: number, molString: string, scaffolds: IColoredScaffold[],
    molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean,
    alignByFirstSubstructure: boolean, details: object = {}, substructureObj?: ISubstruct, renderingOptions?: any,
    engine?: ChemicalEngine): ImageData {
    const fetchMolObj: IMolRenderingInfo =
      this._fetchMol(molString, scaffolds, molRegenerateCoords,
        scaffoldRegenerateCoords, details, alignByFirstSubstructure, engine);
    const molCtx = fetchMolObj.molCtx;
    const substruct = fetchMolObj.substruct;
    //merge row highlight data with substruct object from fetchMolObj
    const newSubstruct: ISubstruct = {};
    if (substructureObj && scaffolds.length) {
      const newAtoms = substruct.atoms ? substruct.atoms.concat(substructureObj.atoms ?? []) :
        substructureObj.atoms ?? [];
      const newAtomsUnique = newAtoms.filter((item, pos) => newAtoms.indexOf(item) === pos);
      const newBonds = substruct.bonds ? substruct.bonds.concat(substructureObj.bonds ?? []) :
        substructureObj.bonds ?? [];
      const newBondsUnique = newBonds.filter((item, pos) => newBonds.indexOf(item) === pos);
      if (substruct.highlightAtomColors) {
        if (substructureObj.highlightAtomColors) {
          for (const key in substruct.highlightAtomColors)
            substructureObj.highlightAtomColors[key] = substruct.highlightAtomColors[key];
        }
      }
      if (substruct.highlightBondColors) {
        if (substructureObj.highlightBondColors) {
          for (const key in substruct.highlightBondColors)
            substructureObj.highlightBondColors[key] = substruct.highlightBondColors[key];
        }
      }
      newSubstruct.atoms = newAtomsUnique;
      newSubstruct.bonds = newBondsUnique;
      newSubstruct.highlightAtomColors = substructureObj.highlightAtomColors;
      newSubstruct.highlightBondColors = substructureObj.highlightBondColors;
    }

    const canvas = this.ensureCanvasSize(width, height);
    const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
    if (molCtx.mol != null) {
      const crux = fetchMolObj.engine === ChemicalEngine.Crux;
      const drawn = drawRdKitMoleculeToOffscreenCanvas(molCtx, width, height, canvas,
        scaffolds.length ? substructureObj ? newSubstruct : substruct : substructureObj ?? null, renderingOptions,
        crux);
      if (!drawn && crux) {
        // Crux refused the drawing: RDKit draws this molecule
        warnCruxFallback('drawing', `RDKit draws what Crux refuses to, e.g. ${molString}`);
        return this._rendererGetOrCreate(width, height, molString, scaffolds, molRegenerateCoords,
          scaffoldRegenerateCoords, alignByFirstSubstructure, details, substructureObj, renderingOptions,
          ChemicalEngine.RDKit);
      }
    } else {
      // draw a crossed rectangle
      ctx.clearRect(0, 0, width, height);
      drawErrorCross(ctx, width, height);
    }

    return ctx.getImageData(0, 0, !Math.floor(width) ? 1 : width, !Math.floor(height) ? 1 : height);
  }

  _fetchRender(
    width: number, height: number, molString: string, scaffolds: IColoredScaffold[],
    molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean,
    alignByFirstSubstructure: boolean, details: object = {}, substructureObj?: ISubstruct, renderOptions?: {[key: string]: any}): ImageData {
    this._syncEngine();
    const name = width + ' || ' + height + ' || ' +
      molString + ' || ' + JSON.stringify(scaffolds) + ' || ' +
      molRegenerateCoords + ' || ' + scaffoldRegenerateCoords + ' || ' +
      ((details as any).haveReferenceSmarts || false).toString() + ' || ' + JSON.stringify(substructureObj) + ' || ' + JSON.stringify(renderOptions);
    return this.rendersCache.getOrCreate(name, (_: any) => this._rendererGetOrCreate(width, height,
      molString, scaffolds, molRegenerateCoords, scaffoldRegenerateCoords,
      alignByFirstSubstructure, details, substructureObj, renderOptions));
  }

  /**
   * Rasterize a fully-prepared molecule string (SMILES, or a laid-out molblock the caller has already
   * aligned / labelled) to a cached `ImageData` at the given device size on a TRANSPARENT background,
   * reusing this renderer's shared `molCache` and `rendersCache`.
   *
   * Background-free on purpose: the cache key is then the molecule and size alone, so a caller that
   * tints per cell composites this raster over its own background instead of baking a per-cell colour
   * into the key — which would make every tinted cell a distinct entry and thrash the shared LRU.
   *
   * For callers that paint molecules onto a canvas of their own (the SAR Matrix grid) rather than
   * going through `render()`. The returned `ImageData` is owned by `rendersCache` — read it, don't
   * mutate it.
   */
  getCachedMolImageData(molString: string, width: number, height: number): ImageData {
    return this._fetchRender(width, height, molString, [], false, false, false, {}, undefined,
      {clearBackground: false});
  }

  _drawMolecule(x: number, y: number, w: number, h: number, onscreenCanvas: HTMLCanvasElement,
    molString: string, scaffolds: IColoredScaffold[],
    molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean, cellStyle: DG.GridCellStyle,
    alignByFirstSubstructure: boolean, details: object = {}, substructureObj?: ISubstruct, renderOptions?: Record<string, any>): void {
    const vertical = cellStyle !== undefined && cellStyle !== null ? cellStyle.textVertical : false;

    if (vertical) {
      h += w;
      w = h - w;
      h -= w;
    }
    if (Math.floor(w) < MIN_MOL_IMAGE_SIZE || Math.floor(h) < MIN_MOL_IMAGE_SIZE)
      return;

    const imageData = this._fetchRender(w, h, molString, scaffolds,
      molRegenerateCoords, scaffoldRegenerateCoords,
      alignByFirstSubstructure, details, substructureObj, renderOptions);

    if (vertical) {
      const ctx = onscreenCanvas.getContext('2d', {willReadFrequently: true})!;
      ctx.save();
      const scl = ctx.getTransform();
      ctx.resetTransform();
      ctx.translate(x + h, y);
      ctx.rotate(Math.PI / 2);
      if (scl.m11 < 1 || scl.m22 < 1)
        ctx.scale(scl.m11, scl.m22);
      const f = new OffscreenCanvas(imageData.width, imageData.height)!;
      f.getContext('2d')?.putImageData(imageData, 0, 0);
      ctx.drawImage(f, 0, 0);
      ctx.restore();
    } else
      onscreenCanvas.getContext('2d', {willReadFrequently: true})!.putImageData(imageData, x, y);
  }

  _initScaffoldString(col: DG.Column, tagName: string): IColoredScaffold[] {
    const scaffoldString = col.getTag(tagName) ?? col.temp?.[tagName];
    if (scaffoldString?.endsWith(this.WHITE_MOLBLOCK_SUFFIX)) {
      const objToRemoveFrom = col.getTag(tagName) ? col.tags : col.temp;
      if (objToRemoveFrom[tagName])
        delete objToRemoveFrom[tagName];
      return [];
    }
    return scaffoldString ? [{molecule: scaffoldString}] : [];
  }

  _initScaffoldArray(col: any, tagName: string, isTempCol?: boolean): IColoredScaffold[] {
    const scaffoldArrStr = !isTempCol ? col.getTag(tagName) : col ? col[tagName] : null;
    const getSortedScaffolds = (): IColoredScaffold[] => {
      let scaffoldArr: IColoredScaffold[] = JSON.parse(scaffoldArrStr);
      if (scaffoldArr.length > 0 && !scaffoldArr[0].hasOwnProperty('molecule'))
        scaffoldArr = scaffoldArr.reduce((acc: any, obj: any) => acc.concat(Object.values(obj)[0]), []);

      const scaffoldArrSorted = scaffoldArr.sort((a: any, b: any) => {
        const getNumAtoms = (molecule: string) => {
          if (molecule && !DG.chem.Sketcher.isEmptyMolfile(molecule)) {
            const mol = this._fetchMol(molecule, [], false, false, {}, false).molCtx.mol;
            if (mol)
              return mol.get_num_atoms();
            return 0;
          }
          return 0;
        };
        a.priority = a.priority ?? 1;
        b.priority = b.priority ?? 1;

        if (a.priority !== b.priority)
          return b.priority - a.priority;


        return getNumAtoms(a.molecule) - getNumAtoms(b.molecule);
      });
      return scaffoldArrSorted;
    };
    if (scaffoldArrStr) {
      const sortedScaffolds = this.sortedScaffoldsCache
        .getOrCreate(scaffoldArrStr, () => getSortedScaffolds());
      return sortedScaffolds;
    }
    return [];
  }

  render(g: any, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    this.picker.ensureGlobalListeners();

    const molString = gridCell.cell.value;
    if (molString == null || molString === '')
      return;

    const r = window.devicePixelRatio;
    x = r * x;
    y = r * y;
    w = r * w;
    h = r * h;

    // value-based drawing (coming from HtmlCellRenderer.renderValue)
    if (gridCell.cell.column == null) {
      this._drawMolecule(x, y, w, h, g.canvas, molString, [], false, false, cellStyle, false, undefined, undefined, {fixedScale: null});
      return;
    }

    const colTemp = gridCell.cell.column.temp;
    const highlightInfo = this.getHighlightTagInfo(colTemp, gridCell);
    // fixed scale enabled by default, if explicitely set to false - then disabled
    const fixedScaleEnabled = gridCell.cell.column.getTag(FIXED_SCALE_TAG) !== 'false';
    const renderOpts = fixedScaleEnabled ? undefined : {fixedScale: null};
    let mhSubstruct: ISubstruct | undefined;
    try {
      const mhData = getMonomerHover();
      if (mhData && mhData.dataFrameId == gridCell.cell.column.dataFrame.id &&
        mhData.gridRowIdx === gridCell.cell.dataFrame.mouseOverRowIdx &&
        mhData.seqColName ===gridCell.cell.column?.getTag(ChemTags.SEQUENCE_SRC_COL)
      )
        mhSubstruct = mhData.getSubstruct();
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
    }

    // TODO: make both filtering scaffold and single highlight scaffold appear
    if (mhSubstruct) {
      this._drawMolecule(x, y, w, h, g.canvas,
        molString, [], false, false, cellStyle, false, undefined, mhSubstruct, renderOpts);
    } else if (highlightInfo.scaffolds && highlightInfo.alignByFirstSubstructure) {
      this._drawMolecule(x, y, w, h, g.canvas,
        molString, highlightInfo.scaffolds, false, false, cellStyle, highlightInfo.alignByFirstSubstructure, undefined, undefined, renderOpts);
    } else
      this.highlightByScaffoldCol(g, x, y, w, h, gridCell, cellStyle, colTemp, molString, highlightInfo.scaffolds, renderOpts);
  }

  getHighlightTagInfo(colTemp: any, gridCell: DG.GridCell): IHighlightTagInfo {
    const col = gridCell.cell.column;
    const highlightScaffold = getSyncTag(col, HIGHLIGHT_BY_SCAFFOLD_COL_SYNC, HIGHLIGHT_BY_SCAFFOLD_COL) === 'true';
    const filter = this._initScaffoldArray(colTemp, FILTER_SCAFFOLD_TAG, true); //expected molBlock
    const align = (this._initScaffoldString(col, ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG) ?? [])
      // temporary concatination to support both tags
      .concat(this._initScaffoldString(col, ALIGN_BY_SCAFFOLD_TAG) ?? []);
    if (!highlightScaffold) {
      for (const s of align)
        s.color = NO_SCAFFOLD_COLOR;
    }
    const highlight = this._initScaffoldArray(col, HIGHLIGHT_BY_SCAFFOLD_TAG);
    const scaffoldTreeHighlight = this._initScaffoldArray(col, SCAFFOLD_TREE_HIGHLIGHT);
    const alignByStructure = !!(filter.length && filter[0].align || align.length);
    const scaffolds = filter.concat(align).concat(scaffoldTreeHighlight).concat(highlight);
    return {scaffolds: scaffolds?.length ? scaffolds : undefined, alignByFirstSubstructure: alignByStructure};
  }

  highlightByScaffoldCol(g: any, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle, colTemp: any, molString: string, highlightScaffolds?: IColoredScaffold[], renderOpts?: Record<string, any>): void {
    let molRegenerateCoords = getSyncTag(gridCell.cell.column, REGENERATE_COORDS_SYNC, REGENERATE_COORDS) === 'true';
    let scaffoldRegenerateCoords = false;
    const df = gridCell.cell.dataFrame;
    let rowScaffoldCol = null;
    let haveParentMol = false;
    let parentMolScaffoldMolString;

    // if given, take the 'scaffold-col' col
    if (colTemp) {
      let rowScaffoldColName = getSyncTag(gridCell.cell.column, SCAFFOLD_COL_SYNC, SCAFFOLD_COL);
      if (!rowScaffoldColName) {
        rowScaffoldColName = colTemp[PARENT_MOL_COL];
        haveParentMol = !!rowScaffoldColName;
      }
      if (rowScaffoldColName) {
        const rowScaffoldColProbe = df.columns.byName(rowScaffoldColName);
        if (rowScaffoldColProbe !== null) {
          if (haveParentMol) {
            const parentMolScaffoldColName = getSyncTag(rowScaffoldColProbe, SCAFFOLD_COL_SYNC, SCAFFOLD_COL);
            if (parentMolScaffoldColName) {
              const idx = gridCell.tableRowIndex; // TODO: supposed to be != null?
              parentMolScaffoldMolString = df.get(parentMolScaffoldColName, idx!);
            }
          }
          scaffoldRegenerateCoords = getSyncTag(rowScaffoldColProbe, REGENERATE_COORDS_SYNC, REGENERATE_COORDS) === 'true';
          molRegenerateCoords = scaffoldRegenerateCoords;
          rowScaffoldCol = rowScaffoldColProbe;
        }
      }
    }

    const idx = gridCell.tableRowIndex; // TODO: supposed to be != null?

    //check for column with per-row ISubstruct objects for highlight
    let substructObj: ISubstruct | undefined = undefined;
    let alignByScaffold = false;
    if (colTemp[ChemTemps.SUBSTRUCT_COL]) {
      const rawSubstructCol = df.columns.byName(colTemp[ChemTemps.SUBSTRUCT_COL]);
      if (rawSubstructCol)
        substructObj = rawSubstructCol.get(idx!);
    } else {
      const [_gridCol, tableCol, _temp] = getGridCellColTemp(gridCell);
      const substructList = (getSubstructProviders(tableCol?.temp) ?? [])
        .map((p) => {
          let res: ISubstruct | undefined;
          try {
            res = p.getSubstruct(gridCell.tableRowIndex);
            if (res?.alignByScaffold) {
              alignByScaffold = true;
              const scaffold: IColoredScaffold = {
                molecule: res?.alignByScaffold,
                highlight: false,
                align: true,
                color: NO_SCAFFOLD_COLOR,
              };
              if (!highlightScaffolds)
                highlightScaffolds = [scaffold];
              else
                highlightScaffolds.push(scaffold);
            }
          } catch (err: any) {
            const [errMsg, errStack] = errInfo(err);
            _package.logger.error(errMsg, undefined, errStack);
          }
          return res;
        });
      substructObj = mergeSubstructs(substructList);
    }

    if (rowScaffoldCol == null || rowScaffoldCol.name === gridCell.cell.column.name) {
      this._drawMolecule(x, y, w, h, g.canvas, molString, highlightScaffolds ?? [],
        molRegenerateCoords, false, cellStyle, alignByScaffold, {}, substructObj, renderOpts);
    } else {
      // drawing with a per-row scaffold
      const scaffoldMolString = df.get(rowScaffoldCol.name, idx!);
      const highlightScaffold = getSyncTag(gridCell.cell.column, HIGHLIGHT_BY_SCAFFOLD_COL_SYNC, HIGHLIGHT_BY_SCAFFOLD_COL) === 'true';
      const details = (haveParentMol ? {
        mappedDummiesAreRGroups: true,
        useCoordGen: false,
        haveReferenceSmarts: true,
        parentMolScaffoldMolString,
      } : {});
      const scaffoldFromColumn: IColoredScaffold[] = highlightScaffold ? [{molecule: scaffoldMolString}] :
        [{molecule: scaffoldMolString, color: NO_SCAFFOLD_COLOR}];
      const totalScaffolds = highlightScaffolds ? scaffoldFromColumn.concat(highlightScaffolds) : scaffoldFromColumn;
      this._drawMolecule(x, y, w, h, g.canvas, molString, totalScaffolds, molRegenerateCoords,
        scaffoldRegenerateCoords, cellStyle, true, details, substructObj, renderOpts);
    }
  }

  // -- Atom picker event hooks ----------------------------------------------
  // Slim shims that delegate to AtomPickerController. State and logic live
  // in `./atom-picker-controller.ts`; the renderer only proxies grid events.

  override onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void {
    this.picker.onMouseLeave(gridCell, e);
  }

  override onKeyDown(gridCell: DG.GridCell, e: KeyboardEvent): void {
    this.picker.onKeyDown(gridCell, e);
  }

  override onClick(gridCell: DG.GridCell<any>, e: MouseEvent): void {
    this.onMouseDown(gridCell, e);
  }

  override onMouseDown(gridCell: DG.GridCell<any>, e: MouseEvent): void {
    this.picker.onMouseDown(gridCell, e);
  }

  override onMouseUp(gridCell: DG.GridCell<any>, e: MouseEvent): void {
    e.preventDefault();
  }

  override onMouseMove(gridCell: DG.GridCell<any>, e: MouseEvent): void {
    if (e.shiftKey || e.ctrlKey || e.metaKey)
      this.onMouseDown(gridCell, e);
  }

  /**
   * Returns atom positions in cell-local CSS pixels, cached by "molString|WxH".
   *
   * Uses the same `IMolContext` as the renderer (via `_fetchMol`) so the SVG
   * layout matches the canvas-rasterized layout exactly — a different mol
   * instance can produce different 2D coords (kekulize / molBlockWedging).
   * Scales by DPR to match `render()`, then divides SVG coords back by DPR.
   */
  private _getCellAtomPositions(molString: string, cssWidth: number, cssHeight: number):
      CellInteractiveInfo | null {
    const dpr = window.devicePixelRatio || 1;
    const w = Math.round(cssWidth * dpr);
    const h = Math.round(cssHeight * dpr);
    const key = molString + '|' + w + 'x' + h;
    const hit = this.atomPositionsCache.get(key);
    if (hit) return hit;

    try {
      const svgString = this._molToSvg(molString, w, h);
      if (!svgString)
        return null;
      // Attach to document so getBBox() works (requires layout context).
      const host = document.createElement('div');
      host.style.position = 'absolute';
      host.style.left = '-9999px';
      host.style.top = '-9999px';
      host.style.width = w + 'px';
      host.style.height = h + 'px';
      host.innerHTML = svgString;
      document.body.appendChild(host);
      try {
        const svgEl = host.querySelector('svg') as SVGSVGElement | null;
        if (!svgEl) return null;
        const extracted = extractAtomPositionsFromSvg(svgEl);
        // Divide SVG coords by DPR to get cell-local CSS pixels.
        const positions = new Map<number, {x: number, y: number}>();
        for (const [idx, p] of extracted.positions.entries())
          positions.set(idx, {x: p.x / dpr, y: p.y / dpr});
        const info: CellInteractiveInfo = {positions, bondAtoms: extracted.bondAtoms};
        this.atomPositionsCache.set(key, info);
        return info;
      } finally {
        if (host.parentNode) host.parentNode.removeChild(host);
      }
    } catch {
      return null;
    }
  }

  _molToSvg(
    molString: string, w: number, h: number, scaffolds: IColoredScaffold[] = [], alignByFirstSubstructure = false,
    engine?: ChemicalEngine,
  ): string | null {
    // Mol is owned by molCache — DO NOT delete it.
    const molRenderingInfo = this._fetchMol(molString, scaffolds, false, false, {}, alignByFirstSubstructure, engine);
    if (molRenderingInfo.engine === ChemicalEngine.Crux) {
      return cruxOrRdKit('drawing', () => this._infoToSvg(molRenderingInfo, w, h, scaffolds),
        () => this._molToSvg(molString, w, h, scaffolds, alignByFirstSubstructure, ChemicalEngine.RDKit));
    }
    return this._infoToSvg(molRenderingInfo, w, h, scaffolds);
  }

  private _infoToSvg(molRenderingInfo: IMolRenderingInfo, w: number, h: number, scaffolds: IColoredScaffold[],
  ): string | null {
    const mol = molRenderingInfo.molCtx.mol;
    if (!mol || !mol.is_valid())
      return null;

    const details: {[k: string]: any} = {};
    for (const k of Object.keys(RDKIT_COMMON_RENDER_OPTS))
      details[k] = RDKIT_COMMON_RENDER_OPTS[k];
    details.width = w;
    details.height = h;
    // Scaffold-derived highlights (substructure filter, align, scaffold tree, highlight-by-scaffold).
    const substruct = scaffolds.length ? molRenderingInfo.substruct : null;
    details.atoms = substruct?.atoms ?? [];
    details.bonds = substruct?.bonds ?? [];
    details.highlightAtomColors = substruct?.highlightAtomColors ?? {};
    details.highlightBondColors = substruct?.highlightBondColors ?? {};
    // Mirror kekulize / molBlockWedging from drawRdKitMoleculeToOffscreenCanvas
    // so the SVG layout matches the canvas layout exactly.
    if (!molRenderingInfo.molCtx.kekulize)
      details.kekulize = false;
    if (molRenderingInfo.molCtx.useMolBlockWedging) {
      details.useMolBlockWedging = true;
      details.wedgeBonds = false;
      details.addChiralHs = false;
    }

    return mol.get_svg_with_highlights(JSON.stringify(details));
  }

  toSvg(gridCell: DG.GridCell, w: number, h: number): string | null {
    const highlightInfo = gridCell.cell.column != null ?
      this.getHighlightTagInfo(gridCell.cell.column.temp, gridCell) :
      undefined;
    return this._molToSvg(gridCell.cell.value, w, h,
      highlightInfo?.scaffolds ?? [], highlightInfo?.alignByFirstSubstructure ?? false);
  }
}

function hasNonZeroZCoords(molfile: string, numAtoms: number): boolean {
  const moveCursorToIdx = (steps: number, symbol: string) => {
    for (let i = 0; i < steps; i++)
      dataBeginIdx = molfile.indexOf(symbol, dataBeginIdx) + 1;
  };
  let dataBeginIdx = 0;
  const headerLinesNum = 4;
  //jump to atoms block
  moveCursorToIdx(headerLinesNum, '\n');
  if (MolfileHandler.isMolfileV2K(molfile)) {
    const zCoordShift = 20;
    const coordDigitsNum = 10;
    for (let i = 0; i < numAtoms; i++) {
      if (parseFloat(molfile.substring(dataBeginIdx + zCoordShift, dataBeginIdx + zCoordShift + coordDigitsNum)))
        return true;
      //go to next row
      dataBeginIdx = molfile.indexOf('\n', dataBeginIdx) + 1;
    }
  } else if (MolfileHandler.isMolfileV3K(molfile)) {
    //jump to atoms block
    const atomCountsLinesNum = 3;
    moveCursorToIdx(atomCountsLinesNum, '\n');
    for (let i = 0; i < numAtoms; i++) {
      //go to z coordinate start
      moveCursorToIdx(7, ' ');
      //get z coordinate end
      const zCoordEnd = molfile.indexOf(' ', dataBeginIdx);
      if (parseFloat(molfile.substring(dataBeginIdx, zCoordEnd)))
        return true;
      //go to next row
      dataBeginIdx = molfile.indexOf('\n', dataBeginIdx) + 1;
    }
  }
  return false;
}

function mergeSubstructs(substructList: (ISubstruct | undefined)[]): ISubstruct | undefined {
  // Filter out undefined entries from providers that don't apply to this row.
  const defined = substructList.filter((s): s is ISubstruct => s !== undefined);
  if (defined.length === 0) return undefined;
  if (defined.length === 1) return defined[0];
  // Delegate to chem-meta: unions atoms/bonds/colors; later providers win on collisions.
  return mergeSubstructsLib(defined);
}
