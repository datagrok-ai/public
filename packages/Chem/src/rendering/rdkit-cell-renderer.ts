/* eslint-disable max-len */
/**
 * RDKit-based molecule cell renderer.
 * */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMonomerHover, getSubstructProviders, mergeSubstructs as mergeSubstructsLib,
  addSubstructProvider, ISubstructProvider}
  from '@datagrok-libraries/chem-meta/src/types';
import {ChemTags, ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

import {
  ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG,
  ALIGN_BY_SCAFFOLD_TAG, FILTER_SCAFFOLD_TAG,
  CHEM_INTERACTIVE_SELECTION_EVENT, CHEM_ATOM_PICKER_LINKED_COL,
  CHEM_MOL3D_HOVER_EVENT,
  FIXED_SCALE_TAG,
  HIGHLIGHT_BY_SCAFFOLD_COL, HIGHLIGHT_BY_SCAFFOLD_COL_SYNC,
  HIGHLIGHT_BY_SCAFFOLD_TAG, MIN_MOL_IMAGE_SIZE, PARENT_MOL_COL,
  REGENERATE_COORDS, REGENERATE_COORDS_SYNC,
  SCAFFOLD_COL, SCAFFOLD_COL_SYNC, SCAFFOLD_TREE_HIGHLIGHT,
  getSyncTag,
} from '../constants';
import {hexToPercentRgb} from '../utils/chem-common';
import {_rdKitModule, drawErrorCross, drawRdKitMoleculeToOffscreenCanvas,
  RDKIT_COMMON_RENDER_OPTS} from '../utils/chem-common-rdkit';
import {IMolContext, getMolSafe} from '../utils/mol-creation_rdkit';
import {getGridCellColTemp} from '../utils/ui-utils';
import {errInfo} from '../utils/err-info';
import {mapAtomIndices2Dto3D, AtomIndexMapping} from '../utils/atom-index-mapper';
import {
  findNearestAtom, computeSelectedBonds, extractAtomPositionsFromSvg,
} from '../utils/chem-atom-picker-utils';
import {handleMol3DHoverEvent, Mol3DHoverRendererDeps} from './mol3d-hover-handler';

import {_package} from '../package';

/** Shape of the cross-package atom selection event. */
interface ChemSelectionEvent {
  column: DG.Column;
  rowIdx: number;
  atoms: number[];
  persistent?: boolean;
  clearAll?: boolean;
  mapping3D?: AtomIndexMapping | null;
  mol3DColumnName?: string;
}

/** Substruct provider with atom-picker metadata fields. */
interface AtomPickerProvider extends ISubstructProvider {
  __atomPicker?: boolean;
  __rowIdx?: number;
  __atoms?: Set<number>;
}

interface IMolRenderingInfo {
  //mol: RDMol | null; // null when molString is invalid?
  molCtx: IMolContext;
  substruct: ISubstruct;
  molString: string;
}

/** Per-cell layout info for atom hit-testing. Cached by "molString|WxH". */
interface CellInteractiveInfo {
  /** Atom index → center in cell-local CSS pixels. */
  positions: Map<number, {x: number, y: number}>;
  /** Bond index → [atomA, atomB] pair from RDKit SVG class annotations. */
  bondAtoms: Map<number, [number, number]>;
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
  alighByFirstSubtruct: boolean,
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
  rdKitModule: RDModule;
  canvasCounter: number;
  molCache: DG.LruCache<String, IMolRenderingInfo> = new DG.LruCache<String, IMolRenderingInfo>();
  rendersCache: DG.LruCache<String, ImageData> = new DG.LruCache<String, ImageData>();
  sortedScaffoldsCache: DG.LruCache<String, IColoredScaffold[]> = new DG.LruCache<String, IColoredScaffold[]>();
  canvasReused: OffscreenCanvas;

  // ---- in-grid atom picker state ------------------------------------------
  /** Per-cell layout info cached by "molString|w|h" — atom positions + bond connectivity. */
  atomPositionsCache: DG.LruCache<string, CellInteractiveInfo> =
    new DG.LruCache<string, CellInteractiveInfo>();
  /** Last atom processed by hover — dedup guard so we don't re-fire on every pixel. */
  private _lastHoveredAtom: {col: string, rowIdx: number, atomIdx: number, erase?: boolean} | null = null;
  /** Transient preview atom (no modifier hover). Cleared on mouse-leave or shift. */
  private _previewAtomIdx: number | null = null;
  /** Row the preview is on. */
  private _previewRowIdx: number = -1;
  /** True when the active preview was set by a 3D→2D hover event rather than a 2D hover.
   *  A 3D-sourced preview must survive 2D mousemoves and is cleared only by the
   *  corresponding 3D "cursor left the atom" event (`atom3DSerial === null`). */
  private _previewFrom3D: boolean = false;

  constructor(rdKitModule: RDModule) {
    super();
    this.rdKitModule = rdKitModule;
    this.canvasCounter = 0;
    this.canvasReused = new OffscreenCanvas(this.defaultWidth, this.defaultHeight);
    this.molCache.onItemEvicted = function(obj: { [_: string]: any }) {
      obj.mol?.delete();
    };
  }

  // -- Interactive atom picker -----------------------------------------------

  private _getProviders(col: DG.Column): AtomPickerProvider[] {
    return (col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as AtomPickerProvider[];
  }

  private _getProviderForRow(col: DG.Column, rowIdx: number): AtomPickerProvider | undefined {
    return this._getProviders(col).find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
  }

  /** Builds a selection event, attaching 3D mapping when a linked Mol3D column exists. */
  private _buildSelectionEvent(
    col: DG.Column, rowIdx: number, atoms: number[], persistent: boolean = true,
  ): ChemSelectionEvent {
    const event: ChemSelectionEvent = {column: col, rowIdx, atoms, persistent};
    try {
      const df = col.dataFrame;
      if (df) {
        const linkedColName = this._getLinkedMol3DColName(col);
        const mol3DCol = linkedColName ? df.col(linkedColName) : null;
        if (mol3DCol && rowIdx >= 0) {
          const smiles2D = col.get(rowIdx);
          const pose3D = mol3DCol.get(rowIdx);
          if (smiles2D && pose3D) {
            const mapping = mapAtomIndices2Dto3D(this.rdKitModule, smiles2D, pose3D);
            if (mapping) {
              event.mapping3D = mapping;
              event.mol3DColumnName = mol3DCol.name;
            }
          }
        }
      }
    } catch {
      /* mapping failed */
    }
    return event;
  }

  /** Clears the render cache so the 2D cell fully redraws.
   *  DG.LruCache has no public `.clear()`, so we zero the backing arrays
   *  via a narrow cast to unblock GC of canvas-backed ImageData entries. */
  private _clearRendersCache(): void {
    const old = this.rendersCache as unknown as {
      K?: unknown[]; V?: unknown[]; items?: Record<string, number>;
      size?: number; head?: number; tail?: number;
    };
    old.K?.fill(undefined);
    old.V?.fill(undefined);
    old.items = {};
    old.size = 0;
    this.rendersCache.onItemEvicted = null;
    this.rendersCache = new DG.LruCache<String, ImageData>();
  }

  private _isSameAtom(colName: string, rowIdx: number, atomIdx: number, erase?: boolean): boolean {
    const last = this._lastHoveredAtom;
    if (!last) return false;
    return last.col === colName &&
      last.rowIdx === rowIdx &&
      last.atomIdx === atomIdx &&
      (erase === undefined || last.erase === erase);
  }

  /** Returns `true` (same atom, caller should bail) or records it and returns `false`. */
  private _trackHoveredAtom(
    colName: string, rowIdx: number, atomIdx: number, erase?: boolean,
  ): boolean {
    if (this._isSameAtom(colName, rowIdx, atomIdx, erase)) return true;
    this._lastHoveredAtom = {col: colName, rowIdx, atomIdx, erase};
    return false;
  }

  /** Returns true when the picker is active for the column (a Mol3D link exists). */
  private _isPickerActive(col: DG.Column): boolean {
    return !!this._getLinkedMol3DColName(col);
  }

  /** Returns the linked Mol3D column name, preferring persistent tags over session temp. */
  private _getLinkedMol3DColName(col: DG.Column): string | null {
    const fromTags = col.tags[CHEM_ATOM_PICKER_LINKED_COL];
    if (fromTags) return fromTags;
    const fromTemp = col.temp?.[CHEM_ATOM_PICKER_LINKED_COL];
    return typeof fromTemp === 'string' ? fromTemp : null;
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

  _fetchMolGetOrCreate(molString: string, scaffolds: IColoredScaffold[], molRegenerateCoords: boolean,
    details: object = {}, alignByFirstSubstr: boolean): IMolRenderingInfo {
    let molCtx: IMolContext;
    let mol = null;
    let substruct: ISubstruct = {};
    if ((details as any).isSubstructure) {
      const mappedDummiesAreRGroups = (details as any).mappedDummiesAreRGroups || false;
      if (molString.includes(' H ') || molString.includes('V3000')) {
        molCtx = getMolSafe(molString, {mergeQueryHs: true, mappedDummiesAreRGroups}, _rdKitModule);
        mol = molCtx.mol;
      } else {
        try {
          if (!molString.includes('\n') && molString.length > 5000)
            throw new Error('Invalid molecule string'); // do not attempt to parse very long SMILES, will cause MOB.
          mol = this.rdKitModule.get_qmol(molString);
          mol.convert_to_aromatic_form();
        } catch (e) {
          if (mol) {
            mol.delete();
            mol = null;
          }
        }
        molCtx = {mol: mol, kekulize: false, isQMol: true, useMolBlockWedging: false};
      }
    } else {
      molCtx = getMolSafe(molString, details, _rdKitModule);
      mol = molCtx.mol;
    }

    if (mol !== null) {
      try {
        let molHasOwnCoords = (mol.has_coords() > 0);
        let molHasRebuiltCoords = false;
        const scaffoldIsMolBlock = scaffolds.length && alignByFirstSubstr ?
          DG.chem.isMolBlock(scaffolds[0].molecule) : null;
        const alignedByFirstSubstr = scaffoldIsMolBlock && alignByFirstSubstr;
        const {haveReferenceSmarts, parentMolScaffoldMolString} = (details as any);
        if (alignedByFirstSubstr) {
          const rdKitScaffoldMolCtx = this._fetchMol(scaffolds[0].molecule,
            parentMolScaffoldMolString ? [{molecule: parentMolScaffoldMolString}] : [], molRegenerateCoords,
            false, {mergeQueryHs: true, isSubstructure: !parentMolScaffoldMolString}, false).molCtx;
          const rdKitScaffoldMol = rdKitScaffoldMolCtx.mol;
          if (rdKitScaffoldMol) {
            rdKitScaffoldMol.normalize_depiction(0);
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
                mol.generate_aligned_coords(rdKitScaffoldMol, JSON.stringify(alignOpts)) :
                mol.get_substruct_match(mol!);
            } catch {
              // exceptions should not be thrown anymore by RDKit, but let's play safe
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
              {mergeQueryHs: true, isSubstructure: true}, false).molCtx.mol;
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
        console.error(
          'In _fetchMolGetOrCreate: RDKit crashed, possibly a malformed molString molecule: `' + molString + '`');
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
    };
  }

  _addAtomsOrBonds(fromAtomsOrBonds: number[], toAtomsOrBonds: number[]): void {
    for (let j = 0; j < fromAtomsOrBonds.length; j++) {
      if (!toAtomsOrBonds?.includes(fromAtomsOrBonds[j]))
        toAtomsOrBonds?.push(fromAtomsOrBonds[j]);
    }
  }

  _fetchMol(molString: string, scaffolds: IColoredScaffold[], molRegenerateCoords: boolean,
    scaffoldRegenerateCoords: boolean, details: object = {}, alignByFirstSubstructure: boolean): IMolRenderingInfo {
    const name = molString + ' || ' + JSON.stringify(scaffolds) + ' || ' +
      molRegenerateCoords + ' || ' + scaffoldRegenerateCoords + ' || ' +
      ((details as any).isSubstructure || false).toString() +
      ((details as any).haveReferenceSmarts || false).toString();
    return this.molCache.getOrCreate(name, (_: any) =>
      this._fetchMolGetOrCreate(molString, scaffolds, molRegenerateCoords, details, alignByFirstSubstructure));
  }

  _rendererGetOrCreate(
    width: number, height: number, molString: string, scaffolds: IColoredScaffold[],
    molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean,
    alignByFirstSubstructure: boolean, details: object = {}, substructureObj?: ISubstruct, renderingOptions?: any): ImageData {
    const fetchMolObj: IMolRenderingInfo =
      this._fetchMol(molString, scaffolds, molRegenerateCoords,
        scaffoldRegenerateCoords, details, alignByFirstSubstructure);
    const rdKitMolCtx = fetchMolObj.molCtx;
    const rdKitMol = rdKitMolCtx.mol;//fetchMolObj.mol;
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

    const canvas = this.ensureCanvasSize(width, height);//new OffscreenCanvas(width, height);
    const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
    this.canvasCounter++;
    if (rdKitMol != null) {
      drawRdKitMoleculeToOffscreenCanvas(rdKitMolCtx, width, height, canvas,
        scaffolds.length ? substructureObj ? newSubstruct : substruct : substructureObj ?? null, renderingOptions);
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
    const name = width + ' || ' + height + ' || ' +
      molString + ' || ' + JSON.stringify(scaffolds) + ' || ' +
      molRegenerateCoords + ' || ' + scaffoldRegenerateCoords + ' || ' +
      ((details as any).haveReferenceSmarts || false).toString() + ' || ' + JSON.stringify(substructureObj) + ' || ' + JSON.stringify(renderOptions);

    return this.rendersCache.getOrCreate(name, (_: any) => this._rendererGetOrCreate(width, height,
      molString, scaffolds, molRegenerateCoords, scaffoldRegenerateCoords,
      alignByFirstSubstructure, details, substructureObj, renderOptions));
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
    } else {
      //const image = offscreenCanvas.getContext('2d')!.getImageData(0, 0, w, h);
      onscreenCanvas.getContext('2d', {willReadFrequently: true})!.putImageData(imageData, x, y);
    }
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

  // Registered once on first render; grid rebuilds get re-attached automatically.
  private _globalListenersAttached = false;
  private _wiredCanvases: WeakSet<HTMLCanvasElement> = new WeakSet();

  render(g: any, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    if (!this._globalListenersAttached) {
      this._globalListenersAttached = true;
      // Subscribe once to the 3D→2D bridge fired by BiostructureViewer's Molstar viewer.
      grok.events.onCustomEvent(CHEM_MOL3D_HOVER_EVENT).subscribe(
        (args: unknown) => this._onMol3DHoverEvent(args));
    }
    // Wire shift-drag once per canvas instance (idempotent via _wiredCanvases).
    this._attachGridCanvasShiftListener(gridCell);

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
    } else if (highlightInfo.scaffolds && highlightInfo.alighByFirstSubtruct) {
      this._drawMolecule(x, y, w, h, g.canvas,
        molString, highlightInfo.scaffolds, false, false, cellStyle, highlightInfo.alighByFirstSubtruct, undefined, undefined, renderOpts);
    } else
      this.highlightByScaffoldCol(g, x, y, w, h, gridCell, cellStyle, colTemp, molString, highlightInfo.scaffolds, renderOpts);

    // Pre-warm the atom-positions cache on idle so hover responds instantly.
    if (this._isPickerActive(gridCell.cell.column)) {
      const cssW = w / r;
      const cssH = h / r;
      const dpr = window.devicePixelRatio || 1;
      const cacheKey = molString + '|' + Math.round(cssW * dpr) + 'x' + Math.round(cssH * dpr);
      if (!this.atomPositionsCache.has(cacheKey)) {
        (typeof requestIdleCallback === 'function' ? requestIdleCallback : setTimeout)(
          () => { this._getCellAtomPositions(molString, cssW, cssH); },
        );
      }
    }
  }

  getHighlightTagInfo(colTemp: any, gridCell: DG.GridCell): IHighlightTagInfo {
    const filter = this._initScaffoldArray(colTemp, FILTER_SCAFFOLD_TAG, true); //expected molBlock
    const align = (this._initScaffoldString(gridCell.cell.column, ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG) ?? [])
      // temporary concatination to support both tags
      .concat(this._initScaffoldString(gridCell.cell.column, ALIGN_BY_SCAFFOLD_TAG) ?? []);
    const highlight = this._initScaffoldArray(gridCell.cell.column, HIGHLIGHT_BY_SCAFFOLD_TAG);
    const scaffoldTreeHighlight = this._initScaffoldArray(gridCell.cell.column, SCAFFOLD_TREE_HIGHLIGHT);
    const alignByStructure = !!(filter.length && filter[0].align || align.length);
    const scaffolds = filter.concat(align).concat(scaffoldTreeHighlight).concat(highlight);
    return {scaffolds: scaffolds?.length ? scaffolds : undefined, alighByFirstSubtruct: alignByStructure};
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

  // -- Hover-paint atom picker ------------------------------------------------
  // Shift+hover paints atoms; Ctrl+Shift erases. No click needed — grid focus
  // stays on the current cell so the context panel remains open. Escape clears.
  // Wired via DG.GridCellRenderer overrides (no document-wide listeners).

  /** Resolves the atom under the cursor. Coordinates are grid-canvas-relative
   *  (e.offsetX/Y) and translated into cell-local space via gridCell.bounds. */
  private _hitTestAtomInCell(gridCell: DG.GridCell, e: MouseEvent): {
    cellInfo: CellInteractiveInfo; nearest: number;
  } | null {
    if (gridCell.tableRowIndex == null || gridCell.tableRowIndex < 0) return null;
    const molString: string = gridCell.cell.value;
    if (!molString || DG.chem.Sketcher.isEmptyMolfile(molString)) return null;

    const bounds = gridCell.bounds;
    const pointerCellX = e.offsetX - bounds.left;
    const pointerCellY = e.offsetY - bounds.top;
    if (pointerCellX < 0 || pointerCellY < 0 ||
        pointerCellX > bounds.width || pointerCellY > bounds.height) return null;

    const cellInfo = this._getCellAtomPositions(molString, bounds.width, bounds.height);
    if (!cellInfo || cellInfo.positions.size === 0) return null;

    const nearest = findNearestAtom(cellInfo.positions, pointerCellX, pointerCellY);
    if (nearest === null) return null;

    return {cellInfo, nearest};
  }

  /**
   * No-modifier hover handler — highlights ONLY the single atom under the
   * cursor (preview mode). Moving to a different atom highlights that one
   * instead. Moving off all atoms clears the preview.
   *
   * `onMouseLeave` handles preview-clear when the cursor leaves the cell
   * toward a non-molecule column.
   */
  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const col = gridCell.tableColumn;
    if (!col || col.semType !== DG.SEMTYPE.MOLECULE || !this._isPickerActive(col))
      return;

    const hit = this._hitTestAtomInCell(gridCell, e);
    const isErase = e.shiftKey && e.ctrlKey;
    const isPaint = e.shiftKey && !e.ctrlKey;

    // -- Shift / Ctrl+Shift: persistent paint/erase mode --
    if (isPaint || isErase) {
      this._previewAtomIdx = null;
      if (!hit) return;
      const rowIdx = gridCell.tableRowIndex!;
      if (this._trackHoveredAtom(col.name, rowIdx, hit.nearest, isErase)) return;
      if (isErase)
        this._removeAtomFromRow(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
      else
        this._addAtomToRow(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
      this._clearRendersCache();
      gridCell.grid.invalidate();
      return;
    }

    // -- No modifier: preview mode --
    if (!hit) {
      // A 3D-sourced preview must survive 2D mousemoves; it is cleared only by
      // the corresponding 3D "cursor left the atom" event (atom3DSerial === null).
      if (this._previewAtomIdx !== null && !this._previewFrom3D)
        this._removePreviewAtom();
      this._lastHoveredAtom = null;
      return;
    }

    const rowIdx = gridCell.tableRowIndex!;
    if (this._trackHoveredAtom(col.name, rowIdx, hit.nearest)) return;

    this._previewFrom3D = false; // 2D hover takes ownership
    this._setPreviewAtom(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
    gridCell.grid.invalidate();
  }

  /** Clears the 2D-sourced preview when the cursor leaves a molecule cell. */
  override onMouseLeave(gridCell: DG.GridCell, _e: MouseEvent): void {
    const col = gridCell.tableColumn;
    if (!col || col.semType !== DG.SEMTYPE.MOLECULE || !this._isPickerActive(col))
      return;
    if (this._previewAtomIdx !== null && !this._previewFrom3D)
      this._removePreviewAtom();
    this._lastHoveredAtom = null;
  }

  // Cast is safe: Mol3DHoverRendererDeps mirrors the private members by name.
  private _onMol3DHoverEvent(args: unknown): void {
    handleMol3DHoverEvent(this as unknown as Mol3DHoverRendererDeps, args);
  }

  /**
   * Attaches capture-phase listeners scoped to the grid's own DOM, idempotent
   * per canvas instance (tracked via `_wiredCanvases`):
   * - mousemove on the canvas (gated on `e.shiftKey`): handles shift-drag
   *   paint/erase — Datagrok captures native row-selection at Dart level and
   *   bypasses the renderer's `onMouseMove` for shift-moves.
   * - keydown on `grid.root` (Escape): keyboard events go to the focused
   *   element and bubble up; `grid.root` catches Escape from any descendant,
   *   while the renderer's `onKeyDown` only fires when a molecule cell is focused.
   */
  private _attachGridCanvasShiftListener(gridCell: DG.GridCell): void {
    let canvas: HTMLCanvasElement | undefined;
    try { canvas = gridCell.grid.overlay ?? gridCell.grid.canvas; } catch { return; }
    if (!canvas || this._wiredCanvases.has(canvas)) return;
    this._wiredCanvases.add(canvas);
    const grid = gridCell.grid;
    canvas.addEventListener('mousemove', (e: MouseEvent) => {
      if (e.shiftKey) this._onShiftDragMouseMove(e);
    }, true);
    let gridRoot: HTMLElement | null = null;
    try { gridRoot = grid.root; } catch { /* fall through */ }
    if (gridRoot) {
      gridRoot.addEventListener('keydown', (e: KeyboardEvent) => {
        if (e.key === 'Escape') this._clearPickerSelection(grid);
      }, true);
    }
  }

  /** Escape when a molecule cell is focused — Datagrok routes the keydown here
   *  at Dart level and it does NOT bubble to grid.root. The grid-root listener
   *  in `_attachGridCanvasShiftListener` covers non-molecule focused cells. */
  override onKeyDown(gridCell: DG.GridCell, e: KeyboardEvent): void {
    if (e.key === 'Escape')
      this._clearPickerSelection(gridCell.grid);
  }

  /** Clears all picker providers and broadcasts a persistent clear-all event
   *  so 3D viewers drop their cached overpaint. Called from both Escape paths. */
  private _clearPickerSelection(grid: DG.Grid): void {
    const df = grid?.dataFrame;
    if (!df) return;
    const molCols = df.columns.toList().filter(
      (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE && this._isPickerActive(c));
    if (molCols.length === 0) return;
    for (const col of molCols) {
      col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = this._getProviders(col)
        .filter((p) => !p.__atomPicker);
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
        column: col, rowIdx: -1, atoms: [], persistent: true, clearAll: true,
      });
    }
    this._lastHoveredAtom = null;
    this._previewAtomIdx = null;
    this._previewFrom3D = false;
    this._clearRendersCache();
    grid.invalidate();
  }

  /** Shift-drag handler from the canvas capture listener. Re-does the grid
   *  hit-test from page coords because Datagrok captures shift-mousemove
   *  at the Dart level and does not forward it to the renderer's onMouseMove. */
  private _onShiftDragMouseMove(e: MouseEvent): void {
    const isErase = e.shiftKey && e.ctrlKey;
    const isPaint = e.shiftKey && !e.ctrlKey;
    if (!isPaint && !isErase) return;

    const grid = grok.shell.tv?.grid;
    if (!grid) return;
    let gridRoot: HTMLElement | null = null;
    try { gridRoot = grid.root; } catch { return; }
    if (!gridRoot) return;

    const gridRect = gridRoot.getBoundingClientRect();
    const localX = e.clientX - gridRect.left;
    const localY = e.clientY - gridRect.top;
    if (localX < 0 || localY < 0 || localX > gridRect.width || localY > gridRect.height) return;

    let gridCell: DG.GridCell | null = null;
    try { gridCell = grid.hitTest(localX, localY); } catch { return; }
    const col = gridCell?.tableColumn;
    if (!gridCell || !col || col.semType !== DG.SEMTYPE.MOLECULE) return;
    if (!this._isPickerActive(col)) return;

    const molString: string = gridCell.cell.value;
    if (!molString || DG.chem.Sketcher.isEmptyMolfile(molString)) return;

    const bounds = gridCell.bounds;
    const pointerCellX = e.clientX - (gridRect.left + bounds.left);
    const pointerCellY = e.clientY - (gridRect.top + bounds.top);
    if (pointerCellX < 0 || pointerCellY < 0 ||
        pointerCellX > bounds.width || pointerCellY > bounds.height) return;

    const cellInfo = this._getCellAtomPositions(molString, bounds.width, bounds.height);
    if (!cellInfo || cellInfo.positions.size === 0) return;

    const nearest = findNearestAtom(cellInfo.positions, pointerCellX, pointerCellY);
    if (nearest === null) return;

    this._previewAtomIdx = null;
    const rowIdx = gridCell.tableRowIndex!;
    if (this._trackHoveredAtom(col.name, rowIdx, nearest, isErase)) return;

    if (isErase)
      this._removeAtomFromRow(col, rowIdx, nearest, cellInfo.bondAtoms);
    else
      this._addAtomToRow(col, rowIdx, nearest, cellInfo.bondAtoms);

    this._clearRendersCache();
    grid.invalidate();
  }

  private _getPersistentAtoms(col: DG.Column, rowIdx: number): Set<number> {
    const prov = this._getProviderForRow(col, rowIdx);
    return new Set<number>(prov?.__atoms ?? []);
  }

  private _addAtomToRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    const current = this._getPersistentAtoms(col, rowIdx);
    if (current.has(atomIdx)) return; // already selected
    current.add(atomIdx);
    this._updateRowSelection(col, rowIdx, [...current], {add: true}, bondAtoms);
  }

  private _removeAtomFromRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    const prov = this._getProviderForRow(col, rowIdx);
    if (!prov?.__atoms?.has(atomIdx)) return;

    prov.__atoms.delete(atomIdx);

    if (prov.__atoms.size === 0) {
      col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = this._getProviders(col).filter(
        (p) => !p.__atomPicker || p.__rowIdx !== rowIdx);
    } else {
      // Rebuild the provider with the reduced atom set.
      const pickColor: number[] = [1.0, 1.0, 0.0, 1.0];
      const atomsArr = [...prov.__atoms];
      const highlightAtomColors: {[k: number]: number[]} = {};
      for (const a of atomsArr) highlightAtomColors[a] = pickColor;
      const {bondsArr, highlightBondColors} =
        computeSelectedBonds(prov.__atoms, bondAtoms, pickColor);

      prov.getSubstruct = (ridx) => {
        if (!this._isPickerActive(col)) return undefined;
        if (ridx !== rowIdx) return undefined;
        return {atoms: atomsArr, bonds: bondsArr, highlightAtomColors, highlightBondColors};
      };
    }

    this._clearRendersCache();
    const atoms = prov.__atoms.size > 0 ? [...prov.__atoms] : [];
    grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT,
      this._buildSelectionEvent(col, rowIdx, atoms));
  }

  /** Shows persistent atoms + the hovered atom as a non-persistent preview.
   *  Restores `__atoms` to the persistent-only set afterward so the preview
   *  atom doesn't leak into storage. */
  private _setPreviewAtom(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    this._previewAtomIdx = atomIdx;
    this._previewRowIdx = rowIdx;
    const persistent = this._getPersistentAtoms(col, rowIdx);
    const combined = new Set(persistent);
    combined.add(atomIdx);
    this._updateRowSelection(col, rowIdx, [...combined], {add: true}, bondAtoms, true, false);
    const prov = this._getProviderForRow(col, rowIdx);
    if (prov) prov.__atoms = persistent;
  }

  /** Restores display to the persistent selection, dropping the preview atom. */
  private _removePreviewAtom(): void {
    const prevAtom = this._previewAtomIdx;
    const prevRow = this._previewRowIdx;
    this._previewAtomIdx = null;
    this._previewFrom3D = false;
    if (prevAtom === null) return;

    const grid = grok.shell.tv?.grid;
    if (!grid) return;
    const df = grid.dataFrame;
    if (!df) return;
    const molCol = df.columns.toList().find(
      (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE && this._isPickerActive(c));
    if (!molCol) return;

    const persistent = this._getPersistentAtoms(molCol, prevRow);
    if (persistent.size > 0) {
      // Restore just the persistent atoms in 2D and 3D.
      // Fire event as non-persistent so the cache keeps the stable Shift-paint set.
      const cellInfo = this._getCellAtomPositions(
        molCol.get(prevRow), 100, 100);
      const bondAtoms = cellInfo?.bondAtoms ?? new Map();
      this._updateRowSelection(molCol, prevRow, [...persistent], {add: true}, bondAtoms, true, false);
    } else {
      // No persistent atoms — clear everything.
      molCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = this._getProviders(molCol)
        .filter((p) => !p.__atomPicker || p.__rowIdx !== prevRow);
      // Fire event to clear 3D highlights too (non-persistent).
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
        column: molCol, rowIdx: prevRow, atoms: [], persistent: false,
      });
    }
    grid.invalidate();
  }


  /**
   * Writes an atom-picker ISubstructProvider for `rowIdx`, replacing the
   * previous one for that row. Other rows' providers are preserved.
   */
  private _updateRowSelection(col: DG.Column, rowIdx: number, boxed: number[],
    modifiers: {add: boolean}, bondAtoms: Map<number, [number, number]>,
    fire3DEvent: boolean = true, persistent: boolean = true): void {
    const prior = this._getProviderForRow(col, rowIdx);
    const current = new Set<number>(prior?.__atoms ?? []);

    // No modifier → replace existing selection with boxed atoms.
    // Shift (add) → keep existing and union with boxed atoms.
    if (!modifiers.add) current.clear();
    for (const a of boxed) current.add(a);

    // Yellow (#FFFF00): avoids collision with Datagrok's green row-selection tint.
    const pickColor: number[] = [1.0, 1.0, 0.0, 1.0];
    const atomsArr = [...current];
    const highlightAtomColors: {[k: number]: number[]} = {};
    for (const a of atomsArr) highlightAtomColors[a] = pickColor;
    const {bondsArr, highlightBondColors} =
      computeSelectedBonds(current, bondAtoms, pickColor);

    const provider: AtomPickerProvider = {
      getSubstruct: (ridx) => {
        if (!this._isPickerActive(col)) return undefined;
        if (ridx !== rowIdx) return undefined;
        return {atoms: atomsArr, bonds: bondsArr, highlightAtomColors, highlightBondColors};
      },
    };
    provider.__atomPicker = true;
    provider.__rowIdx = rowIdx;
    provider.__atoms = current;

    // Replace only the prior provider for this row; keep picks on other rows.
    col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = this._getProviders(col).filter(
      (p) => !p.__atomPicker || p.__rowIdx !== rowIdx);
    addSubstructProvider(col.temp, provider);

    if (fire3DEvent) {
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT,
        this._buildSelectionEvent(col, rowIdx, atomsArr, persistent));
    }
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
      // Mol is owned by molCache — DO NOT delete it.
      const molRenderingInfo = this._fetchMol(molString, [], false, false, {}, false);
      const mol = molRenderingInfo.molCtx.mol;
      if (!mol || !mol.is_valid()) return null;

      const details: {[k: string]: any} = {};
      for (const k of Object.keys(RDKIT_COMMON_RENDER_OPTS))
        details[k] = RDKIT_COMMON_RENDER_OPTS[k];
      details.width = w;
      details.height = h;
      details.atoms = [];
      details.bonds = [];
      details.highlightAtomColors = {};
      details.highlightBondColors = {};
      // Mirror kekulize / molBlockWedging from drawRdKitMoleculeToOffscreenCanvas
      // so the SVG layout matches the canvas layout exactly.
      if (!molRenderingInfo.molCtx.kekulize) details.kekulize = false;
      if (molRenderingInfo.molCtx.useMolBlockWedging) {
        details.useMolBlockWedging = true;
        details.wedgeBonds = false;
        details.addChiralHs = false;
      }

      const svgString = mol.get_svg_with_highlights(JSON.stringify(details));
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
