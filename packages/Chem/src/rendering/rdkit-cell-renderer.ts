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
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

import {
  ALIGN_BY_SCAFFOLD_LAYOUT_PERSISTED_TAG,
  ALIGN_BY_SCAFFOLD_TAG, FILTER_SCAFFOLD_TAG,
  CHEM_INTERACTIVE_SELECTION_EVENT,
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
import {mapAtomIndices2Dto3D} from '../utils/atom-index-mapper';

import {_package} from '../package';

interface IMolRenderingInfo {
  //mol: RDMol | null; // null when molString is invalid?
  molCtx: IMolContext;
  substruct: ISubstruct;
  molString: string;
}

/** Layout info needed to hit-test a molecule cell interactively. Computed
 *  once per (molString, width, height) by rendering the molecule to SVG at
 *  the cell's dimensions and parsing atom + bond class annotations. */
interface CellInteractiveInfo {
  /** Atom index → center position in cell-local CSS pixels. */
  positions: Map<number, {x: number, y: number}>;
  /** Bond index → pair of atom indices that the bond connects. Extracted
   *  from RDKit SVG bond paths with classes like `bond-K atom-A atom-B`. */
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
    // Forward per-atom/per-bond labels and annotation notes from a temp
    // substructure into the main one so providers can layer annotations
    // (e.g. pKa values, R-group labels, contribution scores) on top of the
    // existing highlight pipeline.
    //
    // Cast through `any` because the `atomLabels` / `atomNotes` /
    // `bondNotes` fields were added to `ISubstruct` in this branch's
    // chem-meta change but the published `@datagrok-libraries/chem-meta`
    // package on npm is still the older 1.x without them. CI fetches the
    // npm version, so the type check fails on a fresh `npm install`.
    // Casts erase the type check; runtime is unaffected because the
    // properties are optional and the runtime JS doesn't care about
    // TypeScript types. Once chem-meta is republished with the new
    // fields, these casts can be removed.
    if (tempSubstr) {
      const t = tempSubstr as any;
      const m = mainSubstr as any;
      if (t.atomLabels) {
        m.atomLabels ??= {};
        Object.assign(m.atomLabels, t.atomLabels);
      }
      if (t.atomNotes) {
        m.atomNotes ??= {};
        Object.assign(m.atomNotes, t.atomNotes);
      }
      if (t.bondNotes) {
        m.bondNotes ??= {};
        Object.assign(m.bondNotes, t.bondNotes);
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

  // ---- in-grid box picker state ------------------------------------------
  /** Per-cell interactive layout info, cached by "molString|w|h" so
   *  successive drags on the same cell are O(1). Holds both atom positions
   *  (for box hit-testing) and bond atom-pairs (for deriving which bonds
   *  are "selected" as a consequence of the user's atom selection — a
   *  bond is highlighted iff both of its atoms are in the picked set). */
  atomPositionsCache: DG.LruCache<string, CellInteractiveInfo> =
    new DG.LruCache<string, CellInteractiveInfo>();
  /** Whether the single document-level mousedown listener has been
   *  registered. We use ONE global listener (capture phase) instead of
   *  per-canvas listeners because:
   *  (a) Datagrok's grid stops dispatching `onMouseDown` after a table
   *      switch, and
   *  (b) the canvas from `g.canvas` (render context) is NOT the element
   *      that receives mouse events — it's a rendering-only canvas.
   *  The document listener catches ALL mousedowns, finds the active grid
   *  via `grok.shell.tv.grid`, and hit-tests to determine if a molecule
   *  cell was clicked. */
  private _documentListenerAttached = false;
  /** Drag state shared across the global mousemove/mouseup listeners.
   *  Null when no drag is active.
   *
   *  We use capture-phase document listeners (`addEventListener('mouseup',
   *  handler, true)`) rather than grid renderer lifecycle methods because
   *  Datagrok's grid event pipeline can swallow `onMouseUp` before it reaches
   *  the renderer. Capture-phase fires at document BEFORE the event walks
   *  down to its target, so any `stopPropagation` on the canvas cannot stop
   *  it from firing. */
  /** Tracks the last atom processed by hover so we don't re-fire on
   *  every mousemove pixel within the same atom's radius. */
  private _lastHoveredAtom: {col: string, rowIdx: number, atomIdx: number, erase?: boolean} | null = null;
  /** The currently previewed atom (hover without Alt). Shown alongside
   *  the persistent selection but removed when the cursor moves away. */
  private _previewAtomIdx: number | null = null;
  /** The row the preview is on. */
  private _previewRowIdx: number = -1;

  constructor(rdKitModule: RDModule) {
    super();
    this.rdKitModule = rdKitModule;
    this.canvasCounter = 0;
    this.canvasReused = new OffscreenCanvas(this.defaultWidth, this.defaultHeight);
    this.molCache.onItemEvicted = function(obj: { [_: string]: any }) {
      obj.mol?.delete();
    };
  }

  // -- Interactive atom picker auto-detection --------------------------------

  /** Cache for _isPickerActive to avoid scanning columns on every event. */
  private _pickerActiveCache: {dfId: string, active: boolean} | null = null;

  /** Returns true when the interactive atom picker should be active.
   *  Auto-detects based on whether the DataFrame has associated
   *  Molecule3D (docking poses) or HELM columns — no manual toggle needed. */
  private _isPickerActive(col: DG.Column): boolean {
    const df = col.dataFrame;
    if (!df) return false;
    if (this._pickerActiveCache?.dfId === df.id)
      return this._pickerActiveCache.active;
    const cols = df.columns.toList();
    const active =
      cols.some((c) => c.semType === DG.SEMTYPE.MOLECULE3D) ||
      cols.some((c) => c.semType === (DG.SEMTYPE as any).MACROMOLECULE && (c.meta as any)?.units === 'helm');
    this._pickerActiveCache = {dfId: df.id, active};
    return active;
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

  render(g: any, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    // Register document-level listeners for the hover-paint atom picker.
    // Alt+hover: "paint" atoms by moving the cursor over them while
    // holding Alt. No click needed, so the grid never changes the
    // current cell and the context panel (e.g. AutoDock) stays open.
    if (!this._documentListenerAttached) {
      this._documentListenerAttached = true;
      document.addEventListener('mousemove', (e: MouseEvent) => {
        this._onDocumentMouseMove(e);
      });
      document.addEventListener('mousedown', (e: MouseEvent) => {
        this._onDocumentMouseDown(e);
      }, true);
      document.addEventListener('keydown', (e: KeyboardEvent) => {
        this._onDocumentKeyDown(e);
      });
    }

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

    // Pre-warm the atom-positions cache so hover highlighting responds
    // instantly instead of blocking on first mousemove. _getCellAtomPositions
    // takes CSS dimensions and internally scales by DPR to build the key.
    if (this._isPickerActive(gridCell.cell.column)) {
      const cssW = w / r;
      const cssH = h / r;
      // Mirror the cache-key logic from _getCellAtomPositions.
      const dpr = window.devicePixelRatio || 1;
      const cacheKey = molString + '|' + Math.round(cssW * dpr) + 'x' + Math.round(cssH * dpr);
      if (!this.atomPositionsCache.has(cacheKey)) {
        const mol = molString;
        const cw = cssW;
        const ch = cssH;
        (typeof requestIdleCallback === 'function' ? requestIdleCallback : setTimeout)(
          () => {this._getCellAtomPositions(mol, cw, ch);},
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

  // ========================================================================
  // In-grid box picker: mousedown on a molecule cell starts a drag; the box
  // overlay is drawn with an HTML div positioned over the grid canvas; on
  // mouseup the atoms enclosed by the box are added to that row's pick set
  // via an ISubstructProvider on the column's temp store, and the grid is
  // invalidated to repaint the cell with the new highlights.
  // ========================================================================

  /** Removes every `.chem-grid-box-picker` overlay element from the document.
   *  Safety net for leftover overlays. */
  private _sweepOverlays(): void {
    const orphans = document.querySelectorAll('.chem-grid-box-picker');
    for (let i = 0; i < orphans.length; i++) {
      const el = orphans[i];
      if (el.parentNode) el.parentNode.removeChild(el);
    }
  }

  /** Global mousedown handler (capture phase on document). Finds the active
   *  grid via `grok.shell.tv.grid`, hit-tests to determine which cell was
   *  clicked, and starts a drag if it's a molecule cell with the picker
   *  enabled. This replaces both the renderer's `onMouseDown` lifecycle
   *  (which stops dispatching after a table switch) and per-canvas
   *  listeners (which failed because the rendering canvas is not the
   *  interactive canvas). */
  // -- Hover-paint atom picker ------------------------------------------------
  // Hold Alt and move the mouse over atoms in a SMILES cell to "paint" them
  // into the selection. No click is needed, so the grid never changes the
  // current cell and the context panel (e.g. AutoDock Molstar viewer) stays
  // open. Alt+click on an already-selected atom removes it (toggle).

  /** Resolves the molecule cell and atom under the cursor. Returns null
   *  if the cursor isn't over a valid picker-active molecule cell atom. */
  private _hitTestAtom(e: MouseEvent): {
    grid: DG.Grid; gridCell: DG.GridCell; cellInfo: CellInteractiveInfo;
    nearest: number; pointerCellX: number; pointerCellY: number;
  } | null {
    const grid = grok.shell.tv?.grid;
    if (!grid) return null;
    let gridRoot: HTMLElement | null = null;
    try {gridRoot = grid.root;} catch {return null;}
    if (!gridRoot) return null;

    const gridRect = gridRoot.getBoundingClientRect();
    const localX = e.clientX - gridRect.left;
    const localY = e.clientY - gridRect.top;
    if (localX < 0 || localY < 0 ||
        localX > gridRect.width || localY > gridRect.height) return null;

    let gridCell: DG.GridCell | null = null;
    try {gridCell = grid.hitTest(localX, localY);} catch {return null;}
    if (!gridCell?.tableColumn || gridCell.tableRowIndex == null ||
        gridCell.tableRowIndex < 0) return null;
    if (gridCell.tableColumn.semType !== DG.SEMTYPE.MOLECULE) return null;
    if (!this._isPickerActive(gridCell.tableColumn)) return null;

    const molString: string = gridCell.cell.value;
    if (!molString || DG.chem.Sketcher.isEmptyMolfile(molString)) return null;

    const bounds = gridCell.bounds;
    const pointerCellX = e.clientX - (gridRect.left + bounds.left);
    const pointerCellY = e.clientY - (gridRect.top + bounds.top);
    if (pointerCellX < 0 || pointerCellY < 0 ||
        pointerCellX > bounds.width || pointerCellY > bounds.height) return null;

    const cellInfo = this._getCellAtomPositions(molString, bounds.width, bounds.height);
    if (!cellInfo || cellInfo.positions.size === 0) return null;

    const nearest = this._findNearestAtom(cellInfo.positions, pointerCellX, pointerCellY);
    if (nearest === null) return null;

    return {grid, gridCell, cellInfo, nearest, pointerCellX, pointerCellY};
  }

  /**
   * Mousemove handler — two modes:
   *
   * **Without Alt (preview):** highlights ONLY the single atom under the
   * cursor. Moving to a different atom highlights that one instead. Moving
   * away from all atoms clears the preview. No persistent selection.
   *
   * **With Alt (paint):** each atom the cursor passes over is ADDED to
   * the persistent selection. The selection accumulates until cleared
   * with Escape.
   */
  private _onDocumentMouseMove(e: MouseEvent): void {
    // Only activate the picker when the current cell is a Molecule3D
    // (docking pose) or HELM column. If the user clicked on SMILES,
    // binding energy, or any other column, hovering does nothing.
    const df = grok.shell.tv?.grid?.dataFrame;
    if (df) {
      const curCol = df.currentCol;
      const is3D = curCol?.semType === DG.SEMTYPE.MOLECULE3D;
      const isHelm = curCol?.semType === (DG.SEMTYPE as any).MACROMOLECULE &&
        (curCol?.meta as any)?.units === 'helm';
      if (!is3D && !isHelm) return;
    }

    const hit = this._hitTestAtom(e);

    if (e.altKey || e.shiftKey) {
      // ---- Alt/Shift+hover: PAINT / ERASE mode ------------------------------
      // Alt+hover: adds atoms to persistent selection.
      // Shift+hover: removes atoms from persistent selection.
      this._previewAtomIdx = null;

      if (!hit) return;
      const col = hit.gridCell.tableColumn!;
      const colName = col.name;
      const rowIdx = hit.gridCell.tableRowIndex!;

      // Skip if same atom AND same mode as last action (avoid re-firing per pixel).
      const eraseMode = e.shiftKey && !e.altKey;
      if (this._lastHoveredAtom &&
          this._lastHoveredAtom.col === colName &&
          this._lastHoveredAtom.rowIdx === rowIdx &&
          this._lastHoveredAtom.atomIdx === hit.nearest &&
          this._lastHoveredAtom.erase === eraseMode) return;
      this._lastHoveredAtom = {col: colName, rowIdx, atomIdx: hit.nearest, erase: eraseMode};

      if (eraseMode) {
        // Shift: remove atom from persistent selection.
        this._removeAtomFromRow(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
      } else {
        // Alt: add atom to persistent selection.
        this._addAtomToRow(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
      }
      // Clear render cache so the 2D cell fully redraws.
      this.rendersCache = new DG.LruCache<String, ImageData>();
      hit.grid.invalidate();
    } else {
      // ---- No Alt: PREVIEW mode (persistent + one hovered atom) ----------
      // Shows all persistently selected atoms PLUS the single hovered atom.
      // Moving away clears the preview but keeps persistent atoms.
      if (!hit) {
        if (this._previewAtomIdx !== null)
          this._removePreviewAtom();

        this._lastHoveredAtom = null;
        return;
      }
      const col = hit.gridCell.tableColumn!;
      const colName = col.name;
      const rowIdx = hit.gridCell.tableRowIndex!;

      // Same atom — skip.
      if (this._lastHoveredAtom &&
          this._lastHoveredAtom.col === colName &&
          this._lastHoveredAtom.rowIdx === rowIdx &&
          this._lastHoveredAtom.atomIdx === hit.nearest) return;
      this._lastHoveredAtom = {col: colName, rowIdx, atomIdx: hit.nearest};

      // Show persistent selection + this one preview atom.
      this._setPreviewAtom(col, rowIdx, hit.nearest, hit.cellInfo.bondAtoms);
      hit.grid.invalidate();
    }
  }

  /** Escape key clears the persistent Alt+hover selection. */
  private _onDocumentKeyDown(e: KeyboardEvent): void {
    if (e.key !== 'Escape') return;
    const grid = grok.shell.tv?.grid;
    if (!grid) return;
    const df = grid.dataFrame;
    if (!df) return;
    const molCol = df.columns.toList().find(
      (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE && this._isPickerActive(c));
    if (!molCol) return;

    // Remove all atom-picker providers from the column.
    type TaggedProvider = ISubstructProvider & {__atomPicker?: boolean};
    const providers = (molCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as TaggedProvider[];
    molCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = providers.filter((p) => !p.__atomPicker);
    this._lastHoveredAtom = null;
    this._previewAtomIdx = null;
    // Clear render cache so the 2D cell redraws without highlights.
    this.rendersCache = new DG.LruCache<String, ImageData>();
    grid.invalidate();

    // Fire persistent clear-all event so 3D cache clears completely
    // — including comparison molecules on other rows.
    grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
      column: molCol, rowIdx: -1, atoms: [], persistent: true, clearAll: true,
    });
  }

  /** Mousedown: kept minimal — only clears stale overlays. */
  private _onDocumentMouseDown(e: MouseEvent): void {
    this._sweepOverlays();
  }

  /** Returns the index of the atom whose center is nearest to the click
   *  point, but only if it falls within `CLICK_TOLERANCE_PX` (CSS pixels).
   *  Returns null if the user clicked too far from any atom. */
  private _findNearestAtom(
    positions: Map<number, {x: number, y: number}>,
    clickX: number, clickY: number,
  ): number | null {
    const CLICK_TOLERANCE_PX = 12;
    const tolSq = CLICK_TOLERANCE_PX * CLICK_TOLERANCE_PX;
    let nearestIdx: number | null = null;
    let nearestDistSq = tolSq;
    for (const [idx, p] of positions.entries()) {
      const dx = p.x - clickX;
      const dy = p.y - clickY;
      const dsq = dx * dx + dy * dy;
      if (dsq < nearestDistSq) {
        nearestDistSq = dsq;
        nearestIdx = idx;
      }
    }
    return nearestIdx;
  }

  /** Returns the persistent (Alt+painted) atoms for a given row. These
   *  are the atoms stored in the provider's `__atoms` set. */
  private _getPersistentAtoms(col: DG.Column, rowIdx: number): Set<number> {
    type TaggedProvider = ISubstructProvider & {__atomPicker?: boolean, __rowIdx?: number, __atoms?: Set<number>};
    const existing = (col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as TaggedProvider[];
    const prior = existing.find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
    return new Set<number>(prior?.__atoms ?? []);
  }

  /** Adds a single atom to the persistent selection (Alt+hover paint mode).
   *  Unlike _toggleAtomInRow, this only ADDS — never removes. */
  private _addAtomToRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    const current = this._getPersistentAtoms(col, rowIdx);
    if (current.has(atomIdx)) return; // already selected
    current.add(atomIdx);
    this._updateRowSelection(col, rowIdx, [...current], {add: true}, bondAtoms);
  }

  /** Removes a single atom from the persistent selection (Shift+hover).
   *  Works the same way as Escape (direct provider manipulation) but
   *  for a single atom instead of clearing all. */
  private _removeAtomFromRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    type TaggedProvider = ISubstructProvider & {__atomPicker?: boolean, __rowIdx?: number, __atoms?: Set<number>};
    const existing = (col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as TaggedProvider[];
    const prov = existing.find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
    if (!prov?.__atoms?.has(atomIdx)) return; // not in selection

    prov.__atoms.delete(atomIdx);

    if (prov.__atoms.size === 0) {
      // Last atom — remove the provider entirely (same as Escape).
      col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = existing.filter(
        (p) => !p.__atomPicker || p.__rowIdx !== rowIdx);
    } else {
      // Rebuild the provider with the reduced atom set.
      const pickColor: number[] = [1.0, 1.0, 0.0, 1.0];
      const atomsArr = [...prov.__atoms];
      const highlightAtomColors: {[k: number]: number[]} = {};
      for (const a of atomsArr) highlightAtomColors[a] = pickColor;
      const {bondsArr, highlightBondColors} =
        RDKitCellRenderer._computeSelectedBonds(prov.__atoms, bondAtoms, pickColor);

      // Replace the getSubstruct callback with updated atoms.
      prov.getSubstruct = (ridx) => {
        if (!this._isPickerActive(col)) return undefined;
        if (ridx !== rowIdx) return undefined;
        return {atoms: atomsArr, bonds: bondsArr, highlightAtomColors, highlightBondColors};
      };
    }

    // Clear render cache + invalidate (same pattern as Escape).
    this.rendersCache = new DG.LruCache<String, ImageData>();

    // Fire persistent event so 3D cache + Molstar update.
    const atoms = prov.__atoms.size > 0 ? [...prov.__atoms] : [];
    const eventArgs: any = {column: col, rowIdx, atoms, persistent: true};
    try {
      const df = col.dataFrame;
      if (df) {
        const mol3DCol = df.columns.toList().find(
          (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE3D);
        if (mol3DCol && rowIdx >= 0) {
          const smiles2D = col.get(rowIdx);
          const pose3D = mol3DCol.get(rowIdx);
          if (smiles2D && pose3D) {
            const mapping = mapAtomIndices2Dto3D(this.rdKitModule, smiles2D, pose3D);
            if (mapping) {
              eventArgs.mapping3D = mapping;
              eventArgs.mol3DColumnName = mol3DCol.name;
            }
          }
        }
      }
    } catch {/* mapping failed */}
    grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, eventArgs);
  }

  /** Sets the preview to show persistent atoms PLUS one hovered atom.
   *  The hovered atom is tracked separately so it can be removed when
   *  the cursor moves away without losing the persistent selection.
   *
   *  Important: we call _updateRowSelection to render the combined set,
   *  then RESTORE `__atoms` to the persistent-only set so the preview
   *  atom doesn't leak into the persistent storage. */
  private _setPreviewAtom(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    this._previewAtomIdx = atomIdx;
    this._previewRowIdx = rowIdx;
    const persistent = this._getPersistentAtoms(col, rowIdx);
    const combined = new Set(persistent);
    combined.add(atomIdx);
    // Render the combined set (persistent + preview) in both 2D and 3D.
    // Mark as non-persistent so the 3D cache isn't overwritten.
    this._updateRowSelection(col, rowIdx, [...combined], {add: true}, bondAtoms, true, false);
    // Restore __atoms to persistent-only so _getPersistentAtoms stays correct.
    type TaggedProvider = ISubstructProvider & {__atomPicker?: boolean, __rowIdx?: number, __atoms?: Set<number>};
    const providers = (col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as TaggedProvider[];
    const prov = providers.find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
    if (prov) prov.__atoms = persistent;
  }

  /** Removes just the preview atom, restoring the display to only the
   *  persistent (Alt+painted) selection. Also re-broadcasts the persistent
   *  atoms after a delay so Molstar picks them up after any viewer rebuild
   *  triggered by hover-row changes. */
  private _removePreviewAtom(): void {
    const prevAtom = this._previewAtomIdx;
    const prevRow = this._previewRowIdx;
    this._previewAtomIdx = null;
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
      // Fire event as non-persistent so the cache keeps the stable Alt set.
      const cellInfo = this._getCellAtomPositions(
        molCol.get(prevRow), 100, 100);
      const bondAtoms = cellInfo?.bondAtoms ?? new Map();
      this._updateRowSelection(molCol, prevRow, [...persistent], {add: true}, bondAtoms, true, false);
    } else {
      // No persistent atoms — clear everything.
      type TaggedProvider = ISubstructProvider & {__atomPicker?: boolean};
      const providers = (molCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as TaggedProvider[];
      molCol.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = providers.filter(
        (p) => !p.__atomPicker || (p as any).__rowIdx !== prevRow);
      // Fire event to clear 3D highlights too (non-persistent).
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
        column: molCol, rowIdx: prevRow, atoms: [], persistent: false,
      });
    }
    grid.invalidate();
  }

  /** Toggles a single atom in a row's atom-picker selection. Adds the atom
   *  if it's not in the current selection, removes it if it is. Used by
   *  Alt+click for one-atom-at-a-time fine-tuning. Also re-derives the
   *  highlighted bonds from the new atom set using `bondAtoms`. */
  private _toggleAtomInRow(col: DG.Column, rowIdx: number, atomIdx: number,
    bondAtoms: Map<number, [number, number]>): void {
    type TaggedProvider = ISubstructProvider & {__atomPicker?: boolean, __rowIdx?: number, __atoms?: Set<number>};
    const existing = (col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as TaggedProvider[];
    const prior = existing.find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
    const current = new Set<number>(prior?.__atoms ?? []);
    if (current.has(atomIdx)) current.delete(atomIdx);
    else current.add(atomIdx);
    // Build the fresh provider (same shape as _updateRowSelection writes).
    const pickColor: number[] = [1.0, 1.0, 0.0, 1.0]; // intense yellow #FFFF00
    const atomsArr = [...current];
    const highlightAtomColors: {[k: number]: number[]} = {};
    for (const a of atomsArr) highlightAtomColors[a] = pickColor;
    const {bondsArr, highlightBondColors} =
      RDKitCellRenderer._computeSelectedBonds(current, bondAtoms, pickColor);
    const provider: TaggedProvider = {
      getSubstruct: (ridx) => {
        if (!this._isPickerActive(col)) return undefined;
        if (ridx !== rowIdx) return undefined;
        return {atoms: atomsArr, bonds: bondsArr, highlightAtomColors, highlightBondColors};
      },
    };
    provider.__atomPicker = true;
    provider.__rowIdx = rowIdx;
    provider.__atoms = current;
    col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = existing.filter(
      (p) => !p.__atomPicker || p.__rowIdx !== rowIdx);
    addSubstructProvider(col.temp, provider);

    // Notify cross-package listeners (e.g. BiostructureViewer's 3D Molstar
    // viewer) about the updated selection so they can mirror the highlights.
    // Include a 3D atom-index mapping if a MOLECULE3D column exists in the
    // same dataframe (e.g. AutoDock poses). The mapping is computed via
    // RDKit substructure match between the 2D SMILES and the 3D pose, so
    // BiostructureViewer doesn't need to import RDKit or the mapper.
    const eventArgs: any = {column: col, rowIdx, atoms: atomsArr};
    try {
      const df = col.dataFrame;
      if (df) {
        const mol3DCol = df.columns.toList().find(
          (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE3D);
        if (mol3DCol && rowIdx >= 0) {
          const smiles2D = col.get(rowIdx);
          const pose3D = mol3DCol.get(rowIdx);
          if (smiles2D && pose3D) {
            const mapping = mapAtomIndices2Dto3D(this.rdKitModule, smiles2D, pose3D);
            if (mapping) {
              eventArgs.mapping3D = mapping;
              eventArgs.mol3DColumnName = mol3DCol.name;
            }
          }
        }
      }
    } catch {/* mapping failed — 3D viewers will fall back to direct index */}
    grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, eventArgs);
  }

  /** Given a set of selected atom indices and a bondAtoms map, returns the
   *  bond indices where BOTH endpoint atoms are in the selection, and a
   *  matching highlightBondColors map with every selected bond painted the
   *  same color. Shared by _updateRowSelection and _toggleAtomInRow. */
  private static _computeSelectedBonds(
    atoms: Set<number>,
    bondAtoms: Map<number, [number, number]>,
    color: number[],
  ): {bondsArr: number[], highlightBondColors: {[k: number]: number[]}} {
    const bondsArr: number[] = [];
    const highlightBondColors: {[k: number]: number[]} = {};
    for (const [bondIdx, [a1, a2]] of bondAtoms.entries()) {
      if (atoms.has(a1) && atoms.has(a2)) {
        bondsArr.push(bondIdx);
        highlightBondColors[bondIdx] = color;
      }
    }
    return {bondsArr, highlightBondColors};
  }

  /**
   * Merges the freshly boxed atoms into the row's current atom-picker selection
   * (replace / add / remove depending on modifiers) and writes a tagged
   * ISubstructProvider back onto the column's temp store. Providers are keyed
   * by row index, so picks on other rows survive untouched.
   *
   * The provider closes over `col` and checks `_isPickerActive()` at call
   * time, so the highlights auto-hide when the DataFrame no longer has
   * associated 3D/HELM columns. The picks themselves are preserved in the
   * provider's __atoms set.
   */
  private _updateRowSelection(col: DG.Column, rowIdx: number, boxed: number[],
    modifiers: {add: boolean}, bondAtoms: Map<number, [number, number]>,
    fire3DEvent: boolean = true, persistent: boolean = true): void {
    type TaggedProvider = ISubstructProvider & {__atomPicker?: boolean, __rowIdx?: number, __atoms?: Set<number>};
    const existing = (col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] ?? []) as TaggedProvider[];
    const prior = existing.find((p) => p.__atomPicker && p.__rowIdx === rowIdx);
    const current = new Set<number>(prior?.__atoms ?? []);

    // No modifier → replace existing selection with boxed atoms.
    // Alt (add) → keep existing and union with boxed atoms.
    if (!modifiers.add) current.clear();
    for (const a of boxed) current.add(a);

    // Build the fresh provider. Intense yellow was chosen over magenta or
    // green because green collides with Datagrok's row-selection tint
    // (hard to tell "row is selected" apart from "atoms are picked") and
    // pure yellow reads unambiguously as "highlighted region" on the
    // white cell background — classic text-highlighter look.
    const pickColor: number[] = [1.0, 1.0, 0.0, 1.0]; // intense yellow #FFFF00
    const atomsArr = [...current];
    const highlightAtomColors: {[k: number]: number[]} = {};
    for (const a of atomsArr) highlightAtomColors[a] = pickColor;

    // Derive the bonds whose BOTH endpoints are in the current selection.
    // Those get highlighted in the same magenta as the atoms, making the
    // selected region read as a connected sub-molecule rather than a
    // scatter of disconnected atoms.
    const {bondsArr, highlightBondColors} =
      RDKitCellRenderer._computeSelectedBonds(current, bondAtoms, pickColor);

    const provider: TaggedProvider = {
      getSubstruct: (ridx) => {
        if (!this._isPickerActive(col)) return undefined;
        if (ridx !== rowIdx) return undefined;
        return {atoms: atomsArr, bonds: bondsArr, highlightAtomColors, highlightBondColors};
      },
    };
    provider.__atomPicker = true;
    provider.__rowIdx = rowIdx;
    provider.__atoms = current;

    // Replace only the prior provider for this row; keep picks on other rows
    col.temp[ChemTemps.SUBSTRUCT_PROVIDERS] = existing.filter(
      (p) => !p.__atomPicker || p.__rowIdx !== rowIdx);
    addSubstructProvider(col.temp, provider);

    // Only fire the 3D event for persistent (Alt+hover) selections,
    // not for preview-only highlights which are 2D-only and transient.
    if (fire3DEvent) {
      const eventArgs2: any = {column: col, rowIdx, atoms: atomsArr, persistent};
      try {
        const df = col.dataFrame;
        if (df) {
          const mol3DCol = df.columns.toList().find(
            (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE3D);
          if (mol3DCol && rowIdx >= 0) {
            const smiles2D = col.get(rowIdx);
            const pose3D = mol3DCol.get(rowIdx);
            if (smiles2D && pose3D) {
              const mapping = mapAtomIndices2Dto3D(this.rdKitModule, smiles2D, pose3D);
              if (mapping) {
                eventArgs2.mapping3D = mapping;
                eventArgs2.mol3DColumnName = mol3DCol.name;
              }
            }
          }
        }
      } catch {/* mapping failed */}
      grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, eventArgs2);
    }
  }

  /**
   * Compute atom positions for a cell in cell-local CSS pixels.
   *
   * To guarantee the SVG layout matches the canvas-rasterized layout
   * pixel-for-pixel, we use the SAME `IMolContext` the cell renderer uses
   * (via `_fetchMol`). Calling `get_mol` directly here would create a
   * different mol instance whose 2D coordinates may have been generated
   * with different options (kekulization, mol-block wedging, etc.) — that
   * was the source of the position drift the user reported.
   *
   * The cell renderer's `render()` multiplies cell CSS width/height by
   * `devicePixelRatio` before drawing, and uses `RDKIT_COMMON_RENDER_OPTS`
   * (which includes `fixedScale: 0.07`). We mirror the same DPR scaling
   * and options here, then divide the extracted SVG coordinates back by
   * DPR to convert them to cell-local CSS pixels — the space the pointer
   * events live in.
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
      // Use the SAME cached molctx the renderer uses, so atom 2D coords
      // match exactly. The mol is owned by molCache — DO NOT delete it.
      const molRenderingInfo = this._fetchMol(molString, [], false, false, {}, false);
      const mol = molRenderingInfo.molCtx.mol;
      if (!mol || !mol.is_valid()) return null;

      // Mirror the canvas render path's options exactly. RDKit uses the
      // same internal layout for both `draw_to_canvas_with_highlights`
      // and `get_svg_with_highlights` when given identical opts on the
      // same mol object.
      const details: {[k: string]: any} = {};
      for (const k of Object.keys(RDKIT_COMMON_RENDER_OPTS))
        details[k] = RDKIT_COMMON_RENDER_OPTS[k];
      details.width = w;
      details.height = h;
      details.atoms = [];
      details.bonds = [];
      details.highlightAtomColors = {};
      details.highlightBondColors = {};
      // Mirror the kekulize / molBlockWedging behaviour from
      // drawRdKitMoleculeToOffscreenCanvas so the SVG layout cannot drift
      // from the canvas layout because of these options.
      if (!molRenderingInfo.molCtx.kekulize) details.kekulize = false;
      if (molRenderingInfo.molCtx.useMolBlockWedging) {
        details.useMolBlockWedging = true;
        details.wedgeBonds = false;
        details.addChiralHs = false;
      }

      const svgString = mol.get_svg_with_highlights(JSON.stringify(details));
      // innerHTML on a real div in the main document so layout/getBBox work.
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
        const extracted = this._extractAtomPositionsFromSvg(svgEl);
        // Convert from DPR-scaled SVG coordinates back to cell-local CSS
        // pixels so they can be compared against mouse event coordinates.
        const positions = new Map<number, {x: number, y: number}>();
        for (const [idx, p] of extracted.positions.entries())
          positions.set(idx, {x: p.x / dpr, y: p.y / dpr});
        // Bond pairs are atom-index pairs — no coordinate scaling needed.
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

  private _extractAtomPositionsFromSvg(svgEl: SVGSVGElement):
      {positions: Map<number, {x: number, y: number}>, bondAtoms: Map<number, [number, number]>} {
    const positions = new Map<number, {x: number, y: number}>();
    /** bond index → the two atom indices this bond connects. Populated
     *  below from `class="bond-K atom-A atom-B"` on RDKit's bond paths. */
    const bondAtoms = new Map<number, [number, number]>();

    // Text-labeled heteroatoms
    const texts = svgEl.querySelectorAll('text[class*="atom-"]');
    for (let i = 0; i < texts.length; i++) {
      const t = texts[i];
      const cls = t.getAttribute('class') || '';
      const m = /(?:^|\s)atom-(\d+)/.exec(cls);
      if (!m) continue;
      const idx = parseInt(m[1], 10);
      try {
        const bb = (t as unknown as SVGGraphicsElement).getBBox();
        positions.set(idx, {x: bb.x + bb.width / 2, y: bb.y + bb.height / 2});
      } catch {/* ignore */}
    }

    // Carbons — average bond endpoints extracted from <path> d attributes.
    // Also populate bondAtoms as a side-effect so we don't need to re-walk
    // the path list when the picker wants bond connectivity.
    const bondEnds = new Map<number, Array<{x: number, y: number}>>();
    const bondPaths = svgEl.querySelectorAll('path[class*="bond-"]');
    for (let i = 0; i < bondPaths.length; i++) {
      const p = bondPaths[i];
      const cls = p.getAttribute('class') || '';
      const atomIds: number[] = [];
      const atomMatches = cls.match(/atom-(\d+)/g) || [];
      for (const mm of atomMatches) {
        const n = parseInt(mm.slice(5), 10);
        if (!Number.isNaN(n)) atomIds.push(n);
      }
      if (atomIds.length !== 2) continue;

      // Record the bond-index → atom-pair mapping. RDKit emits multiple
      // path segments per bond for aromatic/double/triple bonds (one
      // segment per line), all sharing the same bond-K class — but
      // `bondAtoms.set` is idempotent for the same key, so overwriting is
      // harmless.
      const bondMatch = /(?:^|\s)bond-(\d+)/.exec(cls);
      if (bondMatch) {
        const bondIdx = parseInt(bondMatch[1], 10);
        bondAtoms.set(bondIdx, [atomIds[0], atomIds[1]]);
      }

      const d = p.getAttribute('d') || '';
      const coords: Array<{x: number, y: number}> = [];
      const re = /[ML]\s*([\-\d.]+)[\s,]+([\-\d.]+)/g;
      let mm: RegExpExecArray | null;
      while ((mm = re.exec(d)) !== null)
        coords.push({x: parseFloat(mm[1]), y: parseFloat(mm[2])});
      if (coords.length < 2) continue;
      if (!bondEnds.has(atomIds[0])) bondEnds.set(atomIds[0], []);
      if (!bondEnds.has(atomIds[1])) bondEnds.set(atomIds[1], []);
      bondEnds.get(atomIds[0])!.push(coords[0]);
      bondEnds.get(atomIds[1])!.push(coords[coords.length - 1]);
    }
    for (const [idx, pts] of bondEnds.entries()) {
      if (positions.has(idx)) continue;
      const cx = pts.reduce((s, p) => s + p.x, 0) / pts.length;
      const cy = pts.reduce((s, p) => s + p.y, 0) / pts.length;
      positions.set(idx, {x: cx, y: cy});
    }

    return {positions, bondAtoms};
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
  // Multiple active providers — delegate to the chem-meta merge helper, which
  // unions atoms / bonds / highlightAtomColors / highlightBondColors and the
  // atomLabels / atomNotes / bondNotes annotation fields. Later providers
  // override earlier ones on color collisions, matching the layered-highlight
  // semantics used elsewhere in the renderer.
  return mergeSubstructsLib(defined);
}
