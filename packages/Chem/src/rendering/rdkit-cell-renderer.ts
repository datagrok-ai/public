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
  CHEM_ATOM_PICKER_TAG, CHEM_INTERACTIVE_SELECTION_EVENT,
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
  private _drag: {
    grid: DG.Grid;
    tableCol: DG.Column;
    rowIdx: number;
    // Cell bounds expressed in viewport (client) CSS pixels so we can draw
    // the overlay with position: fixed and clientX/Y directly.
    cellClientLeft: number;
    cellClientTop: number;
    cellWidth: number;
    cellHeight: number;
    // Mouse position at mousedown, in cell-local CSS pixels.
    startX: number;
    startY: number;
    curX: number;
    curY: number;
    /** `add` = Alt held (additive — adds boxed atoms to existing selection
     *  rather than replacing it).
     *
     *  Alt is the only modifier we can use: Datagrok's grid claims Shift+
     *  Click for row range-select and Ctrl/Cmd+Click for row multi-select
     *  (standard spreadsheet conventions), so neither modifier ever
     *  dispatches `onMouseDown` to a cell renderer. Alt is not part of the
     *  spreadsheet selection vocabulary, so the grid lets it through. */
    modifiers: {add: boolean};
    positions: Map<number, {x: number, y: number}>;
    bondAtoms: Map<number, [number, number]>;
    overlay: HTMLDivElement;
    onMove: (e: MouseEvent) => void;
    onUp: (e: MouseEvent) => void;
  } | null = null;

  constructor(rdKitModule: RDModule) {
    super();
    this.rdKitModule = rdKitModule;
    this.canvasCounter = 0;
    this.canvasReused = new OffscreenCanvas(this.defaultWidth, this.defaultHeight);
    this.molCache.onItemEvicted = function(obj: { [_: string]: any }) {
      obj.mol?.delete();
    };
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
    // Register a single document-level mousedown listener (once, on the
    // first render). This replaces per-canvas listeners which failed
    // because g.canvas is a rendering canvas, not the interactive one.
    if (!this._documentListenerAttached) {
      this._documentListenerAttached = true;
      document.addEventListener('mousedown', (e: MouseEvent) => {
        this._onDocumentMouseDown(e);
      }, true); // capture phase — fires before anything can stopPropagation
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

  /** Tears down an active drag: removes the document listeners and any
   *  overlay(s) from the DOM and clears `_drag`. Safe to call when no drag
   *  is active.
   *
   *  We sweep ALL `.chem-grid-box-picker` elements from the document, not
   *  just the one in `_drag`, to defend against orphan overlays that might
   *  exist if Datagrok's grid dispatches a second `onMouseDown` during the
   *  lifecycle of the first drag (which we observed in practice — the
   *  symptom was a box that persisted after mouseup). */
  private _endDrag(): void {
    const d = this._drag;
    if (d) {
      document.removeEventListener('mousemove', d.onMove, true);
      document.removeEventListener('mouseup', d.onUp, true);
    }
    this._sweepOverlays();
    this._drag = null;
  }

  /** Removes every `.chem-grid-box-picker` element from the document. Called
   *  as a safety net from `_endDrag` and from `onMouseDown` before a new
   *  overlay is created. */
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
  private _onDocumentMouseDown(e: MouseEvent): void {
    this._endDrag();
    this._sweepOverlays();

    if (e.button !== 0) return;
    if ((e.buttons & 1) === 0) return;

    // Get the current table view's grid. This is the grid the user is
    // interacting with, regardless of which table was opened first.
    const grid = grok.shell.tv?.grid;
    if (!grid) return;

    // Hit-test. We need grid-canvas-relative coordinates. `e.target` is
    // the actual element that received the event — compute the offset
    // from there. For a document listener, e.target is usually the grid's
    // interactive overlay element.
    const targetEl = e.target as HTMLElement;
    if (!targetEl) return;
    const targetRect = targetEl.getBoundingClientRect();
    const localX = e.clientX - targetRect.left;
    const localY = e.clientY - targetRect.top;

    let gridCell: DG.GridCell | null = null;
    try {gridCell = grid.hitTest(localX, localY);} catch {return;} // hitTest can throw if the grid is in a transitional state
    if (!gridCell || !gridCell.tableColumn || gridCell.tableRowIndex == null ||
        gridCell.tableRowIndex < 0) return;
    if (gridCell.tableColumn.semType !== DG.SEMTYPE.MOLECULE) return;

    const pickerTag =
      (gridCell.tableColumn.tags as any)[CHEM_ATOM_PICKER_TAG] ??
      gridCell.tableColumn.getTag(CHEM_ATOM_PICKER_TAG);
    if (pickerTag !== 'true') return;

    const molString: string = gridCell.cell.value;
    if (!molString || DG.chem.Sketcher.isEmptyMolfile(molString)) return;

    // Cell position in viewport (client) CSS pixels.
    const bounds = gridCell.bounds;
    const canvasRect = targetRect;
    const cellClientLeft = canvasRect.left + bounds.left;
    const cellClientTop = canvasRect.top + bounds.top;

    const pointerCellX = e.clientX - cellClientLeft;
    const pointerCellY = e.clientY - cellClientTop;
    if (pointerCellX < 0 || pointerCellY < 0 ||
        pointerCellX > bounds.width || pointerCellY > bounds.height) return;

    // Lazily compute the cell's interactive layout info (atom positions
    // and bond atom-pairs), keyed by molString and DPR-scaled dimensions
    // so the cache lines up with what the canvas draws.
    const cellInfo = this._getCellAtomPositions(molString, bounds.width, bounds.height);
    if (!cellInfo || cellInfo.positions.size === 0) return;
    const positions = cellInfo.positions;
    const bondAtoms = cellInfo.bondAtoms;

    // Fixed-position overlay anchored at the mouse-down point.
    const overlay = document.createElement('div');
    overlay.className = 'chem-grid-box-picker';
    Object.assign(overlay.style, {
      position: 'fixed',
      left: e.clientX + 'px',
      top: e.clientY + 'px',
      width: '0px',
      height: '0px',
      border: '1px dashed #4a8',
      background: 'rgba(102, 204, 102, 0.15)',
      pointerEvents: 'none',
      boxSizing: 'border-box',
      zIndex: '10000',
    } as Partial<CSSStyleDeclaration>);
    document.body.appendChild(overlay);

    // Capture-phase document listeners for drag tracking.
    const onMove = (mv: MouseEvent) => this._handleDragMove(mv);
    const onUp = (mu: MouseEvent) => this._handleDragUp(mu);
    document.addEventListener('mousemove', onMove, true);
    document.addEventListener('mouseup', onUp, true);

    this._drag = {
      grid: grid,
      tableCol: gridCell.tableColumn,
      rowIdx: gridCell.tableRowIndex,
      cellClientLeft,
      cellClientTop,
      cellWidth: bounds.width,
      cellHeight: bounds.height,
      startX: pointerCellX,
      startY: pointerCellY,
      curX: pointerCellX,
      curY: pointerCellY,
      modifiers: {add: e.altKey},
      positions,
      bondAtoms,
      overlay,
      onMove,
      onUp,
    };
    e.preventDefault();
    e.stopPropagation();
  }

  private _handleDragMove(e: MouseEvent): void {
    const d = this._drag;
    if (!d) return;

    // Track Alt state continuously throughout the drag. Reading `e.altKey`
    // only at mousedown or only at mouseup is unreliable: at mousedown the
    // user may not yet be holding Alt (they pressed it mid-drag); at
    // mouseup the browser sometimes loses modifier state if a synthetic
    // event is dispatched in the middle. The most reliable signal is
    // whatever was tracked on the last mousemove sample.
    d.modifiers.add = e.altKey;

    // cell-local coordinates, clamped so the rectangle never extends past
    // the molecule cell where the drag started
    const cellX = Math.max(0, Math.min(d.cellWidth, e.clientX - d.cellClientLeft));
    const cellY = Math.max(0, Math.min(d.cellHeight, e.clientY - d.cellClientTop));
    d.curX = cellX;
    d.curY = cellY;
    const vpLeft = d.cellClientLeft + Math.min(d.startX, d.curX);
    const vpTop = d.cellClientTop + Math.min(d.startY, d.curY);
    d.overlay.style.left = vpLeft + 'px';
    d.overlay.style.top = vpTop + 'px';
    d.overlay.style.width = Math.abs(d.curX - d.startX) + 'px';
    d.overlay.style.height = Math.abs(d.curY - d.startY) + 'px';

    // Visual feedback: tint the overlay yellow when Alt is held mid-drag,
    // so the user can see the picker has recognised the modifier and will
    // add (rather than replace) the existing selection. Yellow matches
    // the pick color so the preview color is the color the atoms will
    // take once released. A darker gold border keeps the dashed outline
    // visible against the white cell background — a pure #FFFF00 border
    // on white disappears.
    if (e.altKey) {
      // additive — gold border + translucent yellow fill, matching picks
      d.overlay.style.borderColor = '#c9a400';
      d.overlay.style.background = 'rgba(255, 255, 0, 0.25)';
    } else {
      // replace — default green (distinct from yellow so drag vs.
      // applied state is immediately legible)
      d.overlay.style.borderColor = '#4a8';
      d.overlay.style.background = 'rgba(102, 204, 102, 0.15)';
    }
  }

  private _handleDragUp(e: MouseEvent): void {
    const d = this._drag;
    if (!d) return;

    // Capture data before _endDrag clears state + removes the overlay.
    const x0 = Math.min(d.startX, d.curX);
    const x1 = Math.max(d.startX, d.curX);
    const y0 = Math.min(d.startY, d.curY);
    const y1 = Math.max(d.startY, d.curY);
    const moved = Math.hypot(x1 - x0, y1 - y0) > 3;
    const positions = d.positions;
    const bondAtoms = d.bondAtoms;
    const tableCol = d.tableCol;
    const rowIdx = d.rowIdx;
    const grid = d.grid;
    // Combine the latest tracked-during-mousemove modifier state with the
    // mouseup event's own state. We OR them so Alt held at any recent
    // moment wins — defends against the browser losing modifier state on
    // the mouseup event itself.
    const modifiers = {add: e.altKey || d.modifiers.add};
    const startX = d.startX;
    const startY = d.startY;

    this._endDrag();

    // Click (no significant movement). If Alt was held and the click was
    // close to a specific atom, toggle that single atom — adds it if not
    // already in the selection, removes it if it is. This is the
    // fine-tuning complement to Alt+drag.
    if (!moved) {
      if (modifiers.add) {
        const nearest = this._findNearestAtom(positions, startX, startY);
        if (nearest !== null) {
          this._toggleAtomInRow(tableCol, rowIdx, nearest, bondAtoms);
          grid.invalidate();
        }
      }
      return;
    }

    const boxed: number[] = [];
    for (const [idx, p] of positions.entries()) {
      if (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1)
        boxed.push(idx);
    }
    this._updateRowSelection(tableCol, rowIdx, boxed, modifiers, bondAtoms);
    grid.invalidate();
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
        // Same tag check as _updateRowSelection so disabling the picker
        // hides Alt+click toggled atoms too.
        const tagVal = col.tags[CHEM_ATOM_PICKER_TAG] ?? col.getTag(CHEM_ATOM_PICKER_TAG);
        if (tagVal === 'false') return undefined;
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
    grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
      column: col, rowIdx, atoms: atomsArr,
    });
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
   * The provider closes over `col` and checks `CHEM_ATOM_PICKER_TAG` at call
   * time, so disabling the picker via the "Atom picker" column-property
   * checkbox instantly hides all existing highlights for that column (the
   * provider returns `undefined` while disabled). Re-enabling brings them
   * back — the picks themselves are preserved in the provider's __atoms set.
   */
  private _updateRowSelection(col: DG.Column, rowIdx: number, boxed: number[],
    modifiers: {add: boolean}, bondAtoms: Map<number, [number, number]>): void {
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
        // Check the column tag at render time so disabling the atom picker
        // via the column property checkbox instantly hides the highlights.
        // Reading both the tag proxy and getTag() so either write path
        // (col.tags[..] = v / col.setTag(..)) is honored.
        const tagVal = col.tags[CHEM_ATOM_PICKER_TAG] ?? col.getTag(CHEM_ATOM_PICKER_TAG);
        if (tagVal === 'false') return undefined;
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

    // Notify cross-package listeners (same event as _updateRowSelection).
    grok.events.fireCustomEvent(CHEM_INTERACTIVE_SELECTION_EVENT, {
      column: col, rowIdx, atoms: atomsArr,
    });
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
