/**
 * RDKit-based molecule cell renderer.
 * */

import * as DG from 'datagrok-api/dg';
import {_rdKitModule, drawErrorCross, drawRdKitMoleculeToOffscreenCanvas} from '../utils/chem-common-rdkit';
import {IMolContext, getMolSafe} from '../utils/mol-creation_rdkit';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import { ALIGN_BY_SCAFFOLD_TAG, FILTER_SCAFFOLD_TAG, HIGHLIGHT_BY_SCAFFOLD_TAG, SCAFFOLD_TREE_HIGHLIGHT } from '../constants';
import { hexToPercentRgb } from '../utils/chem-common';

export interface ISubstruct {
  atoms?: number[],
  bonds?: number[],
  highlightAtomColors?: {[key: number]: number[] | null},
  highlightBondColors?: {[key: number]: number[] | null}
}

interface IMolRenderingInfo {
  //mol: RDMol | null; // null when molString is invalid?
  molCtx: IMolContext;
  substruct: ISubstruct;
  molString: string;
}

export interface IColoredScaffold {
  molecule: string,
  color?: string,
  isSuperstructure?: string
}

export interface IHighlightTagInfo {
  scaffolds?: IColoredScaffold[],
  alighByFirstSubtruct: boolean,
}

export function _addColorsToBondsAndAtoms(mainSubstr: ISubstruct, color?: string, tempSubstr?: ISubstruct): void {
  if (color) {
    const colorArr = hexToPercentRgb(color);
    const substrToTakeAtomsFrom = tempSubstr ?? mainSubstr;
    if (substrToTakeAtomsFrom.atoms) {
      for (let j = 0; j < substrToTakeAtomsFrom.atoms.length; j++) {
        mainSubstr.highlightAtomColors ??= {};
        mainSubstr.highlightAtomColors[substrToTakeAtomsFrom.atoms[j]] = colorArr;
      };
    }
    if (substrToTakeAtomsFrom.bonds) {
      for (let j = 0; j < substrToTakeAtomsFrom.bonds.length; j++) {
        mainSubstr.highlightBondColors ??= {};
        mainSubstr.highlightBondColors[substrToTakeAtomsFrom.bonds[j]] = colorArr;
      };
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
  canvasReused: OffscreenCanvas;

  constructor(rdKitModule: RDModule) {
    super();
    this.rdKitModule = rdKitModule;
    this.canvasCounter = 0;
    this.canvasReused = new OffscreenCanvas(this.defaultWidth, this.defaultHeight);
    this.molCache.onItemEvicted = function(obj: {[_ : string]: any}) {
      obj.mol?.delete();
    };
  }

  ensureCanvasSize(w: number, h: number) : OffscreenCanvas {
    if (this.canvasReused.width < w || this.canvasReused.height < h)
      this.canvasReused = new OffscreenCanvas(Math.max(this.defaultWidth, w), Math.max(this.defaultHeight, h));
    return this.canvasReused;
  }

  get name(): string {return 'RDKit cell renderer';}
  get cellType(): DG.SemType {return DG.SEMTYPE.MOLECULE;}
  get defaultWidth(): number {return 200;}
  get defaultHeight(): number {return 100;}

  getDefaultSize(): {width: number, height: number} {
    return {width: this.defaultWidth, height: this.defaultHeight};
  }

  _fetchMolGetOrCreate(molString: string, scaffolds: IColoredScaffold[], molRegenerateCoords: boolean,
    details: object = {}, alignByFirstSubstr: boolean): IMolRenderingInfo {
    let molCtx: IMolContext;
    let mol = null;
    let substruct: ISubstruct = {};
    if ((details as any).isSubstructure) {
      if (molString.includes(' H ') || molString.includes('V3000')) {
        molCtx = getMolSafe(molString, { mergeQueryHs: true }, _rdKitModule);
        mol = molCtx.mol;
      } else {
        try {
          mol = this.rdKitModule.get_qmol(molString);
          mol.convert_to_aromatic_form();
        } catch (e) {
          if (mol) {
            mol.delete();
            mol = null;
          }
        }
        molCtx = { mol: mol, kekulize: false, isQMol: true, useMolBlockWedging: false };
      }
    } else {
      molCtx = getMolSafe(molString, details, _rdKitModule);
      mol = molCtx.mol;
    }

    if (mol !== null) {
      try {
        let molHasOwnCoords = (mol.has_coords() > 0);
        const scaffoldIsMolBlock = scaffolds.length ? DG.chem.isMolBlock(scaffolds[0].molecule) : null;
        const alignedByFirstSubstr = scaffoldIsMolBlock && alignByFirstSubstr;
        if (alignedByFirstSubstr) {
          const rdKitScaffoldMolCtx = this._fetchMol(scaffolds[0].molecule, [], molRegenerateCoords, false,
            { mergeQueryHs: true, isSubstructure: true }, false).molCtx;
          const rdKitScaffoldMol = rdKitScaffoldMolCtx.mol;
          if (rdKitScaffoldMol) {
            rdKitScaffoldMol.normalize_depiction(0);
            if (molHasOwnCoords)
              mol.normalize_depiction(0);
            let substructString = '';
            try {
              substructString = !scaffolds[0].isSuperstructure ? mol.generate_aligned_coords(rdKitScaffoldMol, JSON.stringify({
                useCoordGen: true,
                allowRGroups: true,
                acceptFailure: false,
                alignOnly: molHasOwnCoords,
              })) : mol.get_substruct_match(mol!);
            } catch {
              // exceptions should not be thrown anymore by RDKit, but let's play safe
            }
            if (substructString === '') {
              substruct = {};
              if (molHasOwnCoords)
                mol.straighten_depiction(true);
            } else
              substruct = JSON.parse(substructString);
              if (scaffolds[0].color) {
                if (!substruct.atoms)
                  substruct.atoms = [];
                if (!substruct.bonds)
                  substruct.bonds = [];
                _addColorsToBondsAndAtoms(substruct, scaffolds[0].color);
              }
          }
        }
        for (let i = alignedByFirstSubstr ? 1 : 0; i < scaffolds.length; i++) {
          const substructMol = this._fetchMol(scaffolds[i].molecule, [], molRegenerateCoords, false,
            { mergeQueryHs: true, isSubstructure: true }, false).molCtx.mol;
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
              };
            }
          }
        }
        if (mol.has_coords() === 0 || molRegenerateCoords) {
          mol.set_new_coords(molRegenerateCoords);
          molHasOwnCoords = false;
        }
        if (!scaffoldIsMolBlock) {
          mol.normalize_depiction(molHasOwnCoords ? 0 : 1);
          mol.straighten_depiction(molHasOwnCoords);
        } else if (!molHasOwnCoords)
          mol.normalize_depiction(0);

        molCtx.useMolBlockWedging = (mol.has_coords() === 2);
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
    };
  }

  _fetchMol(molString: string, scaffolds: IColoredScaffold[], molRegenerateCoords: boolean,
    scaffoldRegenerateCoords: boolean, details: object = {}, multiScaffold: boolean): IMolRenderingInfo {    
    const name = molString + ' || ' + JSON.stringify(scaffolds) + ' || ' +
      molRegenerateCoords + ' || ' + scaffoldRegenerateCoords +
      (Object.keys(details).length ? ' || ' + JSON.stringify(details) : '');
    return this.molCache.getOrCreate(name, (_: any) =>
      this._fetchMolGetOrCreate(molString, scaffolds, molRegenerateCoords, details, multiScaffold));
  }

  _rendererGetOrCreate(
    width: number, height: number, molString: string, scaffolds: IColoredScaffold[],
    highlightScaffold: boolean, molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean, multiScaffold: boolean): ImageData {
    const fetchMolObj : IMolRenderingInfo =
      this._fetchMol(molString, scaffolds, molRegenerateCoords, scaffoldRegenerateCoords, {}, multiScaffold);
    const rdKitMolCtx = fetchMolObj.molCtx;
    const rdKitMol = rdKitMolCtx.mol;//fetchMolObj.mol;
    const substruct = fetchMolObj.substruct;

    const canvas = this.ensureCanvasSize(width, height);//new OffscreenCanvas(width, height);
    const ctx = canvas.getContext('2d', {willReadFrequently: true})!;
    this.canvasCounter++;
    if (rdKitMol != null)
      drawRdKitMoleculeToOffscreenCanvas(rdKitMolCtx, width, height, canvas, highlightScaffold ? substruct : null);
    else {
      // draw a crossed rectangle
      ctx.clearRect(0, 0, width, height);
      drawErrorCross(ctx, width, height);
    }

    return ctx.getImageData(0, 0, !Math.floor(width) ? 1 : width, !Math.floor(height) ? 1 : height);
  }

  _fetchRender(
    width: number, height: number, molString: string, scaffolds: IColoredScaffold[],
    highlightScaffold: boolean, molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean, multiScaffold: boolean): ImageData {
    const name = width + ' || ' + height + ' || ' +
      molString + ' || ' + JSON.stringify(scaffolds) + ' || ' + highlightScaffold + ' || ' +
      molRegenerateCoords + ' || ' + scaffoldRegenerateCoords;

    return this.rendersCache.getOrCreate(name, (_: any) => this._rendererGetOrCreate(width, height,
      molString, scaffolds, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords, multiScaffold));
  }

  _drawMolecule(x: number, y: number, w: number, h: number, onscreenCanvas: HTMLCanvasElement,
    molString: string, scaffolds: IColoredScaffold[], highlightScaffold: boolean,
    molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean, cellStyle: DG.GridCellStyle, multiScaffold: boolean): void {
    const vertical = cellStyle !== undefined && cellStyle !== null ? cellStyle.textVertical : false;

    if (vertical) {
      h += w;
      w = h - w;
      h -= w;
    }
    const imageData = this._fetchRender(w, h, molString, scaffolds,
      highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords, multiScaffold);

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

  _initScaffoldString(colTemp: any, tagName: string): IColoredScaffold[] {
    let scaffoldString = colTemp ? colTemp[tagName] : null;
    if (scaffoldString?.endsWith(this.WHITE_MOLBLOCK_SUFFIX)) {
      if (colTemp[tagName])
        delete colTemp[tagName];
      return [];
    }
    return scaffoldString ? [{molecule: scaffoldString}] : [];
  }

  _initScaffoldArray(col: any, tagName: string, isTempCol?: boolean): IColoredScaffold[] {
    const scaffoldArrStr = !isTempCol ? col.getTag(tagName) : col ? col[tagName] : null;
    if (scaffoldArrStr) {
      const scaffoldArr: IColoredScaffold[] = JSON.parse(scaffoldArrStr);    
      const scaffoldArrFinal: IColoredScaffold[] = [];
      scaffoldArr.forEach((it) => {
        if (!it.molecule.endsWith(this.WHITE_MOLBLOCK_SUFFIX))
          scaffoldArrFinal.push(it);
      });
      if(!scaffoldArrFinal.length)
        col.setTag(tagName, '');
      return scaffoldArrFinal;
    }
    return [];
  }

  render(g: any, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    const molString = gridCell.cell.value;
    if (molString == null || molString === '')
      return;

    const r = window.devicePixelRatio;
    x = r * x; y = r * y;
    w = r * w; h = r * h;

    // value-based drawing (coming from HtmlCellRenderer.renderValue)
    if (gridCell.cell.column == null) {
      this._drawMolecule(x, y, w, h, g.canvas, molString, [], false, false, false, cellStyle, false);
      return;
    }

    const colTemp = gridCell.cell.column.temp;
    const highlightInfo = this.getHighlightTagInfo(colTemp, gridCell);
    
    // TODO: make both filtering scaffold and single highlight scaffold appear

    if (highlightInfo.scaffolds) {
      this._drawMolecule(x, y, w, h, g.canvas,
        molString, highlightInfo.scaffolds, true, false, false, cellStyle, highlightInfo.alighByFirstSubtruct);
    } else
      this.highlightByScaffoldCol(g, x, y, w, h, gridCell, cellStyle, colTemp, molString);
  }

  getHighlightTagInfo(colTemp: any, gridCell: DG.GridCell): IHighlightTagInfo {
    const filter = this._initScaffoldString(colTemp, FILTER_SCAFFOLD_TAG); //expected molBlock
    const align = this._initScaffoldString(colTemp, ALIGN_BY_SCAFFOLD_TAG);
    const highlight = this._initScaffoldArray(gridCell.cell.column, HIGHLIGHT_BY_SCAFFOLD_TAG);
    const scaffoldTreeHighlight = this._initScaffoldArray(gridCell.cell.column, SCAFFOLD_TREE_HIGHLIGHT);
    const alignByStructure = !!(filter.length || align.length || scaffoldTreeHighlight.length);
    const scaffolds = filter.length ? filter : align.length ? align : scaffoldTreeHighlight.length ? scaffoldTreeHighlight : highlight;
    return {scaffolds: scaffolds?.length ? scaffolds : undefined, alighByFirstSubtruct: alignByStructure};
  }

  highlightByScaffoldCol(g: any, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle, colTemp: any, molString: string): void {
    let molRegenerateCoords = colTemp && colTemp['regenerate-coords'] === 'true';
    let scaffoldRegenerateCoords = false;
    const df = gridCell.cell.dataFrame;
    let rowScaffoldCol = null;

    // if given, take the 'scaffold-col' col
    if (colTemp && colTemp['scaffold-col']) {
      const rowScaffoldColName = colTemp['scaffold-col'];
      const rowScaffoldColProbe = df.columns.byName(rowScaffoldColName);
      if (rowScaffoldColProbe !== null) {
        const scaffoldColTemp = rowScaffoldColProbe.temp;
        scaffoldRegenerateCoords = scaffoldColTemp && scaffoldColTemp['regenerate-coords'] === 'true';
        molRegenerateCoords = scaffoldRegenerateCoords;
        rowScaffoldCol = rowScaffoldColProbe;
      }
    }

    if (rowScaffoldCol == null || rowScaffoldCol.name === gridCell.cell.column.name) {
      // regular drawing
      this._drawMolecule(x, y, w, h, g.canvas, molString, [], false, molRegenerateCoords, false, cellStyle, false);
    } else {
      // drawing with a per-row scaffold
      const idx = gridCell.tableRowIndex; // TODO: supposed to be != null?
      const scaffoldMolString = df.get(rowScaffoldCol.name, idx!);
      const highlightScaffold = colTemp && colTemp['highlight-scaffold'] === 'true';
      this._drawMolecule(x, y, w, h, g.canvas,
        molString, [{molecule: scaffoldMolString}], highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords, cellStyle, false);
    }
  }
}
