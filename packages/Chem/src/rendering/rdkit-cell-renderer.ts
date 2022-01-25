/**
 * RDKit-based molecule cell renderer.
 * */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {convertToRDKit} from '../analysis/r-group-analysis';
import {drawRdKitMoleculeToOffscreenCanvas} from '../utils/chem-common-rdkit';
import {RDModule, RDMol} from '../rdkit-api';
import {isMolBlock} from '../utils/chem-utils';

export class GridCellRendererProxy extends DG.GridCellRenderer {
  renderer: DG.GridCellRenderer;
  _cellType: string;

  constructor(renderer: DG.GridCellRenderer, cellType: string) {
    super();
    this.renderer = renderer;
    this._cellType = cellType;
  }

  get defaultWidth(): number | null { return this.renderer.defaultWidth;  }
  get defaultHeight(): number | null { return this.renderer.defaultHeight; }

  get name(): string { return this.renderer.name; }
  get cellType(): string { return this._cellType; }

  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    this.renderer.render(g, x, y, w, h, gridCell, cellStyle);
  }


  renderInternal(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
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
  molCache: DG.LruCache = new DG.LruCache();
  rendersCache: DG.LruCache = new DG.LruCache();

  constructor(rdKitModule: RDModule) {
    super();
    this.rdKitModule = rdKitModule;
    this.canvasCounter = 0;

    this.molCache.onItemEvicted = function(obj: {[_ : string]: any} | null) {
      obj!.mol?.delete();
      obj!.mol = null;
      obj!.substruct = null;
    };

    this.rendersCache.onItemEvicted = function(obj: {[_ : string]: any} | null) {
      obj!.canvas = null;
    };
  }

  get name() {return 'RDKit cell renderer';}
  get cellType() {return DG.SEMTYPE.MOLECULE;}
  get defaultWidth() {return 200;}
  get defaultHeight() {return 100;}

  _fetchMolGetOrCreate(molString: string, scaffoldMolString: string, molRegenerateCoords: boolean) {
    let mol: RDMol | null = null;
    let substructJson = '{}';
    try {
      mol = this.rdKitModule.get_mol(molString);
    } catch (e) {
      try {
        mol = this.rdKitModule.get_mol(molString, '{"kekulize":false}');
      } catch (e2) {
        console.error(
          'Chem | In _fetchMolGetOrCreate: RDKit .get_mol crashes on a molString: `' + molString + '`');
        mol = null;
      }
    }
    if (mol) {
      try {
        if (mol.is_valid()) {
          const scaffoldIsMolBlock = isMolBlock(scaffoldMolString);
          if (scaffoldIsMolBlock) {
            const rdKitScaffoldMol = this._fetchMol(scaffoldMolString, '', molRegenerateCoords, false).mol;
            if (rdKitScaffoldMol && rdKitScaffoldMol.is_valid()) {
              substructJson = mol.generate_aligned_coords(rdKitScaffoldMol, true, true, false);
              if (substructJson === '')
                substructJson = '{}';
            }
          } else if (molRegenerateCoords) {
            const molBlock = mol.get_new_coords(true);
            mol.delete();
            mol = this.rdKitModule.get_mol(molBlock);
          }
          if (!scaffoldIsMolBlock || molRegenerateCoords) {
            mol!.normalize_2d_molblock();
            mol!.straighten_2d_layout();
          }
        }
        if (!mol!.is_valid()) {
          console.error(
            'Chem | In _fetchMolGetOrCreate: RDKit mol is invalid on a molString molecule: `' + molString + '`');
          mol!.delete();
        }
      } catch (e) {
        console.error(
          'In _fetchMolGetOrCreate: RDKit crashed, possibly a malformed molString molecule: `' + molString + '`');
      }
    }
    return {mol: mol, substruct: JSON.parse(substructJson), molString: molString, scaffoldMolString: scaffoldMolString};
  }

  _fetchMol(
    molString: string, scaffoldMolString: string,
    molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean) {
    const name = molString + ' || ' + scaffoldMolString + ' || ' +
      molRegenerateCoords + ' || ' + scaffoldRegenerateCoords;
    return this.molCache.getOrCreate(name, (_: any) =>
      this._fetchMolGetOrCreate(molString, scaffoldMolString, molRegenerateCoords));
  }

  _rendererGetOrCreate(
    width: number, height: number, molString: string, scaffoldMolString: string,
    highlightScaffold: boolean, molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean) {
    const fetchMolObj = this._fetchMol(molString, scaffoldMolString, molRegenerateCoords, scaffoldRegenerateCoords);
    const rdKitMol = fetchMolObj.mol;
    const substruct = fetchMolObj.substruct;

    const canvasId = '_canvas-rdkit-' + this.canvasCounter;
    const canvas = new OffscreenCanvas(width, height);
    this.canvasCounter++;
    if (rdKitMol != null)
      drawRdKitMoleculeToOffscreenCanvas(rdKitMol, width, height, canvas, highlightScaffold ? substruct : null);
    else {
      // draw a crossed rectangle
      const ctx = canvas.getContext('2d');
      ctx!.lineWidth = 1;
      ctx!.strokeStyle = '#EFEFEF';
      ctx!.beginPath();
      ctx!.moveTo(0, 0);
      ctx!.lineTo(width, height);
      ctx!.stroke();
      ctx!.beginPath();
      ctx!.moveTo(width, 0);
      ctx!.lineTo(0, height);
      ctx!.stroke();
    }
    return {canvas: canvas};
  }

  _fetchRender(
    width: number, height: number, molString: string, scaffoldMolString: string,
    highlightScaffold: boolean, molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean) {
    const name = width + ' || ' + height + ' || ' +
      molString + ' || ' + scaffoldMolString + ' || ' + highlightScaffold + ' || ' +
      molRegenerateCoords + ' || ' + scaffoldRegenerateCoords;
    return this.rendersCache.getOrCreate(name, (_: any) =>
      this._rendererGetOrCreate(width, height,
        molString, scaffoldMolString, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords));
  }

  _drawMolecule(x: number, y: number, w: number, h: number, onscreenCanvas: HTMLCanvasElement,
    molString: string, scaffoldMolString: string, highlightScaffold: boolean,
    molRegenerateCoords: boolean, scaffoldRegenerateCoords: boolean) {
    const r = window.devicePixelRatio;
    x = r * x; y = r * y;
    w = r * w; h = r * h;
    const renderObj = this._fetchRender(w, h,
      molString, scaffoldMolString, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords);
    const offscreenCanvas = renderObj.canvas;
    const image = offscreenCanvas.getContext('2d').getImageData(0, 0, w, h);
    const context = onscreenCanvas.getContext('2d');
    context!.putImageData(image, x, y);
  }

  _initScaffoldString(colTemp: any, tagName: string) {
    let scaffoldString = colTemp ? colTemp[tagName] : null;
    if (scaffoldString?.endsWith(this.WHITE_MOLBLOCK_SUFFIX)) {
      scaffoldString = null;
      if (colTemp[tagName])
        delete colTemp[tagName];
    }
    return scaffoldString;
  }

  render(g: any, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: any) {
    const molString = convertToRDKit(gridCell.cell.value);
    if (molString == null || molString === '')
      return;


    // value-based drawing (coming from HtmlCellRenderer.renderValue)
    if (gridCell.cell.column == null) {
      this._drawMolecule(x, y, w, h, g.canvas, molString, '', false, false, false);
      return;
    }

    const colTemp = gridCell.cell.column.temp;

    const singleScaffoldHighlightMolString = this._initScaffoldString(colTemp, 'chem-scaffold');
    const singleScaffoldFilterMolString = this._initScaffoldString(colTemp, 'chem-scaffold-filter'); // expected a molBlock
    const singleScaffoldMolString = singleScaffoldFilterMolString ?? singleScaffoldHighlightMolString;
    // TODO: make both filtering scaffold and single highlight scaffold appear

    if (singleScaffoldMolString) {
      this._drawMolecule(x, y, w, h, g.canvas,
        molString, singleScaffoldMolString, true, false, false);
    } else {
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
        this._drawMolecule(x, y, w, h, g.canvas, molString, '', false, molRegenerateCoords, false);
      } else {
        // drawing with a per-row scaffold
        const idx = gridCell.tableRowIndex; // TODO: supposed to be != null?
        const scaffoldMolString = df.get(rowScaffoldCol.name, idx!);
        const highlightScaffold = colTemp && colTemp['highlight-scaffold'] === 'true';
        this._drawMolecule(x, y, w, h, g.canvas,
          molString, scaffoldMolString, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords);
      }
    }
  }
}
