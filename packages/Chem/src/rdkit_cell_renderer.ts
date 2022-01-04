/**
 * RDKit-based molecule cell renderer.
 * */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {convertToRDKit} from './analysis/r_group';

export class RDKitCellRenderer extends DG.GridCellRenderer {
  readonly WHITE_MOLBLOCK_SUFFIX = `
  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
`;
  rdKitModule: any;
  canvasCounter: number;
  molCache: DG.LruCache = new DG.LruCache();
  rendersCache: DG.LruCache = new DG.LruCache();

  constructor(rdKitModule: any) {
    super();
    this.rdKitModule = rdKitModule;
    this.canvasCounter = 0;

    this.molCache.onItemEvicted = function(obj: {[_ : string]: any} | null) {
      obj!.mol?.delete();
      obj!.mol = null;
      obj!.substruct = null;
      obj = null; // ? GC definitely delete
    };

    this.rendersCache.onItemEvicted = function(obj: {[_ : string]: any} | null) {
      obj!.canvas = null;
      obj = null; // ? GC definitely delete
    };
  }

  get name() {
    return 'RDKit cell renderer';
  }
  get cellType() {
    return DG.SEMTYPE.MOLECULE;
  }
  get defaultWidth() {
    return 200;
  }
  get defaultHeight() {
    return 100;
  }

  _isMolBlock(molString: string) {
    return molString.includes('M  END');
  }

  _fetchMolGetOrCreate(molString: string, scaffoldMolString: string, molRegenerateCoords: boolean) {
    let mol = null;
    let substructJson = '{}';

    try {
      mol = this.rdKitModule.get_mol(molString);
    } catch (e) {
      try {
        mol = this.rdKitModule.get_mol(molString, '{"kekulize":false}');
      } catch (e2) {
        console.error(
          'In _fetchMolGetOrCreate: RDKit .get_mol crashes on a molString: `' + molString + '`');
        mol = null;
      }
    }
    if (mol) {
      try {
        if (mol.is_valid()) {
          const scaffoldIsMolBlock = this._isMolBlock(scaffoldMolString);
          if (scaffoldIsMolBlock) {
            const rdkitScaffoldMol = this._fetchMol(scaffoldMolString, '', molRegenerateCoords, false).mol;
            if (rdkitScaffoldMol && rdkitScaffoldMol.is_valid()) {
              substructJson = mol.generate_aligned_coords(rdkitScaffoldMol, true, true, false);
              if (substructJson === '') {
                substructJson = '{}';
              }
            }
          } else if (molRegenerateCoords) {
            const molBlock = mol.get_new_coords(true);
            mol.delete();
            mol = this.rdKitModule.get_mol(molBlock);
          }
          if (!scaffoldIsMolBlock || molRegenerateCoords) {
            mol.normalize_2d_molblock();
            mol.straighten_2d_layout();
          }
        }
        if (!mol.is_valid()) {
          console.error(
            'In _fetchMolGetOrCreate: RDKit mol is invalid on a molString molecule: `' + molString + '`');
          mol.delete();
          mol = null;
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
    const rdkitMol = fetchMolObj.mol;
    const substruct = fetchMolObj.substruct;

    const canvasId = '_canvas-rdkit-' + this.canvasCounter;
    const canvas = new OffscreenCanvas(width, height);
    this.canvasCounter++;
    if (rdkitMol != null) {
      this._drawMoleculeToCanvas(rdkitMol, width, height, canvas, substruct, highlightScaffold);
    } else {
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

  _drawMoleculeToCanvas(
    rdkitMol: any, w: number, h: number, canvas: OffscreenCanvas,
    substruct: Object, highlightScaffold: boolean) {
    const opts = {
      'clearBackground': false,
      'offsetx': 0, 'offsety': 0,
      'width': Math.floor(w),
      'height': Math.floor(h),
      'bondLineWidth': 1,
      'fixedScale': 0.07,
      'minFontSize': 9,
      'highlightBondWidthMultiplier': 12,
      'dummyIsotopeLabels': false,
      'atomColourPalette': {
        16: [0.498, 0.247, 0.0],
        9: [0.0, 0.498, 0.0],
        17: [0.0, 0.498, 0.0],
      },
    };
    if (highlightScaffold) {
      Object.assign(opts, substruct);
    }
    rdkitMol.draw_to_canvas_with_highlights(canvas, JSON.stringify(opts));
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
      if (colTemp[tagName]) {
        delete colTemp[tagName];
      }
    }
    return scaffoldString;
  }

  render(g: any, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: any) {
    const molString = convertToRDKit(gridCell.cell.value);
    if (molString == null || molString === '') {
      return;
    }

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
