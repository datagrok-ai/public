/**
 * RDKit-based molecule cell renderer.
 * */
class RDKitCellRenderer extends DG.GridCellRenderer {

  constructor() {

    super();
    this.canvasCounter = 0;
    this.molCache = new DG.LruCache();
    this.molCache.onItemEvicted = function (obj) {
      obj.mol.delete();
      obj.mol = null;
      obj.substruct = null;
      obj = null;
    };
    this.rendersCache = new DG.LruCache();
    this.rendersCache.onItemEvicted = function (obj) {
      obj.canvas = null;
      obj = null;
    }

  }

  get name() { return 'RDKit cell renderer'; }
  get cellType() { return DG.SEMTYPE.MOLECULE; }
  get defaultWidth() { return 200; }
  get defaultHeight() { return 100; }

  _isMolBlock(molString) {

    return molString.includes('M  END');

  }

  _fetchMolGetOrCreate(molString, scaffoldMolString, molRegenerateCoords) {

    let mol = null;
    let substructJson = "{}";

    try {
      mol = rdKitModule.get_mol(molString);
    } catch (e) {
      console.error(
        "In _fetchMolGetOrCreate: RDKit .get_mol crashes on a molString: `" + molString + "`");
      mol = null;
    }
    try {
      if (mol.is_valid()) {
        if (this._isMolBlock(scaffoldMolString)) {
          let rdkitScaffoldMol = this._fetchMol(scaffoldMolString, "", molRegenerateCoords, false).mol;
          substructJson = mol.generate_aligned_coords(rdkitScaffoldMol, true, true);
        }
        if (molRegenerateCoords || (substructJson === "" && !this._isMolBlock(molString))) {
          substructJson = "{}";
          let molBlock = mol.get_new_coords(true);
          mol.delete();
          mol = rdKitModule.get_mol(molBlock);
        }
      }
      if (!mol.is_valid()) {
        console.error(
          "In _fetchMolGetOrCreate: RDKit mol is invalid on a molString molecule: `" + molString + "`");
        mol = null;
      }
    } catch (e) {
      console.error(
        "In _fetchMolGetOrCreate: RDKit crashed, possibly a malformed molString molecule: `" + molString + "`");
    }
    return { mol: mol, substruct: JSON.parse(substructJson) };
  }

  _fetchMol(molString, scaffoldMolString, molRegenerateCoords, scaffoldRegenerateCoords) {
    const name = molString + " || " + scaffoldMolString + " || "
      + molRegenerateCoords + " || " + scaffoldRegenerateCoords;
    return this.molCache.getOrCreate(name, (s) =>
      this._fetchMolGetOrCreate(molString, scaffoldMolString, molRegenerateCoords));
  }
  
  _rendererGetOrCreate(width, height, molString, scaffoldMolString, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords) {

    let fetchMolObj = this._fetchMol(molString, scaffoldMolString, molRegenerateCoords, scaffoldRegenerateCoords);
    let rdkitMol = fetchMolObj.mol;
    const substruct = fetchMolObj.substruct;

    const canvasId = '_canvas-rdkit-' + this.canvasCounter;
    let canvas = new OffscreenCanvas(width, height);
    this.canvasCounter++;
    if (rdkitMol != null) {
      this._drawMoleculeToCanvas(rdkitMol, width, height, canvas, substruct, highlightScaffold);
    } else {
      let ctx = canvas.getContext("2d");
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
    return {canvas: canvas};

  }

  _fetchRender(width, height, molString, scaffoldMolString, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords) {

    const name = width + " || " + height + " || "
      + molString + " || " + scaffoldMolString  + " || " + highlightScaffold + " || "
      + molRegenerateCoords + " || " + scaffoldRegenerateCoords;
    return this.rendersCache.getOrCreate(name, (s) =>
      this._rendererGetOrCreate(width, height,
        molString, scaffoldMolString, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords));

  }

  _drawMoleculeToCanvas(rdkitMol, w, h, canvas, substruct, highlightScaffold) {
    let opts = {
      "clearBackground": false,
      "offsetx": 0, "offsety": 0,
      "width": Math.floor(w),
      "height": Math.floor(h),
      "bondLineWidth": 1,
      "minFontSize": 11,
      'highlightBondWidthMultiplier': 12 
    };
    if (highlightScaffold) {
      Object.assign(opts, substruct);
    }
    rdkitMol.draw_to_canvas_with_highlights(canvas, JSON.stringify(opts));
  }

  _drawMolecule(x, y, w, h, onscreenCanvas,
                molString, scaffoldMolString, highlightScaffold,
                molRegenerateCoords, scaffoldRegenerateCoords) {

    const r = window.devicePixelRatio;
    x = r * x; y = r * y;
    w = r * w; h = r * h;
    const renderObj = this._fetchRender(w, h,
      molString, scaffoldMolString, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords);
    let offscreenCanvas = renderObj.canvas;
    let image = offscreenCanvas.getContext('2d').getImageData(0, 0, w, h);
    let context = onscreenCanvas.getContext('2d');
    context.putImageData(image, x, y);
  }

  render(g, x, y, w, h, gridCell, cellStyle) {

    let molString = gridCell.cell.value;
    if (molString == null || molString === '')
      return;

    const colTags = gridCell.tableColumn.tags;
    let singleScaffoldMolString = colTags && colTags['chem-scaffold'];
    
    if (singleScaffoldMolString) {

      this._drawMolecule(x, y, w, h, g.canvas,
        molString, singleScaffoldMolString, true, false, false);

    } else {

      let molRegenerateCoords = colTags && colTags['regenerate-coords'] === 'true';
      let scaffoldRegenerateCoords = false;
      let df = gridCell.grid.dataFrame;
      let rowScaffoldCol = null;

      // if given, take the 'scaffold-col' col
      if (colTags && colTags['scaffold-col']) {
        let rowScaffoldColName = colTags['scaffold-col'];
        let rowScaffoldColProbe = df.columns.byName(rowScaffoldColName);
        if (rowScaffoldColProbe !== null) {
          const scaffoldColTags = rowScaffoldColProbe.tags;
          scaffoldRegenerateCoords = scaffoldColTags && scaffoldColTags['regenerate-coords'] === 'true';
          molRegenerateCoords = scaffoldRegenerateCoords;
          rowScaffoldCol = rowScaffoldColProbe;
        }
      }

      if (rowScaffoldCol == null || rowScaffoldCol.name === gridCell.tableColumn.name) {
        // regular drawing
        this._drawMolecule(x, y, w, h, g.canvas, molString, "", false, molRegenerateCoords, false);
      } else {
        // drawing with a per-row scaffold
        let idx = gridCell.tableRowIndex;
        let scaffoldMolString = df.get(rowScaffoldCol.name, idx);
        let highlightScaffold = colTags && colTags['highlight-scaffold'] === 'true';
        this._drawMolecule(x, y, w, h, g.canvas,
          molString, scaffoldMolString, highlightScaffold, molRegenerateCoords, scaffoldRegenerateCoords);
      }
    }
  }
}