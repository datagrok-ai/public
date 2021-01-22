/**
 * RDKit-based molecule cell renderer.
 * */
class RDKitCellRenderer extends DG.GridCellRenderer {

  constructor() {

    super();
    this.canvasCounter = 0;
    this.emptyMol = rdKitModule.get_mol("");
    this.molCache = new DG.LruCache();
    this.molCache.onItemEvicted = function (mol) {
      mol.delete();
      mol = null;
    };
    this.rendersCache = new DG.LruCache();
    this.rendersCache.onItemEvicted = function (obj) {
      obj.canvas = null;
      obj.canvasId = null;
      obj = null;
    }

  }

  get name() { return 'RDKit cell renderer'; }
  get cellType() { return DG.SEMTYPE.MOLECULE; }
  get defaultWidth() { return 200; }
  get defaultHeight() { return 100; }

  _molIsInMolBlock(molString, rdkitMol) {

    const smilesMolString = rdkitMol.get_smiles();
    if (smilesMolString === molString)
      return false;
    const cxsmilesMolString = rdkitMol.get_cxsmiles();
    if (cxsmilesMolString === molString)
      return false;
    const inchiMolString = rdkitMol.get_inchi();
    if (inchiMolString === molString)
      return false;
    return true;

  }

  _fetchMol(molString, scaffoldMolString = "") {
    const name = molString + " || " + scaffoldMolString;
    const emptyMol = this.emptyMol;
    return this.molCache.getOrCreate(name, function (s) {
      let mol = emptyMol;
      try {
        mol = rdKitModule.get_mol(molString);
        if (!mol.is_valid()) {
          mol = emptyMol;
        }
      } catch (e) {
        console.error(
          "Possibly a malformed molecule (rendering, no scaffold): `" + s + "`");
      }
      return mol;
    });
  }

  _fetchRender(width, height, molString, scaffoldMolString = "") {

    const name = width + " || " + height + " || " +  molString + " || " + scaffoldMolString;
    return this.rendersCache.getOrCreate(name, (s) => {

      let rdkitMol = this._fetchMol(molString, scaffoldMolString);
      let rdkitScaffoldMol = this._fetchMol(scaffoldMolString);

      if (scaffoldMolString !== "") {
        try {
          if (this._molIsInMolBlock(scaffoldMolString, rdkitScaffoldMol)) {
            const substructJson = rdkitMol.get_substruct_match(rdkitScaffoldMol);
            if (substructJson !== '{}') {
              rdkitMol.generate_aligned_coords(rdkitScaffoldMol, true);
            }
          }
        } catch (e) {
          console.error(
            "Possibly a malformed molecule (rendering, scaffolds): `" + s + "`");
        }
      }
      const canvasId = '_canvas-rdkit-' + this.canvasCounter;
      let canvas = window.document.createElement('canvas');
      canvas.setAttribute('id', canvasId);
      this.canvasCounter++;
      this._drawMoleculeToCanvas(rdkitMol, width, height, canvas);

      return { canvas: canvas, canvasId: canvasId };
    });
  }

  _drawMoleculeToCanvas(rdkitMol, w, h, canvas) {
    const opts = {
      "clearBackground": false,
      "offsetx": 0, "offsety": 0,
      "width": Math.floor(w),
      "height": Math.floor(h),
      "bondLineWidth": 1,
      "minFontSize": 11
    }
    canvas.width = w;
    canvas.height = h;
    rdkitMol.draw_to_canvas_with_highlights(canvas, JSON.stringify(opts));
  }

  _drawMolecule(x, y, w, h, onscreenCanvas, molString, scaffoldMolString = "") {

    const renderObj = this._fetchRender(w, h, molString, scaffoldMolString);
    let offscreenCanvas = renderObj.canvas;
    let image = offscreenCanvas.getContext('2d').getImageData(0, 0, w, h);
    let context = onscreenCanvas.getContext('2d');
    context.putImageData(image, x, y);

  }

  render(g, x, y, w, h, gridCell, cellStyle) {

    let molString = gridCell.cell.value;
    if (molString == null || molString === '')
      return;

    const molCol = gridCell.tableColumn.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    let singleScaffoldMolString = molCol ? molCol.tags['chem-scaffold'] : null;

    if (singleScaffoldMolString) {

      this._drawMolecule(x, y, w, h, g.canvas, molString, singleScaffoldMolString);

    } else {

      let df = gridCell.tableColumn.dataFrame;
      const rowScaffoldCol = (() => {

        // if given, take the 'scaffold-col' col
        let colTags = gridCell.tableColumn.tags;
        if (colTags && colTags['scaffold-col']) {
          let rowScaffoldColName = colTags['scaffold-col'];
          let rowScaffoldColProbe = df.columns.byName(rowScaffoldColName);
          if (rowScaffoldColProbe !== null) {
            return rowScaffoldColProbe;
          }
        }
        return null;

      })();

      if (rowScaffoldCol == null || rowScaffoldCol.name === gridCell.tableColumn.name) {
        // regular drawing
        this._drawMolecule(x, y, w, h, g.canvas, molString);
      } else {
        // drawing with a per-row scaffold
        let idx = gridCell.tableRowIndex;
        let scaffoldMolString = df.get(rowScaffoldCol.name, idx);
        this._drawMolecule(x, y, w, h, g.canvas, molString, scaffoldMolString);
      }
    }
  }
}