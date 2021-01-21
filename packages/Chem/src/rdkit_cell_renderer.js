/**
 * RDKit-based molecule cell renderer.
 * */
class RDKitCellRenderer extends DG.GridCellRenderer {

  constructor() {
    super();

    this.canvasCounter = 0;
    this.molCache = new DG.LruCache();
    this.molCache.onItemEvicted = function (mol) {
      mol?.delete();
      mol = null;
    };
    this.rendersCache = new DG.LruCache();
    this.rendersCache.onItemEvicted = function (obj) {
      let element = document.getElementById(obj.canvasId);
      element?.parentNode?.removeChild(element);
      element = null;
      obj.canvasId = null;
      obj.canvas = null;
      obj = null;
    }
  }

  get name() { return 'RDKit cell renderer'; }
  get cellType() { return DG.SEMTYPE.MOLECULE; }
  get defaultWidth() { return 200; }
  get defaultHeight() { return 100; }

  render(g, x, y, w, h, gridCell, cellStyle) {

    const emptyMol = rdKitModule.get_mol("");
    let canvasCounter = this.canvasCounter;
    let molString = gridCell.cell.value;
    if (molString == null || molString === '')
      return;

    let molCache = this.molCache;
    let rendersCache = this.rendersCache;

    const fetchMol = function (molString) {
      return molCache.getOrCreate(molString, (s) => {
        let mol = emptyMol;
        try {
          mol = rdKitModule.get_mol(s);
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

    const drawMoleculeToCanvas = function (rdkitMol, w, h, canvas) {
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

    const fetchRender = function (width, height, molString, scaffoldMolString = "") {

      const name = width + " || " + height + " || " +  molString + " || " + scaffoldMolString;
      return rendersCache.getOrCreate(name, (s) => {

        let rdkitMol = fetchMol(molString);
        let rdkitScaffoldMol = fetchMol(scaffoldMolString);

        if (scaffoldMolString !== "") {
          try {
            if (molIsInMolBlock(scaffoldMolString, rdkitScaffoldMol)) {
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
        const canvasId = '_canvas-rdkit-' + canvasCounter;
        let canvas = window.document.createElement('canvas');
        canvas.setAttribute('id', canvasId);
        canvasCounter++;
        drawMoleculeToCanvas(rdkitMol, w, h, canvas);
        return { canvas: canvas, canvasId: canvasId };
      });
    }

    const drawMolecule = function (molString, scaffoldMolString = "") {

      const renderObj = fetchRender(w, h, molString, scaffoldMolString);
      let offscreenCanvas = renderObj.canvas;
      let image = offscreenCanvas.getContext('2d').getImageData(0, 0, w, h);
      let context = g.canvas.getContext('2d');
      context.putImageData(image, x, y);

    }

    let molIsInMolBlock = function (molString, rdkitMol) {
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

    const molCol = gridCell.tableColumn.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    let singleScaffoldMolString = molCol ? molCol.tags['chem-scaffold'] : null;

    if (singleScaffoldMolString) {

      drawMolecule(molString, singleScaffoldMolString);

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
        drawMolecule(molString);
      } else {
        let idx = gridCell.tableRowIndex;
        let scaffoldMolString = df.get(rowScaffoldCol.name, idx);
        drawMolecule(molString, scaffoldMolString);
      }
    }
  }
}