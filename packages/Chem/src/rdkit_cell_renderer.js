/**
 * RDKit-based molecule cell renderer.
 * */
class RDKitCellRenderer extends DG.GridCellRenderer {

  constructor() {
    super();
    this.molCache = new DG.LruCache();
    this.molCache.onItemEvicted = function (mol) {
      mol.delete();
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

  render(g, x, y, w, h, gridCell, cellStyle) {

    let value = gridCell.cell.value;
    if (value == null || value === '')
      return;

    let molCache = this.molCache;

    const fetchMolNoScaffold = function (molString) {
      return molCache.getOrCreate(molString, (s) => {
        try {
          let mol = rdKitModule.get_mol(s);
          return mol.is_valid() ? mol : rdKitModule.get_mol("");
        } catch (e) {
          return rdKitModule.get_mol("");
        }
      });
    }

    const fetchMolOnScaffold = function (molString, mol, scaffoldMolString, scaffoldMol) {
      return molCache.getOrCreate(molString + " || " + scaffoldMolString, (s) => {
        try {
          if (molIsInMolBlock(scaffoldMolString, scaffoldMol)) {
            const substructJson = mol.get_substruct_match(scaffoldMol);
            if (substructJson !== '{}') {
              let newMol = rdKitModule.get_mol(molString);
              newMol.generate_aligned_coords(scaffoldMol, true);
              return newMol;
            }
          }
          return mol;
        } catch (e) {
          return mol;
        }
      });
    }

    let mol = fetchMolNoScaffold(value);

    if (!mol.is_valid())
      return;

    let drawMolecule = function (rdkitMol) {
      //rdkitMol.draw_to_canvas_with_offset(g.canvas, x, -y, w, h);
      const opts = {
        "clearBackground": false,
        "offsetx": Math.floor(x),
        "offsety": -Math.floor(y),
        "width": Math.floor(w),
        "height": Math.floor(h),
        "bondLineWidth": 1,
        "minFontSize": 11
      }
      rdkitMol.draw_to_canvas_with_highlights(g.canvas, JSON.stringify(opts));
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

    const drawMoleculeWithScaffold = function (scaffoldMolString, rdkitMol, rdkitMolSmiles, scaffoldCache) {

      let scaffoldMol = fetchMolNoScaffold(scaffoldMolString);
      const molWithScaffold = fetchMolOnScaffold(rdkitMolSmiles, rdkitMol, scaffoldMolString, scaffoldMol);
      drawMolecule(molWithScaffold);

    }

    const molCol = gridCell.tableColumn.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    let singleScaffoldMolString = molCol ? molCol.tags['chem-scaffold'] : null;

    if (singleScaffoldMolString) {
      drawMoleculeWithScaffold(singleScaffoldMolString, mol, value, this.molScaffoldCache);
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

        // TODO: deprecate
        // otherwise, find the 'row-scaffold' col
        for (let j = 0; j < df.columns.length; ++j) {
          let col = df.columns.byIndex(j);
          let tags = col.tags;
          if (tags && tags['row-scaffold']) {
            return col;
          }
        }

        return null;

      })();

      if (rowScaffoldCol == null || rowScaffoldCol.name === gridCell.tableColumn.name) {
        // regular drawing
        drawMolecule(mol);
      } else {
        let idx = gridCell.tableRowIndex;
        let scaffoldMolString = df.get(rowScaffoldCol.name, idx);
        drawMoleculeWithScaffold(scaffoldMolString, mol, value, this.molScaffoldCache);
      }
    }
  }
}