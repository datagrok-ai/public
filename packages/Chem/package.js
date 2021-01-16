var rdKitModule = null;
var rdKitWorkerProxy = null;

//name: Chem
class ChemPackage extends DG.Package {

  constructor() {
    super();
    this.initialized = false;
  }

  /** Guaranteed to be executed exactly once before the execution of any function below */
  async init() {
    if (!initialized) {
      this.name = "Chem";
      rdKitModule = await initRDKitModule();
      console.log('RDKit (package) initialized');
      rdKitModule.prefer_coordgen(false);
      rdKitWorkerProxy = new RdKitWorkerProxy(this.webRoot);
      await rdKitWorkerProxy.moduleInit();
      this.STORAGE_NAME = 'rdkit_descriptors';
      this.KEY = 'selected';
      this.rdKitRendererCache = new DG.LruCache();
      this.rdKitRendererCache.onItemEvicted = (mol) => {
        mol.delete();
      };
      this.initialized = true;
    }
  }

  //name: SubstructureFilter
  //description: RDKit-based substructure filter
  //tags: filter
  //output: filter result
  substructureFilter() {
    return new SubstructureFilter();
  }

  _svgDiv(mol) {
    let root = ui.div();
    root.innerHTML = mol.get_svg();
    return root;
  }

  //name: getCLogP
  //input: string smiles {semType: Molecule}
  //output: double cLogP
  getCLogP(smiles) {
    let mol = rdKitModule.get_mol(smiles);
    return JSON.parse(mol.get_descriptors()).CrippenClogP;
  }

  //name: RDKitInfo
  //tags: panel, widgets
  //input: string smiles {semType: Molecule}
  //output: widget result
  rdkitInfoPanel(smiles) {
    let mol = rdKitModule.get_mol(smiles);
    return new DG.Widget(ui.divV([
      this._svgDiv(mol),
      ui.divText(`${this.getCLogP(smiles)}`)
    ]));
  }

  //name: rdkitCellRenderer
  //tags: cellRenderer, cellRenderer-Molecule
  //meta-cell-renderer-sem-type: Molecule
  //output: grid_cell_renderer result
  async rdkitCellRenderer() {
    let props = DG.toJs(await this.getProperties());
    if (props.Renderer && props.Renderer === 'RDKit') {
      return new RDKitCellRenderer();
    }
  }

  //name: similarityScoring
  //input: column molStringsColumn
  //input: string molString
  //input: bool sorted
  //output: dataframe result
  similarityScoring(molStringsColumn, molString, sorted) {
    try {
      let result = chemSimilarityScoring(molStringsColumn, molString, {'sorted': sorted});
      if (result == null) {
        return DG.DataFrame.create();
      }
      return (sorted ? result : DG.DataFrame.fromColumns([result]));
    } catch (e) {
      console.error("In similarityScoring: " + e.toString());
      throw e;
    }
  }

  //name: getSimilarities
  //input: column molStringsColumn
  //input: string molString
  //output: dataframe result
  getSimilarities(molStringsColumn, molString) {
    try {
      let result = chemGetSimilarities(molStringsColumn, molString);
      // TODO: get rid of a wrapping DataFrame and be able to return Columns
      return result ? DG.DataFrame.fromColumns([result]) : DG.DataFrame.create();
    } catch (e) {
      console.error("In getSimilarities: " + e.toString());
      throw e;
    }
  }

  //name: findSimilar
  //input: column molStringsColumn
  //input: string molString
  //input: int limit
  //input: int cutoff
  //output: dataframe result
  findSimilar(molStringsColumn, molString, aLimit, aCutoff) {
    try {
      let result = chemFindSimilar(molStringsColumn, molString, {limit: aLimit, cutoff: aCutoff});
      return result ? result : DG.DataFrame.create();
    } catch (e) {
      console.error("In getSimilarities: " + e.toString());
      throw e;
    }
  }

  //name: substructureSearch
  //input: column molStringsColumn
  //input: string molString
  //input: bool substructLibrary
  //output: column result
  async substructureSearch(molStringsColumn, molString, substructLibrary) {

    try {
      let result =
        substructLibrary ?
          await chemSubstructureSearchLibrary(molStringsColumn, molString) :
          chemSubstructureSearchGraph(molStringsColumn, molString);
      return DG.Column.fromList('object', 'bitset', [result]);
    } catch (e) {
      console.error("In substructureSearch: " + e.toString());
      throw e;
    }
  }

  //name: molColumnPropertyPanel
  //input: column molColumn
  //tags: panel
  //output: widget result
  molColumnPropertyPanel(molColumn) {
    return getMolColumnPropertyPanel(molColumn);
  }

  //tags: app
  descriptorsApp(context) {
    let defaultSmiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';
    let sketcherValue = defaultSmiles;

    let windows = grok.shell.windows;
    windows.showToolbox = false;
    windows.showHelp = false;
    windows.showProperties = false;

    let table = DG.DataFrame.create();
    table.name = 'Descriptors';
    let view = grok.shell.addTableView(table);

    let dsDiv = ui.divV([], 'grok-prop-panel');
    dsDiv.appendChild(this.descriptorsWidget(defaultSmiles).root);

    let sketcher = grok.chem.sketcher((smiles, molfile) => {
      sketcherValue = smiles;
      ChemPackage.removeChildren(dsDiv);
      dsDiv.appendChild(this.descriptorsWidget(smiles).root);
    }, defaultSmiles);
    let addButton = ui.bigButton('ADD', async () => {
      this.getSelected().then(selected => {
        grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'smiles', selected).then(t => {
          let columnNames = table.columns.names();
          if ((table.columns.length !== selected.length + 1) || selected.some(s => !columnNames.includes(s))) {
            table = DG.DataFrame.create();
            table.name = 'Descriptors';
            view.dataFrame = table;
            for (let col of t.columns.toList())
              table.columns.addNew(col.name, col.type);
          }
          table.rows.addNew(t.columns.toList().map(c => c.get(0)));
        });
      });
    });
    addButton.style.marginTop = '12px';
    let skDiv = ui.divV([sketcher, addButton], 'grok-prop-panel,dlg-sketcher,pure-form');

    let skNode = view.dockManager.dock(skDiv, DG.DOCK_TYPE.RIGHT, null, 'Sketcher', 0.25);
    view.dockManager.dock(dsDiv, DG.DOCK_TYPE.DOWN, skNode, 'Descriptors', 0.5);

    grok.events.onViewRemoved.subscribe((v) => {
      if (v.name === view.name) {
        windows.showToolbox = true;
        windows.showHelp = true;
        windows.showProperties = true;
      }
    });
  }

  //name: Chem Descriptors
  //tags: panel, widgets
  //input: string smiles { semType: Molecule }
  //output: widget result
  descriptorsWidget(smiles) {
    let widget = new DG.Widget(ui.div());
    let result = ui.div();
    let selectButton = ui.bigButton('SELECT', async () => {
      ChemPackage.openDescriptorsDialog(await this.getSelected(), async (selected) => {
        await grok.dapi.userDataStorage.postValue(this.STORAGE_NAME, this.KEY, JSON.stringify(selected));
        update();
      });
    });
    selectButton.style.marginTop = '20px';

    let update = () => {
      ChemPackage.removeChildren(result);
      result.appendChild(ui.loader());
      this.getSelected().then(selected => {
        grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${smiles}`), 'smiles', selected).then(table => {
          ChemPackage.removeChildren(result);
          let map = {};
          for (let descriptor of selected)
            map[descriptor] = table.col(descriptor).get(0);
          result.appendChild(ui.tableFromMap(map));
        });
      });
    }

    widget.root.appendChild(result);
    widget.root.appendChild(selectButton);

    update();

    return widget;
  }

  //description: Get selected descriptors
  async getSelected() {
    let str = await grok.dapi.userDataStorage.getValue(this.STORAGE_NAME, this.KEY);
    let selected = (str != null && str !== '') ? JSON.parse(str) : [];
    if (selected.length === 0) {
      selected = (await grok.chem.descriptorsTree())['Lipinski']['descriptors'].slice(0, 3).map(p => p['name']);
      await grok.dapi.userDataStorage.postValue(this.STORAGE_NAME, this.KEY, JSON.stringify(selected));
    }
    return selected;
  }

  //description: Open descriptors selection dialog
  static openDescriptorsDialog(selected, onOK) {
    grok.chem.descriptorsTree().then(descriptors => {
      let tree = ui.tree();
      tree.root.style.maxHeight = '400px';

      let groups = {};
      let items = [];

      for (let groupName in descriptors) {
        let group = tree.group(groupName, null, false);
        group.enableCheckBox();
        groups[groupName] = group;

        for (let descriptor of descriptors[groupName]['descriptors']) {
          let item = group.item(descriptor['name'], descriptor);
          item.enableCheckBox(selected.includes(descriptor['name']));
          items.push(item);
        }
      }

      let clear = ui.button('NONE', () => {
        for (let g in groups) groups[g].checked = false;
        for (let i of items) i.checked = false;
      });

      ui.dialog('Chem Descriptors')
        .add(clear)
        .add(tree.root)
        .onOK(() => onOK(items.filter(i => i.checked).map(i => i.value['name'])))
        .show();
    });
  }

  renderMolRdKitCanvasCache(molString, canvas, x, y, w, h) {

    let mol = this.rdKitRendererCache.getOrCreate(molString, (s) => {
      try {
        return rdKitModule.get_mol(s);
      } catch (e) {
        return rdKitModule.get_mol("");
      }
    });

    const opts = {
      "clearBackground": false,
      "offsetx": Math.floor(x),
      "offsety": -Math.floor(y),
      "width": Math.floor(w),
      "height": Math.floor(h)
    };
    mol.draw_to_canvas_with_highlights(canvas, JSON.stringify(opts));

  }

  get rdKitModule() {

    return rdKitModule;

  }

  //description: Removes all children from node
  static removeChildren(node) {
    while (node.firstChild)
      node.removeChild(node.firstChild);
  }
}