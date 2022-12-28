import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

export let _package = new DG.Package();


export async function chemblSubstructureSearch(molecule: string): Promise<DG.DataFrame> {
    try {
        const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
        const smarts = mol.get_smarts();
        mol?.delete();
        let df = await grok.data.query(`${_package.name}:patternSubstructureSearch`, {'pattern': smarts, 'maxRows': 100});
        if (df == null)
            return DG.DataFrame.create();
        return df;
    } catch (e: any) {
        console.error('In SubstructureSearch: ' + e.toString());
        throw e;
    }
}

export async function chemblSimilaritySearch(molecule: string): Promise<DG.DataFrame> {
    try {
        const mol = (await grok.functions.call('Chem:getRdKitModule')).get_mol(molecule);
        const smiles = mol.get_smiles();
        mol?.delete();
        let df = await grok.data.query(`${_package.name}:patternSimilaritySearch`, {'pattern': smiles, 'maxRows': 100});
        if (df == null)
            return DG.DataFrame.create();
        return df;
    } catch (e: any) {
        console.error('In SimilaritySearch: ' + e.toString());
        throw e;
    }

}

//name: chemblSearchWidgetLocalDb
//tags: widgets
//input: string mol {semType: Molecule}
//input: string searchType
//output: widget result
export function chemblSearchWidgetLocalDb(mol: string, substructure: boolean = false): DG.Widget {
    const headerHost = ui.divH([]);
    const compsHost = ui.divH([ui.loader(), headerHost]);
    const panel = ui.divV([compsHost]);
    const searchFunc = substructure ?
      async () => chemblSubstructureSearch(mol) :
      async () => chemblSimilaritySearch(mol);

    searchFunc().then((t: any) => {
          compsHost.removeChild(compsHost.firstChild!);
          if (t == null) {

              compsHost.appendChild(ui.divText('No matches'));
              return;
          }
          t.col('smiles').semType = 'Molecule';
          t.col('smiles').setTag('cell.renderer', 'Molecule');


          const grid = t.plot.grid();
          compsHost.appendChild(grid.root);
          headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
              t.name = `"DrugBank Similarity Search"`;
              grok.shell.addTableView(t);
          }, 'Open compounds as table'));
          compsHost.style.overflowY = 'auto';
      }
    )
      .catch((err: any) => {
          if (compsHost.children.length > 0) {
              compsHost.removeChild(compsHost.firstChild!);
          }
          const div = ui.divText('No matches');
          ui.tooltip.bind(div, `${err}`);
          compsHost.appendChild(div);
      });
    return new DG.Widget(panel);
}

//name: Chembl Substructure Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function chemblSubstructureSearchPanel(mol: string): DG.Widget {
    return mol ? chemblSearchWidgetLocalDb(mol, true) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: Chembl Similarity Search
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function chemblSimilaritySearchPanel(mol: string): DG.Widget {
    return mol ? chemblSearchWidgetLocalDb(mol) : new DG.Widget(ui.divText('SMILES is empty'));
}

//name: ChemblBrowser
//tags: app
export async function Browser() {
  // Filter inputs
  let molecule = ui.moleculeInput('Substructure', 'C1CCCCC1');
  let subName = ui.stringInput('Subname', '');
  let molregno = ui.stringInput('Molregno', '');
  let ro5Violation = ui.choiceInput('RO5 Violations', 'All', ['All', '0', '1', '2', '3', '4']);
  let maxPhase = ui.choiceInput('Max Phase', 'All', ['All', '0', '1', '2', '3', '4']);
  let molecule_type = ui.choiceInput('Molecule type', 'All', ['Protein', 'Oligonucleotide', 'Unknown', 'Antibody', 'Oligosaccharide', 'Unclassified', 'Enzyme', 'Cell', 'All']);
  let clear = ui.button([ui.iconFA('trash-alt'), 'Clear filters'], () => clearFilters());

  let controlPanel = ui.form([
      molecule,
      subName,
      molregno,
      ro5Violation,
      maxPhase,
      molecule_type
  ]);
  $(controlPanel).addClass('ui-form-condensed');

  // Filter handlers
  molecule.onChanged(() => update());
  subName.onChanged(() => update());
  molregno.onChanged(() => findByMolregno());
  ro5Violation.onChanged(() => update());
  maxPhase.onChanged(() => update());
  molecule_type.onChanged(() => update());

  let v: DG.View;
  let r: DG.Viewer;

  async function initView() {
      let parser = document.createElement('a');
      //@ts-ignore
      parser.href = window.location;
      let pathSegments = parser.href.split("?");

      // if we came to app by link with molregno in URL
      if (pathSegments.length === 2 && pathSegments[1].includes("molregno")) {
          let parsedMolregno = parseInt(pathSegments[1].split("=")[1]);
          molregno.value = parsedMolregno.toString();
          let data = await grok.data.query(`${_package.name}:cbFindByMolregno`, {'molregno': parsedMolregno});
          data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
          v = grok.shell.newView('Chembl Browser');
          r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
          v.append(ui.divV([r]));
          v.toolbox = ui.div(controlPanel);
          v.box = true;

      }
      // if we came to app by link with set of parameters in URL
      else if (pathSegments.length > 2 && pathSegments[1].includes("substructure")) {

          let parsedSubstructure = pathSegments[1].split("=")[1];
          let parsedSubname = pathSegments[2].split("=")[1];
          let parsedRo5 = parseInt(pathSegments[3].split("=")[1]);
          let parsedMaxPhase = parseInt(pathSegments[4].split("=")[1]);
          let parsedMoleculeType = pathSegments[5].split("=")[1];

          let queryParameters = {
              'substructure': parsedSubstructure,
              'subname': parsedSubname,
              'num_ro5_violations': parsedRo5,
              'max_phase': parsedMaxPhase,
              'molecule_type': parsedMoleculeType
          };

          molecule.value = parsedSubstructure;
          subName.value = parsedSubname;
          ro5Violation.value = parsedRo5.toString();
          maxPhase.value = parsedMaxPhase.toString();
          molecule_type.value = parsedMoleculeType;

          let data = await grok.data.query(`${_package.name}:cbChemblBrowserQuery`, queryParameters);
          data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
          v = grok.shell.newView('Chembl Browser');
          r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
          v.append(ui.divV([r]));
          v.toolbox = ui.div(controlPanel);
          v.box = true;
      } else {
          let data = await grok.data.query(`${_package.name}:cbAllChemblStructures`, {});
          data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
          v = grok.shell.newView('Chembl Browser');
          r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
          v.append(ui.divV([r]));
          v.toolbox = ui.div(controlPanel);
          v.box = true;
      }
  }

  async function update() {
      let ro5 = ro5Violation.value == 'All' ? -1 : parseInt(ro5Violation.value!);
      let max_phase = maxPhase.value == 'All' ? -1 : parseInt(maxPhase.value!);
      let queryParameters = {
          'substructure': molecule.value,
          'subname': subName.value,
          'num_ro5_violations': ro5,
          'max_phase': max_phase,
          'molecule_type': molecule_type.value
      };
      let query = await grok.data.query(`${_package.name}:cbChemblBrowserQuery`, queryParameters);
      if (query.rowCount == 0) {
          grok.shell.info("No results for this filter were found")
      } else {
          query.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
          v.root.children[0].remove();
          r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
          r.setOptions({look:look});
          console.log(r.getOptions());
          v.toolbox = ui.div(controlPanel);
          v.append(ui.divV([r]));
          v.box = true;
          v.path = '';
          v.path = `/apps/Chemblbrowser?substructure=${molecule.value}?subname=${subName.value}?num_ro5_violations=${ro5}?max_phase=${max_phase}?molecule_type=${molecule_type.value}`;
      }
  }

  async function findByMolregno() {
      let query = await grok.data.query(`${_package.name}:cbFindByMolregno`, {'molregno': parseInt(molregno.value)});
      if (query.rowCount == 0) {
          grok.shell.info("No results for this filter were found")
      } else {
          query.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
          v.root.children[0].remove();
          r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
          v.toolbox = ui.div(controlPanel);
          v.append(ui.divV([r]));
          v.box = true;
          v.path = '/apps/Chemblbrowser?molregno=' + molregno.value;
      }
  }

  async function clearFilters() {
      molecule.value = '';
      subName.value = '';
      ro5Violation.value = 'All';
      maxPhase.value = 'All';
      molecule_type.value = 'All';
      molregno.value = '';
      v.path = '';
  }

  await initView();

  const look = {
    "sketchState": {
        "#type": "SketchState",
        "elementStates": [
            {
                "left": 142,
                "top": 94,
                "width": 100,
                "height": 20,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByMolregno",
                    "column": "molregno",
                    "format": null
                }
            },
            {
                "left": 142,
                "top": 118,
                "width": 100,
                "height": 20,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByMolregno",
                    "column": "compound_name",
                    "format": null
                }
            },
            {
                "left": 10,
                "top": 190,
                "width": 112,
                "height": 20,
                "type": "html",
                "viewerSettings": {
                    "markup": "<label class=\"d4-sketch-column-name\">max_phase</label>"
                }
            },
            {
                "left": 142,
                "top": 214,
                "width": 100,
                "height": 20,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByMolregno",
                    "column": "chembl_id",
                    "format": null
                }
            },
            {
                "left": 10,
                "top": 118,
                "width": 112,
                "height": 20,
                "type": "html",
                "viewerSettings": {
                    "markup": "<label class=\"d4-sketch-column-name\">compound_name</label>"
                }
            },
            {
                "left": 142,
                "top": 142,
                "width": 100,
                "height": 20,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByMolregno",
                    "column": "num_ro5_violations",
                    "format": null
                }
            },
            {
                "left": 10,
                "top": 142,
                "width": 112,
                "height": 20,
                "type": "html",
                "viewerSettings": {
                    "markup": "<label class=\"d4-sketch-column-name\">num_ro5_violations</label>"
                }
            },
            {
                "left": 10,
                "top": 94,
                "width": 112,
                "height": 20,
                "type": "html",
                "viewerSettings": {
                    "markup": "<label class=\"d4-sketch-column-name\">molregno</label>"
                }
            },
            {
                "left": 142,
                "top": 166,
                "width": 100,
                "height": 20,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByMolregno",
                    "column": "synonyms",
                    "format": null
                }
            },
            {
                "left": 10,
                "top": 166,
                "width": 112,
                "height": 20,
                "type": "html",
                "viewerSettings": {
                    "markup": "<label class=\"d4-sketch-column-name\">synonyms</label>"
                }
            },
            {
                "left": 10,
                "top": 214,
                "width": 112,
                "height": 20,
                "type": "html",
                "viewerSettings": {
                    "markup": "<label class=\"d4-sketch-column-name\">chembl_id</label>"
                }
            },
            {
                "left": 142,
                "top": 190,
                "width": 100,
                "height": 20,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByMolregno",
                    "column": "max_phase",
                    "format": null
                }
            },
            {
                "left": 7,
                "top": 20,
                "width": 235,
                "height": 70,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByMolregno",
                    "column": "canonical_smiles",
                    "format": null
                }
            }
        ]
    }
  };
}
