/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {data} from "datagrok-api/grok";
import {inputs, stringInput} from "datagrok-api/ui";

export let _package = new DG.Package();
let packageName = 'Chemblapi'

//name: test
//input: string s
export function test(s) {
    grok.shell.info(_package.webRoot);
}

let r = null;
let v = null;


//name: Browser
//tags: app
export async function Browser() {

    // Filter inputs
    let molecule = ui.moleculeInput('Substructure', 'C1CCCCC1');
    let subName = ui.stringInput('Subname', '');
    let chemblid = ui.stringInput('Chemblid', '');
    let ro5Violation = ui.choiceInput('RO5 Violations', 'All', ['All', '0', '1', '2', '3', '4']);
    let maxPhase = ui.choiceInput('Max Phase', 'All', ['All', '0', '1', '2', '3', '4']);
    let molecule_type = ui.choiceInput('Molecule type', 'All', ['Protein', 'Oligonucleotide', 'Unknown', 'Antibody', 'Oligosaccharide', 'Unclassified', 'Enzyme', 'Cell', 'All']);
    let clear = ui.button([ui.iconFA('trash-alt'), 'Clear filters'], () => clearFilters());

    let controlPanel = ui.form([
        molecule,
        subName,
        chemblid,
        ro5Violation,
        maxPhase,
        molecule_type,
        clear
    ]);
    $(controlPanel).addClass('ui-form-condensed'); //TODO: use ui.smallForm instead

    molecule.onChanged(() => update());
    subName.onChanged(() => update());
    chemblid.onChanged(() => findByChemblid());
    ro5Violation.onChanged(() => update());
    maxPhase.onChanged(() => update());
    molecule_type.onChanged(() => update());


    async function initView() {
        let queryParameters = {};
        let parser = document.createElement('a');
        parser.href = window.location;
        let pathSegments = parser.href.split("?");


        // if we came to app by link with chemblid in URL
        if (pathSegments.length === 2 && pathSegments[1].includes("chemblid")) {
            let parsedChemblid = parseInt(pathSegments[1].split("=")[1]);
            chemblid.value = parsedChemblid;
            let data = await grok.data.query(`${packageName}:MoleculeJson`, {'molecule_chembl_id__exact': `CHEMBL${parsedChemblid}`});
            console.error(data)
            data.col('molecules/molecule_structures/canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
            v = grok.shell.newView('Chembl Browser');
            r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);

            v.append(r);
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

            queryParameters = {
                'smile': parsedSubstructure,
                'subname': parsedSubname,
                'num_ro5_violations': parsedRo5,
                'max_phase': parsedMaxPhase,
                'molecule_type': parsedMoleculeType
            };

            molecule.value = parsedSubstructure;
            subName.value = parsedSubname;
            ro5Violation.value = parsedRo5;
            maxPhase = parsedMaxPhase;
            molecule_type.value = parsedMoleculeType;



        }
        else {
            queryParameters = {}
        }
        await displayQuerry(queryParameters);
    }

    async function displayQuerry(queryParameters) {
        let data = await grok.data.query(`${packageName}:MoleculeJson`, queryParameters);
        data.col('molecules/molecule_structures/canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
        v = grok.shell.newView('Chembl Browser');
        r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
        v.append(ui.divV([r.root]));
        v.toolbox = ui.div(controlPanel);
        v.box = true;
    }

    async function update() {
        let ro5 = ro5Violation.value == 'All' ? -1 : parseInt(ro5Violation.value);
        let max_phase = maxPhase.value == 'All' ? -1 : parseInt(maxPhase.value);
        let queryParameters = {
            'subname': subName.value,
            'num_ro5_violations': ro5,
            'max_phase': max_phase,
            'molecule_type': molecule_type.value
        };
        let query = null;
        try {
            query = (molecule.value === '') ?
                await grok.data.query(`${packageName}:MoleculeJson`) :
                await grok.data.query(`${packageName}:SubstructureSmile`,
                    {'smile': molecule.value});
        } catch (e) {

            console.error(e);


        }


        if (query.rowCount == 0) {
            grok.shell.info("No results for this filter were found")
        } else {
            query.col('molecules/molecule_structures/canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
            v.root.children[0].remove();
            r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
            r.setOptions({look: look});
            r.box = true
            console.log(r.getOptions());
            v.toolbox = ui.div(controlPanel);
            v.append(r);
            v.box = true;
            v.path = '';
            v.path = `/apps/Chemblapibrowser?substructure=${molecule.value}?subname=${subName.value}?num_ro5_violations=${ro5}?max_phase=${max_phase}?molecule_type=${molecule_type.value}`;
        }
    }

    async function findByChemblid() {


        let query = await grok.data.query(`${packageName}:MoleculeJson`, {'molecule_chembl_id__exact': `CHEMBL${chemblid.value}`});
        if (query.rowCount == 0) {
            grok.shell.info("No results for this filter were found")
        } else {
            query.col('molecules/molecule_structures/canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
            v.root.children[0].remove();
            r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
            r.box = true
            v.toolbox = ui.div(controlPanel);
            v.append(r);
            v.box = true;
            v.path = '/apps/Chemblapibrowser?chemblid=' + chemblid.value;
        }
    }

    async function clearFilters() {
        molecule.value = '';
        subName.value = '';
        ro5Violation.value = 'All';
        maxPhase.value = 'All';
        molecule_type.value = 'All';
        chemblid.value = null;
        v.path = '';
    }

    await initView();

}


const look = {
    "sketchState": {
        "#type": "SketchState",
        "elementStates": [

            {
                "left": 142,
                "top": 118,
                "width": 100,
                "height": 20,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByChemblid",
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
                    "table": "FindByChemblid",
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
                    "table": "FindByChemblid",
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
                    "markup": "<label class=\"d4-sketch-column-name\">chemblid</label>"
                }
            },
            {
                "left": 142,
                "top": 166,
                "width": 100,
                "height": 20,
                "type": "field",
                "viewerSettings": {
                    "table": "FindByChemblid",
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
                    "table": "FindByChemblid",
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
                    "table": "FindByChemblid",
                    "column": "molecules/molecule_structures/canonical_smiles",
                    "format": null
                }
            }
        ]
    }
};


//name: ChemblSubstructureSearch
//input: string mol {semType: Molecule}
//output: dataframe dbResult
export async function
chemblSubstructureSearch(mol) {
    try {
        let df = await grok.data.query('Chemblapi:SubstructureSmile', {'smile': mol})

        if (df == null) {
            return DG.DataFrame.create();
        }
        return df;
    } catch (e) {
        console.error("In SubstructureSearch: " + e.toString());
        throw e;
    }
}

//name: ChemblSimilaritySearch
//input: string molecule {semType: Molecule}
//output: dataframe drbResult
export async function chemblSimilaritySearch(molecule) {
    try {
        let df = await grok.data.query('Chemblapi:SimilaritySmileScore', {'smile': molecule, 'score': 40})
        console.warn(df)
        if (df == null) {
            return DG.DataFrame.create();
        }
        return df;

    } catch (e) {
        console.error(e);
        return null;
    }

}

//name: Chembl Search Widget
//tags: widgets
//input: string mol {semType: Molecule}
//input: string searchType
//output: widget result
export function ChemblSearchWidget(mol, searchType) {
    let headerHost = ui.divH([]);
    let compsHost = ui.divH([ui.loader(), headerHost]);
    let panel = ui.divV([compsHost]);
    let search = {
        'similarity': async () => chemblSimilaritySearch(mol),
        'substructure': async () => chemblSubstructureSearch(mol)
    }
    if (!searchType in search) {

        throw "DrugBankSearch: No such search type" + searchType;
    }

    search[searchType]().then(t => {
            compsHost.removeChild(compsHost.firstChild);
            if (t == null) {

                compsHost.appendChild(ui.divText('No matches'));
                return;
            }
            compsHost.appendChild(DG.Viewer.grid(t).root)


            headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
                t.name = `"DrugBank Similarity Search"`;
                grok.shell.addTableView(t);
            }, 'Open compounds as table'));
            compsHost.style.overflowY = 'auto';
        }
    )
        .catch(err => {
            if (compsHost.children.length > 0) {
                compsHost.removeChild(compsHost.firstChild);
            }
            let div = ui.divText('No matches');
            ui.tooltip.bind(div, `${err}`);
            compsHost.appendChild(div);
        });
    return new DG.Widget(panel);
}

//name Chembl Get by Id
//input string id
//output dataframe
export async function getById(id){
    if (!id.toLowerCase().startsWith("chembl")){
        id = "CHEMBL"+ id
    }
    try {
        return await grok.data.query(`${packageName}:MoleculeJson`, {'molecule_chembl_id__exact': id});
    }
    catch (e){
        console.error(e);
        return null;
    }

}

//name: Chembl Substructure Search Widget
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function ChemblSubstructureSearchPanel(mol) {
    return ChemblSearchWidget(mol, 'substructure');
}

//name: Chembl Similarity Search Widget
//tags: panel, widgets
//input: string mol {semType: Molecule}
//output: widget result
//condition: true
export function ChemblSimilaritySearchPanel(mol) {
    return ChemblSearchWidget(mol, 'similarity');
}
