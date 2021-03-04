/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();
let packageName = 'Chemblbrowser'

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
    let molecule = ui.moleculeInput('', 'C1CCCCC1');
    let subName = ui.stringInput('', 'aspirin');
    let molregno = ui.intInput('', 42);
    let ro5Violation = ui.choiceInput('', '0', ['0','1','2','3','4'] );
    let molecule_types = ['Protein','Oligonucleotide','Unknown','Antibody','Oligosaccharide','Unclassified','Enzyme','Cell',''];
    let molecule_type = ui.choiceInput('', '', molecule_types );
    let showSynonyms = ui.divH([ui.boolInput('', false, () => update("ShowSynonyms", {}))]);
    $(showSynonyms).css('align-items','center');
    let showChemblID = ui.divH([ui.boolInput('', false, () => update("ShowChemblID", {}))]);
    $(showChemblID).css('align-items','center');
    let clear = ui.button('Clear filters', () => clearFilters());



    let controlPanel = ui.divV([
        ui.panel([ui.label('Find by molecular substructure'), molecule]),
        ui.panel([ui.label('Find by subname'), subName]),
        ui.panel([ui.label('Find by molregno'), molregno]),
        ui.panel([ui.label('Find by RO5 Violations'), ro5Violation]),
        ui.panel([ui.label('Filter by molecule type'), molecule_type]),
        ui.panel([ui.label('Show synonyms'), showSynonyms]),
        ui.panel([ui.label('Show ChemblID'), showChemblID]),
        ui.panel([clear])
    ]);

    // Filter handlers
    molecule.onChanged(() => update("FindBySubstructure", {'sub': molecule.stringValue}));
    subName.onChanged(() => update("FindByName", {'sub': subName.stringValue}));
    molregno.onChanged(() => update("FindByMolregno", {'molregno': molregno.value}));
    ro5Violation.onChanged(() => update("FindByRO5", {'num_ro5_violations': parseInt(ro5Violation.value)}));
    molecule_type.onChanged(() => update("FilterByMoleculeType", {'molecule_type': molecule_type.value}));


    async function initView() {
        let data = await grok.data.query(`${packageName}:allChemblStructures`, {});
        data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
        v = grok.shell.newView('Chembl Browser');
        r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
        v.append(ui.divV([r]));
        v.toolbox = ui.div(controlPanel);
        v.box = true;
    }

    async function update(queryName, queryParameters) {
        let query = await grok.data.query(`${packageName}:${queryName}`, queryParameters);
        query.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
        v.root.children[0].remove();
        r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
        v.toolbox = ui.div(controlPanel);
        v.append(ui.divV([r]));
        v.box = true;
    }

    async function clearFilters() {
        molecule.value = '';
        subName.value = '';
        ro5Violation.value = '';
        molregno.value = null;
    }

    initView();

}





