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
    let molecule = ui.moleculeInput('Substructure', 'C1CCCCC1');
    let subName = ui.stringInput('Subname', 'aspirin');
    let molregno = ui.intInput('Molregno', 42);
    let ro5Violation = ui.choiceInput('RO5 Violations', '0', ['0','1','2','3','4'] );
    let maxPhase = ui.choiceInput('Max Phase', '0', ['0','1','2','3','4'] );
    let molecule_types = ['Protein','Oligonucleotide','Unknown','Antibody','Oligosaccharide','Unclassified','Enzyme','Cell'];
    let molecule_type = ui.choiceInput('Molecule type', '', molecule_types );
    let showSynonyms = ui.boolInput('Show synonyms', false, () => update("ShowSynonyms", {}));
    let showChemblID = ui.boolInput('Show ChemblID', false, () => update("ShowChemblID", {}));
    let clear = ui.button([ui.iconFA('trash-alt'), 'Clear filters'], () => clearFilters());

    let controlPanel = ui.form([
        molecule,
        subName,
        molregno,
        ro5Violation,
        maxPhase,
        molecule_type,
        showSynonyms,
        showChemblID,
        clear
    ]);
    $(controlPanel).addClass('ui-form-condensed'); //TODO: use ui.smallForm instead

    // Filter handlers
    molecule.onChanged(() => update("FindBySubstructure", {'sub': molecule.stringValue}));
    subName.onChanged(() => update("FindByName", {'sub': subName.stringValue}));
    molregno.onChanged(() => update("FindByMolregno", {'molregno': molregno.value}));
    ro5Violation.onChanged(() => update("FindByRO5", {'num_ro5_violations': parseInt(ro5Violation.value)}));
    maxPhase.onChanged(() => update("FilterByMaxPhase", {'max_phase': parseInt(maxPhase.value)}));
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
        maxPhase.value = '';
        molecule_type.value = '';
        molregno.value = null;
    }

    await initView();

}





