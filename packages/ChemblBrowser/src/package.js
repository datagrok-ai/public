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
    let subName = ui.stringInput('Subname', '');
    let molregno = ui.stringInput('Molregno', '');
    let ro5Violation = ui.choiceInput('RO5 Violations', 'All', ['All','0','1','2','3','4'] );
    let maxPhase = ui.choiceInput('Max Phase', 'All', ['All', '0','1','2','3','4'] );
    let molecule_type = ui.choiceInput('Molecule type', 'All', ['Protein','Oligonucleotide','Unknown','Antibody','Oligosaccharide','Unclassified','Enzyme','Cell', 'All']);
    let clear = ui.button([ui.iconFA('trash-alt'), 'Clear filters'], () => clearFilters());

    let controlPanel = ui.form([
        molecule,
        subName,
        molregno,
        ro5Violation,
        maxPhase,
        molecule_type,
        clear
    ]);
    $(controlPanel).addClass('ui-form-condensed'); //TODO: use ui.smallForm instead

    // Filter handlers
    molecule.onChanged(() => update());
    subName.onChanged(() => update());
    molregno.onChanged(() => findByMolregno());
    ro5Violation.onChanged(() => update());
    maxPhase.onChanged(() => update());
    molecule_type.onChanged(() => update());


    async function initView() {

        let parser = document.createElement('a');
        parser.href = window.location;
        let pathSegments = parser.pathname.split('?');
        if (pathSegments.length>2) {console.log(pathSegments)};

        let data = await grok.data.query(`${packageName}:allChemblStructures`, {});
        data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
        v = grok.shell.newView('Chembl Browser');
        r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
        v.append(ui.divV([r]));
        v.toolbox = ui.div(controlPanel);
        v.box = true;
    }

    async function update() {
        let ro5 = ro5Violation.value == 'All' ? -1 : parseInt(ro5Violation.value);
        let max_phase = ro5Violation.value == 'All' ? -1 : parseInt(maxPhase.value);
        let queryParameters = {'substructure': molecule.value, 'subname': subName.value,'num_ro5_violations': ro5,'max_phase': max_phase, 'molecule_type': molecule_type.value};
        let query = await grok.data.query(`${packageName}:ChemblBrowserQuery`, queryParameters);
        if (query.rowCount == 0) {
            grok.shell.info("No results for this filter were found")
        } else {
            query.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
            v.root.children[0].remove();
            r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
            v.toolbox = ui.div(controlPanel);
            v.append(ui.divV([r]));
            v.box = true;

            v.path = '';
            v.path = `/apps/Chemblbrowser/substructure=${molecule.value}?subname=${subName.value}?num_ro5_violations=${ro5}?max_phase=${max_phase}?molecule_type=${molecule_type.value}`;
        }
    }

    async function findByMolregno() {
        let query = await grok.data.query(`${packageName}:FindByMolregno`, {'molregno': parseInt(molregno.value)});
        if (query.rowCount == 0) {
            grok.shell.info("No results for this filter were found")
        } else {
            query.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
            v.root.children[0].remove();
            r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, query);
            v.toolbox = ui.div(controlPanel);
            v.append(ui.divV([r]));
            v.box = true;
            v.path = '/apps/Chemblbrowser/?molregno=' + molregno.value;
        }
    }

    async function clearFilters() {
        molecule.value = '';
        subName.value = '';
        ro5Violation.value = 'All';
        maxPhase.value = 'All';
        molecule_type.value = 'All';
        molregno.value = null;
        v.path = '';
    }

    await initView();

}





