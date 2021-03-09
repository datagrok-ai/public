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
        let pathSegments = parser.href.split("?");


        // if we came to app by link with molregno in URL
        if (pathSegments.length == 2 && pathSegments[1].includes("molregno")) {
            let parsedMolregno = parseInt(pathSegments[1].split("=")[1]);
            molregno.value = parsedMolregno;
            let data = await grok.data.query(`${packageName}:FindByMolregno`, {'molregno': parsedMolregno});
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

            console.log(queryParameters);

            molecule.value = parsedSubstructure;
            subName.value = parsedSubname;
            ro5Violation.value = parsedRo5;
            maxPhase = parsedMaxPhase;
            molecule_type.value = parsedMoleculeType;

            let data = await grok.data.query(`${packageName}:ChemblBrowserQuery`, queryParameters);
            data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
            v = grok.shell.newView('Chembl Browser');
            r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
            v.append(ui.divV([r]));
            v.toolbox = ui.div(controlPanel);
            v.box = true;
        } else {
            let data = await grok.data.query(`${packageName}:allChemblStructures`, {});
            data.col('canonical_smiles').semType = DG.SEMTYPE.MOLECULE;
            v = grok.shell.newView('Chembl Browser');
            r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);
            v.append(ui.divV([r]));
            v.toolbox = ui.div(controlPanel);
            v.box = true;
        }
    }

    async function update() {
        let ro5 = ro5Violation.value == 'All' ? -1 : parseInt(ro5Violation.value);
        let max_phase = ro5Violation.value == 'All' ? -1 : parseInt(maxPhase.value);
        let queryParameters = {
            'substructure': molecule.value,
            'subname': subName.value,
            'num_ro5_violations': ro5,
            'max_phase': max_phase,
            'molecule_type': molecule_type.value
        };
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
            v.path = `/apps/Chemblbrowser?substructure=${molecule.value}?subname=${subName.value}?num_ro5_violations=${ro5}?max_phase=${max_phase}?molecule_type=${molecule_type.value}`;
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
            v.path = '/apps/Chemblbrowser?molregno=' + molregno.value;
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





