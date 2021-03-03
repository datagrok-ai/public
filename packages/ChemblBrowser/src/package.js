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
//name: Browser
//tags: app
export async function Browser() {

    async function FindBySubstructure() {
        let sub = molecule.stringValue;
        let query = await grok.data.query(`${packageName}:FindBySubstructure`, {'sub': sub});
        grok.shell.addTableView(query);

    }

    async function FindByName() {
        let query = await grok.data.query(`${packageName}:FindByName`, {'sub': subName.stringValue});
        grok.shell.addTableView(query);

    }


    async function FindByMolregno() {
        let query = await grok.data.query(`${packageName}:FindByMolregno`, {'sub': molregno.value});
        grok.shell.addTableView(query);
    }



    let molecule = ui.moleculeInput('', '');
    let subName = ui.stringInput('Find by name', '');
    let molregno = ui.intInput('Find by molregno');

    molecule.onChanged(() => FindBySubstructure());
    subName.onChanged(() => FindByName());
    molregno.onChanged(() => FindByMolregno());



    let inputs = [molecule, subName,molregno];
    let controlPanel = ui.divV([ui.h2("Find By Substructure"), ui.inputs(inputs)]);

    let data = await grok.data.query(`${packageName}:allChemblStructures`, {});
    let v = grok.shell.newView('Chembl Browser');
    let r = DG.Viewer.fromType(DG.VIEWER.TILE_VIEWER, data);

    v.append(ui.divV([controlPanel, r]));


    v.box = true;





}





