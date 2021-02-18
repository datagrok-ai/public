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

  let v = grok.shell.newView('demo: ChemblBrowser');
    let proteinClassification = ui.divV([
        ui.divH([
            ui.divText('Protein Classification'),

            ui.bigButton('RUN', () => {
                grok.data.query(`${packageName}:proteinClassification`, {}).then(t => grok.shell.addTableView(t));
            }),

            ui.bigButton('RUN in cards', async () => {
             let data = await grok.data.query(`${packageName}:proteinClassification`, {});
                grok.shell.newView( "Protein Classification", [ ui.virtualView(data.rowCount, (i) => ui.card(ui.tableFromMap({
                pref_name: data.row(i).pref_name,
                definition: data.row(i).definition,
                class_level: data.row(i).class_level,
            }))).root]
                )})])]);



    let allChemblStructures = ui.divV([
        ui.divH([
            ui.divText('Chembl Structures'),
            ui.bigButton('RUN', () => {
                grok.data.query(`${packageName}:allChemblStructures`, {}).then(t => grok.shell.addTableView(t));
            }),

            ui.bigButton('RUN in cards', async () => {
                let data = await grok.data.query(`${packageName}:allChemblStructures`, {});
                grok.shell.newView( "Chembl Structures", [ ui.virtualView(data.rowCount, (i) => ui.card(ui.tableFromMap({
                        canonical_smiles: data.row(i).canonical_smiles,
                        molregno: data.row(i).molregno,
                    }))).root]
                )}),
        ])]);

    let target = ui.stringInput('Target', 'CHEMBL1827');

    let compoundActivityDetailsForTarget = ui.divV([
     ui.divH([
         ui.divText('Compound activity details for @target'),
         target,
         ui.bigButton('RUN', () => {
             grok.data.query(`${packageName}:compoundActivityDetailsForTarget`, {'target':target.stringValue}).then(t => grok.shell.addTableView(t));
         })
     ])]);


    let selectiveForTarget = ui.stringInput(' Selective For Target', 'CHEMBL301');
    let OverTarget = ui.stringInput('Over Target', 'CHEMBL4036');

    let compoundsSelectiveToOneTargetOverSecond = ui.divV([
        ui.divH([
            ui.divText('Compounds which are selective to one target over a second target'),
            selectiveForTarget, OverTarget,
            ui.bigButton('RUN', () => {
                grok.data.query(`${packageName}:compoundsSelectiveToOneTargetOverSecond`, {'selectiveFor':selectiveForTarget.stringValue, 'over': OverTarget.stringValue}).then(t => grok.shell.addTableView(t));
            })
        ])]);





    v.append(ui.divV([
        proteinClassification,  allChemblStructures, compoundActivityDetailsForTarget, compoundsSelectiveToOneTargetOverSecond]));


}



