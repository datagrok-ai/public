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
  let query = ui.choiceInput('query', 'allChemblStructures', ['allChemblStructures', 'allChemblIdsWithInchiKeys','proteinClassification']);
  let container = ui.div();
  v.append(container);

  let inputs = [query]
  container.appendChild(ui.inputs(inputs));

  container.appendChild(ui.bigButton('RUN', () => {
    grok.data.query(`${packageName}:${query.caption}`, {}).then(t => grok.shell.addTableView(t));
  }));


}
