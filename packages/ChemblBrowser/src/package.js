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


        async function molregnoView(molregno) {
            let data = grok.data.query(`${packageName}:CompoundProperties`, {'molregno': parseInt(molregno)})
            let view = grok.shell.addTableView(data);
        }

        async function initView()  {

        let data = await grok.data.query(`${packageName}:allChemblStructures`, {});
        let container = ui.divH([],'card-container');
        for (let i=0; i<data.rowCount; i++){
            let h2 = ui.h2(`molregno: ${data.col('molregno').get(i)}`);
            let h3 = ui.h3(`standard_inchi: ${data.col('standard_inchi').get(i)}`);
            h2.addEventListener('click', e => molregnoView(h2.innerText))
            h2.addEventListener('dblclick', e => grok.shell.info(h3.innerText))
            container.append(ui.card(ui.div([
                h2, grok.chem.svgMol(data.col('canonical_smiles').get(i),150, 100)]
            )))
        }
        grok.shell.newView( "demo: ChemblBrowser", [container])
        $('.card-container').css('flex-wrap','wrap');
    }

    initView();
}





