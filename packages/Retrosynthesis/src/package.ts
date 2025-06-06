/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/aizynthfinder.css';
import {getUserConfigsFromDocker, syncConfig} from './config-utils';
import {CONFIGS_PATH, DEMO_MOLECULE} from './const';
import {updateRetrosynthesisWidget} from './utils';

export const _package = new DG.Package();


//name: Chemistry | Retrosynthesis
//tags: panel, chem, widgets
//meta.allowAddAsColumn: false
//condition: true
//input: string smiles { semType: Molecule }
//output: widget result
export function retroSynthesisPath(molecule: string): DG.Widget {
  const w = new DG.Widget(ui.div('', 'retrosynthesis-widget-div'));
  updateRetrosynthesisWidget(molecule, w);
  return w;
}


//name: retrosynthesisTopMenu
export function retrosynthesisTopMenu(): void {
  (grok.shell.v as DG.TableView).addViewer('Retrosynthesis Viewer');
}

//name: GetUserConfigs
//output: list<string> configs
export async function getUserConfigs(): Promise<string[]> {
  return getUserConfigsFromDocker();
}

//name: Retrosynthesis Demo
//description: Generate retrosynthesis paths
//meta.demoPath: Cheminformatics | Retrosynthesis
export async function retrosynthesisDemo(): Promise<void> {
  await grok.functions.call('Chem:initChemAutostart');
  const view = DG.View.create();
  view.name = 'Retrosynthesis Demo';

  const sketcher: DG.chem.Sketcher = new DG.chem.Sketcher();
  sketcher.setSmiles('COc1ccc2c(c1)c(CC(=O)N3CCCC3C(=O)Oc4ccc(C)cc4OC)c(C)n2C(=O)c5ccc(Cl)cc5')
  const retrosynthesisDiv = ui.div('', 'retrosynthesis-demo');

  const container = ui.divH([
    sketcher.root,
    retrosynthesisDiv,
  ], {style: {height: '100%'}});

  view.append(container);

  let demoInited = false;
  sketcher.onChanged.subscribe(async () => {
    const smiles = sketcher.getSmiles();
    if (smiles) {
      try {
        ui.empty(retrosynthesisDiv);
        ui.setUpdateIndicator(retrosynthesisDiv, true, 'Calculating retrosyntehsis paths...');
        const widget = await retroSynthesisPath(!demoInited ? DEMO_MOLECULE : smiles);
        demoInited = true;
        retrosynthesisDiv.append(widget.root);
        ui.setUpdateIndicator(retrosynthesisDiv, false);
      } catch (e) {
        grok.shell.error('Invalid or empty molecule');
      }
    }
  });
  grok.shell.addPreview(view);
}


