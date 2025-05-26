/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/aizynthfinder.css';
import {createPathsTreeTabs, isFragment, TAB_ID} from './utils';
import {ReactionData, Tree} from './aizynth-api';
import {DEMO_DATA, SAMPLE_TREE} from './mock-data';
import {configIcon, getUserConfigsFromDocker, KEY, STORAGE_NAME, syncConfig} from './config-utils';
import { CONFIGS_PATH } from './const';

export const _package = new DG.Package();
const DEMO_MOLECULE = 'demo_molecule';

//name: CalculateRetroSynthesisPaths
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * 1
//input: string molecule = "O=C1Nc2ccccc2C(C2CCCCC2)=NC1" { semType: Molecule }
//input: string configDir
//output: string paths
export async function calculateRetroSynthesisPaths(molecule: string, configDir?: string): Promise<string> {
  let progressBar: DG.TaskBarProgressIndicator | null = null;
  try {
    const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
    if (configDir) {
      const syncStartTime = performance.now();
      progressBar = DG.TaskBarProgressIndicator.create(`Config synchronization in progress...`);
      await syncConfig(`${CONFIGS_PATH}/${configDir}`);
      console.log(`Synchronized ${CONFIGS_PATH}/${configDir} in ${performance.now() - syncStartTime} ms`);
      progressBar.close();
    }
    const startTime = performance.now();
    progressBar = DG.TaskBarProgressIndicator.create(`Calculating retrosynthesis paths...`);
    const currentUser = await grok.dapi.users.current();
    const userId = currentUser.id;
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/aizynthfind', {
      method: 'POST',
      body: JSON.stringify({smiles: molecule, user_id: userId, config_dir: configDir !== '' ? `${CONFIGS_PATH}/${configDir}` : ''}),
      headers: {'Content-Type': 'application/json'},
    });
    console.log(`Request to aizynthfinder finished in ${performance.now() - startTime} ms`);
    const resJson = await response.json();
    if (!resJson['success'])
      throw new Error(`Error occured during paths generation: ${resJson['error']}`);
    return resJson['result'];
  } finally {
    progressBar?.close();
  }
}


//name: Chemistry | Retrosynthesis
//tags: panel, chem, widgets
//condition: true
//input: string smiles { semType: Molecule }
//output: widget result
export async function retroSynthesisPath(molecule: string): Promise<DG.Widget> {
  if (molecule !== DEMO_MOLECULE) {
    if (!molecule || DG.chem.Sketcher.isEmptyMolfile(molecule))
      return new DG.Widget(ui.divText('Molecule is empty'));
    if (DG.chem.isSmarts(molecule) || isFragment(molecule))
      return new DG.Widget(ui.divText('Not applicable for smarts or moleculer fragments'));

    //check molecule is valid and convert to smiles
    try {
      molecule = DG.chem.convert(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
    } catch {
      return new DG.Widget(ui.divText('Molecule is possibly malformed'));
    }
  }

  const configName = grok.userSettings.getValue(STORAGE_NAME, KEY);

  let paths: Tree[] = [];
  if (molecule !== DEMO_MOLECULE) {
    try {
      const result = await grok.functions.call('Retrosynthesis:calculateRetroSynthesisPaths',
        {molecule: molecule, configDir: configName ?? ''});
      const reactionData: ReactionData = JSON.parse(result);
      paths = reactionData?.data?.[0]?.trees;
    } catch (e: any) {
      return new DG.Widget(ui.divText(e?.message ?? e));
    }
  } else
    paths = DEMO_DATA;

  try {
    if (paths.length) {
      const pathsObjects: {[key: string]:{smiles: string, type: string}[]} = {};
      let currentTreeObjId: string | null = null;
      const tabControl = createPathsTreeTabs(paths, pathsObjects, false);
      const settings = configIcon();
      const addToWorkSpace = ui.icons.add(() => {
        if (currentTreeObjId) {
          const df = DG.DataFrame.fromObjects(pathsObjects[currentTreeObjId]);
          if (df) {
            df.name = tabControl.currentPane.name;
            grok.shell.addTableView(df);
          }
        }
      }, 'Add paths to workspace as dataframe');
      addToWorkSpace.classList.add('retrosynthesis-add-to-workspace-icon');
      const iconsDiv = ui.divH([addToWorkSpace, settings], 'retrosynthesis-icons-div');
      const updateCurrentPane = () => {
        tabControl.currentPane.content.append(iconsDiv);
        currentTreeObjId = tabControl.currentPane.content.getAttribute(TAB_ID);
      };
      tabControl.onTabChanged.subscribe(() => {
        updateCurrentPane();
      });
      updateCurrentPane();
      const w = new DG.Widget(tabControl.root);
      //workaround to make tree visible in undocked panel
      ui.tools.waitForElementInDom(w.root).then(() => {
        if (w.root.closest('.dialog-floating')) {
          const accPanel = w.root.closest('.panel-content') as HTMLElement;
          if (accPanel) {
            const sub = ui.onSizeChanged(accPanel).subscribe((_) => {
              const waitParentEl = w.root.parentElement;
              if (waitParentEl?.classList.contains('grok-wait'))
                waitParentEl.style.height = `100%`;
              sub.unsubscribe();
            });
          }
        }
      });
      return w;
    } else
      return new DG.Widget(ui.divText('No paths found for the molecule'));
  } catch {
    return new DG.Widget(ui.divText('Error processing retrosynthesis data'));
  }
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


