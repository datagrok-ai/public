/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/aizynthfinder.css';
import {AiZynthFinderViewer} from './aizynthfinder-viewer';
import {createPathsTreeTabs, isFragment} from './utils';
import {ReactionData, Tree} from './aizynth-api';
import {SAMPLE_TREE} from './mock-data';

export const _package = new DG.Package();

const STORAGE_NAME = 'retrosynthesis';
const KEY = 'config';

//name: CalculateRetroSynthesisPaths
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: string molecule = "O=C1Nc2ccccc2C(C2CCCCC2)=NC1" { semType: Molecule }
//input: string configName
//output: string paths
export async function calculateRetroSynthesisPaths(molecule: string, configName?: string): Promise<string> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const startTime = performance.now();
  const currentUser = await grok.dapi.users.current();
  const userId = currentUser.id;
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/aizynthfind', {
    method: 'POST',
    body: JSON.stringify({smiles: molecule, user_id: userId, config_name: configName}),
    headers: {'Content-Type': 'application/json'},
  });
  console.log(`Request to aizynthfinder finished in ${performance.now() - startTime} ms`);
  const resJson = await response.json();
  if (!resJson['success'])
    throw new Error('Error occured during paths generation');

  return resJson['result'];
}


export async function addUserDefinedConfig(file: DG.FileInfo): Promise<void> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const fileContent = await file.readAsString();
  const currentUser = await grok.dapi.users.current();
  const userId = currentUser.id;

  console.log(JSON.stringify({config: fileContent, user_id: userId, config_name: file.name}));
  try {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/add_user_config', {
      method: 'POST',
      body: JSON.stringify({config: fileContent, user_id: userId, config_name: file.name}),
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to add configuration: ${errorText}`);
    }

    grok.shell.info('Configuration added successfully');
  } catch (error) {
    grok.shell.error(`Error adding configuration: ${error}`);
    throw error;
  }
}


//name: Chemistry | Retrosynthesis
//tags: panel, chem, widgets
//condition: true
//input: string smiles { semType: Molecule }
//output: widget result
export async function retroSynthesisPath(molecule: string): Promise<DG.Widget> {
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

  const configName = grok.userSettings.getValue(STORAGE_NAME, KEY);

  // Call retrosynthesis function
  // let result: string;
  // try {
  //   result = await grok.functions.call('Retrosynthesis:calculateRetroSynthesisPaths',
  //     {molecule: molecule, configName: configName});
  // } catch (e: any) {
  //   return new DG.Widget(ui.divText(e));
  // }

  // Parse and process reaction data
  try {
   // const reactionData: ReactionData = JSON.parse(result);
   // const paths: Tree[] = reactionData?.data?.[0]?.trees;
    const paths = SAMPLE_TREE;
    if (paths.length) {
      const w = new DG.Widget(createPathsTreeTabs(paths, false).root);
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
      let userConfigs: string[] | null = null;
      const settings = ui.icons.settings(async () => {
        const configFileInput = ui.input.file('File', {nullable: true});
        const addConfigIcon = ui.icons.add(async () => {
          if (configFileInput.value) {
            const currentConfig = configChoice.value;
            if (userConfigs!.includes(configFileInput.value.name)) {
              grok.shell.error(`Config ${configFileInput.value.name} already exists`);
              return;
            }
            await addUserDefinedConfig(configFileInput.value);
            userConfigs = userConfigs!.concat([configFileInput.value.name]);
            ui.empty(choicesDiv);
            configChoice = ui.input.choice('Current config', {
              nullable: false,
              value: currentConfig,
              items: ['default'].concat(userConfigs),
            });
            choicesDiv.append(configChoice.root);
          }
          configFileInput.value = null;
        });
        addConfigIcon.classList.add('retrosynthesis-add-config-icon');
        if (!userConfigs)
          userConfigs = await getUserConfigs();
        let configChoice = ui.input.choice('Current config', {
          nullable: false,
          value: configName ?? 'default',
          items: ['default'].concat(userConfigs),
        });
        const choicesDiv = ui.div(configChoice.root);
        const settingsDiv = ui.divV([
          choicesDiv,
          ui.divH([configFileInput.root, addConfigIcon]),
        ]);
        const dlg = ui.dialog('Settings')
          .add(settingsDiv)
          .onOK(() => {
            if (configFileInput.value)
              grok.userSettings.add(STORAGE_NAME, KEY, configChoice.value!);
            grok.shell.info(`Current config saved`);
          });
        dlg.root.classList.add('retrosynthesis-settings-dlg');
        dlg.show();
      });
      settings.classList.add('retrosynthesis-settings-icon');
      w.root.append(settings);
      return w;
    } else
      return new DG.Widget(ui.divText('No paths found for the molecule'));
  } catch {
    return new DG.Widget(ui.divText('Error processing retrosynthesis data'));
  }
}


//name: Retrosynthesis Viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/chem-similarity-search-viewer.svg
export function retrosynthesisViewer(): AiZynthFinderViewer {
  return new AiZynthFinderViewer();
}

//name: retrosynthesisTopMenu
export function retrosynthesisTopMenu(): void {
  (grok.shell.v as DG.TableView).addViewer('Retrosynthesis Viewer');
}

//name: GetUserConfigs
//output: list<string> configs
export async function getUserConfigs(): Promise<string[]> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const currentUser = await grok.dapi.users.current();
  const userId = currentUser.id;

  try {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, `/get_user_configs`, {
      method: 'POST',
      body: JSON.stringify({user_id: userId}),
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to get configuration files: ${errorText}`);
    }

    const configs = await response.json();
    return configs.configs ?? [];
  } catch (error) {
    grok.shell.error(`Error getting configuration files: ${error}`);
    throw error;
  }
}

