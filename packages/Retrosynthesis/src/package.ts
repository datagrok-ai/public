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


//name: CalculateRetroSynthesisPaths
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: string molecule = "O=C1Nc2ccccc2C(C2CCCCC2)=NC1" { semType: Molecule }
//output: string paths
export async function calculateRetroSynthesisPaths(molecule: string): Promise<string> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const startTime = performance.now();
  const currentUser = await grok.dapi.users.current();
  const userId = currentUser.id;
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/aizynthfind', {
    method: 'POST',
    body: JSON.stringify({smiles: molecule, id: userId}),
    headers: {'Content-Type': 'application/json'},
  });
  console.log(`Request to aizynthfinder finished in ${performance.now() - startTime} ms`);
  const resJson = await response.json();
  if (!resJson['success'])
    throw new Error('Error occured during paths generation');

  return resJson['result'];
}


export async function setUserDefinedConfig(file: DG.FileInfo): Promise<void> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const fileContent = await file.readAsString();
  const currentUser = await grok.dapi.users.current();
  const userId = currentUser.id;

  try {
    const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/set_config', {
      method: 'POST',
      body: JSON.stringify({config: fileContent, id: userId}),
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to set configuration: ${errorText}`);
    }

    grok.shell.info('Configuration updated successfully');
  } catch (error) {
    grok.shell.error(`Error setting configuration: ${error}`);
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

  // Call retrosynthesis function
  let result: string;
  try {
    result = await grok.functions.call('Retrosynthesis:calculateRetroSynthesisPaths',
      {molecule: molecule});
  } catch (e: any) {
    return new DG.Widget(ui.divText(e));
  }

  // Parse and process reaction data
  try {
    const reactionData: ReactionData = JSON.parse(result);
    const paths: Tree[] = reactionData?.data?.[0]?.trees;
    //const paths = SAMPLE_TREE;
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
      const settings = ui.icons.settings(() => {
        const configFileInput = ui.input.file('Config', {nullable: true});
        w.root.append(configFileInput.root);
        const dlg = ui.dialog('Settings')
          .add(configFileInput.root)
          .onOK(() => {
            if (configFileInput.value)
              setUserDefinedConfig(configFileInput.value);
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

