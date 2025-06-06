/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/aizynthfinder.css';
import {createPathsTreeTabs, isFragment, TAB_ID} from './tree-creation-utils';
import {ReactionData, Tree} from './aizynth-api';
import {DEMO_DATA, SAMPLE_TREE} from './mock-data';
import {configIcon, KEY, STORAGE_NAME} from './config-utils';
import {DEMO_MOLECULE} from './const';

export async function updateRetrosynthesisWidget(molecule: string, w: DG.Widget) {
  ui.empty(w.root);
  ui.setUpdateIndicator(w.root, true, 'Calculating paths...');
  const updateWidgetRoot = (el: HTMLElement) => {
    ui.setUpdateIndicator(w.root, false);
    w.root.append(el);
  };
  if (molecule !== DEMO_MOLECULE) {
    if (!molecule || DG.chem.Sketcher.isEmptyMolfile(molecule)) {
      updateWidgetRoot(ui.divText('Molecule is empty'));
      return;
    }
    if (DG.chem.isSmarts(molecule) || isFragment(molecule)) {
      updateWidgetRoot(ui.divText('Not applicable for smarts or moleculer fragments'));
      return;
    }

    //check molecule is valid and convert to smiles
    try {
      molecule = DG.chem.convert(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
    } catch {
      updateWidgetRoot(ui.divText('Molecule is possibly malformed'));
      return;
    }
  }

  const configName = grok.userSettings.getValue(STORAGE_NAME, KEY);

  let paths: Tree[] = [];
  if (molecule !== DEMO_MOLECULE) {
    try {

      const res = await grok.functions.call('Retrosynthesis:run_aizynthfind',
        {molecule: molecule, config: configName ?? ''});
      paths = JSON.parse(res);
    } catch (e: any) {
      updateWidgetRoot(ui.divText(e?.message ?? e));
      return;
    }
  } else
    paths = DEMO_DATA;

  try {
    if (paths.length) {
      const pathsObjects: {[key: string]:{smiles: string, type: string}[]} = {};
      let currentTreeObjId: string | null = null;
      const tabControl = createPathsTreeTabs(paths, pathsObjects, false);
      const settings = configIcon(molecule, w);
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
      updateWidgetRoot(tabControl.root);
    } else
      updateWidgetRoot(ui.divText('No paths found for the molecule'));
  } catch {
    updateWidgetRoot(ui.divText('Error processing retrosynthesis data'));
  }
}
