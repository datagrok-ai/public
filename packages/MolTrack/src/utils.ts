import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let openedView: DG.ViewBase | null = null;

export function createRegistrationNode(treeNode: DG.TreeViewGroup) {
  openedView?.close();
  openedView = DG.View.create();
  openedView.name = 'Register Entities';

  const datasetInput = ui.input.file('Dataset');
  const scopeChoices = ['compounds', 'batches', 'assays', 'assay_runs', 'assay_results'];
  const scopeInput = ui.input.choice('Scope', {value: scopeChoices[0], items: scopeChoices});
  const noteInput = ui.input.string('Note', {value: ''});
  const gridDiv = ui.div('', 'registration-results-grid');
  const registerButton = ui.bigButton('REGISTER', async () => {
  });

  registerButton.classList.add('registration-run-button');

  const addToWorkspaceButton = ui.icons.add(() => {
  }, 'Add registration results to workspace');

  openedView.setRibbonPanels([[addToWorkspaceButton]]);
  openedView.root.append(ui.divV([
    datasetInput.root,
    scopeInput.root,
    noteInput.root,
    registerButton,
    gridDiv
  ], { style: { height: '100%', gap: '8px' } }));

  grok.shell.addPreview(openedView);
}