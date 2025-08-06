import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ErrorHandling, Scope } from './constants'
import '../css/moltrack.css';

let openedView: DG.ViewBase | null = null;

export async function getMolTrackContainer() {
    return await grok.dapi.docker.dockerContainers.filter('moltrack').first();
}

//input: dynamic treeNode
export function createRegistrationNode(treeNode: DG.TreeViewGroup) {
    openedView?.close();
    openedView = DG.View.create();
    openedView.name = 'Register Entities';

    const datasetInput = ui.input.file('Dataset');

    const scopeChoices = Object.values(Scope);
    const scopeInput = ui.input.choice('Scope', { value: scopeChoices[0], items: scopeChoices });

    const errorHandlingChoices = Object.values(ErrorHandling);
    const errorHandlingInput = ui.input.choice('Error handling', { value: errorHandlingChoices[0], items: errorHandlingChoices });

    const mappingInput = ui.input.textArea('Mapping', { value: '' });
    const gridDiv = ui.div('','moltrack-register-res-div');
    const registerButton = ui.bigButton('REGISTER', async () => {
        //Register
        ui.setUpdateIndicator(gridDiv, true);

        const file = datasetInput.value;
        if (!file) {
            grok.shell.warning('Please upload a dataset.');
            return;
        }
        const result = await grok.functions.call('Moltrack:registerBulk', {
            csv_file: file,
            scope: scopeInput.value,
            mapping: mappingInput.value,
            errorHandling: errorHandlingInput.value,
        });
        ui.empty(gridDiv);
        if (result instanceof DG.DataFrame)
            gridDiv.appendChild(result.plot.grid().root);

        ui.setUpdateIndicator(gridDiv, false);
    });

    registerButton.classList.add('moltrack-run-register-button');

    const addToWorkspaceButton = ui.icons.add(() => {
    }, 'Add registration results to workspace');

    openedView.setRibbonPanels([[addToWorkspaceButton]]);
    openedView.root.append(ui.divV([
        datasetInput.root,
        scopeInput.root,
        mappingInput.root,
        errorHandlingInput.root,
        registerButton,
        gridDiv
    ], { style: { height: '100%', gap: '8px', width: '100%' } }));

    grok.shell.addPreview(openedView);
}
