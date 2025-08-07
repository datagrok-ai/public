import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ErrorHandling, Scope } from './constants'
import '../css/moltrack.css';
import { delay } from '@datagrok-libraries/utils/src/test';

let openedView: DG.ViewBase | null = null;

export async function getMolTrackContainer() {
    return await grok.dapi.docker.dockerContainers.filter('name = "moltrack"').first();
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
    const gridDiv = ui.div('', 'moltrack-register-res-div');
    let df: DG.DataFrame | null = null;
    const registerButton = ui.bigButton('REGISTER', async () => {
        ui.setUpdateIndicator(gridDiv, true);
        const file = datasetInput.value;
        if (file && !file.id) {
            grok.shell.warning('Please upload a dataset.');
            ui.setUpdateIndicator(gridDiv, false);
            return;
        }
        df = await grok.functions.call('Moltrack:registerBulk', {
            csv_file: file,
            scope: scopeInput.value,
            mapping: mappingInput.value,
            errorHandling: errorHandlingInput.value,
        });
        ui.empty(gridDiv);
        if (df)
        {
            df.name = `Register Entities: ${scopeInput.value}`;
            gridDiv.append(df.plot.grid().root);
        }
        ui.setUpdateIndicator(gridDiv, false);
    });
    registerButton.classList.add('moltrack-run-register-button');
    const addToWorkspaceButton = ui.icons.add(() => {
        if (df) {
            const tv = grok.shell.addTablePreview(df);
            adjustIdColumnWidth(tv);
        }
    }, 'Add registration results to workspace');

    openedView.setRibbonPanels([[addToWorkspaceButton]]);
    openedView.root.append(ui.divV([
        datasetInput.root,
        scopeInput.root,
        errorHandlingInput.root,
        mappingInput.root,
        registerButton,
        gridDiv
    ], { style: { height: '100%', gap: '8px', width: '100%' } }));

    grok.shell.addPreview(openedView);
}

async function adjustIdColumnWidth(tv: DG.TableView) {
   await delay(100);
   const idCol = tv.grid.col('id');
   if (idCol)
    idCol.width = 100;
}
