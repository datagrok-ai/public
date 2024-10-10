import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export class ServiceLogsApp extends DG.ViewBase {
    static DEFAULT_LIMIT: number = 1000;

    static createView(parentCall: DG.FuncCall): DG.ViewBase {
        return DG.View.fromViewAsync(async () => {
            const view = DG.View.create();
            view.name = 'Service Logs';
            view.path = '';
            view.parentCall = parentCall;

            let currentTable = DG.DataFrame.create();
            const tableGridRoot: HTMLDivElement = document.createElement('div');
            tableGridRoot.append(DG.Grid.create(currentTable).root);
            const limitInput: DG.InputBase<number | null> = ui.input.int('Limit records', {min: 1, max: 10000, tooltipText: 'Limits number of returned records', value: ServiceLogsApp.DEFAULT_LIMIT, showPlusMinus: false});
            const serviceChoiceInput: DG.ChoiceInput<string> = ui.input.choice<string>('Service name', {items: await grok.dapi.docker.getAvailableServices(), tooltipText: 'Choose service', nullable: false}) as DG.ChoiceInput<string>;
            let runButton: HTMLButtonElement;
            runButton = ui.bigButton('GET LOGS', async () => {
                const progress = DG.TaskBarProgressIndicator.create(`Loading ${serviceChoiceInput.value} logs...`);
                runButton.disabled = true;
                const loader = ui.loader();
                loader.style.top = 'auto';
                try {
                    while(tableGridRoot.firstChild)
                        tableGridRoot.removeChild(tableGridRoot.firstChild);
                    tableGridRoot.append(loader);
                    currentTable = await grok.dapi.docker.getServiceLogs(serviceChoiceInput.value,limitInput.value ??  ServiceLogsApp.DEFAULT_LIMIT);
                    tableGridRoot.removeChild(loader);
                    tableGridRoot.append(DG.Grid.create(currentTable).root);
                } catch (e: any) {
                    tableGridRoot.removeChild(loader);
                    tableGridRoot.append(DG.Grid.create(DG.DataFrame.create()).root);
                    console.error(e);
                    new DG.Balloon().error(e.toString());
                } finally {
                    progress.close();
                    runButton.disabled = false;
                }
            }, 'Executes request to get service logs');

            const form = ui.forms.condensed([serviceChoiceInput, limitInput]);
            ui.forms.addButtons(form, [runButton]);
            form.style.flexBasis = '15%';
            tableGridRoot.style.flexBasis = '85%';
            const content = ui.divV([
                form,
                tableGridRoot
            ]);
            content.style.flexBasis = '100%';
            view.root.append(content);
            view.setRibbonPanels([[ui.iconFA('plus', () => grok.shell.addTable(currentTable),
                'Add table to workspace')]]);
            return view;
        });
    }
}
