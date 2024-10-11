import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export class ServiceLogsApp extends DG.ViewBase {
    static readonly DEFAULT_LIMIT: number = 1000;
    static readonly TIMESTAMP_COLUMN_NAME = 'timestamp';
    static readonly MESSAGE_COLUMN_NAME = 'message';

    private readonly tableGridRoot: HTMLDivElement;
    private readonly limitInput: DG.InputBase<number | null>;
    private readonly serviceChoiceInput: DG.ChoiceInput<string | null>;
    private readonly runButton: HTMLButtonElement;
    private currentTable: DG.DataFrame | undefined;
    private isInit: boolean = false;

    constructor(parentCall: DG.FuncCall) {
        super();
        this.parentCall = parentCall;
        this.name = 'Service Logs';
        this.tableGridRoot = document.createElement('div');
        this.limitInput = ui.input.int('Limit records', {min: 1, max: 10000, tooltipText: 'Limits number of returned records', value: ServiceLogsApp.DEFAULT_LIMIT, showPlusMinus: false});
        this.serviceChoiceInput = ui.input.choice<string | null>('Service name', {items: [], tooltipText: 'Choose service', nullable: false});
        this.runButton = ui.bigButton('GET LOGS', async () => this.getLogs(), 'Executes request to get service logs');
        const form = ui.forms.condensed([this.serviceChoiceInput, this.limitInput]);
        ui.forms.addButtons(form, [this.runButton]);
        form.style.flexBasis = '15%';
        this.tableGridRoot.style.flexBasis = '85%';
        const content = ui.divV([
            form,
            this.tableGridRoot
        ]);
        content.style.flexBasis = '100%';
        this.root.append(content);
        this.setRibbonPanels([[ui.iconFA('plus', () => grok.shell.addTable(this.currentTable ?? DG.DataFrame.create()),
            'Add table to workspace')]]);
    }

    async init(): Promise<void> {
        this.serviceChoiceInput.items = await grok.dapi.docker.getAvailableServices();
        this.isInit = true;
        await this.getLogs();
    }

    private async getLogs(): Promise<void> {
        if (!this.isInit)
            throw 'ServiceLogsApp wasn\'t init';
        if (!this.serviceChoiceInput.value) {
            new DG.Balloon().error('Service name should be chosen');
            return;
        }
        const progress = DG.TaskBarProgressIndicator.create(`Loading ${this.serviceChoiceInput.value} logs...`);
        this.runButton.disabled = true;
        this.serviceChoiceInput.enabled = false;
        const loader = ui.loader();
        loader.style.top = 'auto';
        try {
            while(this.tableGridRoot.firstChild)
                this.tableGridRoot.removeChild(this.tableGridRoot.firstChild);
            this.tableGridRoot.append(loader);
            this.currentTable = await grok.dapi.docker.getServiceLogs(this.serviceChoiceInput.value!,this.limitInput.value ?? ServiceLogsApp.DEFAULT_LIMIT);
            this.currentTable.getCol(ServiceLogsApp.MESSAGE_COLUMN_NAME).semType = 'Text';
        } catch (e: any) {
            console.error(e);
            new DG.Balloon().error(e.toString());
        } finally {
            this.tableGridRoot.removeChild(loader);
            const grid = DG.Grid.create(this.currentTable ?? DG.DataFrame.create());
            if (grid.table.rowCount > 0) {
                grid.sort([ServiceLogsApp.TIMESTAMP_COLUMN_NAME], [false]);
                grid.col(ServiceLogsApp.TIMESTAMP_COLUMN_NAME)!.width = 175;
                grid.col(ServiceLogsApp.MESSAGE_COLUMN_NAME)!.width = 800;
            }
            this.tableGridRoot.append(grid.root);
            progress.close();
            this.runButton.disabled = false;
            this.serviceChoiceInput.enabled = true;
        }
    }
}
