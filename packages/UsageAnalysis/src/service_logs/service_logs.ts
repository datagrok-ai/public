import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as rx from "rxjs";

export class ServiceLogsApp extends DG.ViewBase {
    static readonly DEFAULT_LIMIT: number = 1000;
    static readonly TIMESTAMP_COLUMN_NAME = 'timestamp';
    static readonly MESSAGE_COLUMN_NAME = 'message';

    private readonly tableGridRoot: HTMLDivElement;
    private readonly limitInput: DG.InputBase<number | null>;
    private readonly serviceChoiceInput: DG.ChoiceInput<string | null>;
    private currentTable: DG.DataFrame | undefined;

    constructor(parentCall: DG.FuncCall, services: string[], service?: string, limit?: number) {
        super();
        this.parentCall = parentCall;
        this.name = 'Service Logs';
        this.tableGridRoot = document.createElement('div');
        this.updateGrid();
        this.limitInput = ui.input.int('Limit', {min: 1, max: 100000, tooltipText: 'Limits number of returned records', showPlusMinus: false, value: limit});
        (this.limitInput.input as HTMLInputElement).placeholder = `${ServiceLogsApp.DEFAULT_LIMIT}`;
        this.serviceChoiceInput = ui.input.choice<string | null>('Service', {items: services, tooltipText: 'Choose service',
            nullable: false, value: service ? service : services.length > 0 ? services[0] : null});
        rx.merge(this.limitInput.onInput, this.serviceChoiceInput.onInput).subscribe(async (_) => {
            if (this.limitInput.input.classList.contains('d4-invalid'))
                return;
            this.updatePath();
            await this.getLogs();
        });
        const form = ui.forms.condensed([this.serviceChoiceInput, this.limitInput]);
        form.style.flexBasis = '10%';
        this.tableGridRoot.style.flexBasis = '90%';
        const content = ui.divV([form, this.tableGridRoot]);
        content.style.flexBasis = '100%';
        this.root.append(content);
        this.setRibbonPanels([[ui.iconFA('plus', () => grok.shell.addTable(this.currentTable ?? DG.DataFrame.create()),
            'Add table to workspace')]]);
        setTimeout(() => this.updatePath(), 300);
    }

     async getLogs(): Promise<void> {
        if (!this.serviceChoiceInput.value) {
            new DG.Balloon().error('Service name should be chosen');
            return;
        }
        const progress = DG.TaskBarProgressIndicator.create(`Loading ${this.serviceChoiceInput.value} logs...`);
        this.serviceChoiceInput.enabled = false;
        this.limitInput.enabled = false;
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
            this.updateGrid();
            progress.close();
            this.serviceChoiceInput.enabled = true;
            this.limitInput.enabled = true;
        }
    }

    private updateGrid(): void {
        const grid = DG.Grid.create(this.currentTable ?? DG.DataFrame.create());
        if (grid.table.rowCount > 0) {
            grid.sort([ServiceLogsApp.TIMESTAMP_COLUMN_NAME], [false]);
            grid.col(ServiceLogsApp.TIMESTAMP_COLUMN_NAME)!.width = 175;
            grid.col(ServiceLogsApp.MESSAGE_COLUMN_NAME)!.width = 800;
        }
        this.tableGridRoot.append(grid.root);
    }

     private updatePath(): void {
        let segments: string[] = window.location.href.split('/');
        if (segments[segments.length - 1] !== 'service-logs')
            segments = segments.slice(0, segments.length - 1);
        const currentService = this.serviceChoiceInput.value;
        if (currentService)
            segments.push(currentService!);
        let url: string = segments.join('/') + `?limit=${this.limitInput.value ?? ServiceLogsApp.DEFAULT_LIMIT}`;
        window.history.pushState(null, `${currentService}`, url);
    }
}
