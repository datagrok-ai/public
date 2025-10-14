import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export class ServiceLogsApp extends DG.ViewBase {
    static readonly DEFAULT_LIMIT: number = 1000;
    static readonly APP_NAME = 'Service Logs';
    private readonly logsRoot: HTMLDivElement;
    private readonly logsRootContent: HTMLDivElement;
    private readonly limitInput: DG.InputBase<number | null>;
    private readonly serviceChoiceInput: DG.ChoiceInput<string | null>;
    private currentLogs: string = '';

    constructor(parentCall: DG.FuncCall, services: string[], service?: string, limit?: number, preview: boolean = false) {
        super();
        this.parentCall = parentCall;
        this.name = preview && service ? service : ServiceLogsApp.APP_NAME;
        this.logsRoot = document.createElement('div');
        this.logsRootContent = document.createElement('div');
        this.logsRoot.classList.add('div-textarea', 'ui-report-textarea', 'ui-log-textarea');
        this.logsRoot.style.display = 'flex';
        this.logsRoot.style.flexDirection = 'column';
        this.logsRoot.style.minHeight = preview && service ? '95%' : '90%';
        this.logsRoot.style.overflowY = 'auto';
        this.logsRootContent.style.display = 'flex';
        this.logsRootContent.style.flexDirection = 'column-reverse';
        this.logsRoot.append(this.logsRootContent);
        this.limitInput = ui.input.int('Limit', {min: 1, max: 100000, tooltipText: 'Limits number of returned records', showPlusMinus: false, value: limit});
        (this.limitInput.input as HTMLInputElement).placeholder = `${ServiceLogsApp.DEFAULT_LIMIT}`;
        this.serviceChoiceInput = ui.input.choice<string | null>('Service', {items: services, tooltipText: 'Choose service',
            nullable: false, value: service ? service : services.length > 0 ? services[0] : null});
        if (preview && service)
            this.serviceChoiceInput.root.style.display = 'none';
        this.subs.push(
            this.serviceChoiceInput.onInput.subscribe(async (_) => {
                if (this.limitInput.input.classList.contains('d4-invalid'))
                    return;
                this.updatePath();
                await this.getLogs();
            })
        );
        const form = ui.forms.condensed([this.serviceChoiceInput, this.limitInput]);
        const content = ui.divV([form, this.logsRoot]);
        this.root.append(content);
        const downloadIcon = ui.iconFA('arrow-to-bottom', () => DG.Utils.download(`${this.serviceChoiceInput.value ?? ''}.logs`, new Blob([this.currentLogs], {
                type: 'text/plain'}),  'text/plain'),
            'Download logs');
        const refreshIcon = ui.iconFA('sync', async () => await this.getLogs(), 'Refresh logs');
        this.setRibbonPanels([[downloadIcon, refreshIcon]]);
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
            while(this.logsRootContent.firstChild)
                this.logsRootContent.removeChild(this.logsRootContent.firstChild);
            this.logsRoot.append(loader);
            this.currentLogs = await grok.dapi
                .docker.getServiceLogs(this.serviceChoiceInput.value!,this.limitInput.value ?? ServiceLogsApp.DEFAULT_LIMIT) ?? '';
        } catch (e: any) {
            console.error(e);
            new DG.Balloon().error(e.toString());
        } finally {
            this.updateLogArea();
            this.logsRoot.removeChild(loader);
            progress.close();
            this.serviceChoiceInput.enabled = true;
            this.limitInput.enabled = true;
        }
    }

    private updateLogArea(): void {
        for (const log of this.currentLogs.split('\n')) {
            this.logsRootContent.append(ui.span([log]));
        }
    }

     private updatePath(): void {
        window.history.pushState(null, `${this.serviceChoiceInput.value}`, this.getAppPath());
    }

    public getAppPath(): string {
        const uri = new URL(window.location.href);
        let segments: string[] = (uri.origin + uri.pathname).split('/');
        if (segments[segments.length - 1] !== 'service-logs')
            segments = segments.slice(0, segments.length - 1);
        const currentService = this.serviceChoiceInput.value;
        if (currentService) {
            if (this.serviceChoiceInput.root.style.display == 'none')
                segments.push('apps', 'usage', 'service-logs');
            segments.push(currentService!);
        }
        return segments.join('/') + `?limit=${this.limitInput.value ?? ServiceLogsApp.DEFAULT_LIMIT}`;
    }
}
