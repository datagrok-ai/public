import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';

/** Where the events come from. See `core/docs/LOG_ARCHITECTURE.md` for the two cloud tiers. */
export const SOURCE_PLATFORM = 'Platform events';
export const SOURCE_CLOUDWATCH = 'CloudWatch';
export const SOURCE_ARCHIVE = 'Archive';

const RANGES: {[key: string]: number} = {
  'Last 15 minutes': 15,
  'Last hour': 60,
  'Last 4 hours': 240,
  'Last day': 24 * 60,
  'Last week': 7 * 24 * 60,
};

/**
 * Explorer over the three log tiers: this instance's own events, the CloudWatch
 * hot tier (30 days, queryable), and the write-once S3 archive (cold, listed by
 * key). The two cloud tiers are admin-only server-side — the archive spans the
 * whole AWS estate, so it is other tenants' data on a tenant instance.
 */
export class CloudLogsApp extends DG.ViewBase {
  static readonly APP_NAME = 'Cloud Logs';
  static readonly DEFAULT_LIMIT = 1000;

  private readonly sourceInput: DG.ChoiceInput<string | null>;
  private readonly connectionInput: DG.ChoiceInput<string | null>;
  private readonly groupInput: DG.ChoiceInput<string | null>;
  private readonly prefixInput: DG.InputBase<string>;
  private readonly rangeInput: DG.ChoiceInput<string | null>;
  private readonly filterInput: DG.InputBase<string>;
  private readonly limitInput: DG.InputBase<number | null>;

  private readonly objectsHost = ui.box();
  private readonly eventsHost = ui.box();
  private connections: DG.DataConnection[] = [];
  private events = DG.DataFrame.create();

  constructor(parentCall: DG.FuncCall) {
    super();
    this.parentCall = parentCall;
    this.name = CloudLogsApp.APP_NAME;

    this.sourceInput = ui.input.choice<string | null>('Source', {
      items: [SOURCE_PLATFORM, SOURCE_CLOUDWATCH, SOURCE_ARCHIVE],
      value: SOURCE_CLOUDWATCH, nullable: false,
      onValueChanged: () => this.onSourceChanged(),
    });
    this.connectionInput = ui.input.choice<string | null>('Connection', {
      items: [], nullable: true,
      tooltipText: 'S3 or AWS connection. Leave empty to use the instance role.',
      onValueChanged: () => this.loadGroups(),
    });
    this.groupInput = ui.input.choice<string | null>('Group', {items: [], nullable: true});
    this.prefixInput = ui.input.string('Prefix', {
      value: 'cw/',
      tooltipText: 'Key prefix in the archive bucket. cw/ holds the CloudWatch deliveries ' +
        '(cloudwatch/ once the CloudFormation template is applied), alb/ and vpc-flow/ the ' +
        'rest. Narrow by date: cw/2026/08/26/',
    });
    this.rangeInput = ui.input.choice<string | null>('Time', {
      items: Object.keys(RANGES), value: 'Last hour', nullable: false,
    });
    this.filterInput = ui.input.string('Filter', {
      tooltipText: 'CloudWatch filter pattern, not a regular expression: a bare word matches a ' +
        'substring, ?a ?b is OR, -x excludes',
    });
    this.limitInput = ui.input.int('Limit', {
      min: 1, max: 100000, showPlusMinus: false, value: CloudLogsApp.DEFAULT_LIMIT,
    });

    const form = ui.forms.condensed([this.sourceInput, this.connectionInput, this.groupInput,
      this.prefixInput, this.rangeInput, this.filterInput, this.limitInput]);
    this.root.append(ui.splitV([
      ui.box(form, {style: {maxHeight: '190px'}}),
      ui.splitH([this.objectsHost, this.eventsHost], null, true),
    ]));

    this.setRibbonPanels([[
      ui.icons.sync(() => this.run(), 'Run'),
      ui.iconFA('arrow-to-bottom', () => this.download(), 'Download as CSV'),
    ]]);
    this.onSourceChanged();
    this.init();
  }

  private get source(): string {
    return this.sourceInput.value ?? SOURCE_CLOUDWATCH;
  }

  private async init(): Promise<void> {
    this.connections = (await grok.dapi.connections.list())
      .filter((c) => c.dataSource === 'S3');
    this.connectionInput.items = this.connections.map((c) => c.nqName);
    await this.loadGroups();
  }

  /** Only the inputs the current source actually uses stay on the form. */
  private onSourceChanged(): void {
    const cloud = this.source !== SOURCE_PLATFORM;
    const archive = this.source === SOURCE_ARCHIVE;
    this.connectionInput.root.style.display = cloud ? '' : 'none';
    this.groupInput.root.style.display = this.source === SOURCE_CLOUDWATCH ? '' : 'none';
    this.prefixInput.root.style.display = archive ? '' : 'none';
    this.filterInput.root.style.display = archive ? 'none' : '';
    this.rangeInput.root.style.display = archive ? 'none' : '';
    this.objectsHost.style.display = archive ? '' : 'none';
  }

  private connectionId(): string {
    const name = this.connectionInput.value;
    return this.connections.find((c) => c.nqName === name)?.id ?? '';
  }

  private async loadGroups(): Promise<void> {
    if (this.source !== SOURCE_CLOUDWATCH)
      return;
    try {
      const groups = await grok.dapi.log.getCloudLogGroups({connection: this.connectionId()});
      this.groupInput.items = groups;
      this.groupInput.value ??= groups[0];
    } catch (e: any) {
      this.groupInput.items = [];
      grok.shell.warning(`Could not list log groups: ${e.message ?? e}`);
    }
  }

  async run(): Promise<void> {
    const progress = DG.TaskBarProgressIndicator.create(`Loading ${this.source} logs...`);
    try {
      if (this.source === SOURCE_ARCHIVE)
        await this.showArchiveObjects();
      else
        this.showEvents(this.source === SOURCE_CLOUDWATCH
          ? await this.queryCloudWatch() : await this.queryPlatform());
    } catch (e: any) {
      grok.shell.error(e.message ?? String(e));
    } finally {
      progress.close();
    }
  }

  private range(): {start: dayjs.Dayjs, end: dayjs.Dayjs} {
    const end = dayjs();
    return {start: end.subtract(RANGES[this.rangeInput.value!] ?? 60, 'minute'), end: end};
  }

  private queryCloudWatch(): Promise<DG.DataFrame> {
    if (!this.groupInput.value)
      throw new Error('Choose a log group');
    const {start, end} = this.range();
    return grok.dapi.log.getCloudLogEvents(this.groupInput.value, start, end, {
      connection: this.connectionId(),
      filter: this.filterInput.value,
      limit: this.limitInput.value ?? CloudLogsApp.DEFAULT_LIMIT,
    });
  }

  /** Projected onto the cloud sources' columns so switching source keeps the grid. */
  private async queryPlatform(): Promise<DG.DataFrame> {
    const {start, end} = this.range();
    const events = await grok.dapi.log.where({start: start, end: end})
      .by(this.limitInput.value ?? CloudLogsApp.DEFAULT_LIMIT).list();
    const text = this.filterInput.value?.toLowerCase();
    const rows = events
      .map((e) => ({
        time: e.eventTime,
        source: e.eventType?.friendlyName ?? '',
        message: e.description ?? e.name ?? '',
      }))
      .filter((r) => !text || r.message.toLowerCase().includes(text));
    return DG.DataFrame.fromColumns([
      DG.Column.dateTime('time', rows.length).init((i) => rows[i].time),
      DG.Column.fromStrings('source', rows.map((r) => r.source)),
      DG.Column.fromStrings('stream', rows.map(() => '')),
      DG.Column.fromStrings('message', rows.map((r) => r.message)),
    ]);
  }

  /** The archive is two-step: list keys, then decode the one you pick. */
  private async showArchiveObjects(): Promise<void> {
    const connection = this.connectionId();
    if (!connection)
      throw new Error('Choose the connection that points at the archive bucket');
    const objects = await grok.dapi.log.getArchiveObjects(connection, {
      prefix: this.prefixInput.value,
      limit: this.limitInput.value ?? CloudLogsApp.DEFAULT_LIMIT,
    });
    ui.empty(this.objectsHost);
    if (objects.rowCount === 0) {
      this.objectsHost.append(ui.divText('No objects under this prefix'));
      return;
    }
    objects.onCurrentRowChanged.subscribe(async () => {
      const key = objects.get('key', objects.currentRowIdx);
      const progress = DG.TaskBarProgressIndicator.create(`Decoding ${key}...`);
      try {
        this.showEvents(await grok.dapi.log.getArchiveEvents(connection, key));
      } catch (e: any) {
        grok.shell.error(e.message ?? String(e));
      } finally {
        progress.close();
      }
    });
    this.objectsHost.append(DG.Viewer.grid(objects, {showRowHeader: false}).root);
  }

  private showEvents(events: DG.DataFrame): void {
    this.events = events;
    ui.empty(this.eventsHost);
    if (events.rowCount === 0) {
      this.eventsHost.append(ui.divText('No events'));
      return;
    }
    const grid = DG.Viewer.grid(events, {showRowHeader: false, allowBlockSelection: false});
    grid.onCellPrepare((gc) => {
      if (gc.gridColumn.name === 'message')
        gc.style.font = '13px monospace';
    });
    this.eventsHost.append(grid.root);
  }

  private download(): void {
    if (this.events.rowCount === 0) {
      grok.shell.warning('Nothing to download');
      return;
    }
    DG.Utils.download(`${this.source.toLowerCase().replace(' ', '-')}.csv`,
      this.events.toCsv(), 'text/csv');
  }
}
