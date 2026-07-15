import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaToolbox} from '../ua-toolbox';
import {UaView} from './ua';

export const VEX_INDEX_URL = 'https://data.datagrok.ai/vex/index.json';

export interface VexImage {
  repo: string;
  category: string;
  tag: string;
  description?: string;
  scanned?: string;
  critical: number;
  high: number;
  medium: number;
  low: number;
  total: number;
  vex: {json: string, csv: string, html: string};
}

const COLOR_CRITICAL = 0xFF7B2D24;
const COLOR_HIGH = 0xFFFF5A00;
const COLOR_MEDIUM = 0xFFFFA500;

export function vexImagesToDataFrame(rows: VexImage[]): DG.DataFrame {
  const int32 = (name: string, get: (r: VexImage) => number) =>
    DG.Column.fromInt32Array(name, Int32Array.from(rows.map((r) => get(r) ?? 0)));
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('channel', rows.map((r) => r.tag === 'bleeding-edge' ? 'bleeding-edge' : 'release')),
    DG.Column.fromStrings('category', rows.map((r) => r.category ?? '')),
    DG.Column.fromStrings('service', rows.map((r) => r.repo)),
    DG.Column.fromStrings('version', rows.map((r) => r.tag)),
    int32('critical', (r) => r.critical),
    int32('high', (r) => r.high),
    int32('medium', (r) => r.medium),
    int32('low', (r) => r.low),
    int32('total', (r) => r.total),
    DG.Column.fromStrings('report', rows.map((r) => r.vex?.html ?? '')),
    DG.Column.fromStrings('description', rows.map((r) => r.description ?? '')),
  ]);
  df.name = 'Vulnerabilities';
  return df;
}

export async function fetchVexImages(): Promise<VexImage[]> {
  const resp = await grok.dapi.fetchProxy(VEX_INDEX_URL);
  const index = await resp.json();
  return [...(index.services ?? []), ...(index.bleeding_edge ?? [])];
}

export class VulnerabilitiesView extends UaView {
  private images: VexImage[] = [];
  private summaryHost!: HTMLElement;
  private detailsHost!: HTMLElement;
  private detailsHeader!: HTMLElement;
  private asOfLabel!: HTMLElement;
  private cards = new Map<string, HTMLElement>();

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Vulnerabilities';
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics';
    this.asOfLabel = ui.divText('as of —', 'ua-metrics-asof');
    this.summaryHost = ui.div([ui.loader()], 'ua-metrics-grid-host');
    this.detailsHeader = ui.divText('Select an image to see its CVEs', 'ua-metrics-panel-title');
    this.detailsHost = ui.div([], 'ua-metrics-grid-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');

    const header = ui.divH([
      this.asOfLabel,
      ui.div([], {style: {flex: '1'}}),
      ui.link('Security docs', 'https://datagrok.ai/help/datagrok/solutions/teams/it/security'),
      ui.link('Full VEX index', 'https://data.datagrok.ai/vex/index.html'),
      refreshBtn,
    ], 'ua-metrics-header');

    const cardsRow = ui.divH(
      ['critical', 'high', 'medium', 'low'].map((s) => this.makeCard(s)), 'ua-metrics-cards-row');

    const summaryPanel = ui.div([
      ui.divH([ui.divText('Scanned images (VEX)', 'ua-metrics-panel-title')], 'ua-metrics-panel-header'),
      this.summaryHost,
    ], 'ua-metrics-panel');
    const detailsPanel = ui.div([
      ui.divH([this.detailsHeader], 'ua-metrics-panel-header'),
      this.detailsHost,
    ], 'ua-metrics-panel');

    this.root.append(ui.div([header, cardsRow, summaryPanel, detailsPanel], 'ua-metrics-root'));
    await this.refresh();
  }

  private makeCard(severity: string): HTMLElement {
    const value = ui.divText('—', 'ua-metrics-card-value');
    const root = ui.div([
      ui.divText(severity[0].toUpperCase() + severity.slice(1), 'ua-metrics-card-title'),
      value,
      ui.divText('released images', 'ua-metrics-card-sub'),
    ], 'ua-metrics-card');
    this.cards.set(severity, value);
    return root;
  }

  private async refresh(): Promise<void> {
    this.summaryHost.innerHTML = '';
    this.summaryHost.append(ui.loader());
    try {
      const resp = await grok.dapi.fetchProxy(VEX_INDEX_URL);
      const index = await resp.json();
      const services: VexImage[] = index.services ?? [];
      const bleedingEdge: VexImage[] = index.bleeding_edge ?? [];
      this.images = [...services, ...bleedingEdge];
      this.asOfLabel.textContent = `as of ${index.generated ?? '—'}`;

      for (const s of ['critical', 'high', 'medium', 'low']) {
        const total = services.reduce((acc, r) => acc + ((r as any)[s] ?? 0), 0);
        this.cards.get(s)!.textContent = `${total}`;
      }

      this.summaryHost.innerHTML = '';
      if (this.images.length === 0) {
        this.summaryHost.append(ui.divText('No published VEX reports found.', 'ua-metrics-degraded'));
        return;
      }
      this.summaryHost.append(this.buildSummaryGrid());
    } catch (e) {
      this.summaryHost.innerHTML = '';
      this.summaryHost.append(ui.divText(`Failed to load ${VEX_INDEX_URL}: ${e}`, 'ua-metrics-degraded'));
    }
  }

  private buildSummaryGrid(): HTMLElement {
    const rows = this.images;
    const df = vexImagesToDataFrame(rows);
    df.onCurrentRowChanged.subscribe(() => {
      if (df.currentRowIdx >= 0)
        this.loadDetails(rows[df.currentRowIdx]);
    });

    const grid = DG.Viewer.grid(df, {
      'showColumnGridlines': false,
      'allowBlockSelection': false,
      'showCurrentCellOutline': false,
    });
    grid.sort(['critical', 'high', 'medium'], [false, false, false]);
    const descCol = grid.col('description');
    if (descCol)
      descCol.width = 320;
    grid.onCellPrepare((gc) => VulnerabilitiesView.colorSeverityCell(gc));
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    return grid.root;
  }

  private static colorSeverityCell(gc: DG.GridCell): void {
    const v = gc.cell?.value as number;
    if (v == null || v === 0)
      return;
    if (gc.gridColumn.name === 'critical')
      gc.style.textColor = COLOR_CRITICAL;
    else if (gc.gridColumn.name === 'high')
      gc.style.textColor = COLOR_HIGH;
    else if (gc.gridColumn.name === 'medium')
      gc.style.textColor = COLOR_MEDIUM;
  }

  private async loadDetails(image: VexImage): Promise<void> {
    this.detailsHeader.innerHTML = '';
    this.detailsHeader.append(
      ui.divText(`${image.repo}:${image.tag}`, 'ua-metrics-panel-title'),
      ui.div([], {style: {width: '12px'}}),
      ui.link('HTML report', image.vex.html),
      ui.div([], {style: {width: '8px'}}),
      ui.link('OpenVEX', image.vex.json),
      ui.div([], {style: {width: '8px'}}),
      ui.link('CSV', image.vex.csv),
    );
    this.detailsHost.innerHTML = '';
    this.detailsHost.append(ui.loader());
    try {
      const resp = await grok.dapi.fetchProxy(image.vex.csv);
      const csv = resp.ok ? await resp.text() : '';
      this.detailsHost.innerHTML = '';
      if (csv.trim().split('\n').length <= 1) {
        this.detailsHost.append(ui.divText('No known CVEs in this image.', 'ua-metrics-degraded'));
        return;
      }
      const df = DG.DataFrame.fromCsv(csv);
      df.name = `${image.repo}:${image.tag} CVEs`;
      const grid = DG.Viewer.grid(df, {'showColumnGridlines': false});
      grid.root.style.width = '100%';
      grid.root.style.height = '100%';
      this.detailsHost.append(grid.root);
    } catch (e) {
      this.detailsHost.innerHTML = '';
      this.detailsHost.append(ui.divText(`Failed to load ${image.vex.csv}: ${e}`, 'ua-metrics-degraded'));
    }
  }
}
