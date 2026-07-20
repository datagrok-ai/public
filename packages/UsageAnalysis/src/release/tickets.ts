import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {UaToolbox} from '../ua-toolbox';
import {defaultNextVersion, JIRA_BROWSE, ReleaseContext, openInWorkspaceIcon} from './data';

interface JiraIssue { key: string; fields: Record<string, any>; }

export const MAIN_LABEL = 'Main';

/** Builds the JiraConnect filter object for the next-release fixVersion (equality-only JQL builder). */
export function fixVersionFilter(version: string): Record<string, string> {
  return {project: 'GROK', fixVersion: `"${version}"`};
}

export async function fetchReleaseTickets(version: string): Promise<JiraIssue[]> {
  return await grok.functions.call('JiraConnect:GetJiraTicketsByFilter', {filter: fixVersionFilter(version)});
}

/** A ticket still needs attention unless it's Done or Won't Fix (by either status or resolution). */
export function isActionable(issue: JiraIssue): boolean {
  const status = issue.fields.status?.name ?? '';
  const res = issue.fields.resolution?.name ?? '';
  return status !== 'Done' && status !== 'Won\'t Fix' && res !== 'Done' && res !== 'Won\'t Fix';
}

/** True when the ticket carries the "Main" label (case-insensitive). */
export function hasMainLabel(issue: JiraIssue): boolean {
  const target = MAIN_LABEL.toLowerCase();
  return (issue.fields.labels || []).some((l: string) => (l ?? '').toLowerCase() === target);
}

// Grid cell backgrounds (ARGB — grid cells take numeric colors, not design tokens).
function statusColor(s: string): number {
  const t = (s || '').toLowerCase();
  if (t.includes('done') || t.includes('closed') || t.includes('resolved')) return 0xFFE6F4EA; // green
  if (t.includes('progress') || t.includes('review') || t.includes('testing')) return 0xFFFFF3E0; // amber
  return 0xFFF3F3F3; // to-do / open / backlog — grey
}
function priorityColor(p: string): number {
  const t = (p || '').toLowerCase();
  if (t.includes('blocker') || t.includes('highest') || t.includes('critical')) return 0xFFFDE7E7; // red
  if (t === 'high') return 0xFFFCE7D6; // orange
  if (t.includes('medium')) return 0xFFFCF7DA; // yellow
  return 0xFFF3F3F3; // low / lowest — grey
}
function resolutionColor(r: string): number {
  if (r === 'Done') return 0xFFE6F4EA;
  if (r === 'Won\'t Fix') return 0xFFF0F0F0;
  return 0xFFFFF3E0; // unresolved / other — amber (still actionable)
}

export function ticketsToDataFrame(issues: JiraIssue[]): DG.DataFrame {
  const str = (get: (r: JiraIssue) => any) => issues.map((r) => `${get(r) ?? ''}`);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('ticket', issues.map((r) => `<a href="${JIRA_BROWSE}${r.key}" target="_blank">${r.key}</a>`)),
    DG.Column.fromStrings('summary', str((r) => r.fields.summary)),
    DG.Column.fromStrings('status', str((r) => r.fields.status?.name)),
    DG.Column.fromStrings('priority', str((r) => r.fields.priority?.name)),
    DG.Column.fromStrings('assignee', str((r) => r.fields.assignee?.displayName)),
    DG.Column.fromStrings('resolution', issues.map((r) => r.fields.resolution?.name ?? 'Unresolved')),
    DG.Column.fromStrings('type', str((r) => r.fields.issuetype?.name)),
  ]);
  df.name = 'Release tickets';
  df.getCol('ticket').setTag('cell.renderer', 'html');
  return df;
}

export class ReleaseTicketsView extends UaView {
  private versionInput!: DG.InputBase;
  private titleLabel!: HTMLElement;
  private mainTitle!: HTMLElement;
  private otherTitle!: HTMLElement;
  private mainHost!: HTMLElement;
  private otherHost!: HTMLElement;
  private ticketsDf: DG.DataFrame | null = null;

  constructor(private ctx: ReleaseContext, uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Tickets';
    this.ctx.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics ua-metrics-fill';
    this.versionInput = ui.input.string('Fix version',
      {value: defaultNextVersion(), onValueChanged: () => this.refresh()});
    this.titleLabel = ui.divText('Tickets marked for the next release', 'ua-metrics-asof');
    this.mainTitle = ui.divText('MAIN', 'ua-metrics-panel-title');
    this.otherTitle = ui.divText('Other', 'ua-metrics-panel-title');
    this.mainHost = ui.div([], 'ua-metrics-grid-host');
    this.otherHost = ui.div([], 'ua-metrics-grid-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([this.versionInput.root, this.titleLabel,
      ui.div([], {style: {flex: '1'}}), refreshBtn], 'ua-metrics-header');

    // MAIN-labelled tickets are shown in their own panel, above the rest.
    const mainPanel = ui.div([ui.divH([this.mainTitle], 'ua-metrics-panel-header'), this.mainHost], 'ua-metrics-panel');
    const otherPanel = ui.div([ui.divH([this.otherTitle, openInWorkspaceIcon(() => this.ticketsDf)],
      'ua-metrics-panel-header'), this.otherHost], 'ua-metrics-panel ua-metrics-panel-grow');
    this.root.append(ui.div([header, mainPanel, otherPanel], 'ua-metrics-root'));
    await this.refresh();
  }

  private async refresh(): Promise<void> {
    const version = (this.versionInput.value ?? '').trim();
    this.mainHost.innerHTML = '';
    this.otherHost.innerHTML = '';
    this.otherHost.append(ui.loader());
    if (!version) {
      this.otherHost.innerHTML = '';
      this.otherHost.append(ui.divText('Enter the next-release fix version.', 'ua-metrics-degraded'));
      return;
    }
    try {
      const issues = await fetchReleaseTickets(version);
      this.otherHost.innerHTML = '';
      if (!issues || issues.length === 0) {
        this.otherHost.append(ui.divText(
          `No GROK tickets with fixVersion "${version}" (or JiraConnect credentials are not configured on this server).`,
          'ua-metrics-degraded'));
        this.titleLabel.textContent = `fixVersion "${version}" — 0 tickets`;
        return;
      }
      const main = issues.filter(hasMainLabel);
      const other = issues.filter((r) => !hasMainLabel(r));
      const open = issues.filter(isActionable).length;
      const mainOpen = main.filter(isActionable).length;
      this.titleLabel.textContent =
        `fixVersion "${version}" — ${issues.length} tickets, ${open} open (${mainOpen} MAIN)`;
      this.mainTitle.textContent = `MAIN — ${mainOpen} open of ${main.length}`;
      this.otherTitle.textContent = `Other — ${open - mainOpen} open of ${other.length}`;
      this.renderGrid(this.mainHost, main, 'No MAIN-labelled tickets.');
      this.ticketsDf = this.renderGrid(this.otherHost, other, 'No other tickets.') ?? this.ticketsDf;
    } catch (e) {
      this.otherHost.innerHTML = '';
      this.otherHost.append(ui.divText(
        `Could not reach JiraConnect (${e}). Install the JiraConnect package and configure Jira credentials.`,
        'ua-metrics-degraded'));
    }
  }

  /** Builds a color-coded tickets grid into `host`; returns the dataframe (or null if empty). */
  private renderGrid(host: HTMLElement, issues: JiraIssue[], emptyMsg: string): DG.DataFrame | null {
    host.innerHTML = '';
    if (issues.length === 0) {
      host.append(ui.divText(emptyMsg, 'ua-metrics-degraded'));
      return null;
    }
    const df = ticketsToDataFrame(issues);
    const grid = DG.Viewer.grid(df, {showColumnGridlines: false, allowBlockSelection: false});
    grid.onCellPrepare((gc) => {
      if (!gc.isTableCell)
        return;
      const v = gc.cell.value as string;
      if (gc.gridColumn.name === 'status')
        gc.style.backColor = statusColor(v);
      else if (gc.gridColumn.name === 'priority')
        gc.style.backColor = priorityColor(v);
      else if (gc.gridColumn.name === 'resolution')
        gc.style.backColor = resolutionColor(v);
    });
    grid.sort(['priority', 'status'], [true, true]);
    const summaryCol = grid.col('summary');
    if (summaryCol)
      summaryCol.width = 360;
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    host.append(grid.root);
    return df;
  }
}
