import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from '../tabs/ua';
import {UaToolbox} from '../ua-toolbox';
import {defaultNextVersion, JIRA_BROWSE, ReleaseContext, openInWorkspaceIcon} from './data';

interface JiraIssue { key: string; fields: Record<string, any>; }

/** Builds the JiraConnect filter object for the next-release fixVersion (equality-only JQL builder). */
export function fixVersionFilter(version: string): Record<string, string> {
  return {project: 'GROK', fixVersion: `"${version}"`};
}

export async function fetchReleaseTickets(version: string): Promise<JiraIssue[]> {
  return await grok.functions.call('JiraConnect:GetJiraTicketsByFilter', {filter: fixVersionFilter(version)});
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
  private gridHost!: HTMLElement;
  private titleLabel!: HTMLElement;
  private ticketsDf: DG.DataFrame | null = null;

  constructor(private ctx: ReleaseContext, uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Tickets';
    this.ctx.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics ua-metrics-fill';
    this.versionInput = ui.input.string('Fix version', {value: defaultNextVersion(), onValueChanged: () => this.refresh()});
    this.titleLabel = ui.divText('Tickets marked for the next release', 'ua-metrics-panel-title');
    this.gridHost = ui.div([], 'ua-metrics-grid-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh());
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([
      this.versionInput.root, ui.div([], {style: {flex: '1'}}), refreshBtn,
    ], 'ua-metrics-header');

    const panel = ui.div([ui.divH([this.titleLabel, openInWorkspaceIcon(() => this.ticketsDf)], 'ua-metrics-panel-header'),
      this.gridHost], 'ua-metrics-panel ua-metrics-panel-grow');
    this.root.append(ui.div([header, panel], 'ua-metrics-root'));
    await this.refresh();
  }

  private async refresh(): Promise<void> {
    const version = (this.versionInput.value ?? '').trim();
    this.gridHost.innerHTML = '';
    this.gridHost.append(ui.loader());
    if (!version) {
      this.gridHost.innerHTML = '';
      this.gridHost.append(ui.divText('Enter the next-release fix version.', 'ua-metrics-degraded'));
      return;
    }
    try {
      const issues = await fetchReleaseTickets(version);
      this.gridHost.innerHTML = '';
      if (!issues || issues.length === 0) {
        this.gridHost.append(ui.divText(
          `No GROK tickets with fixVersion "${version}" (or JiraConnect credentials are not configured on this server).`,
          'ua-metrics-degraded'));
        this.titleLabel.textContent = `fixVersion "${version}" — 0 tickets`;
        return;
      }
      const df = ticketsToDataFrame(issues);
      this.ticketsDf = df;
      const open = issues.filter((r) => (r.fields.resolution?.name ?? 'Unresolved') === 'Unresolved').length;
      this.titleLabel.textContent = `fixVersion "${version}" — ${issues.length} tickets, ${open} unresolved`;
      const grid = DG.Viewer.grid(df, {showColumnGridlines: false, allowBlockSelection: false});
      grid.sort(['status', 'priority'], [true, true]);
      const summaryCol = grid.col('summary');
      if (summaryCol)
        summaryCol.width = 360;
      grid.root.style.width = '100%';
      grid.root.style.height = '100%';
      this.gridHost.append(grid.root);
    } catch (e) {
      this.gridHost.innerHTML = '';
      this.gridHost.append(ui.divText(
        `Could not reach JiraConnect (${e}). Install the JiraConnect package and configure Jira credentials.`,
        'ua-metrics-degraded'));
    }
  }
}
