import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {debounceTime} from 'rxjs/operators';

import {UaView} from '../tabs/ua';
import {UaToolbox} from '../ua-toolbox';
import {defaultNextVersion, ReleaseContext, openInWorkspaceIcon} from './data';
import {GitHubIssue, isOpen} from '../github/issue';
import {issuesGrid} from '../github/issues-grid';

const MILESTONE_TYPING_MS = 500;
const MAIN_LABEL = 'main';

export function isMain(issue: GitHubIssue): boolean {
  return issue.labels.some((l) => l.toLowerCase() === MAIN_LABEL);
}

export interface TicketCounts {
  total: number;
  open: number;
  main: number;
  mainOpen: number;
}

export function countTickets(issues: GitHubIssue[]): TicketCounts {
  const main = issues.filter(isMain);
  return {total: issues.length, open: issues.filter(isOpen).length, main: main.length,
    mainOpen: main.filter(isOpen).length};
}

export class ReleaseTicketsView extends UaView {
  milestoneInput!: DG.InputBase<string>;
  titleLabel!: HTMLElement;
  mainTitle!: HTMLElement;
  otherTitle!: HTMLElement;
  mainHost!: HTMLElement;
  otherHost!: HTMLElement;
  otherDf: DG.DataFrame | null = null;
  run = 0;

  constructor(private ctx: ReleaseContext, uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Tickets';
    this.ctx.refresh.subscribe(() => { if (this.initialized) this.refresh(); });
  }

  async initViewers(): Promise<void> {
    this.root.className = 'grok-view ui-box ua-metrics ua-metrics-fill';
    this.milestoneInput = ui.input.string('Milestone', {value: defaultNextVersion()});
    this.subs.push(this.milestoneInput.onChanged.pipe(debounceTime(MILESTONE_TYPING_MS))
      .subscribe(() => this.refresh()));
    this.titleLabel = ui.divText('Issues of core and public in the release milestone', 'ua-metrics-asof');
    this.mainTitle = ui.divText('MAIN', 'ua-metrics-panel-title');
    this.otherTitle = ui.divText('Other', 'ua-metrics-panel-title');
    this.mainHost = ui.div([], 'ua-metrics-grid-host');
    this.otherHost = ui.div([], 'ua-metrics-grid-host');

    const refreshBtn = ui.bigButton('Refresh', () => this.refresh(true));
    refreshBtn.prepend(ui.iconFA('sync-alt'), ' ');
    refreshBtn.classList.add('ua-metrics-btn-secondary');
    const header = ui.divH([this.milestoneInput.root, this.titleLabel,
      ui.div([], {style: {flex: '1'}}), refreshBtn], 'ua-metrics-header');

    const mainPanel = ui.div([ui.divH([this.mainTitle], 'ua-metrics-panel-header'), this.mainHost], 'ua-metrics-panel');
    const otherPanel = ui.div([ui.divH([this.otherTitle, openInWorkspaceIcon(() => this.otherDf)],
      'ua-metrics-panel-header'), this.otherHost], 'ua-metrics-panel ua-metrics-panel-grow');
    this.root.append(ui.div([header, mainPanel, otherPanel], 'ua-metrics-root'));
    await this.refresh();
  }

  async refresh(force = false): Promise<void> {
    const milestone = this.milestoneInput.value.trim();
    const run = ++this.run;
    this.mainHost.replaceChildren();
    this.otherHost.replaceChildren(ui.loader());
    this.otherDf = null;
    if (!milestone) {
      this.otherHost.replaceChildren(ui.divText('Enter the release milestone.', 'ua-metrics-degraded'));
      return;
    }
    let issues: GitHubIssue[];
    try {
      issues = await this.ctx.tickets(milestone, force);
    } catch (e) {
      if (run === this.run) {
        this.titleLabel.textContent = `Milestone "${milestone}"`;
        this.otherHost.replaceChildren(ui.divText(`Could not read GitHub issues: ${e instanceof Error ? e.message : e}`,
          'ua-metrics-degraded'));
      }
      return;
    }
    if (run === this.run)
      this.render(milestone, issues);
  }

  render(milestone: string, issues: GitHubIssue[]): void {
    const c = countTickets(issues);
    this.titleLabel.textContent = `Milestone "${milestone}" — ${c.total} issues, ${c.open} open (${c.mainOpen} MAIN)`;
    this.mainTitle.textContent = `MAIN — ${c.mainOpen} open of ${c.main}`;
    this.otherTitle.textContent = `Other — ${c.open - c.mainOpen} open of ${c.total - c.main}`;
    this.renderGrid(this.mainHost, issues.filter(isMain), 'No MAIN-labelled issues.');
    this.otherDf = this.renderGrid(this.otherHost, issues.filter((i) => !isMain(i)),
      c.total === 0 ? `No issues in milestone "${milestone}" of core or public.` : 'No other issues.');
  }

  renderGrid(host: HTMLElement, issues: GitHubIssue[], emptyMessage: string): DG.DataFrame | null {
    if (issues.length === 0) {
      host.replaceChildren(ui.divText(emptyMessage, 'ua-metrics-degraded'));
      return null;
    }
    const grid = issuesGrid(issues, 'Release tickets');
    host.replaceChildren(grid.root);
    return grid.dataFrame;
  }
}
