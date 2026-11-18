import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { _package, _properties } from '../package';
import {findIssueRefs, GitHubIssue, IssueRef, isOpen} from '../github/issue';
import {GitHubIssuesClient} from '../github/issues-client';
import {issuesGrid} from '../github/issues-grid';

const MANUAL_TICKETS_FILE = 'manual-tickets.json';

class Priority {
  static BLOCKER: string = 'BLOCKER';
  static INFO: string = 'INFO';
  static ERROR: string = 'ERROR';
  static RESOLVED: string = 'RESOLVED';

  static codeMapping: { [key: string]: number } = {
    [Priority.BLOCKER]: 3,
    [Priority.ERROR]: 2,
    [Priority.INFO]: 1,
    [Priority.RESOLVED]: 0,
  };

  static getMaxPriority(priorities: string[]): string {
    if (!priorities || priorities.length === 0)
      return '';

    return priorities.reduce((maxPriority, currentPriority) => {
      const maxCode = Priority.codeMapping[maxPriority] || 0;
      const currentCode = Priority.codeMapping[currentPriority] || 0;
      return currentCode > maxCode ? currentPriority : maxPriority;
    });
  }
}

function verdictPriority(issue: GitHubIssue, platformVersion: string | undefined): string {
  if (!isOpen(issue))
    return Priority.RESOLVED;
  if (issue.priority === 'Blocker')
    return Priority.BLOCKER;
  const thisVersion = !!platformVersion && !!issue.milestone?.startsWith(platformVersion);
  if (issue.type === 'Bug' && issue.priority !== 'Low' && thisVersion)
    return Priority.ERROR;
  return Priority.INFO;
}

class Verdict {
  priority: string;
  category: string;
  widget: HTMLElement;
  constructor(priority: string, category: string, widget: HTMLElement) {
    this.priority = priority;
    this.widget = widget;
    this.category = category;
  }
}

function getIconForVerdict(priority: string): HTMLElement {
  let icon: HTMLElement | undefined;
  if (priority == Priority.BLOCKER) {
    icon = ui.iconFA('exclamation-triangle', null, 'Is considered a blocker for the release');
    icon.style.color = DG.Color.toRgb(DG.Color.red);
  } else if (priority == Priority.ERROR) {
    icon = ui.iconFA('exclamation-circle', null, 'Errors to fix');
    icon.style.color = DG.Color.toRgb(DG.Color.red);
  } else if (priority == Priority.INFO) {
    icon = ui.iconFA('info-circle', null, 'Good to know, doesn\'t require immediate action');
    icon.style.color = DG.Color.toRgb(DG.Color.blue);
  } if (priority == Priority.RESOLVED) {
    icon = ui.iconFA('info-circle', null, 'Already resolved');
    icon.style.color = DG.Color.toRgb(DG.Color.green);
  }
  return icon!;
}

export class TestDashboardWidget extends DG.JsViewer {
  static isOpen: boolean = false;
  ticketRefs = new Map<string, IssueRef>();
  githubIssues: GitHubIssuesClient;

  constructor(githubIssues: GitHubIssuesClient) {
    super();
    this.githubIssues = githubIssues;
  }

  onFrameAttached(dataFrame: DG.DataFrame): void {
    this.root.childNodes.forEach((node, _idx, _parent) => node.remove());
    this.root.appendChild(ui.wait(async () => {
      let tableNames: RegExp[] = [/Benchmark/];
      let promises: Promise<void>[] = tableNames.map((name) => new Promise((resolve, reject) => {
        const checkCondition = () => {
          try {
            if (grok.shell.tables.find((df) => df.name.match(name) != null)) {
              resolve();
            } else {
              setTimeout(checkCondition, 1000);
            }
          } catch (error) {
            reject(error);
          }
        };
        checkCondition();
      }));
      await Promise.all(promises);
      let d = ui.div();
      d.style.padding = '15px';
      let verdicts: Verdict[] = [];
      verdicts.push(...this.verdictsOnCorrectness(dataFrame));
      verdicts.push(...this.verdictsOnPerformance());
      verdicts.push(...(await this.verdictsForTickets()));
      verdicts.push(...this.getBuildInfos(dataFrame));
      let status = Priority.getMaxPriority(verdicts.map((v) => v.priority));
      if (status == Priority.BLOCKER)
        d.append(ui.h1('Platform is not release-ready ❌❌❌', { style: { color: DG.Color.toRgb(DG.Color.red), marginBottom: '5px' } }));
      else
        d.append(ui.h1('Green light! 🚀🚀🚀', { style: { color: DG.Color.toRgb(DG.Color.darkGreen) } }));
      
      verdicts.sort((a, b) => {
        let priorityA = Priority.codeMapping[a.priority];
        let priorityB = Priority.codeMapping[b.priority];
        if (priorityA != priorityB)
          return priorityB - priorityA;
        return a.category == b.category ? 0 :
              (a.category < b.category ? 1 : -1);
      });

      let list: HTMLElement = this.composeResultList(verdicts);
      d.append(list);
      return d;
    }));
  }

  verdictsOnCorrectness(df: DG.DataFrame): Verdict[] {
    let verdicts: Verdict[] = [];
    if (df.col('jira') != null)
      this.collectTickets(df.col('jira')!);
    if (df.col('ignore?') != null)
      verdicts.push(...this.ignoredTests(df));
    verdicts.push(...this.majorPackages(df));
    return verdicts;
  }
  verdictsOnPerformance(): Verdict[] {
    let verdicts: Verdict[] = [];
    let df = grok.shell.tables.find((df) => df.name.match(/Benchmark/));
    if (df == null)
      return verdicts;
    verdicts.push(...this.performanceDowngrades(df));
    return verdicts;
  }
  collectTickets(ticketCol: DG.Column<string>): void {
    for (let i = 0; i < ticketCol.length; i++)
      this.watch(findIssueRefs(ticketCol.getString(i)));
  }

  watch(refs: IssueRef[]): void {
    for (const ref of refs)
      this.ticketRefs.set(ref.key, ref);
  }
  majorPackages(df: DG.DataFrame): Verdict[] {
    try {
      let packageList = [...(new Set(df.col('test')?.categories.map(e => e.split(':')[0]) ?? []))];
      let testColumn: DG.Column = df.col('test')!;
      let failingColumn: DG.Column<boolean> = df.col('needs_attention')!;
      let brokenPackages: Map<string, string[]> = new Map();
      let checkedPackages = new Set();
      for (let i = 0; i < packageList.length; i++)
        for (let row = 0; row < testColumn.length; row++)
          if (testColumn.getString(row)?.startsWith(packageList[i])) {
            checkedPackages.add(packageList[i]);
            if (failingColumn.get(row)) {
              if (!brokenPackages.has(packageList[i]))
                brokenPackages.set(packageList[i], []);
              brokenPackages.get(packageList[i])!.push(testColumn.getString(row).slice(packageList[i].length));
            }
          }

      let verdicts: Verdict[] = [];
      for (const packageName of brokenPackages.keys()) {
        let packageFilter= '';
        let widget: DG.Widget = DG.Widget.fromRoot(ui.span([packageName + ' has broken tests']));
        widget.root.onclick = (ev: MouseEvent) => {
          grok.shell.tv.grid.scrollToCell('test', 0);
          if (packageFilter.length === 0){
            df.filter.init((row) => testColumn.getString(row)?.startsWith(packageName));
            packageFilter = packageName;
          }
          else {
            df.filter.init((row) => true);
            packageFilter = '';
          }
          df.selection.init((row) => testColumn.getString(row)?.startsWith(packageName) && (failingColumn.get(row) ?? false));
          grok.shell.tv.grid.sort(['needs_attention'], [false]);
        };
        ui.tooltip.bind(widget.root, 'Click to filter');
        verdicts.push(new Verdict(Priority.BLOCKER, 'Critical package failure', widget.root));
      }

      if (!packageList.every((p) => checkedPackages.has(p))) {
        for (const packageName of packageList.filter((p) => !checkedPackages.has(p))) {
          let widget: DG.Widget = DG.Widget.fromRoot(ui.span([packageName + ' tests weren\'t run']));
          verdicts.push(new Verdict(Priority.BLOCKER, 'Missing required packages tests', widget.root));
        }
      }

      return verdicts;
    } catch (x) {
      return [new Verdict(Priority.BLOCKER, 'Unhandled exception', ui.div([
        ui.span(['Failed to verify major packages']),
        ui.span([x])
      ]))];
    }
  }
  getBuildInfos(df: DG.DataFrame): Verdict[] {
    try {
      let builds = ['1', '2', '3', '4', '5'];
      builds.sort((a, b) => parseInt(a) - parseInt(b));
      let res = ui.table(builds, (build, idx) => {
          return [ui.span([build]),
            ui.tableFromMap({
              'Build': df.col(build)?.tags?.get('build'),              
              'Instance': df.col(build)?.tags?.get('instance'),
              'Commit': df.col(build)?.tags?.get('build_commit'),
              'Build date': df.col(build)?.tags?.get('build_date')})
          ];
        }
      );
      return [new Verdict(Priority.INFO, 'Builds info', res)];
    } catch (x) {
      return [new Verdict(Priority.ERROR, 'Unhandled exception', ui.div([
        ui.span(['Failed to get build infos']),
        ui.span([x])
      ]))];
    }
  }
  ignoredTests(df: DG.DataFrame): Verdict[] {
    let verdicts: Verdict[] = [];
    let testColumn: DG.Column = df.col('test')!;
    let ignoreColumn: DG.Column = df.col('ignore?')!;
    let widgets: DG.Widget[] = [];
    for (let i = 0; i < df.rowCount; i++) {
      if (ignoreColumn.get(i)) {
        widgets.push(DG.Widget.fromRoot(ui.div([testColumn.getString(i)],{style: {marginBottom: '5px'}})));
        let test = testColumn.getString(i);
        widgets[widgets.length - 1].root.onclick = (ev: MouseEvent) => {
          let df = grok.shell.tv.grid.dataFrame;
          let rows = df.rows;
          let idx = -1;
          for(let j = 0; j < df.rowCount && idx === -1; j++)
            if (rows.get(j).get('test') ===test)
              idx = j;
          
          if(idx >= 0)
            grok.shell.tv.grid.scrollToCell('test', idx);
          
        }
      }
    }
    for (let i = 0; i < widgets.length; i++)
      verdicts.push(new Verdict(Priority.INFO, `Ignored tests - ${widgets.length} items`, widgets[i].root));
    return verdicts;
  }

  performanceDowngrades(df: DG.DataFrame): Verdict[] {
    function getLastDuration(df: DG.DataFrame, index: number) {
      let duration = 0;
      for(let i = 1; i <= 5; i++) {
        if (!df.getCol(`${i} avg(duration)`).isNone(index)) {
          duration = df.get(`${i} avg(duration)`, index);
          break;
        }
      }
      return duration;
    }

    try {
      let verdicts: Verdict[] = [];
      let verdictGrid: DG.Grid;
      let performanceObjects= [];
      for (let i = 0; i < df.rowCount; i++) {
        if (df.col('has_suspicious')!.get(i)) {
          let lastDuration: number = getLastDuration(df, i) ?? 0;
          let minDuration: number = df.col('min')!.get(i) ?? 0;
          let percentage: number = Number(((lastDuration / minDuration)).toFixed(2));
          if (percentage > 0.20)
            performanceObjects.push({ 'increased(times)': percentage, test: df.col('test')!.getString(i), from: minDuration, to:  lastDuration});
        }
      }
      let verdictDF = DG.DataFrame.fromObjects(performanceObjects);
     
      if (verdictDF) {
        verdictGrid = (DG.Grid.create(verdictDF) )
        let col = verdictGrid.col('increased(times)');
        if(col)
          col.width = 50;
        let button = ui.button('OPEN IN THE WORKSPACE', ()=> grok.shell.addTableView(verdictDF?.clone()));
        button.classList.add('performance-downgrade-button');
        verdicts.push(new Verdict(Priority.ERROR, 'Performance downgrade', ui.divV([button, verdictGrid.root])));
      }
      return verdicts;
    } catch (x) {
      return [new Verdict(Priority.ERROR, 'Unhandled exception', ui.div([
        ui.span(['Failed to gather performance downgrades']),
        ui.span([x])
      ]))];
    }
  }
  
  async readManualTickets(): Promise<string[]> {
    return await _package.files.exists(MANUAL_TICKETS_FILE) ?
      JSON.parse(await _package.files.readAsText(MANUAL_TICKETS_FILE)) : [];
  }

  async addManualTickets(refs: IssueRef[]): Promise<void> {
    const stored = await this.readManualTickets();
    const added = refs.map((r) => r.key).filter((key) => !stored.includes(key));
    if (added.length > 0)
      await _package.files.writeAsText(MANUAL_TICKETS_FILE, JSON.stringify([...stored, ...added]));
  }

  async verdictsForTickets(): Promise<Verdict[]> {
    try {
      this.watch(findIssueRefs((await this.readManualTickets()).join(' ')));
      const found = await this.githubIssues.findIssues([...this.ticketRefs.values()]);
      const issues = new Map([...found.values()].map((i) => [i.key, i]));
      const groups = new Map<string, GitHubIssue[]>();
      for (const issue of issues.values()) {
        const priority = verdictPriority(issue, _properties['Platform version']);
        groups.set(priority, [...(groups.get(priority) ?? []), issue]);
      }
      return [...groups].map(([priority, group]) => {
        const grid = issuesGrid(group, `GitHub issues [${priority}]`);
        grid.props.rowHeight = 32;
        grid.autoSize(1000, 500);
        return new Verdict(priority, `GitHub issues [${priority}] - ${group.length} items`, grid.root);
      });
    } catch (e) {
      return [new Verdict(Priority.ERROR, 'GitHub issues', ui.div([
        ui.span(['Failed to read GitHub issues: ']),
        ui.span([e instanceof Error ? e.message : `${e}`]),
      ]))];
    }
  }

  composeResultList(verdicts: Verdict[]): HTMLElement {
    let res: DG.Accordion = ui.accordion('');
    let lastAccordion: HTMLDivElement | null;
    for (let i = 0; i < verdicts.length; i++) {
      if (i == 0 || verdicts[i - 1].priority != verdicts[i].priority || verdicts[i - 1].category != verdicts[i].category) {
        let accordionToShow = ui.div([]);
        let pane = res.addPane(verdicts[i].category, () => accordionToShow);
        let icon = getIconForVerdict(verdicts[i].priority);
        icon.style.marginRight = '4px';
        pane.root.firstChild!.insertBefore(icon, pane.root.firstChild!.firstChild!);
        lastAccordion = accordionToShow;
      }
      lastAccordion!.appendChild(ui.div([verdicts[i].widget]));
    }
    res.addPane('Actions', () => ui.actionLink('Watch a GitHub issue', () => {
      const dialog = ui.dialog('Watch issue');
      const ticketInput = ui.input.textArea('issue', {tooltipText: 'core#456, public#12 or an issue URL'});
      dialog.add(ticketInput.root);
      dialog.onOK(async () => {
        const refs = findIssueRefs(ticketInput.value);
        if (refs.length > 0)
          await this.addManualTickets(refs);
        else
          grok.shell.warning(`Couldn't recognize an issue in ${ticketInput.value}`);
      });
      dialog.show();
    }));
    return res.root;
  }

  close(): void {
    TestDashboardWidget.isOpen = false;
    super.close();
  }
}