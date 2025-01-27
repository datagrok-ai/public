import * as DG from 'datagrok-api/dg';
import { IndexPredicate } from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

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

class Verdict {
  priority: string;
  widget: HTMLElement;
  constructor(priority: string, widget: HTMLElement) {
    this.priority = priority;
    this.widget = widget;
  }
}

function getIconForVerdict(priority: string): HTMLElement {
  let icon: HTMLElement | undefined;
  if (priority == Priority.BLOCKER) {
    icon = ui.iconFA('exclamation-triangle');
    icon.style.color = DG.Color.toRgb(DG.Color.red);
  } else if (priority == Priority.ERROR) {
    icon = ui.iconFA('exclamation-circle');
    icon.style.color = DG.Color.toRgb(DG.Color.red);
  } else if (priority == Priority.INFO) {
    icon = ui.iconFA('info-circle');
    icon.style.color = DG.Color.toRgb(DG.Color.blue);
  } if (priority == Priority.RESOLVED) {
    icon = ui.iconFA('info-circle');
    icon.style.color = DG.Color.toRgb(DG.Color.green);
  }
  return icon!;
}

export class TestDashboardWidget extends DG.JsViewer {
  static isOpen: boolean = false;
  constructor() {
    super();
  }

  tickets: Set<string> = new Set();

  onFrameAttached(dataFrame: DG.DataFrame): void {
    this.root.appendChild(ui.wait(async () => {
      let tableNames: RegExp[] = [/Benchmark.*Dashboard/, /Test.*Track.*Dashboard/];
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
      let verdicts: Verdict[] = [];
      verdicts.push(...this.verdictsOnCorrectness(dataFrame));
      verdicts.push(...this.verdictsOnPerformance());
      verdicts.push(...this.verdictsOnTestTrack());
      verdicts.push(...(await this.unaddressedTests(dataFrame)));
      verdicts.push(...(await this.verdictsForTickets()));
      let status = Priority.getMaxPriority(verdicts.map((v) => v.priority));
      if (status == Priority.BLOCKER)
        d.append(ui.h1('Release is not okay to publish ❌❌❌', { style: { color: DG.Color.toRgb(DG.Color.red) } }));
      else
        d.append(ui.h1('Green light! 🚀🚀🚀', { style: { color: DG.Color.toRgb(DG.Color.darkGreen) } }));
      
      verdicts.sort((a, b) => Priority.codeMapping[b.priority] - Priority.codeMapping[a.priority]);

      let list = ui.list(verdicts.map((verdict) => {
        return ui.div([getIconForVerdict(verdict.priority), verdict.widget]);
      }));
      d.append(list);
      d.style.maxWidth = '300px';
      return d;
    }));
  }

  verdictsOnCorrectness(df: DG.DataFrame): Verdict[] {
    let verdicts: Verdict[] = [];
    if (df.col('jira') != null)
      this.jiraTickets(df.col('jira')!);
    verdicts.push(...this.majorPackages(df));
    return verdicts;
  }
  verdictsOnPerformance(): Verdict[] {
    let verdicts: Verdict[] = [];
    let df = grok.shell.tables.find((df) => df.name.match(/Benchmark.*Dashboard/));
    if (df == null)
      return verdicts;
    verdicts.push(...this.performanceDowngrades(df));
    return verdicts;
  }
  verdictsOnTestTrack(): Verdict[] {
    let verdicts: Verdict[] = [];
    let df = grok.shell.tables.find((df) => df.name.match(/Test.*Track.*Dashboard/));
    if (df == null)
      return verdicts;
    if (df.col('jira') != null)
      this.jiraTickets(df.col('jira')!);
    verdicts.push(...this.testTrackDowngrades(df));
    return verdicts;
  }
  jiraTickets(jiraCol: DG.Column<string>): Verdict[] {
    try {
      for (var i = 0; i < jiraCol.length; i++)
        if ((jiraCol.getString(i)?.length ?? 0) > 0)
          jiraCol.getString(i).split(',').forEach((ticket, _) => {
            let match = ticket.match(/GROK\-\d+/);
            console.log(ticket);
            console.log(match);
            if (match != null)
              this.tickets.add(match[0]);
          });
      return [];
    } catch (x) {
      return [new Verdict(Priority.ERROR, ui.div([
        ui.span(['Failed to check jira tickets']),
        ui.span([x])
      ]))];
    }
  }
  majorPackages(df: DG.DataFrame): Verdict[] {
    try {
      let packageList = ['datlas', 'ddt', 'ddtx', 'dinq', 'ApiTests', 'ApiSamples', 'Bio', 'Chem']
      let testColumn: DG.Column = df.col('test')!;
      let failingColumn: DG.Column<boolean> = df.col('failing')!;
      let brokenPackages: Map<string, string[]> = new Map();
      for (var i = 0; i < packageList.length; i++)
        for (var row = 0; row < testColumn.length; row++)
          if (testColumn.getString(row)?.startsWith(packageList[i]))
            if (failingColumn.get(row)) {
              if (!brokenPackages.has(packageList[i]))
                brokenPackages.set(packageList[i], []);
              brokenPackages.get(packageList[i])!.push(testColumn.getString(row).slice(packageList[i].length));
            }

      let verdicts: Verdict[] = [];
      for (const packageName of brokenPackages.keys())
        verdicts.push(new Verdict(Priority.BLOCKER, ui.span(['package ' + packageName + ' has broken tests'])));
      return verdicts;
    } catch (x) {
      return [new Verdict(Priority.BLOCKER, ui.div([
        ui.span(['Failed to verify major packages']),
        ui.span([x])
      ]))];
    }
  }
  async unaddressedTests(df: DG.DataFrame): Promise<Verdict[]> {
    try {
      let verdicts: Verdict[] = [];
      let failingColumn: DG.Column<boolean> = df.col('failing')!;
      let jiraCol = df.col('jira')!;
      let testColumn: DG.Column = df.col('test')!;
      let ownerColumn: DG.Column = df.col('owner')!;

      let unaddressedTests: { [owner: string]: string[] } = {};
      let predicate: IndexPredicate = (row) => {
        return failingColumn.get(row)! && (!jiraCol.getString(row) || jiraCol.getString(row).trim() === '');
      };

      for (let i = 0; i < df.rowCount; i++) {
        if (predicate(i)) {
          let owner = ownerColumn.getString(i) || 'unknown';
          if (!unaddressedTests[owner])
            unaddressedTests[owner] = [];
          unaddressedTests[owner].push(testColumn.getString(i));
        }
      }

      if (Object.keys(unaddressedTests).length > 0) {
        let summary = ui.div([ui.span(['Owners with failing tests:\n'])]);
        for (const [owner, tests] of Object.entries(unaddressedTests)) {
          const login = owner.match(/\w+@datagrok.ai/) ? owner.match(/\w+@datagrok.ai/)?.[0].split('@')[0] : owner
          let userIcon : HTMLElement | undefined;
          if (login == null) {
            userIcon = ui.span([owner]);
          } else {
            const user = await grok.dapi.users.filter(`login="${login}"`).first();
            if (user == undefined)
              userIcon = ui.span([login]);
            else
              userIcon = ui.render(user.toMarkup());
          }
          summary.appendChild(ui.div([userIcon, ui.span([`: ${tests.length} tests`])]));
        }
        verdicts.push(new Verdict(Priority.INFO, ui.div([
          ui.span([summary]),
          ui.button('Copy to Clipboard', () => {
            let slackMessage = '';
            for (const [owner, tests] of Object.entries(unaddressedTests)) {
              const slackOwner = owner.match(/\w+@datagrok.ai/) ? '@' + owner.match(/\w+@datagrok.ai/)?.[0].split('@')[0] : owner;
              slackMessage += `${slackOwner}: ` + tests.length + '\n\n';
            }
      
            navigator.clipboard.writeText(slackMessage).then(() => {
              grok.shell.info('Summary copied to clipboard');
            });
          })
        ])));
      }
      return verdicts;
    } catch (x) {
      return [new Verdict(Priority.ERROR, ui.div([
        ui.span(['Failed to gather unaddressed tests']),
        ui.span([x])
      ]))];
    }
  }
  performanceDowngrades(df: DG.DataFrame): Verdict[] {
    try {
      let verdicts: Verdict[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        if (df.col('has_suspicious')!.get(i)) {
          let lastDuration: number = df.col('1 avg(duration)')!.get(i) ?? 0;
          let minDuration: number = df.col('min')!.get(i) ?? 0;
          let maxDuration: number = df.col('max')!.get(i) ?? 0;
          if (lastDuration - minDuration > maxDuration - lastDuration)
            verdicts.push(new Verdict(Priority.ERROR, ui.span([df.col('test')!.getString(i) + ` is down ${(maxDuration / minDuration - 1) * 100}% in performance`])));
        }
      }
      return verdicts;
    } catch (x) {
      return [new Verdict(Priority.ERROR, ui.div([
        ui.span(['Failed to gather performance downgrades']),
        ui.span([x])
      ]))];
    }
  }
  testTrackDowngrades(df: DG.DataFrame): Verdict[] {
    try {
      let verdicts: Verdict[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        if ((df.col('1')?.getString(i) ?? 'passed') != 'passed') {
          verdicts.push(new Verdict(Priority.ERROR, ui.span([df.col('test')!.getString(i) + ` is failed`])));
        }
      }
      return verdicts;
    } catch (x) {
      return [new Verdict(Priority.ERROR, ui.div([
        ui.span(['Failed to gather test track downgrades']),
        ui.span([x])
      ]))];
    }
  }
  async verdictsForTickets(): Promise<Verdict[]> {
    try {
      let verdicts: Verdict[] = [];
      let jiraCol: DG.Column<string> = DG.Column.fromStrings('jira', [...this.tickets]);
      let severityCol: DG.Column<string> = await grok.functions.call('JiraConnect:getJiraField', {
        ticketColumn: jiraCol,
        field: 'priority:name'
      });
      let statusCol: DG.Column<string> = await grok.functions.call('JiraConnect:getJiraField', {
        ticketColumn: jiraCol,
        field: 'status:name'
      });
      let issueTypeCol: DG.Column<string> = await grok.functions.call('JiraConnect:getJiraField', {
        ticketColumn: jiraCol,
        field: 'issueType:name'
      });
      for (var i = 0; i < jiraCol.length; i++) {
        let priority: string = Priority.INFO;
        if (severityCol.getString(i).startsWith('High'))
          priority = Priority.BLOCKER;
        else if (statusCol.getString(i) == 'Done')
          priority = Priority.RESOLVED;
        else if (issueTypeCol.getString(i) == 'Bug') {
          if (severityCol.getString(i).startsWith('Low'))
            priority = Priority.INFO;
          else
            priority = Priority.ERROR;
        }
        let ticket: string = jiraCol.getString(i);
        verdicts.push(new Verdict(priority, ui.renderCard(DG.SemanticValue.fromValueType(ticket, 'JIRA Ticket'))));//., false)));
      }
      return verdicts;
    } catch (x) {
      return [new Verdict(Priority.ERROR, ui.div([
        ui.span(['Failed to gather jira tickets info']),
        ui.span([x])
      ]))];
    }
  }
  close(): void {
    TestDashboardWidget.isOpen = false;
    super.close();
  }
}