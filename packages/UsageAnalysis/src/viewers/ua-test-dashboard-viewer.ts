import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

class Priority {
  static BLOCKER: string = 'BLOCKER';
  static INFO: string = 'INFO';
  static ERROR: string = 'ERROR';

  static codeMapping: { [key: string]: number } = {
    [Priority.BLOCKER]: 3,
    [Priority.ERROR]: 2,
    [Priority.INFO]: 1
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

export class TestDashboardWidget extends DG.ViewBase {
  static isOpen: boolean = false;
  constructor() {
    super();
    this.name = 'Test Analysis';
    let verdicts: Verdict[] = [];
    verdicts.push(...this.verdictsOnCorrectness());
    let status = Priority.getMaxPriority(verdicts.map((v) => v.priority));
    if (status == Priority.BLOCKER)
      this.root.append(ui.h1('Release is not okay to publish ❌❌❌', { style: { color: DG.Color.toRgb(DG.Color.red) } }));
    else
      this.root.append(ui.h1('Green light! 🚀🚀🚀', { style: { color: DG.Color.toRgb(DG.Color.darkGreen) } }));
    let list = ui.list(verdicts.map((verdict) => verdict.widget));
    this.root.append(list);
    this.root.style.maxWidth = '300px';
    TestDashboardWidget.isOpen = true;
  }

  verdictsOnCorrectness(): Verdict[] {
    let verdicts: Verdict[] = [];
    let df = Array.from(grok.shell.tableViews).find((tv) => tv.name.search(/'Tests.*Dashboard'/))!.dataFrame;
    verdicts.push(...this.jiraTickets(df));
    verdicts.push(...this.majorPackages(df));
    return verdicts;
  }
  jiraTickets(df: DG.DataFrame): Verdict[] {
    let verdicts: Verdict[] = [];
    let tickets: Set<string> = new Set();
    let jiraCol = df.col('jira')!;
    let severityCol = df.col('severity')!;
    for (var i = 0; i < df.rowCount; i++)
      if ((severityCol.getString(i).length ?? 0) > 0)
        jiraCol.getString(i).split(',').forEach((ticket, _) => tickets.add(ticket));
    for (const ticket of tickets)
      verdicts.push(new Verdict(Priority.INFO, ui.link(ticket, 'https://reddata.atlassian.net/jira/software/c/projects/GROK/issues/' + ticket)));
    return verdicts;
  }
  majorPackages(df: DG.DataFrame): Verdict[] {
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
  }
  close(): void {
    TestDashboardWidget.isOpen = false;
    super.close();
  }
}


export function testDashboardWidget() {
  if (TestDashboardWidget.isOpen)
    return;
  grok.shell.addView(new TestDashboardWidget(), DG.DOCK_TYPE.LEFT, 300);
}