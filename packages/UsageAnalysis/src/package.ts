import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UsageWidget} from './widgets/usage-widget';
import {PackageUsageWidget} from './widgets/package-usage-widget';
import '../css/usage_analysis.css';
import '../css/test_track.css';
import {ViewHandler} from './view-handler';
import {TestTrack} from './test-track/app';
import {ReportsWidget} from './widgets/reports-widget';
import {ReportingApp} from './reporting/reporting_app';
import {TestAnalysisManager} from './test-analysis/test-analysis-manager';
import {getDate} from './utils';
import {ServiceLogsApp} from './service_logs/service_logs';
import {TestGridCellHandler} from './handlers/test-grid-cell-handler';
import {initTestStickyMeta} from './test-analysis/sticky-meta-initialization';
import {TestDashboardWidget} from './viewers/ua-test-dashboard-viewer';

export const _package = new DG.Package();
export let _properties: any;
export * from './package.g';

let AVAILABLE_SERVICES: string[];


const highestPriorityIndex = 3;
const priorityLevels : Record<string, number> = {
  'Lowest': 0,
  'Low': 1,
  'Medium': 2,
  'Blocker': highestPriorityIndex,
};
export class PackageFunctions {
  @grok.decorators.init()
  static _initUA(): void {
    _properties = _package.settings;
    DG.ObjectHandler.register(new TestGridCellHandler());
    initTestStickyMeta();
  }


  @grok.decorators.func({
    'meta': {
      'url': '/tests/list',
    },
  })
  static async TestsList(): Promise<DG.DataFrame| undefined> {
    const pacakageTests = await TestAnalysisManager.collectPackageTests();
    const packageTestsListMapped = pacakageTests.map((elem) => {
      return {'name': 'test-package ' + elem.packageName + ': ' + elem.test.category + ': ' + elem.test.name};
    });
    const manualTest = await TestAnalysisManager.collectManualTestNames();
    const manualTestsListMapped = manualTest.map((elem) => {
      return {'name': 'test-manual ' + elem};
    });
    const resultTestsList = manualTestsListMapped.concat(packageTestsListMapped);
    return DG.DataFrame.fromObjects(resultTestsList);
  }


  @grok.decorators.func({
    'meta': {
      'url': '/tests/joinedlist',
    },
  })
  static async TestsListJoined(): Promise<DG.DataFrame| undefined> {
    const pacakageTests = await TestAnalysisManager.collectPackageTests();
    const packageTestsListMapped = pacakageTests.map((elem) => {
      return {'type': 'package ', 'test': elem.packageName + ': ' + elem.test.category + ': ' + elem.test.name};
    });
    const manualTest = await TestAnalysisManager.collectManualTestNames();
    const manualTestsListMapped = manualTest.map((elem) => {
      return {'type': 'manual ', 'test': 'Test Track: ' + elem};
    });
    const resultTestsList = DG.DataFrame.fromObjects(manualTestsListMapped.concat(packageTestsListMapped));

    // const builds: DG.DataFrame = await grok.functions.call('UsageAnalysis:Builds');
    // const id = builds.get('name', 0);

    // const tests = await grok.functions.call('UsageAnalysis:getTestStatusesAcordingDF', { 'buildId': id, 'testslist': resultTestsList });
    grok.shell.addTableView(resultTestsList!);
    return resultTestsList;
  }

  @grok.decorators.func()
  static async TestAnalysisReportForCurrentDay(
    @grok.decorators.param({'type': 'datetime'}) date: any) : Promise<DG.DataFrame> {
    const tests = await TestAnalysisManager.collectPackageTests();
    const testsListMapped = tests.map((elem) => {
      return {'name': 'test-package ' + elem.packageName + ': ' + elem.test.category + ': ' + elem.test.name};
    });
    const testRuns = await grok.functions.call('UsageAnalysis:getServerStartTestResults', {'date': getDate(new Date(date)), 'testslist': DG.DataFrame.fromObjects(testsListMapped)});
    return testRuns;
  }

  @grok.decorators.app({
    'url': '/',
    'browsePath': 'Admin',
    'name': 'Usage Analysis',
  })
  static usageAnalysisApp(
    @grok.decorators.param({'options': {'optional': true, 'meta.url': true}}) path?: string,
    @grok.decorators.param({'options': {'optional': true}}) date?: string,
    @grok.decorators.param({'options': {'optional': true}}) groups?: string,
    @grok.decorators.param({'options': {'optional': true}}) packages?: string,
    @grok.decorators.param({'options': {'optional': true}}) tags?: string,
    @grok.decorators.param({'options': {'optional': true}}) categories?: string,
    @grok.decorators.param({'options': {'optional': true}}) projects?: string): DG.ViewBase | null {
    const handler = new ViewHandler();
    handler.view.parentCall = grok.functions.getCurrentCall();
    handler.init(date, groups, packages, tags, categories, projects, path);
    return handler.view;
  }

  @grok.decorators.func({ meta: { vectorFunc: 'true' } })
  static async getTicketsVerdict(
      @grok.decorators.param({ type: 'column<string>' }) ticketColumn: DG.Column,
      @grok.decorators.param({ type: 'column<string>' }) resultColumn: DG.Column
  ): Promise<void> {

    const n = ticketColumn.length;
    if (n != resultColumn.length)
      throw new Error('Ticket column and result column should have the same length.');

    const ticketRegex = /GROK-\d*/g;
    const issueIdToIdx = new Map<number, Set<string>>();
    const issueIdsOrKeys = new Set<string>();

    for (let i = 0; i < n; i++) {
      const cellValue = ticketColumn.get(i);
      if (!cellValue)
        continue;

      const matches = cellValue.matchAll(ticketRegex);

      for (const match of matches) {
        const ticket = match[0];
        if (issueIdToIdx.has(i))
          issueIdToIdx.get(i)!.add(ticket);
        else
          issueIdToIdx.set(i, new Set<string>([ticket]));
        issueIdsOrKeys.add(ticket);
      }
    }
    const { issues } = await grok.functions.call('JiraConnect:getJiraTicketsBulk',
        { 'issueIdsOrKeys': [...issueIdsOrKeys], 'fields': ['status', 'priority']});
    const resultIssuesInfo = new Map(issues.map((issue: any) => [issue.key, issue]));

    for (let i = 0; i < n; i++) {
      if (!issueIdToIdx.has(i))
        continue;

      const resultStatuses: { status: string; severity: string }[] = [...issueIdToIdx.get(i)!].map((k) => {
        const issueData: any = resultIssuesInfo.get(k);
        return issueData ? { status: issueData.fields.status.name, severity: issueData.fields.priority.name }
            : null;
      }).filter((s): s is { status: string; severity: string } => s !== null);

      if (resultStatuses.length === 0)
        continue;

      let verdict: string;
      if (resultStatuses.some((e) => e.status !== 'Done')) {
        const priority = this.getHighestPriorityLevel(resultStatuses as any);
        if (resultStatuses.some((e) => e.status === 'Done'))
          verdict = `Partially Fixed (${priority})`;
        else
          verdict = `Wasn't Fixed (${priority})`;
      }
      else
        verdict = 'Fixed';

      resultColumn.set(i, verdict);
    }
  }

  // @grok.decorators.func({ meta: { vectorFunc: 'true' } })
  // static async getTicketsVerdict(
  //     @grok.decorators.param({ type: 'column<string>' }) ticketColumn: DG.Column,
  //     @grok.decorators.param({ type: 'column<string>' }) resultColumn: DG.Column,
  //     @grok.decorators.param({ type: 'object' }) progress: DG.ProgressIndicator,
  // ): Promise<void> {
  //
  //   const ticketRegex = /GROK-\d*/g;
  //   const ticketsMap = new Map<string, { status: string; severity: string }>();
  //   const n = ticketColumn.length;
  //
  //   for (let i = 0; i < n; i++) {
  //     const cellValue = ticketColumn.get(i);
  //     if (!cellValue)
  //       continue;
  //
  //     const matches = cellValue.matchAll(ticketRegex);
  //     const resultStatuses: { status: string; severity: string }[] = [];
  //
  //     for (const match of matches) {
  //       const ticket = match[0];
  //       let info = ticketsMap.get(ticket);
  //
  //       if (!info) {
  //         const status = await grok.functions.call('JiraConnect:issueData', { issueKey: ticket });
  //         if (progress.canceled)
  //           return;
  //         if (!status)
  //           continue;
  //
  //         info = { status: status.fields.status.name, severity: status.fields.priority.name };
  //         ticketsMap.set(ticket, info);
  //       }
  //
  //       resultStatuses.push(info);
  //     }
  //
  //     if (resultStatuses.length === 0)
  //       continue;
  //
  //     let verdict: string;
  //     if (resultStatuses.some((e) => e.status !== 'Done')) {
  //       const priority = this.getHighestPriorityLevel(resultStatuses as any);
  //       if (resultStatuses.some((e) => e.status === 'Done'))
  //         verdict = `Partially Fixed (${priority})`;
  //       else
  //         verdict = `Wasn't Fixed (${priority})`;
  //     } else {
  //       verdict = 'Fixed';
  //     }
  //
  //     resultColumn.set(i, verdict);
  //   }
  // }


  // @grok.decorators.func({meta: {vectorFunc: 'true'}})
  // static async getTicketsVerdict(@grok.decorators.param({type: 'column<string>'}) ticketColumn: DG.Column): Promise<DG.Column | undefined> {
  //   const ticketsMap = new Map<string, {status: string, severity: string}>();
  //   const ticketRegex = /GROK-\d*/g;
  //   const ticketColumnList = ticketColumn.toList();
  //   for (let i = 0; i < ticketColumnList.length; i++) {
  //     if (ticketColumnList[i]) {
  //       const matches = ticketColumnList[i].matchAll(ticketRegex);
  //       for (const match of matches) {
  //         const ticket: string = match[0];
  //         const status = (await grok.functions.call('JiraConnect:issueData', {issueKey: ticket}));
  //         if (status && !ticketsMap.has(ticket))
  //           ticketsMap.set(ticket, {status: status.fields.status.name, severity: status.fields.priority.name});
  //       }
  //     }
  //   }
  //
  //   const resultCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, `${ticketColumn} ticket verdict`, ticketColumn.length);
  //   for (let i = 0; i < ticketColumnList.length; i++) {
  //     if (ticketColumnList[i]) {
  //       const matches = ticketColumnList[i].matchAll(ticketRegex);
  //       const resultStatuses = [];
  //       let status = undefined;
  //       for (const match of matches) {
  //         if (ticketsMap.has(match[0]))
  //           resultStatuses.push(ticketsMap.get(match[0]));
  //       }
  //
  //       if (resultStatuses.length > 0) {
  //         if (resultStatuses.some((e)=> e?.status !== 'Done')) {
  //           const priority = this.getHighestPriorityLevel(resultStatuses as any);
  //           if (resultStatuses.some((e)=> e?.status === 'Done'))
  //             status = `Partially Fixed (${priority})`;
  //           else
  //             status = `Wasn\'t Fixed (${priority})`;
  //         } else
  //           status = 'Fixed';
  //         resultCol.set(i, status);
  //       }
  //     }
  //   }
  //
  //   return resultCol;
  // }

  static getHighestPriorityLevel(ticketsMap: {status: string, severity: string}[]) {
    let highestPriority = 'Lowest';
    let highestPriorityIndex = 0;
    for (const [key, value] of ticketsMap.entries()) {
      const index = priorityLevels[value.severity];
      if (index > highestPriorityIndex) {
        highestPriority = (value.severity);
        highestPriorityIndex = index;
      }
      if (index >= highestPriorityIndex)
        break;
    }
    return highestPriority;
  }

  @grok.decorators.app({
    'url': '/tests/manager',
    'browsePath': 'Admin',
    'name': 'Test Track',
  })
  static testTrackApp(): void {
    if (!grok.shell.dockManager.findNode(TestTrack.getInstance().root))
      TestTrack.getInstance().init();
    else
      TestTrack.getInstance().reopen();
  }


  @grok.decorators.app({
    'url': '/reports',
    'browsePath': 'Admin',
    'name': 'Reports',
  })
  static async reportsApp(
    @grok.decorators.param({'options': {'optional': true, 'meta.url': true}}) path?: string): Promise<DG.ViewBase> {
    const parent = grok.functions.getCurrentCall();
    const app = new ReportingApp(parent);
    app.init(path).catch((e) => console.log(e));
    return app.view;
  }


  @grok.decorators.app({
    'url': '/service-logs',
    'browsePath': 'Admin',
    'name': 'Service Logs',
  })
  static serviceLogsApp(
    @grok.decorators.param({'options': {'optional': true, 'meta.url': true}}) path?: string,
    @grok.decorators.param({'type': 'map', 'options': {'optional': true}}) params?: any,
    @grok.decorators.param({'type': 'int', 'options': {'optional': true}}) limit?: number): DG.ViewBase {
    if (path && path.startsWith('/'))
      path = path.slice(1);
    const view = DG.View.fromViewAsync(async () => {
      const currentCall = grok.functions.getCurrentCall();
      AVAILABLE_SERVICES ??= await grok.dapi.docker.getAvailableServices();
      const app = new ServiceLogsApp(currentCall, AVAILABLE_SERVICES, path, limit);
      if (AVAILABLE_SERVICES.length > 0)
        app.getLogs().then((_) => {});
      //@ts-ignore
      return app as DG.View;
    });
    view.name = ServiceLogsApp.APP_NAME;
    return view;
  }


  @grok.decorators.appTreeBrowser({app: 'Service Logs'})
  static async serviceLogsAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
    loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
    const loaderItem = treeNode.item(loaderDiv);
    try {
      AVAILABLE_SERVICES ??= await grok.dapi.docker.getAvailableServices();
      const services = ['datagrok', 'grok_pipe', 'grok-pipe', 'rabbitmq', 'jkg', 'spawner', 'grok-connect', 'grok_connect'];

      const getItem = (service: string) => {
        let icon: HTMLElement;
        if (services.some((s) => service.includes(s)))
          icon = ui.image('/images/entities/grok.png', 10, 10);

        else
          icon = ui.iconFA('server');
        icon.style.marginRight = '3px';
        const span = ui.span([icon, ' ', service]);
        span.style.display = 'flex';
        span.style.alignItems = 'center';
        return span;
      };

      let currentView: DG.View;
      for (const service of AVAILABLE_SERVICES) {
        const node = treeNode.item(getItem(service));
        node.onSelected.subscribe(async (_) => {
          currentView?.close();
          currentView= DG.View.fromViewAsync(async () => {
            const app = new ServiceLogsApp(grok.functions.getCurrentCall(), AVAILABLE_SERVICES, service, ServiceLogsApp.DEFAULT_LIMIT, true);
            await app.getLogs();
            //@ts-ignore
            return app as DG.View;
          });
          currentView.name = service;
          grok.shell.addPreview(currentView);
        });
      }
    } finally {
      loaderItem.remove();
    }
  }

  @grok.decorators.appTreeBrowser({app: 'Reports'})
  static async reportsAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    await treeNode.group('Reports', null, false).loadSources(grok.dapi.reports.by(10));
    await treeNode.group('Rules', null, false).loadSources(grok.dapi.rules.include('actions,actions.assignee').by(10));
  }

  @grok.decorators.dashboard({
    'meta': {
      'canView': 'Developers,Administrators',
    },
    'name': 'Usage',
    'test': 'usageWidget()',
  })
  static usageWidget(): DG.Widget {
    return new UsageWidget();
  }

  @grok.decorators.dashboard({
    'meta': {'canView': 'Developers,Administrators'},
    'name': 'Reports',
    'test': 'reportsWidget()',
  })
  static reportsWidget(): DG.Widget {
    return new ReportsWidget();
  }


  @grok.decorators.func({})
  static packageUsageWidget(
    @grok.decorators.param({'type': 'object'}) pckg: DG.Package): DG.Widget {
    return new PackageUsageWidget(pckg);
  }


  @grok.decorators.func({
    'meta': {showInGallery: 'false', role: 'viewer'},
    'outputs': [
      {
        type: 'viewer',
        name: 'result',
      },
    ],
  })
  static testDashboardsViewer(): TestDashboardWidget {
    return new TestDashboardWidget();
  }

  @grok.decorators.autostart()
  static describeCurrentObj(): void {
    grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      const ent = acc.context;
      if (ent != null && ent.constructor.name === 'Package') {
        const pane = acc.addPane('Usage', () => ui.wait(async () => {
          let widget: HTMLElement;
          try {
            widget = PackageFunctions.packageUsageWidget(ent).root;
          } catch (e) {
            widget = ui.divText('Error on loading', {style: {color: 'var(--failure)'}});
          }
          return widget;
        }));
        const UAlink = ui.link('', async () => {
          grok.shell.v.path = `/apps/UsageAnalysis/Packages?date=this%20week&users=${(await grok.dapi.groups.getGroupsLookup('All users'))[0].id}&packages=${ent.name}`;
          grok.functions.eval('UsageAnalysis:usageAnalysisApp()');
        }, 'Open Usage Analysis');
        UAlink.style.marginLeft = '3px';
        const header = pane.root.querySelector('.d4-accordion-pane-header') as HTMLElement;
        header.appendChild(UAlink);
      }
    });
  }


  @grok.decorators.func({
    'name': 'Create JIRA ticket',
    'description': 'Creates JIRA ticket using current error log',
  })
  static createJiraTicket() {
    grok.data.query('JiraCreateIssue', {
      'createRequest': JSON.stringify({
        'fields': {
          'project': {
            'key': 'GROK',
          },
          'summary': 'test',
          'description': '',
          'issuetype': {
            'name': 'Bug',
          },
        },
      }),
      'updateHistory': false,
    }).then((t) => {
      grok.shell.info('Created');
      console.log(t);
    });
  }
}
