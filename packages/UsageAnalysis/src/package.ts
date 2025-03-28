import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UsageWidget} from './widgets/usage-widget';
import {PackageUsageWidget} from './widgets/package-usage-widget';
import '../css/usage_analysis.css';
import '../css/test_track.css';
import {ViewHandler} from './view-handler';
import {TestTrack} from './test-track/app';
import {ReportsWidget} from "./widgets/reports-widget";
import {ReportingApp} from "./reporting/reporting_app";
import {TestAnalysisManager} from './test-analysis/test-analysis-manager';
import {getDate} from './utils';
import {ServiceLogsApp} from "./service_logs/service_logs";
import {TestGridCellHandler} from './handlers/test-grid-cell-handler';
import {initTestStickyMeta} from './test-analysis/sticky-meta-initialization';
import {TestDashboardWidget} from './viewers/ua-test-dashboard-viewer';

export const _package = new DG.Package();
export let _properties: any;


//tags: init
export function _initUA(): void {
  _properties = _package.settings;
  DG.ObjectHandler.register(new TestGridCellHandler());
  initTestStickyMeta();
}

//name: TestsList 
//meta.url: /tests/list
//output: dataframe df
export async function TestsList(): Promise<DG.DataFrame| undefined> { 
  const pacakageTests = await TestAnalysisManager.collectPackageTests();
  const packageTestsListMapped = pacakageTests.map((elem) => {
    return { 'name':  "test-package " + elem.packageName + ": " + elem.test.category + ": " + elem.test.name };
  });
  const manualTest = await TestAnalysisManager.collectManualTestNames();
  const manualTestsListMapped = manualTest.map((elem) => {
    return { 'name':  "test-manual " + elem };
  });
  const resultTestsList = manualTestsListMapped.concat(packageTestsListMapped);
  return DG.DataFrame.fromObjects(resultTestsList);
}

//name: TestsListJoined 
//meta.url: /tests/joinedlist
//output: dataframe df
export async function TestsListJoined(): Promise<DG.DataFrame| undefined> { 
  
  const pacakageTests = await TestAnalysisManager.collectPackageTests();
  const packageTestsListMapped = pacakageTests.map((elem) => {
    return { 'type':  "package ", 'test': elem.packageName + ": " + elem.test.category + ": " + elem.test.name };
  });
  const manualTest = await TestAnalysisManager.collectManualTestNames();
  const manualTestsListMapped = manualTest.map((elem) => {
    return { 'type':  "manual ", 'test': 'Test Track: ' + elem };
  });
  const resultTestsList = DG.DataFrame.fromObjects(manualTestsListMapped.concat(packageTestsListMapped));

  // const builds: DG.DataFrame = await grok.functions.call('UsageAnalysis:Builds'); 
  // const id = builds.get('name', 0); 

  // const tests = await grok.functions.call('UsageAnalysis:getTestStatusesAcordingDF', { 'buildId': id, 'testslist': resultTestsList });
  grok.shell.addTableView(resultTestsList!);
  return resultTestsList;
}


//name: TestAnalysisReportForCurrentDay
//input: datetime date 
//output: dataframe df
export async function TestAnalysisReportForCurrentDay(date: any) {
  const tests = await TestAnalysisManager.collectPackageTests();
  const testsListMapped = tests.map((elem) => {
    return { 'name':  "test-package " + elem.packageName + ": " + elem.test.category + ": " + elem.test.name };
  });
  const testRuns = await grok.functions.call('UsageAnalysis:getServerStartTestResults', { 'date': getDate(new Date(date)), 'testslist': DG.DataFrame.fromObjects(testsListMapped) });
  return testRuns;
}


//name: Usage Analysis
//tags: app
//meta.url: /
//meta.browsePath: Admin
//input: string path {isOptional: true; meta.url: true}
//input: string date {isOptional: true}
//input: string groups {isOptional: true}
//input: string packages {isOptional: true}
//input: string tags {isOptional: true}
//input: string categories {isOptional: true}
//input: string projects {isOptional: true}
//input: map params {isOptional: true}
//output: view v
export async function usageAnalysisApp(path?: string, date?: string, groups?: string, packages?: string, tags?: string, categories?: string, projects?: string): Promise<DG.ViewBase | null> {
  const handler = new ViewHandler();
  await handler.init(date, groups, packages, tags, categories, projects, path);
  return handler.view;
}

//name: Test Track
//tags: app
//meta.url: /tests/manager
//meta.browsePath: Admin
export function testTrackApp(): void {
  if (!grok.shell.dockManager.findNode(TestTrack.getInstance().root))
    TestTrack.getInstance().init();
  else
    TestTrack.getInstance().reopen(); 
}

//name: Reports
//tags: app
//meta.url: /reports
//meta.browsePath: Admin
//input: string path {isOptional: true; meta.url: true}
//input: map params {isOptional: true}
//output: view v
export async function reportsApp(path?: string): Promise<DG.ViewBase> {
  const parent = grok.functions.getCurrentCall();
  const app = new ReportingApp(parent);
  app.init(path).catch((e) => console.log(e));
  return app.empty;
}


let AVAILABLE_SERVICES: string[];

//name: Service Logs
//tags: app
//meta.url: /service-logs
//meta.browsePath: Admin
//input: string path {isOptional: true; meta.url: true}
//input: map params {isOptional: true}
//input: int limit {isOptional: true}
//output: view v
export function serviceLogsApp(path?: string, params?: any, limit?: number): DG.ViewBase {
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

//input: dynamic treeNode
//input: view browseView
export async function serviceLogsAppTreeBrowser(treeNode: DG.TreeViewGroup, browseView: any) {
  const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  const loaderItem = treeNode.item(loaderDiv);
  try {
    AVAILABLE_SERVICES ??= await grok.dapi.docker.getAvailableServices();
    const services = ['datagrok', 'grok_pipe', 'grok-pipe', 'rabbitmq', 'jkg', 'spawner', 'grok-connect', 'grok_connect'];

    const getItem = (service: string) => {
      let icon: HTMLElement;
      if (services.some((s) => service.includes(s))) {
        icon = ui.image('/images/entities/grok.png', 10, 10);
      }
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

//input: dynamic treeNode
//input: view browseView
export async function reportsAppTreeBrowser(treeNode: DG.TreeViewGroup, browseView: any) {
  await treeNode.group('Reports', null, false).loadSources(grok.dapi.reports.by(10));
  await treeNode.group('Rules', null, false).loadSources(grok.dapi.rules.include('actions,actions.assignee').by(10));
}

//name: Usage
//meta.canView: Developers,Administrators
//output: widget result
//tags: dashboard
//test: usageWidget()
export function usageWidget(): DG.Widget {
  return new UsageWidget();
}

//name: Reports
//meta.canView: Developers,Administrators
//output: widget result
//tags: dashboard
//test: reportsWidget()
export function reportsWidget(): DG.Widget {
  return new ReportsWidget();
}

//name: packageUsageWidget
//input: object package
//output: widget result
export function packageUsageWidget(pack: DG.Package): DG.Widget {
  return new PackageUsageWidget(pack);
}

//name: testDashboardsViewer
//tags: viewer
//output: viewer result
export function testDashboardsViewer(): TestDashboardWidget {
  return new TestDashboardWidget();
}

//tags: autostart
export function describeCurrentObj(): void {
  grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
    const ent = acc.context;
    if (ent != null && ent.constructor.name === 'Package') {
      const pane = acc.addPane('Usage', () => ui.wait(async () => {
        let widget: HTMLElement;
        try {
          widget = packageUsageWidget(ent).root;
        } catch (e) {
          widget = ui.divText('Error on loading', { style: { color: 'var(--failure)' } });
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

//name: Create JIRA ticket
//description: Creates JIRA ticket using current error log  
export function createJiraTicket(msg: string){ 
  grok.data.query('JiraCreateIssue', {
    'createRequest': JSON.stringify({
      'fields': {
        'project': {
          'key': 'GROK',
        },
        'summary': 'test',
        'description':'',
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
