import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { UsageWidget } from './widgets/usage-widget';
import { PackageUsageWidget } from './widgets/package-usage-widget';
import '../css/usage_analysis.css';
import '../css/test_track.css';
import { ViewHandler } from './view-handler';
import { TestTrack } from './test-track/app';
import { ReportsWidget } from "./widgets/reports-widget";
import { ReportingApp } from "./reporting/reporting_app";
import { TestAnalysisManager } from './test-analysis/test-analysis-manager'; 
import { getDate } from './utils';
import dayjs from "dayjs";
import {ServiceLogsApp} from "./service_logs/service_logs";
import { TestGridCellHandler } from './test-grid-cell-handler';

export const _package = new DG.Package();


//tags: init
export function _initUA(): void {
  DG.ObjectHandler.register(new TestGridCellHandler());
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
    return { 'type':  "manual ", 'test': 'Unknown: ' + elem };
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
//input: map params {isOptional: true}
//output: view v
export async function usageAnalysisApp(path?: string, date?: string, groups?: string, packages?: string): Promise<DG.ViewBase | null> {
  const handler = new ViewHandler();
  await handler.init(date, groups, packages, path);
  return handler.view;
}

//name: Test Track
//tags: app
//meta.url: /tests/manager
//meta.browsePath: Admin
//input: string path {isOptional: true; meta.url: true}
//input: map params {isOptional: true}
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
  await app.init(path);
  return app.view!;
}

//name: Service Logs
//tags: app
//meta.url: /service-logs
//meta.browsePath: Admin
//input: string path {isOptional: true; meta.url: true}
//input: map params {isOptional: true}
//input: int limit {isOptional: true}
//output: view v
export async function serviceLogsApp(path?: string, params?: any, limit?: number): Promise<DG.ViewBase> {
  const currentCall = grok.functions.getCurrentCall();
  const services = await grok.dapi.docker.getAvailableServices();
  const app = new ServiceLogsApp(currentCall, services, path, limit);
  if (services.length > 0)
    app.getLogs().then((_) => {});
  return app;
}

//input: dynamic treeNode
//input: view browseView
export async function reportsAppTreeBrowser(treeNode: DG.TreeViewGroup, browseView: DG.BrowseView) {
  await treeNode.group('Reports', null, false).loadSources(grok.dapi.reports.by(10));
  await treeNode.group('Rules', null, false).loadSources(grok.dapi.rules.include('actions,actions.assignee').by(10));
}

//name: Usage
//output: widget result
//tags: dashboard
//test: usageWidget()
export async function usageWidget(): Promise<DG.Widget | null> {
  return await hasAccess() ? new UsageWidget() : null;
}

//name: Reports
//output: widget result
//tags: dashboard
//test: reportsWidget()
export async function reportsWidget(): Promise<DG.Widget | null> {
  return await hasAccess() ? new ReportsWidget() : null;
}

//name: packageUsageWidget
//input: object package
//output: widget result
export function packageUsageWidget(pack: DG.Package): DG.Widget {
  return new PackageUsageWidget(pack);
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

async function hasAccess(): Promise<boolean> {
  const userGroup = await grok.dapi.groups.find(DG.User.current().group.id);
  return userGroup.memberships.some((g) => g.friendlyName === 'Developers' || g.friendlyName === 'Administrators')
      || userGroup.adminMemberships.some((g) => g.friendlyName === 'Developers' || g.friendlyName === 'Administrators');
}
