import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: _initUA
//tags: init
export function _initUA() : void {
  PackageFunctions._initUA();
}

//name: TestsList
//output: dataframe result
//meta.url: /tests/list
export async function TestsList() : Promise<any> {
  return PackageFunctions.TestsList();
}

//name: TestsListJoined
//output: dataframe result
//meta.url: /tests/joinedlist
export async function TestsListJoined() : Promise<any> {
  return PackageFunctions.TestsListJoined();
}

//name: TestAnalysisReportForCurrentDay
//input: datetime date 
//output: dataframe result
export async function TestAnalysisReportForCurrentDay(date: any) : Promise<any> {
  return PackageFunctions.TestAnalysisReportForCurrentDay(date);
}

//name: Usage Analysis
//tags: app
//input: string path { optional: true; meta.url: true }
//input: string date { optional: true }
//input: string groups { optional: true }
//input: string packages { optional: true }
//input: string tags { optional: true }
//input: string categories { optional: true }
//input: string projects { optional: true }
//output: view result
//meta.url: /
//meta.browsePath: Admin
export async function usageAnalysisApp(path?: string, date?: string, groups?: string, packages?: string, tags?: string, categories?: string, projects?: string) : Promise<any> {
  return PackageFunctions.usageAnalysisApp(path, date, groups, packages, tags, categories, projects);
}

//name: Test Track
//tags: app
//meta.url: /tests/manager
//meta.browsePath: Admin
export function testTrackApp() : void {
  PackageFunctions.testTrackApp();
}

//name: Reports
//tags: app
//input: string path { optional: true; meta.url: true }
//output: view result
//meta.url: /reports
//meta.browsePath: Admin
export async function reportsApp(path?: string) : Promise<any> {
  return PackageFunctions.reportsApp(path);
}

//name: Service Logs
//tags: app
//input: string path { optional: true; meta.url: true }
//input: map params { optional: true }
//input: int limit { optional: true }
//output: view result
//meta.url: /service-logs
//meta.browsePath: Admin
export function serviceLogsApp(path?: string, params?: any, limit?: number) : any {
  return PackageFunctions.serviceLogsApp(path, params, limit);
}

//name: serviceLogsAppTreeBrowser
//input: dynamic treeNode 
//meta.role: appTreeBrowser
export async function serviceLogsAppTreeBrowser(treeNode: any) : Promise<void> {
  PackageFunctions.serviceLogsAppTreeBrowser(treeNode);
}

//name: reportsAppTreeBrowser
//input: dynamic treeNode 
export async function reportsAppTreeBrowser(treeNode: any) : Promise<void> {
  PackageFunctions.reportsAppTreeBrowser(treeNode);
}

//name: Usage
//tags: dashboard
//output: widget result
//meta.canView: Developers,Administrators
//test: usageWidget()
export function usageWidget() : any {
  return PackageFunctions.usageWidget();
}

//name: Reports
//tags: dashboard
//output: widget result
//meta.canView: Developers,Administrators
//test: reportsWidget()
export function reportsWidget() : any {
  return PackageFunctions.reportsWidget();
}

//name: packageUsageWidget
//input: object pckg 
//output: widget result
export function packageUsageWidget(pckg: any) : any {
  return PackageFunctions.packageUsageWidget(pckg);
}

//name: testDashboardsViewer
//tags: viewer
//output: viewer result
export function testDashboardsViewer() : any {
  return PackageFunctions.testDashboardsViewer();
}

//name: describeCurrentObj
//tags: autostart
export function describeCurrentObj() : void {
  PackageFunctions.describeCurrentObj();
}

//name: Create JIRA ticket
//description: Creates JIRA ticket using current error log
export function createJiraTicket() : void {
  PackageFunctions.createJiraTicket();
}
