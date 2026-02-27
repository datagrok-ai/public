import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export function _initUA() : void {
  PackageFunctions._initUA();
}

//output: dataframe result
//meta.url: /tests/list
export async function TestsList() : Promise<any> {
  return await PackageFunctions.TestsList();
}

//output: dataframe result
//meta.url: /tests/joinedlist
export async function TestsListJoined() : Promise<any> {
  return await PackageFunctions.TestsListJoined();
}

//input: datetime date 
//output: dataframe result
export async function TestAnalysisReportForCurrentDay(date: any) : Promise<any> {
  return await PackageFunctions.TestAnalysisReportForCurrentDay(date);
}

//name: Usage Analysis
//input: string path { optional: true; meta.url: true }
//input: string date { optional: true }
//input: string groups { optional: true }
//input: string packages { optional: true }
//input: string tags { optional: true }
//input: string categories { optional: true }
//input: string projects { optional: true }
//output: view result
//meta.role: app
//meta.url: /
//meta.browsePath: Admin
export function usageAnalysisApp(path?: string, date?: string, groups?: string, packages?: string, tags?: string, categories?: string, projects?: string) : any {
  return PackageFunctions.usageAnalysisApp(path, date, groups, packages, tags, categories, projects);
}

//input: column<string> ticketColumn 
//input: column<string> resultColumn 
//meta.vectorFunc: true
export async function getTicketsVerdict(ticketColumn: DG.Column, resultColumn: DG.Column) : Promise<void> {
  await PackageFunctions.getTicketsVerdict(ticketColumn, resultColumn);
}

//name: Test Track
//meta.role: app
//meta.url: /tests/manager
//meta.browsePath: Admin
export function testTrackApp() : void {
  PackageFunctions.testTrackApp();
}

//name: Reports
//input: string path { optional: true; meta.url: true }
//output: view result
//meta.role: app
//meta.url: /reports
//meta.browsePath: Admin
export async function reportsApp(path?: string) : Promise<any> {
  return await PackageFunctions.reportsApp(path);
}

//name: Service Logs
//input: string path { optional: true; meta.url: true }
//input: map params { optional: true }
//input: int limit { optional: true }
//output: view result
//meta.role: app
//meta.url: /service-logs
//meta.browsePath: Admin
export function serviceLogsApp(path?: string, params?: any, limit?: number) : any {
  return PackageFunctions.serviceLogsApp(path, params, limit);
}

//input: dynamic treeNode 
//meta.role: appTreeBrowser
//meta.app: Service Logs
export async function serviceLogsAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.serviceLogsAppTreeBrowser(treeNode);
}

//input: dynamic treeNode 
//meta.role: appTreeBrowser
//meta.app: Reports
export async function reportsAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.reportsAppTreeBrowser(treeNode);
}

//name: Usage
//output: widget result
//meta.canView: Developers,Administrators
//meta.role: dashboard
//test: usageWidget()
export function usageWidget() : any {
  return PackageFunctions.usageWidget();
}

//name: Reports
//output: widget result
//meta.canView: Developers,Administrators
//meta.role: dashboard
//test: reportsWidget()
export function reportsWidget() : any {
  return PackageFunctions.reportsWidget();
}

//input: object pckg 
//output: widget result
export function packageUsageWidget(pckg: any) : any {
  return PackageFunctions.packageUsageWidget(pckg);
}

//output: viewer result
//meta.showInGallery: false
//meta.role: viewer
export function testDashboardsViewer() : any {
  return PackageFunctions.testDashboardsViewer();
}

//meta.role: autostart
export function describeCurrentObj() : void {
  PackageFunctions.describeCurrentObj();
}

//name: Create JIRA ticket
//description: Creates JIRA ticket using current error log
export function createJiraTicket() : void {
  PackageFunctions.createJiraTicket();
}
