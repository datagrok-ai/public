import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: _initUA
//tags: init
export function _initUA() {
  return PackageFunctions._initUA();
}

//name: TestsList
//output: dynamic result
//meta.url: /tests/list
export async function TestsList() {
  return PackageFunctions.TestsList();
}

//name: TestsListJoined
//output: dynamic result
//meta.url: /tests/joinedlist
export async function TestsListJoined() {
  return PackageFunctions.TestsListJoined();
}

//name: TestAnalysisReportForCurrentDay
//input: datetime date 
//output: dynamic result
export async function TestAnalysisReportForCurrentDay(date: any) {
  return PackageFunctions.TestAnalysisReportForCurrentDay(date);
}

//name: aaa
//input: string path { optional: true; meta.url: true }
//output: dynamic result
export async function aaa(path: string) {
  return PackageFunctions.aaa(path);
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
export async function usageAnalysisApp(path: string, date: string, groups: string, packages: string, tags: string, categories: string, projects: string) {
  return PackageFunctions.usageAnalysisApp(path, date, groups, packages, tags, categories, projects);
}

//name: Test Track
//tags: app
//meta.url: /tests/manager
//meta.browsePath: Admin
export function testTrackApp() {
  return PackageFunctions.testTrackApp();
}

//name: Reports
//tags: app
//input: string path { optional: true; meta.url: true }
//output: view result
//meta.url: /reports
//meta.browsePath: Admin
export async function reportsApp(path: string) {
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
export function serviceLogsApp(path: string, params: any, limit: number) {
  return PackageFunctions.serviceLogsApp(path, params, limit);
}

//name: serviceLogsAppTreeBrowser
//input: dynamic treeNode 
//input: dynamic browseView 
//output: dynamic result
export async function serviceLogsAppTreeBrowser(treeNode: any, browseView: any) {
  return PackageFunctions.serviceLogsAppTreeBrowser(treeNode, browseView);
}

//name: reportsAppTreeBrowser
//input: dynamic treeNode 
//input: dynamic browseView 
//output: dynamic result
export async function reportsAppTreeBrowser(treeNode: any, browseView: any) {
  return PackageFunctions.reportsAppTreeBrowser(treeNode, browseView);
}

//name: Usage
//tags: dashboard
//output: widget result
//meta.canView: Developers,Administrators
export function usageWidget() {
  return PackageFunctions.usageWidget();
}

//name: Reports
//tags: dashboard
//output: widget result
//meta.canView: Developers,Administrators
export function reportsWidget() {
  return PackageFunctions.reportsWidget();
}

//name: packageUsageWidget
//input: object pckg 
//output: widget result
export function packageUsageWidget(pckg: any) {
  return PackageFunctions.packageUsageWidget(pckg);
}

//name: testDashboardsViewer
//tags: viewer
//output: viewer result
export function testDashboardsViewer() {
  return PackageFunctions.testDashboardsViewer();
}

//name: describeCurrentObj
//tags: autostart
export function describeCurrentObj() {
  return PackageFunctions.describeCurrentObj();
}

//name: Create JIRA ticket
//description: Creates JIRA ticket using current error log
//output: dynamic result
export function createJiraTicket() {
  return PackageFunctions.createJiraTicket();
}
