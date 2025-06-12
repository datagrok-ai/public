import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Queries {
  export async function allTestRuns(benchmarks: boolean): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:AllTestRuns', { benchmarks });
  }

  export async function uniqueUsersCount(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UniqueUsersCount', { date });
  }

  export async function newUsersCount(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:NewUsersCount', { date });
  }

  export async function sessionsCount(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:SessionsCount', { date });
  }

  export async function viewsCount(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ViewsCount', { date });
  }

  export async function connectionsCount(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ConnectionsCount', { date });
  }

  export async function queriesCount(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:QueriesCount', { date });
  }

  export async function testsCount(date: any): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TestsCount', { date });
  }

  export async function topQueriesUsingDataSource(date: string, data_source: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopQueriesUsingDataSource', { date, data_source });
  }

  export async function topUsersOfQuery(name: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopUsersOfQuery', { name, date });
  }

  export async function topUsersOfConnection(name: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopUsersOfConnection', { name, date });
  }

  export async function queries1(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:Queries1', { date });
  }

  export async function topQueries(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopQueries', { date });
  }

  export async function topConnections(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopConnections', { date });
  }

  export async function topDataSources(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopDataSources', { date });
  }

  export async function entityLinks(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:EntityLinks', {});
  }

  export async function functionInfoByFriendlyName(name: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:FunctionInfoByFriendlyName', { name });
  }

  export async function functionInfoBySource(name: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:FunctionInfoBySource', { name, date });
  }

  export async function topDisabledErrors(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopDisabledErrors', { date });
  }

  export async function topPackageErrors(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopPackageErrors', { date });
  }

  export async function topErrorSources(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopErrorSources', { date });
  }

  export async function eventErrors(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:EventErrors', { date });
  }

  export async function reportsCount(date: string, event_id: string): Promise<number> {
    return await grok.data.query('UsageAnalysis:ReportsCount', { date, event_id });
  }

  export async function sameErrors(date: string, event_id: string): Promise<number> {
    return await grok.data.query('UsageAnalysis:SameErrors', { date, event_id });
  }

  export async function topErrors(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopErrors', { date });
  }

  export async function packageInfo(name: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PackageInfo', { name });
  }

  export async function functionInfoByName(name: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:FunctionInfoByName', { name });
  }

  export async function topUsersOfPackage(name: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopUsersOfPackage', { name, date });
  }

  export async function topFunctionsOfPackage(date: string, name: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopFunctionsOfPackage', { date, name });
  }

  export async function topErrorsOfPackage(date: string, name: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopErrorsOfPackage', { date, name });
  }

  export async function topUsersOfFunction(name: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopUsersOfFunction', { name, date });
  }

  export async function topFunctionsOfSource(name: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopFunctionsOfSource', { name, date });
  }

  export async function topUsersOfSource(name: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopUsersOfSource', { name, date });
  }

  export async function topFunctions(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopFunctions', { date });
  }

  export async function topPackageFunctions(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopPackageFunctions', { date });
  }

  export async function events1(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:Events1', { date });
  }

  export async function topPackages(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopPackages', { date });
  }

  export async function topSources(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopSources', { date });
  }

  export async function eventsSources(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:EventsSources', { date });
  }

  export async function eventsUsersSources(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:EventsUsersSources', { date });
  }

  export async function getUsersInGroups(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:GetUsersInGroups', {});
  }

  export async function eventByErrorMessageAndFriendlyName(errorMessage: string, friendlyName: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:EventByErrorMessageAndFriendlyName', { errorMessage, friendlyName });
  }

  export async function updateEventsIsErrorComment(errorMessage: string, friendlyName: string, isError: boolean, comment: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UpdateEventsIsErrorComment', { errorMessage, friendlyName, isError, comment });
  }

  export async function functionErrors(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:FunctionErrors', { date });
  }

  export async function topFunctionErrors(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopFunctionErrors', { date });
  }

  export async function topFunctionDisabledErrors(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopFunctionDisabledErrors', { date });
  }

  export async function topPackagesByError(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopPackagesByError', { date });
  }

  export async function functionsUsage(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:FunctionsUsage', { date });
  }

  export async function functionsContextPane(time_start: number, time_end: number): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:FunctionsContextPane', { time_start, time_end });
  }

  export async function functionsExecTime(function: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:FunctionsExecTime', { function });
  }

  export async function entitiesTags(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:EntitiesTags', {});
  }

  export async function groups(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:Groups', {});
  }

  export async function userById(id: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UserById', { id });
  }

  export async function userInfoByEmailPanel(email: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UserInfoByEmailPanel', { email });
  }

  export async function logActions(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LogActions', {});
  }

  export async function logActionsSummary(eventTime: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LogActionsSummary', { eventTime });
  }

  export async function logActionsSummaryByHours(eventTime: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LogActionsSummaryByHours', { eventTime });
  }

  export async function logErrorSummary(eventTime: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LogErrorSummary', { eventTime });
  }

  export async function logParameterValues(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LogParameterValues', {});
  }

  export async function logSessions(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LogSessions', {});
  }

  export async function logTail(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LogTail', { date });
  }

  export async function uniqueUsersList(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UniqueUsersList', { date });
  }

  export async function totalUsersAndGroups(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TotalUsersAndGroups', {});
  }

  export async function uniqueUsersOverview(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UniqueUsersOverview', { date });
  }

  export async function packagesUsageOverview(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PackagesUsageOverview', { date });
  }

  export async function packagesUsage(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PackagesUsage', { date });
  }

  export async function packagesContextPaneFunctions(time_start: number, time_end: number): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PackagesContextPaneFunctions', { time_start, time_end });
  }

  export async function packagesContextPaneLogs(time_start: number, time_end: number): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PackagesContextPaneLogs', { time_start, time_end });
  }

  export async function packagesContextPaneAudit(time_start: number, time_end: number): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PackagesContextPaneAudit', { time_start, time_end });
  }

  export async function packagesInstallationTime(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PackagesInstallationTime', { date });
  }

  export async function packagesCategories(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PackagesCategories', {});
  }

  export async function groupCompleter(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:GroupCompleter', { sub });
  }

  export async function userCompleter(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UserCompleter', { sub });
  }

  export async function actionCompleter(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ActionCompleter', { sub });
  }

  export async function queriesCompleter(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:QueriesCompleter', { sub });
  }

  export async function projectCompleter(sub: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ProjectCompleter', { sub });
  }

  export async function actionByUserOnDate(action: string, user: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ActionByUserOnDate', { action, user, date });
  }

  export async function groupMembers(group: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:GroupMembers', { group });
  }

  export async function actionsByUserOnDate(user: string, date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ActionsByUserOnDate', { user, date });
  }

  export async function pgStatStatements(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:PgStatStatements', {});
  }

  export async function uniqueUsersPerProject(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UniqueUsersPerProject', { date });
  }

  export async function userAccessFrequencyPerProject(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UserAccessFrequencyPerProject', { date });
  }

  export async function accessCountPerPeriodPerProject(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:AccessCountPerPeriodPerProject', { date });
  }

  export async function projectsList(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ProjectsList', {});
  }

  export async function userReports(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UserReports', { date });
  }

  export async function reportSameErrors(stackTraceHash: string, errorMessage: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ReportSameErrors', { stackTraceHash, errorMessage });
  }

  export async function reportsTop20(packageOwnerId: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ReportsTop20', { packageOwnerId });
  }

  export async function userReportsSingle(reportNumber: number): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UserReportsSingle', { reportNumber });
  }

  export async function reportsMigration(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ReportsMigration', {});
  }

  export async function reportDataMigration(report_id: string, id: string, screenshot: string, details: string, client_settings: string, server_settings: string, errors: string, client_log: string, server_log: string, console: string, queries_log: string, containers_log: string, images_log: string, services: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ReportDataMigration', { report_id, id, screenshot, details, client_settings, server_settings, errors, client_log, server_log, console, queries_log, containers_log, images_log, services });
  }

  export async function getSystemTableSizes(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:GetSystemTableSizes', {});
  }

  export async function benchmarkAnalysis(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:BenchmarkAnalysis', {});
  }

  export async function benchmarksDashboard(instanceFilter: string, lastBuildsNum: number, showNotRun: boolean, showBenchmarks: boolean, showNotCiCd: boolean): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:BenchmarksDashboard', { instanceFilter, lastBuildsNum, showNotRun, showBenchmarks, showNotCiCd });
  }

  export async function testsDashboard(instanceFilter: string, lastBuildsNum: number, packageFilter: string, showNotRun: boolean, showBenchmarks: boolean, showNotCiCd: boolean, versionFilter: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TestsDashboard', { instanceFilter, lastBuildsNum, packageFilter, showNotRun, showBenchmarks, showNotCiCd, versionFilter });
  }

  export async function manualTests(lastBatchesNum: number): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ManualTests', { lastBatchesNum });
  }

  export async function stressTestsDashboard(lastBuildsNum: number): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:StressTestsDashboard', { lastBuildsNum });
  }

  export async function manualTicketCreation(name: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ManualTicketCreation', { name });
  }

  export async function manualTicketFetch(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ManualTicketFetch', {});
  }

  export async function benchmarks(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:Benchmarks', {});
  }

  export async function stressTests(batch_name: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:StressTests', { batch_name });
  }

  export async function lastBuildsBenchmarksCompare(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LastBuildsBenchmarksCompare', {});
  }

  export async function lastBuildsCompare(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LastBuildsCompare', {});
  }

  export async function lastVersionsCompare(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LastVersionsCompare', {});
  }

  export async function testTrack(batchName: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TestTrack', { batchName });
  }

  export async function lastStatuses(path: string, batchToExclude: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:LastStatuses', { path, batchToExclude });
  }

  export async function testingNames(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TestingNames', {});
  }

  export async function testingName(version: string, uid: string, start: string): Promise<string> {
    return await grok.data.query('UsageAnalysis:TestingName', { version, uid, start });
  }

  export async function builds(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:Builds', {});
  }

  export async function getTestStatusesAcordingDF(buildId: string, testslist: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:GetTestStatusesAcordingDF', { buildId, testslist });
  }

  export async function manualTestingTestStatuses(batchName: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:ManualTestingTestStatuses', { batchName });
  }

  export async function topEventsOfUser(date: string, name: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopEventsOfUser', { date, name });
  }

  export async function uniqueUsers(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UniqueUsers', { date });
  }

  export async function uniqueSessions(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UniqueSessions', { date });
  }

  export async function usage(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:Usage', { date });
  }

  export async function topPackagesByUsers(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopPackagesByUsers', { date });
  }

  export async function topUsers(date: string): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:TopUsers', { date });
  }

  export async function schemaColumns(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:SchemaColumns', {});
  }

  export async function uniqueUsersSummary(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UniqueUsersSummary', {});
  }

  export async function usersEventsSummary(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UsersEventsSummary', {});
  }

  export async function usersErrorsSummary(): Promise<DG.DataFrame> {
    return await grok.data.query('UsageAnalysis:UsersErrorsSummary', {});
  }
}

export namespace Funcs {
  export async function initUA(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:InitUA', {});
  }

  export async function testsList(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:TestsList', {});
  }

  export async function testsListJoined(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:TestsListJoined', {});
  }

  export async function testAnalysisReportForCurrentDay(date: any): Promise<any> {
    return await grok.functions.call('UsageAnalysis:TestAnalysisReportForCurrentDay', { date });
  }

  export async function usageAnalysisApp(path: string, date: string, groups: string, packages: string, tags: string, categories: string, projects: string, params: any): Promise<any> {
    return await grok.functions.call('UsageAnalysis:UsageAnalysisApp', { path, date, groups, packages, tags, categories, projects, params });
  }

  export async function testTrackApp(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:TestTrackApp', {});
  }

  export async function reportsApp(path: string, params: any): Promise<any> {
    return await grok.functions.call('UsageAnalysis:ReportsApp', { path, params });
  }

  export async function serviceLogsApp(path: string, params: any, limit: number): Promise<any> {
    return await grok.functions.call('UsageAnalysis:ServiceLogsApp', { path, params, limit });
  }

  export async function serviceLogsAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('UsageAnalysis:ServiceLogsAppTreeBrowser', { treeNode, browseView });
  }

  export async function reportsAppTreeBrowser(treeNode: any, browseView: DG.View): Promise<any> {
    return await grok.functions.call('UsageAnalysis:ReportsAppTreeBrowser', { treeNode, browseView });
  }

  export async function usageWidget(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:UsageWidget', {});
  }

  export async function reportsWidget(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:ReportsWidget', {});
  }

  export async function packageUsageWidget(package: any): Promise<any> {
    return await grok.functions.call('UsageAnalysis:PackageUsageWidget', { package });
  }

  export async function testDashboardsViewer(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:TestDashboardsViewer', {});
  }

  export async function describeCurrentObj(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:DescribeCurrentObj', {});
  }

  //Creates JIRA ticket using current error log  
  export async function createJiraTicket(): Promise<any> {
    return await grok.functions.call('UsageAnalysis:CreateJiraTicket', {});
  }
}
