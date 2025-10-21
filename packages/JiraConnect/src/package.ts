/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { loadProjectsList, loadIssues, loadIssueData, loadProjectData, loadIssuesBulk } from './api/data';
import {ErrorMessageResponse, JiraIssue, JiraIssuesBulkList, Project} from './api/objects';
import { AuthCreds } from './api/objects';
import { getJiraCreds } from './app/credentials';
import { addJiraDetector } from './detectors';
import { JiraGridCellHandler } from './jira-grid-cell-handler';
import '../css/jira_connect.css';
export * from './package.g';
export let _package = new DG.Package();
export const cache = new DG.LruCache<string, JiraIssue>(100);
export const queried = new Set<string>();


export class PackageFunctions {
  @grok.decorators.func()
  static info() {
    grok.shell.info(_package.webRoot);
  }


  @grok.decorators.autostart()
  static async _init(): Promise<void> {
    const projects = await PackageFunctions.projectsList() ?? [];
    if (projects.length > 0) {
      const projectNames = projects.map(e => e.key);
      addJiraDetector(projectNames, '');
      DG.ObjectHandler.register(new JiraGridCellHandler());
    }
  }


  @grok.decorators.func({outputs: [{name: 'projects', type: 'object'}]})
  static async projectsList(): Promise<Project[] | null> {
    const jiraCreds = await getJiraCreds();
    if (!jiraCreds)
      return null;
    return await loadProjectsList(jiraCreds.host, new AuthCreds(jiraCreds.userName, jiraCreds.authKey));
  }


  @grok.decorators.func({outputs: [{name: 'projects', type: 'object'}]})
  static async projectData(projectKey: string): Promise<Project | null> {
    const jiraCreds = await getJiraCreds();
    if (!jiraCreds)
      return null;
    return await loadProjectData(jiraCreds.host, new AuthCreds(jiraCreds.userName, jiraCreds.authKey), projectKey);
  }


  @grok.decorators.func({outputs: [{name: 'issue', type: 'object'}]})
  static async issueData(issueKey: string): Promise<JiraIssue | null> {
    const jiraCreds = await getJiraCreds();
    if (!jiraCreds)
      return null;
    if (!issueKey || issueKey?.length === 0)
      return null;
    if (cache.has(issueKey))
      return cache.get(issueKey)!;

    queried.add(issueKey);
    const issue = await loadIssueData(jiraCreds.host, new AuthCreds(jiraCreds.userName, jiraCreds.authKey), issueKey);
    if (!issue || (issue as ErrorMessageResponse)?.errorMessages)
      return null;
    queried.delete(issueKey);
    cache.set(issueKey, issue);

    return (issue as JiraIssue);
  }


  @grok.decorators.func({meta: {vectorFunc: 'true'}})
  static async getJiraField( @grok.decorators.param({type:'column<string>'})  ticketColumn: DG.Column,
    field: string): Promise<DG.Column | undefined> {
  
    const jiraCreds = await getJiraCreds();
    if (!jiraCreds)
      return undefined;
    const ticketsMap = new Map<string, string>()
    const keys = field.split(':');
    const ticketColumnList = ticketColumn.toList();
    for (let i = 0; i < ticketColumnList.length; i++)
      if (ticketColumnList[i])
        ticketsMap.set(ticketColumnList[i], '');

    const ticketKeys = Array.from(ticketsMap.keys()).map(key => (key ?? '').trim());
    const totalTickets = ticketKeys.length;
    const chunkSize = 100;
    let index = 0;
    while (index < totalTickets) {
      const keysToLoad = ticketKeys.slice(index, index + chunkSize);
      try {
        const loadedIssues = await loadIssues(jiraCreds.host, new AuthCreds(jiraCreds.userName, jiraCreds.authKey),
          0, chunkSize, undefined, keysToLoad);

        upperFor:
        for (let i = 0; i < (loadedIssues?.issues.length ?? 0); i++) {
          let current: any = loadedIssues?.issues[i].fields;

          for (const key of keys) {
            if (current && key in current)
              current = current[key];
            else {
              ticketsMap.set(ticketColumnList[index + i], '');
              continue upperFor;
            }
          }
          if (Array.isArray(current))
            current = (current ?? ['']).join(', ');

          if (loadedIssues?.issues[i]?.key)
            ticketsMap.set(loadedIssues.issues[i].key, (current ?? '') as string);
        }
      } catch (error) {
        console.error(`Error loading issues for index range ${index} - ${index + chunkSize}`, error);
      }

      index += chunkSize;
    }

    let resultList: string[] = [];
    for (let ticket of ticketColumnList)
      resultList.push(ticketsMap.get((ticket ?? '').trim()) ?? '');

    return DG.Column.fromStrings(field, resultList);
  }


  @grok.decorators.func({outputs: [{name: 'result', type: 'object'}]})
  static async getJiraTicketsByFilter(@grok.decorators.param({type: 'object'}) filter?: object): Promise<JiraIssue[]> {
    const jiraCreds = await getJiraCreds();
    if (!jiraCreds)
      return [];
    let result: JiraIssue[] = [];
    let total = -1;
    const chunkSize = 100;
    let startAt = 0;
    while (true) {
      const loadedIssues = await loadIssues(jiraCreds.host, new AuthCreds(jiraCreds.userName, jiraCreds.authKey),
        startAt, chunkSize, filter, undefined);
      if (loadedIssues != null) {
        total = loadedIssues.total;
        result = [...result, ...loadedIssues.issues]
      }

      if (total <= startAt)
        break;
      startAt += chunkSize;
    }
    return result;
  }

  @grok.decorators.func({outputs: [{name: 'result', type: 'object'}]})
  static async getJiraTicketsBulk(
      @grok.decorators.param({type: 'list'}) issueIdsOrKeys: string[],
      @grok.decorators.param({type: 'list'}) fields?: string[],
      @grok.decorators.param({type: 'list'}) expand?: string[]
  ): Promise<JiraIssuesBulkList> {
    const jiraCreds = await getJiraCreds();
    if (!jiraCreds)
      throw new Error('JiraConnect: Missing credentials');

    const result: JiraIssuesBulkList = {
      expand: expand ?? ['names', 'schema'],
      issues: [],
      issueErrors: []
    };

    const chunkSize = 100;
    let index = 0;

    while (index < issueIdsOrKeys.length) {
      const keysToLoad = issueIdsOrKeys.slice(index, index + chunkSize);
      try {
        const loadedIssues = await loadIssuesBulk(jiraCreds.host, new AuthCreds(jiraCreds.userName, jiraCreds.authKey),
            keysToLoad, expand, fields);
        result.issues.push(...loadedIssues.issues);
        result.issueErrors.push(...loadedIssues.issueErrors);

      } catch (error) {
        console.error(`Error loading issues for index range ${index} - ${index + chunkSize}`, error);
      }
      index += chunkSize;
    }

    return result;
  }
}
