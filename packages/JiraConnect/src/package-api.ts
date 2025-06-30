import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('JiraConnect:Info', {});
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('JiraConnect:Init', {});
  }

  export async function projectsList(): Promise<any> {
    return await grok.functions.call('JiraConnect:ProjectsList', {});
  }

  export async function projectData(projectKey: string): Promise<any> {
    return await grok.functions.call('JiraConnect:ProjectData', { projectKey });
  }

  export async function issueData(issueKey: string): Promise<any> {
    return await grok.functions.call('JiraConnect:IssueData', { issueKey });
  }

  export async function getJiraField(ticketColumn: DG.Column, field: string): Promise<any> {
    return await grok.functions.call('JiraConnect:GetJiraField', { ticketColumn, field });
  }

  export async function getJiraTicketsByFilter(filter: any): Promise<any> {
    return await grok.functions.call('JiraConnect:GetJiraTicketsByFilter', { filter });
  }
}
