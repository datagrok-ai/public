import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Jiraconnect:Info', {});
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('Jiraconnect:Init', {});
  }

  export async function projectsList(): Promise<any> {
    return await grok.functions.call('Jiraconnect:ProjectsList', {});
  }

  export async function projectData(projectKey: string): Promise<any> {
    return await grok.functions.call('Jiraconnect:ProjectData', { projectKey });
  }

  export async function issueData(issueKey: string): Promise<any> {
    return await grok.functions.call('Jiraconnect:IssueData', { issueKey });
  }

  export async function getJiraField(field: string): Promise<any> {
    return await grok.functions.call('Jiraconnect:GetJiraField', { field });
  }

  export async function getJiraTicketsByFilter(filter: any): Promise<any> {
    return await grok.functions.call('Jiraconnect:GetJiraTicketsByFilter', { filter });
  }
}
