import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/jiraconnect:Info', {});
  }

  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/jiraconnect:Init', {});
  }

  export async function projectsList(): Promise<any> {
    return await grok.functions.call('@datagrok/jiraconnect:ProjectsList', {});
  }

  export async function projectData(projectKey: string): Promise<any> {
    return await grok.functions.call('@datagrok/jiraconnect:ProjectData', { projectKey });
  }

  export async function issueData(issueKey: string): Promise<any> {
    return await grok.functions.call('@datagrok/jiraconnect:IssueData', { issueKey });
  }

  export async function getJiraField(field: string): Promise<any> {
    return await grok.functions.call('@datagrok/jiraconnect:GetJiraField', { field });
  }

  export async function getJiraTicketsByFilter(filter: any): Promise<any> {
    return await grok.functions.call('@datagrok/jiraconnect:GetJiraTicketsByFilter', { filter });
  }
}
