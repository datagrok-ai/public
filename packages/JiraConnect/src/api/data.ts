/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Project, JiraIssue, AuthCreds, JiraIssuesList } from './objects';


export async function loadProjectsList(host: string, creds: AuthCreds): Promise<Project[] | null> {
    const url = `https://${host}/rest/api/3/project`;

    let requestOptions = buildRequestOptions(url, creds);
    let projects: Project[] | null = await invokeApiFetch(url, requestOptions);
    return projects;
}

export async function loadProjectData(host: string, creds: AuthCreds, projectKey: string){
    const url = `https://${host}/rest/api/3/project/${projectKey}`;

    let requestOptions = buildRequestOptions(url, creds);
    let project: Project | null = await invokeApiFetch(url, requestOptions);
    return project; 
}

export async function loadIssues(host: string, creds: AuthCreds, index: number = 0, count: number = 50, projectId?: number, keys?: string[]): Promise<JiraIssuesList | null> {
    let url = `https://${host}/rest/api/3/search?startAt=${index}&maxResults=${count}`;
    if (projectId || keys) {
        let trimmedKeys = keys?.map(key => key.trim());
        let keysParams = trimmedKeys?.join(',').trim();
        let projectJQL =  `project+in+(${projectId})`;
        let keysJQL =  `key+in+(${keysParams})`;
        let resultJQL = [];

        if(projectId)
            resultJQL.push(projectJQL);
        if(keysParams)
            resultJQL.push(keysJQL);
        
        url = `${url}&jql=${resultJQL.join('+and+')}`;
    }
    let requestOptions = buildRequestOptions(url, creds);
    let issues: JiraIssuesList | null = await invokeApiFetch(url, requestOptions);
    return issues;
}

export async function loadIssueData(host: string, creds: AuthCreds, issueName: string): Promise<JiraIssue | null> {
    let url = `https://${host}/rest/api/3/issue/${issueName}`;
    let requestOptions = buildRequestOptions(url, creds);
    let issue: JiraIssue | null = await invokeApiFetch(url, requestOptions);
    return issue;
}


function buildRequestOptions(url: string, creds: AuthCreds): RequestInit {
    const headers = new Headers();
    headers.append('Accept', 'application/json');
    headers.append('Authorization', `Basic ${btoa(creds.userName + ":" + creds.authKey)}`);
    headers.append('original-url', url);
    headers.append('original-method', 'GET');

    const requestOptions: RequestInit = {
        method: "GET",
        headers: headers,
        redirect: "follow"
    };
    return requestOptions;
}

async function invokeApiFetch(url: string, options: RequestInit): Promise<any> {
    try {
        const response = await grok.dapi.fetchProxy(url, options);
        return await response.json();;
    } catch (error) {
        console.error('Error fetching projects:', error);
    }
    return null;
}