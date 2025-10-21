/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import {AuthCreds, ErrorMessageResponse, JiraIssue, JiraIssuesBulkList, JiraIssuesList, Project} from './objects';


export async function loadProjectsList(host: string, creds: AuthCreds): Promise<Project[] | null> {
    const url = `https://${host}/rest/api/3/project`;

    let requestOptions = buildRequestOptions(url, creds);
    let projects: Project[] | null = await invokeApiFetch(url, requestOptions);
    return projects;
}

export async function loadProjectData(host: string, creds: AuthCreds, projectKey: string) {
    const url = `https://${host}/rest/api/3/project/${projectKey}`;

    let requestOptions = buildRequestOptions(url, creds);
    let project: Project | null = await invokeApiFetch(url, requestOptions);
    return project;
}

// https://developer.atlassian.com/cloud/jira/platform/rest/v3/api-group-issues/#api-rest-api-3-issue-bulkfetch-post
export async function loadIssuesBulk(host: string, creds: AuthCreds, issueIdsOrKeys: string[], expand?: string[], fields?: string[]): Promise<JiraIssuesBulkList> {
    if (issueIdsOrKeys.length == 0)
        return {expand: expand ?? ['names', 'schema'], issues: [], issueErrors: []};
    if (issueIdsOrKeys.length > 100)
        throw new Error('You can\'t request more than 100 issues');
    const url = `https://${host}/rest/api/3/issue/bulkfetch`;
    const requestOptions = buildRequestOptions(url, creds, 'POST');
    requestOptions.body = JSON.stringify({
        expand: expand ?? [],
        fields: fields ?? [],
        fieldsByKeys: false,
        issueIdsOrKeys: issueIdsOrKeys,
        properties: []
    });
    return await invokeApiFetch(url, requestOptions);
}

export async function loadIssues(host: string, creds: AuthCreds, index: number = 0, count: number = 50, filterObject?: object, keys?: string[]): Promise<JiraIssuesList> {
    let url = `https://${host}/rest/api/3/search?startAt=${index}&maxResults=${count}`;
    
    if (filterObject || keys) {
        let trimmedKeys = keys?.map(key => key.trim());
        let keysParams = trimmedKeys?.join('+OR+key+=+')?.trim();
        let keysJQL = `key+=+${keysParams}`;
        
        let filterData = [];
        for (let filterKey of Object.keys(filterObject ?? {})) {
            filterData.push(`${filterKey}+=+${(filterObject as Record<string, any>)[filterKey]}`);
        }
        let resultJQL = [];
        if (filterObject)
            resultJQL.push(filterData.join('+AND+'));
        if (keysParams)
            resultJQL.push(keysJQL);

        url = `${url}&jql=${resultJQL.join('+and+')}`;
    }
    let requestOptions = buildRequestOptions(url, creds);
    let issues: JiraIssuesList | ErrorMessageResponse = await invokeApiFetch(url, requestOptions);

    if ((issues as ErrorMessageResponse).errorMessages) {
        issues = { maxResults: count, startAt: index, issues: [], total: 0 }
        for (let issue of keys ?? []) {
            let issueData = await loadIssueData(host, creds, issue)
            if (issueData && (issueData as JiraIssue).key) {
                issues.issues.push(issueData as JiraIssue);
                issues.total = issues.total + 1;
            }
        }
    }
    return issues as JiraIssuesList;
}

export async function loadIssueData(host: string, creds: AuthCreds, issueName: string): Promise<JiraIssue | null | ErrorMessageResponse> {
    let url = `https://${host}/rest/api/3/issue/${issueName}`;
    let requestOptions = buildRequestOptions(url, creds);
    return await invokeApiFetch(url, requestOptions);
}


function buildRequestOptions(url: string, creds: AuthCreds, method: string = 'GET'): RequestInit {
    const headers = new Headers();
    headers.append('Accept', 'application/json');
    headers.append('Authorization', `Basic ${btoa(creds.userName + ":" + creds.authKey)}`);
    headers.append('original-url', url);
    headers.append('original-method', 'GET');

    const requestOptions: RequestInit = {
        method: method,
        headers: headers,
        redirect: "follow"
    };
    return requestOptions;
}

async function invokeApiFetch(url: string, options: RequestInit): Promise<any> {
    try {
        const response = await grok.dapi.fetchProxy(url, options);
        return await response.json();
    } catch (error) {
        console.error('Error fetching projects:', error);
    }
    return null;
}