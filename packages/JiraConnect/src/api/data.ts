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

function buildJql(filterObject?: object, keys?: string[]): string {
    const parts: string[] = [];
    if (filterObject) {
        const f = Object.keys(filterObject)
            .map((k) => `${k} = ${(filterObject as Record<string, any>)[k]}`).join(' AND ');
        if (f)
            parts.push(f);
    }
    if (keys && keys.length) {
        const k = keys.map((x) => x.trim()).filter((x) => x).join(', ');
        if (k)
            parts.push(`key in (${k})`);
    }
    return parts.join(' AND ');
}

// Uses the enhanced JQL search (/rest/api/3/search/jql). The legacy /rest/api/3/search was removed by
// Atlassian (returns HTTP 410). Token-paginated: fetches all pages and returns them as a JiraIssuesList.
export async function loadIssues(host: string, creds: AuthCreds, filterObject?: object, keys?: string[],
    count: number = 100, fields?: string[]): Promise<JiraIssuesList> {
    const jql = buildJql(filterObject, keys);
    const fieldsParam = fields && fields.length ? fields.join(',') : '*navigable';
    const issues: JiraIssue[] = [];
    let token: string | undefined = undefined;
    let guard = 0;
    do {
        let url = `https://${host}/rest/api/3/search/jql?maxResults=${count}&fields=${encodeURIComponent(fieldsParam)}`;
        if (jql)
            url += `&jql=${encodeURIComponent(jql)}`;
        if (token)
            url += `&nextPageToken=${encodeURIComponent(token)}`;
        const resp = await invokeApiFetch(url, buildRequestOptions(url, creds));
        if (!resp || (resp as ErrorMessageResponse).errorMessages)
            break;
        if (Array.isArray(resp.issues))
            issues.push(...resp.issues);
        token = resp.isLast ? undefined : resp.nextPageToken;
    } while (token && ++guard < 100);
    return {maxResults: count, startAt: 0, total: issues.length, issues};
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