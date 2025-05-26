export declare namespace Funcs {
    function info(): Promise<any>;
    function init(): Promise<any>;
    function projectsList(): Promise<any>;
    function projectData(projectKey: string): Promise<any>;
    function issueData(issueKey: string): Promise<any>;
    function getJiraField(field: string): Promise<any>;
    function getJiraTicketsByFilter(filter: any): Promise<any>;
}
