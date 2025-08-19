import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() {
  return PackageFunctions.info();
}

//name: _init
//tags: autostart
export async function _init() {
  return PackageFunctions._init();
}

//name: projectsList
//output: object projects
export async function projectsList() {
  return PackageFunctions.projectsList();
}

//name: projectData
//input: string projectKey 
//output: object projects
export async function projectData(projectKey: string) {
  return PackageFunctions.projectData(projectKey);
}

//name: issueData
//input: string issueKey 
//output: object issue
export async function issueData(issueKey: string) {
  return PackageFunctions.issueData(issueKey);
}

//name: getJiraField
//input: column<string> ticketColumn 
//input: string field 
//output: column result
//meta.vectorFunc: true
export async function getJiraField(ticketColumn: DG.Column, field: string) {
  return PackageFunctions.getJiraField(ticketColumn, field);
}

//name: getJiraTicketsByFilter
//input: object filter 
//output: object result
export async function getJiraTicketsByFilter(filter?: any) {
  return PackageFunctions.getJiraTicketsByFilter(filter);
}
