import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//meta.role: autostart
export async function _init() : Promise<void> {
  await PackageFunctions._init();
}

//output: object projects
export async function projectsList() : Promise<any> {
  return await PackageFunctions.projectsList();
}

//input: string projectKey 
//output: object projects
export async function projectData(projectKey: string) : Promise<any> {
  return await PackageFunctions.projectData(projectKey);
}

//input: string issueKey 
//output: object issue
export async function issueData(issueKey: string) : Promise<any> {
  return await PackageFunctions.issueData(issueKey);
}

//input: column<string> ticketColumn 
//input: string field 
//output: column result
//meta.vectorFunc: true
export async function getJiraField(ticketColumn: DG.Column, field: string) : Promise<any> {
  return await PackageFunctions.getJiraField(ticketColumn, field);
}

//input: object filter 
//output: object result
export async function getJiraTicketsByFilter(filter?: any) : Promise<any> {
  return await PackageFunctions.getJiraTicketsByFilter(filter);
}

//input: list issueIdsOrKeys 
//input: list fields 
//input: list expand 
//output: object result
export async function getJiraTicketsBulk(issueIdsOrKeys: string[], fields?: string[], expand?: string[]) : Promise<any> {
  return await PackageFunctions.getJiraTicketsBulk(issueIdsOrKeys, fields, expand);
}
