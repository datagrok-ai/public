import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

export async function addJiraDetector(projects: string[], projectDescription:string = ''){
    DG.SemanticValue.registerRegExpDetector(`Jira Ticket`, `^(${(projects.map(e=>`${e}`)).join('|')})-\\d+$`, `${projectDescription}`)
}