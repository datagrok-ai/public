/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {GritApp} from './grit-app';
import {GritIssueHandler} from './grit-issue-handler';
export * from './package.g';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: autostart
//meta.autostartImmediate: true
export function _initGrit(): void {
  DG.ObjectHandler.register(new GritIssueHandler());
}

//name: gritIssueHandler
//description: Custom rendering for grit.issue domain rows (status/priority badges, timeline)
//tags: objectHandler, objectHandler-grit.issue
//output: object handler
export function gritIssueHandler(): DG.ObjectHandler<DG.DomainRow> {
  return new GritIssueHandler();
}

//name: Grit
//description: GRok Issue Tracker — projects, issues, comments over entity-mapped domain schemas
//tags: app
//output: view result
export async function gritApp(): Promise<DG.ViewBase> {
  return await GritApp.run();
}
