/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {GritApp} from './grit-app';
import {GritIssueHandler} from './grit-issue-handler';
import {gritDb} from './generated/db';
import '../css/grit.css';
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
  registerIssueHandleDetectors();
}

/** Registers a `<KEY>-\d+` detector per existing Grit project so that a real issue
 * handle (e.g. `GRITEST-1`) typed into global search is tagged as a `grit.issue`
 * and resolved by the handler (WO-31, ARCHITECTURE §8.6 discovery). Issue numbers
 * are per-project, so the fixed part is the project's own key, not a literal. */
async function registerIssueHandleDetectors(): Promise<void> {
  try {
    const projects = await gritDb.project.query({});
    for (const p of projects) {
      const key = (p.key ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
      if (key.length > 0)
        DG.SemanticValue.registerRegExpDetector('grit.issue', `${key}-\\d+`, `Grit issue (${p.key})`);
    }
  } catch (e) {
    console.error('Grit issue-handle detectors not registered:', e);
  }
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
