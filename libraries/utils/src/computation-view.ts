/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {FunctionView} from './function-view';

/**
 * Base class for handling Compute models (see https://github.com/datagrok-ai/public/blob/master/help/compute/compute.md).
 * In most cases, the computation is a simple {@link Func}
 * Extend it in cases where a behavior or UI not supported by the {@link FunctionView} is needed.
 *
 * It provides the following functionality out-of-the-box, where each section could be customized:
 * - a structured way to represent input and output parameters: {@link parameters}
 * - generic way to generate UI for inputs, outputs, and interactivity (running the model, etc)
 *   - persisting historical results to the db (via {@link parameters})
 * - export (to Excel and PDF): {@link export}
 * - easy loading of historical runs
 * - routing
 * - entering the real, measured (as opposed to predicted) values manually
 * - notifications for changing inputs, completion of computations, etc: {@link onInputChanged}
 * */
export class ComputationView extends FunctionView {
  jiraUrl: string = '';

  constructor(funccall: DG.FuncCall) {
    super(funccall);
  }

  override async init(funccall: DG.FuncCall) {
    super.init(funccall);

    this.buildRibbonMenu();
  }

  buildRibbonMenu() {
    this.ribbonMenu
      .group('Model')
      .group('Export')
      .items(this.supportedExportFormats, async (format: string) => DG.Utils.download(this.exportFilename(format), await this.export(format)))
      .endGroup()
      .items([
        ...this.jiraUrl ? ['Report a bug']: [],
        ...this.helpUrl ? ['Get help']: [],
      ], (item: string) => {
        switch (item) {
        case 'Report a bug':
          window.open(`https://${this.jiraUrl}/CreateIssueDetails!init.jspa?pid=27326&issuetype=1&priority=3&assignee=AZakiev&summary=1651759426359:+user+bug+report&description=System+error:+${grok.shell.lastError}`);
          break;
        case 'Help':
          window.open(this.helpUrl!);
        }
      })
      .endGroup();
  }
}
