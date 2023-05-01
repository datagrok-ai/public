import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UsageWidget} from './usage-widget';
import '../css/usage_analysis.css';
import {ViewHandler} from './view-handler';

export const _package = new DG.Package();


//name: Usage Analysis
//tags: app
export function usageAnalysisApp(): void {
  ViewHandler.getInstance();
  if (!grok.shell.view(ViewHandler.UAname)) ViewHandler.getInstance().init();
}

//input: dynamic header
//output: widget result
//tags: dashboard
export function usageWidget(header: HTMLDivElement): DG.Widget {
  return new UsageWidget(header);
}

//name: Create JIRA ticket
//description: Creates JIRA ticket using current error log
//tags: panel, widgets
//input: string msg {semType: ErrorMessage}
//output: widget result
//condition: true
export function createJiraTicket(msg:string): DG.Widget {
  const root = ui.div();

  const summary = ui.stringInput('Summary', '');
  const description = ui.stringInput('Description', msg);

  const button = ui.bigButton('CREATE', () => {
    grok.data.query('Vnerozin:JiraCreateIssue', {
      'createRequest': JSON.stringify({
        'fields': {
          'project': {
            'key': 'GROK',
          },
          'summary': summary.value,
          'description': description.value,
          'issuetype': {
            'name': 'Bug',
          },
        },
      }),
      'updateHistory': false,
    }).then((t) => {
      grok.shell.info('Created');
    });
  });
  button.style.marginTop = '12px';

  root.appendChild(ui.inputs([summary, description]));
  root.appendChild(button);

  return new DG.Widget(root);
}
