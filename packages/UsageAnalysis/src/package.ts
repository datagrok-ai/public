import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { UsageAnalysisView } from './usage-view';
import { UsageWidget } from "./usage-widget";
import '../css/usage_analysis.css';

export let _package = new DG.Package();

//name: UsageAnalysis
//tags: app
export function usageAnalysisApp(): void {
  grok.shell.addView(new UsageAnalysisView('Usage Analysis'));
}

//output: widget result
//tags: dashboard
export function usageWidget(): DG.Widget {
  return new UsageWidget();
}

//name: Create JIRA ticket
//description: Creates JIRA ticket using current error log
//tags: panel, widgets
//input: string msg {semType: ErrorMessage}
//output: widget result
//condition: true
export function createJiraTicket(msg:string): DG.Widget {
  let root = ui.div();

  let summary = ui.stringInput('Summary', '');
  let description = ui.stringInput('Description', msg);

  let button = ui.bigButton('CREATE', () => {
    grok.data.query('Vnerozin:JiraCreateIssue', {
      'createRequest': JSON.stringify({
        "fields": {
          "project": {
            "key": "GROK"
          },
          "summary": summary.value,
          "description": description.value,
          "issuetype": {
            "name": "Bug"
          }
        }
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
