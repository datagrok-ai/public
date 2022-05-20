import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { EventsView } from './views/events-view';
import { UsageWidget } from "./usage-widget";
import '../css/usage_analysis.css';
import {UaToolbox} from "./ua-toolbox";
import {ErrorsView} from "./views/errors-view";
import {UsersView} from "./views/users-view";
import {OverviewView} from "./views/overview-view";
import {DataView} from "./views/data-view";
import {ViewHandler as ViewHandler} from "./view-handler";
import {FunctionErrorsView} from "./views/function-errors-view";
import { UaView } from './views/ua-view';

export let _package = new DG.Package();

//name: Usage Analysis
//tags: app
export function usageAnalysisApp(): void {
  ViewHandler.getInstance().init();
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
