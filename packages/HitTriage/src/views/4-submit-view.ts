import {HitTriageBaseView} from "./hit-triage-base-view";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {HitTriageApp} from "../hit-triage-app";

export class SubmitView extends HitTriageBaseView {

  constructor(app: HitTriageApp) {
    super(app);
    this.render();
  }

  render(): void {
    ui.empty(this.root);
    const content = ui.divV([
      ui.h1('Summary'),
      ui.div([ui.tableFromMap(this.template.getSummary())]),
      ui.divH([ui.bigButton('SUBMIT', () => grok.shell.info('Good job!'))])
    ])
    this.root.appendChild(content);
  }
}