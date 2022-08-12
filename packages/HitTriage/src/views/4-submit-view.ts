import {HitTriageBaseView} from "./hit-triage-base-view";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {HitTriageApp} from "../hit-triage-app";
import {_package} from "../package";

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
      ui.divH([ui.bigButton('SUBMIT', () => this.submit())])
    ])
    this.root.appendChild(content);
  }

  async submit(): Promise<any> {
    const d = new Date();
    let time = `${d.getFullYear()}-${d.getMonth()}-${d.getDate()}_${d.getHours()}-${d.getMinutes()}-${d.getSeconds()}`;
    let folder = `${time}_${grok.shell.user.login}`;
    await _package.files.writeAsText(`${folder}/session.json`, JSON.stringify(this.app.template));

    await _package.files.writeAsText(`${folder}/molecules.csv`, this.app.template.enrichedTable!.toCsv());

    grok.shell.info('Submitted successfully.')
  }
}