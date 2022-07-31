import {HitTriageBaseView} from "./hit-triage-base-view";
import * as ui from "datagrok-api/ui";
import {HitTriageApp} from "../hit-triage-app";

export class GetMoleculesView extends HitTriageBaseView {

  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Source Molecules';

    const from = ui.choiceInput('From', 'file', ['file', 'database', 'webservice']);
    const content = ui.divV([
      ui.divH([
          ui.divText('Ingest', {style: {'font-weight': 'bold'}}),
          from.root,
          ui.divText(this.template.sourceDescription)],
        {style: {'display': 'flex', 'align-items': 'center', 'gap': '12px'}}
      ),
      this.template.sourceDataFrame!.plot.grid()
    ])

    this.root.appendChild(content);
  }
}