import {HitTriageBaseView} from "./hit-triage-base-view";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import * as DG from 'datagrok-api/dg';
import {HitTriageApp} from "../hit-triage-app";

export class GetMoleculesView extends HitTriageBaseView {

  hitsGrid: DG.Grid;

  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Hit Triage | Ingest';

    this.hitsGrid = this.template.hitsTable!.plot.grid();

    const from = ui.choiceInput('From', 'file', ['file', 'database', 'webservice']);
    const content = ui.divV([
      ui.divH([
          ui.divText('Hits', {style: {'font-weight': 'bold'}}),
          from.root,
          ui.divText(this.template.hitsQuery)],
        {style: {'display': 'flex', 'align-items': 'center', 'gap': '12px'}}
      ),
      ui.divH([
          ui.divText('Targets', {style: {'font-weight': 'bold'}}),
          ui.choiceInput('From', 'file', ['file', 'database', 'webservice']).root,
          ui.divText(this.template.hitsTargetsQuery)],
        {style: {'display': 'flex', 'align-items': 'center', 'gap': '12px'}}
      ),
      this.hitsGrid
    ])

    this.root.appendChild(content);

    this.template.loadData().then((_) => {
      grok.shell.info('Data loaded.');
      this.hitsGrid.dataFrame = this.template.hitsTable!;
    });
  }
}