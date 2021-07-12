import { ClinicalCaseView } from "../clinical-case-view";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import {InputBase, STATS} from "datagrok-api/dg";

export class AdverseEventsView extends DG.ViewBase {

  constructor() {
    super();

    let typesPlot = study.domains.ae.plot.bar( {
      split: 'AEDECOD',
      style: 'dashboard' }).root;

    let bodySystemsPlot = study.domains.ae.plot.bar( {
      split: 'AEBODSYS',
      style: 'dashboard' }).root;

    let timelinesPlot = study.domains.ae.plot.line({
      x: 'week',
      yColumnNames: ['week'],
      chartTypes: ['Stacked Bar Chart'],
      yAggrTypes: [DG.STATS.TOTAL_COUNT],
      split: 'AESEV',
      style: 'dashboard' }).root;

    let grid = study.domains.ae.plot.grid();

    this.root.appendChild(ui.div([
      ui.divH([
        ui.block25([ui.h2('Types'), typesPlot]),
        ui.block25([ui.h2('Body System'), bodySystemsPlot]),
        ui.block50([ui.h2('Events per Week'), timelinesPlot]),
      ], { style: { width: '100%' } }),
      ui.divH([ grid ], { style: { width: '100%' } })
    ]));
  }
}