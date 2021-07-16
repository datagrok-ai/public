import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";

export class AdverseEventsView extends DG.ViewBase {

  constructor() {
    super();

    function bar(categoryColumn: string) {
      return study.domains.ae.plot.bar( {
        split: categoryColumn,
        style: 'dashboard' }).root
    }

    let typesPlot = bar('AEDECOD');
    let bodySystemsPlot = bar('AEBODSYS');
    let causalityPlot = bar('AEREL');
    let outcomePlot = bar('AEOUT');

    let timelinesPlot = study.domains.ae.plot.line({
      x: 'week',
      yColumnNames: ['week'],
      chartTypes: ['Stacked Bar Chart'],
      yAggrTypes: [DG.STATS.TOTAL_COUNT],
      split: 'AESEV',
      style: 'dashboard' }).root;

    let scatterPlot = study.domains.ae.plot.scatter({
      x: 'AESTDY',
      y: 'USUBJID',
      color: 'AESEV',
      markerDefaultSize: 5,
      legendVisibility: 'Never',
      style: 'dashboard'
    }).root;

    let grid = study.domains.ae.plot.grid();
    study.domains.ae.onCurrentRowChanged.subscribe((_) => {
      grok.shell.o = new ClinRow(study.domains.ae.currentRow);
    });

    this.root.appendChild(ui.div([
      ui.divH([
        ui.block25([ui.h2('Types'), typesPlot]),
        ui.block25([ui.h2('Body System'), bodySystemsPlot]),
        ui.block50([ui.h2('Events per Week'), timelinesPlot]),
      ], { style: { width: '100%' } }),
      ui.divH([
        ui.block25([ui.h2('Causality'), causalityPlot]),
        ui.block25([ui.h2('Outcome'), outcomePlot]),
        ui.block50([ui.h2('All Events'), scatterPlot]),
      ], { style: { width: '100%' } }),
      ui.divH([ grid ], { style: { width: '100%' } })
    ]));
  }
}