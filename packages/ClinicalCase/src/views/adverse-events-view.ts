import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";

export class AdverseEventsView extends DG.ViewBase {

  constructor() {
    super();

    let viewerTitle = {style:{
      'color':'var(--grey-6)',
      'margin':'12px 0px 6px 12px',
      'font-size':'16px',
    }};

    function bar(categoryColumn: string, title:string) {
      let chart = study.domains.ae.plot.bar( {
        split: categoryColumn,
        style: 'dashboard' }).root;
      chart.prepend(ui.divText(title, viewerTitle))
      return chart
    }

    let typesPlot = bar('AEDECOD','Types');
    let bodySystemsPlot = bar('AEBODSYS', 'Body System');
    let causalityPlot = bar('AEREL', 'Causality');
    let outcomePlot = bar('AEOUT', 'Outcome');

    let timelinesPlot = study.domains.ae.plot.line({
      x: 'week',
      yColumnNames: ['week'],
      chartTypes: ['Stacked Bar Chart'],
      yAggrTypes: [DG.STATS.TOTAL_COUNT],
      split: 'AESEV',
      style: 'dashboard' }).root;

    timelinesPlot.prepend(ui.divText('Events per Week', viewerTitle));

    let scatterPlot = study.domains.ae.plot.scatter({
      x: 'AESTDY',
      y: 'USUBJID',
      color: 'AESEV',
      markerDefaultSize: 5,
      legendVisibility: 'Never',
      style: 'dashboard'
    }).root;

    scatterPlot.prepend(ui.divText('All Events', viewerTitle))

    let grid = study.domains.ae.plot.grid();
    study.domains.ae.onCurrentRowChanged.subscribe((_) => {
      grok.shell.o = new ClinRow(study.domains.ae.currentRow);
    });

    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      ui.splitH([timelinesPlot, scatterPlot]),
      ui.splitH([typesPlot,bodySystemsPlot,causalityPlot,outcomePlot]),
      ui.splitH([grid.root])
    ]))
    /*
    this.root.appendChild(ui.block([
        ui.block25([typesPlot]),
        ui.block25([bodySystemsPlot]),
        ui.block50([timelinesPlot]),
        ui.block25([causalityPlot]),
        ui.block25([outcomePlot]),
        ui.block50([scatterPlot]),
        ui.block([grid.root])
    ]));
    */
  }
}