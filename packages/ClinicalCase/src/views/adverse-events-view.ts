import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { addTreatmentArm, getUniqueValues } from '../data-preparation/utils';
import { TREATMENT_ARM } from '../constants';

export class AdverseEventsView extends DG.ViewBase {

  aeWithArm: DG.DataFrame;

  constructor(name) {
    super(name);
    this.name = name;
    this.aeWithArm = addTreatmentArm(study.domains.ae.clone(), study.domains.dm, study.domains.ae.columns.names());

    let viewerTitle = {style:{
      'color':'var(--grey-6)',
      'margin':'12px 0px 6px 12px',
      'font-size':'16px',
    }};


    let typesPlot = this.bar('AEDECOD','Types', viewerTitle, TREATMENT_ARM);
    let bodySystemsPlot = this.bar('AEBODSYS', 'Body system', viewerTitle, TREATMENT_ARM);
    let causalityPlot = this.bar('AEREL', 'Causality', viewerTitle, TREATMENT_ARM);
    let outcomePlot = this.bar('AEOUT', 'Outcome', viewerTitle, TREATMENT_ARM);

    let timelinesPlot = this.aeWithArm.plot.line({
      x: 'week',
      yColumnNames: ['week'],
      chartTypes: ['Stacked Bar Chart'],
      yAggrTypes: [DG.STATS.TOTAL_COUNT],
      split: 'AESEV',
      style: 'dashboard' }).root;

    timelinesPlot.prepend(ui.divText('Events per week', viewerTitle));

    let scatterPlot = this.aeWithArm.plot.scatter({
      x: 'AESTDY',
      y: 'USUBJID',
      color: 'AESEV',
      markerDefaultSize: 5,
      legendVisibility: 'Never',
      style: 'dashboard'
    }).root;

    scatterPlot.prepend(ui.divText('All events', viewerTitle));

    let grid = this.aeWithArm.plot.grid();
    this.aeWithArm.onCurrentRowChanged.subscribe((_) => {
      grok.shell.o = new ClinRow(this.aeWithArm.currentRow);
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

  private bar(categoryColumn: string, title:string, viewerTitle: any, splitColumn: string) {
    let chart = this.aeWithArm.plot.bar( {
      split: categoryColumn,
      stack: splitColumn,
      style: 'dashboard',
      legendPosition: 'Top' }).root;
    chart.prepend(ui.divText(title, viewerTitle))
    return chart
  }
}