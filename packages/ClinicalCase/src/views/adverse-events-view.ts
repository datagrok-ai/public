import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { addDataFromDmDomain, getUniqueValues } from '../data-preparation/utils';
import { TREATMENT_ARM } from '../constants';
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkDomainExists } from './utils';


export class AdverseEventsView extends DG.ViewBase implements ILazyLoading {

  aeWithArm: DG.DataFrame;

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/views_help/adverse_events.md';
    //@ts-ignore
    this.basePath = '/adverse-events';
  }

  loaded: boolean;

  load(): void {
    checkDomainExists(['dm', 'ae'], false, this);
 }

  createView(): void {
    this.aeWithArm = addDataFromDmDomain(study.domains.ae.clone(), study.domains.dm, study.domains.ae.columns.names(), [TREATMENT_ARM]);
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


    let legend = DG.Legend.create(this.aeWithArm.columns.byName(TREATMENT_ARM));

    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      ui.splitH([timelinesPlot, scatterPlot]),
      ui.box(legend.root, { style: { maxHeight: '50px' } }),
      ui.splitH([typesPlot,bodySystemsPlot,causalityPlot,outcomePlot]),
      ui.splitH([grid.root])
    ]))

  }

  private bar(categoryColumn: string, title:string, viewerTitle: any, splitColumn: string) {
    let chart = this.aeWithArm.plot.bar( {
      split: categoryColumn,
      stack: splitColumn,
      style: 'dashboard',
      legendVisibility: 'Never'}).root;
    chart.prepend(ui.divText(title, viewerTitle))
    return chart
  }
}