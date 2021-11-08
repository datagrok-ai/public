import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { addDataFromDmDomain } from '../data-preparation/utils';
import { AE_BODY_SYSTEM, AE_CAUSALITY, AE_DECOD_TERM, AE_OUTCOME, AE_SEVERITY, AE_START_DAY, SUBJECT_ID, TREATMENT_ARM } from '../columns-constants';
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains } from './utils';
import { _package } from '../package';
import { requiredColumnsByView } from '../constants';


export class AdverseEventsView extends DG.ViewBase implements ILazyLoading {

  aeWithArm: DG.DataFrame;

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/adverse_events.md`;
    //@ts-ignore
    this.basePath = '/adverse-events';
  }

  loaded: boolean;

  load(): void {
    checkMissingDomains(requiredColumnsByView[this.name], false, this);
 }

  createView(): void {
    this.aeWithArm = addDataFromDmDomain(study.domains.ae.clone(), study.domains.dm, study.domains.ae.columns.names(), [TREATMENT_ARM]);
    let viewerTitle = {style:{
      'color':'var(--grey-6)',
      'margin':'12px 0px 6px 12px',
      'font-size':'16px',
    }};

    let typesPlot = this.bar(AE_DECOD_TERM,'Types', viewerTitle, TREATMENT_ARM);
    let bodySystemsPlot = this.bar(AE_BODY_SYSTEM, 'Body system', viewerTitle, TREATMENT_ARM);
    let causalityPlot = this.bar(AE_CAUSALITY, 'Causality', viewerTitle, TREATMENT_ARM);
    let outcomePlot = this.bar(AE_OUTCOME, 'Outcome', viewerTitle, TREATMENT_ARM);

    let timelinesPlot = this.aeWithArm.plot.line({
      x: 'week',
      yColumnNames: ['week'],
      chartTypes: ['Stacked Bar Chart'],
      yAggrTypes: [DG.STATS.TOTAL_COUNT],
      split: AE_SEVERITY,
      style: 'dashboard' }).root;

    timelinesPlot.prepend(ui.divText('Events per week', viewerTitle));

    let scatterPlot = this.aeWithArm.plot.scatter({
      x: AE_START_DAY,
      y: SUBJECT_ID,
      color: AE_SEVERITY,
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
    legend.root.style.width = '500px';
    legend.root.style.height = '35px';

    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      ui.splitH([timelinesPlot, scatterPlot]),
      ui.divV([
        ui.div(legend.root, { style: { height: '50px' } }),
        ui.splitH([typesPlot,bodySystemsPlot,causalityPlot,outcomePlot], {style: {width: '100%'}})
      ]),
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