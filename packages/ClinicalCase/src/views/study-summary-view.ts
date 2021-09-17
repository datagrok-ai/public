import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { createKaplanMeierDataframe } from "../data-preparation/data-preparation";
import { createKaplanMeierScatterPlot } from "../custom-scatter-plots/custom-scatter-plots";


export class StudySummaryView extends DG.ViewBase {

 validationView: DG.ViewBase;
 errorsByDomain: any;
 errorsByDomainWithLinks: any;

 constructor(name) {
    super(name);

    this.name = name;

    const errorsMap = this.createErrorsMap();
    this.errorsByDomain = errorsMap.withCount;
    this.errorsByDomainWithLinks = errorsMap.withLinks;
    this.buildView();
  }

  async buildView() {
    let subjsPerDay = study.domains.dm.groupBy(['RFSTDTC']).uniqueCount('USUBJID').aggregate();
    let refStartCol = subjsPerDay.col('RFSTDTC');
    for (let i = 0, rowCount = subjsPerDay.rowCount; i < rowCount; i++) {
      if (refStartCol.isNone(i))
        subjsPerDay.rows.removeAt(i);
    }
    let lc = DG.Viewer.lineChart(subjsPerDay);
    if (refStartCol.type != DG.TYPE.DATE_TIME) {
      await subjsPerDay.columns.addNewCalculated('~RFSTDTC', 'Date(${RFSTDTC}, 1, 1)');
      lc.setOptions({ x: '~RFSTDTC', yColumnNames: ['unique(USUBJID)'] });
    }

    let summary = ui.tableFromMap({
      'subjects': study.subjectsCount,
      'sites': study.sitesCount
    });


    let errorsSummary = ui.tableFromMap(this.errorsByDomainWithLinks);

    let kaplanMeierDataframe = createKaplanMeierDataframe();
    let kaplanMeierPlot = createKaplanMeierScatterPlot(kaplanMeierDataframe, 'SUBJID', 'TIME', 'SURVIVAL', 'GROUP');

    let summaryStyle = {style:{
      'color':'var(--grey-6)',
      'margin-top':'8px',
      'font-size':'16px',
    }};
    let viewerTitle = {style:{
      'color':'var(--grey-6)',
      'margin':'12px 0px 6px 12px',
      'font-size':'16px',
    }};
    
    lc.root.prepend(ui.divText('Enrollment by day', viewerTitle));

    let arm = study.domains.dm.plot.bar( { split: 'arm', style: 'dashboard' });
    arm.root.prepend(ui.divText('Arm treatment', viewerTitle));

    let sex = DG.Viewer.barChart(study.domains.dm, { split: 'sex', style: 'dashboard' });
    sex.root.prepend(ui.divText('Sex', viewerTitle));

    let race = DG.Viewer.barChart(study.domains.dm, { split: 'race', style: 'dashboard' });
    race.root.prepend(ui.divText('Race', viewerTitle));

    let age = DG.Viewer.histogram(study.domains.dm, { value: 'age', style: 'dashboard' });
    age.root.prepend(ui.divText('Age', viewerTitle));

    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      ui.splitH([
        ui.panel([
          ui.divText('Summary', summaryStyle),
          summary
        ]),
        ui.panel([
          ui.divText('Errors', summaryStyle),
          errorsSummary
        ]),
      ], {style:{maxHeight:'105px'}}),
      lc.root,
      ui.splitH([
        arm.root,
        sex.root,
        race.root,
        age.root
      ])
    ]))
  }

  private createErrorsMap() {
    const errorsMap = {};
    const errorsMapWithCount = {};
    const validationSummary = study.validationResults.groupBy([ 'Domain' ]).count().aggregate();
    for (let i = 0; i < validationSummary.rowCount; ++i) {
      const domain = validationSummary.get('Domain', i);
      const errorsCount = validationSummary.get('count', i);
      const link = ui.link(errorsCount, {}, '', {id: domain});
      link.addEventListener('click', (event) => {
        grok.shell.v = this.validationView;
        event.stopPropagation();
      });
      errorsMap[ domain ] = link;
      errorsMapWithCount[ domain ] = errorsCount;
    }
    return {withCount: errorsMapWithCount, withLinks: errorsMap};
  }
}
