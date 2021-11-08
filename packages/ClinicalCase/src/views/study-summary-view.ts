import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { cumulativeEnrollemntByDay } from "../data-preparation/data-preparation";
import { CLINICAL_TRIAL_GOV_FIELDS, requiredColumnsByView } from "../constants";
import { HttpService } from "../services/http.service";
import { ILazyLoading } from "../lazy-loading/lazy-loading";
import { checkMissingDomains } from "./utils";
import { _package } from "../package";
import { AGE, RACE, SEX, STUDY_ID, SUBJECT_ID, SUBJ_REF_STDT, TREATMENT_ARM } from "../columns-constants";


export class StudySummaryView extends DG.ViewBase implements ILazyLoading {

 validationView: DG.ViewBase;
 errorsByDomain: any;
 errorsByDomainWithLinks: any;
 studyId: string;
 httpService = new HttpService();

 constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/summary.md`;
    this.path = '/summary';
  }
  
  loaded: boolean;

  load(): void {
    checkMissingDomains(requiredColumnsByView[this.name], false, this);
  }

  createView(){
    this.studyId = study.domains.dm.get(STUDY_ID, 0);
    const errorsMap = this.createErrorsMap();
    this.errorsByDomain = errorsMap.withCount;
    this.errorsByDomainWithLinks = errorsMap.withLinks;
    this.buildView();
  }

  async buildView() {
    let dateCol = SUBJ_REF_STDT;
    let cumulativeCol = 'CUMULATIVE_ENROLLMENT';
    let subjsPerDay = cumulativeEnrollemntByDay(study.domains.dm, dateCol, SUBJECT_ID, cumulativeCol)
    let refStartCol = subjsPerDay.col(dateCol);
    let lc = DG.Viewer.lineChart(subjsPerDay);
    if (refStartCol.type != DG.TYPE.DATE_TIME) {
      await subjsPerDay.columns.addNewCalculated(`~${dateCol}`, `Date(\${${dateCol}}, 1, 1)`);
      lc.setOptions({ x: `~${dateCol}`, yColumnNames: [cumulativeCol] });
    } else {
      lc.setOptions({ x: `${dateCol}`, yColumnNames: [cumulativeCol] });
    }

    let summary = ui.tableFromMap({
      'subjects': study.subjectsCount,
      'sites': study.sitesCount
    });


    let errorsSummary = ui.tableFromMap(this.errorsByDomainWithLinks);

    let summaryStyle = {style:{
      'color':'var(--grey-6)',
      'margin-top':'8px',
      'margin-left': '5px',
      'font-size':'16px',
    }};

    let summaryLinkStyle = {style:{
      'margin-top':'8px',
      'font-size':'16px',
    }};

    let viewerTitle = {style:{
      'color':'var(--grey-6)',
      'margin':'12px 0px 6px 12px',
      'font-size':'16px',
    }};
    
    lc.root.prepend(ui.divText('Enrollment by day', viewerTitle));

    let arm = DG.Viewer.barChart(study.domains.dm, { split: TREATMENT_ARM, style: 'dashboard', barColor: DG.Color.lightBlue });
    arm.root.prepend(ui.divText('Treatment arm', viewerTitle));

    let sex = DG.Viewer.barChart(study.domains.dm, { split: SEX, style: 'dashboard' });
    sex.root.prepend(ui.divText('Sex', viewerTitle));

    let race = DG.Viewer.barChart(study.domains.dm, { split: RACE, style: 'dashboard' });
    race.root.prepend(ui.divText('Race', viewerTitle));

    let age = DG.Viewer.histogram(study.domains.dm, { value: AGE, style: 'dashboard' });
    age.root.prepend(ui.divText('Age', viewerTitle));

    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      ui.splitH([
        ui.panel([
          ui.divH([ui.link(`${this.studyId}`, ()=>{
            this.getStudyInfoFromClinTrialsGov(this.studyId);
          }, '', summaryLinkStyle), 
          ui.divText('summary', summaryStyle)]),
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

  private async getStudyInfoFromClinTrialsGov(studyId: string){
    //http.getStudyData(studyId);
    let result = await this.httpService.getStudyData('R01NS050536', Object.keys(CLINICAL_TRIAL_GOV_FIELDS));
    ui.dialog({ title: 'Study info' })
              .add(ui.tableFromMap(result))
              .onOK(() => {})
              .show();
  }

}
