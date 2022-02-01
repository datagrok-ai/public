import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { cumulativeEnrollemntByDay } from "../data-preparation/data-preparation";
import { CLINICAL_TRIAL_GOV_FIELDS } from "../constants";
import { CLIN_TRIAL_GOV_SEARCH, HttpService } from "../services/http.service";
import { _package } from "../package";
import { AGE, RACE, SEX, STUDY_ID, SUBJECT_ID, SUBJ_REF_STDT, TREATMENT_ARM } from "../columns-constants";
import { ClinicalCaseViewBase } from "../model/ClinicalCaseViewBase";
import $ from "cash-dom";
import { checkDateFormat } from "../data-preparation/utils";
import { checkColumnsAndCreateViewer, checkRequiredColumns, updateDivInnerHTML } from "./utils";


export class StudySummaryView extends ClinicalCaseViewBase {

  validationView: DG.ViewBase;
  errorsByDomain: any;
  errorsByDomainWithLinks: any;
  studyId: string;
  httpService = new HttpService();
  enrollmentChart = ui.box();
  armChart = ui.box();
  sexChart = ui.box();
  raceChart = ui.box();
  ageChart = ui.box();

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/summary.md`;
    this.path = '/summary';
  }

  createView() {
    this.studyId = study.domains.dm.get(STUDY_ID, 0);
    const errorsMap = this.createErrorsMap();
    if (errorsMap) {
      this.errorsByDomain = errorsMap.withCount;
      this.errorsByDomainWithLinks = errorsMap.withLinks;
    }
    this.buildView();
  }

  async buildView() {

    let viewerTitle = {
      style: {
        'color': 'var(--grey-6)',
        'margin': '12px 0px 6px 12px',
        'font-size': '16px',
      }
    };
    checkColumnsAndCreateViewer(
      study.domains.dm,
      [SUBJ_REF_STDT],
      this.enrollmentChart,
      async () => { await this.createCumulativeEnrollmentChart(viewerTitle) },
      'Cumulative enrollment');

    let summary = ui.tableFromMap({
      'subjects': study.subjectsCount,
      'sites': study.sitesCount
    });


    let errorsSummary = this.errorsByDomainWithLinks ?
      ui.tableFromMap(this.errorsByDomainWithLinks) :
      ui.divText('No errors found', { style: { marginTop: '8px', marginLeft: '5px' } });

    let summaryStyle = {
      style: {
        'color': 'var(--grey-6)',
        'margin-top': '8px',
        'margin-left': '5px',
        'font-size': '16px',
      }
    };

    checkColumnsAndCreateViewer(
      study.domains.dm,
      [TREATMENT_ARM],
      this.armChart, () => {
        let arm = DG.Viewer.barChart(study.domains.dm, { split: TREATMENT_ARM, style: 'dashboard', barColor: DG.Color.lightBlue });
        arm.root.prepend(ui.divText('Treatment arm', viewerTitle));
        updateDivInnerHTML(this.armChart, arm.root);
      },
      'Treatment arm');

    checkColumnsAndCreateViewer(
      study.domains.dm,
      [SEX],
      this.sexChart, () => {
        let sex = DG.Viewer.barChart(study.domains.dm, { split: SEX, style: 'dashboard' });
        sex.root.prepend(ui.divText('Sex', viewerTitle));
        updateDivInnerHTML(this.sexChart, sex.root);
      },
      'Sex');

    checkColumnsAndCreateViewer(
      study.domains.dm,
      [RACE],
      this.raceChart, () => {
        let race = DG.Viewer.barChart(study.domains.dm, { split: RACE, style: 'dashboard' });
        race.root.prepend(ui.divText('Race', viewerTitle));
        updateDivInnerHTML(this.raceChart, race.root);
      },
      'Race');

    checkColumnsAndCreateViewer(
      study.domains.dm,
      [AGE],
      this.ageChart, () => {
        let age = DG.Viewer.histogram(study.domains.dm, { value: AGE, style: 'dashboard' });
        age.root.prepend(ui.divText('Age', viewerTitle));
        updateDivInnerHTML(this.ageChart, age.root);
      },
      'Age');

    this.root.className = 'grok-view ui-box';
    this.root.append(ui.splitV([
      ui.splitH([
        ui.panel([
          ui.divText(`${this.studyId} summary`, summaryStyle),
          summary
        ]),
        ui.panel([
          ui.divText('Errors', summaryStyle),
          errorsSummary
        ]),
      ], { style: { maxHeight: '105px' } }),
      this.enrollmentChart,
      ui.splitH([
        this.armChart,
        this.sexChart,
        this.raceChart,
        this.ageChart
      ])
    ]))
  }

  private async createCumulativeEnrollmentChart(viewerTitle: any) {
    let dateCol = SUBJ_REF_STDT;
    let cumulativeCol = 'CUMULATIVE_ENROLLMENT';
    let incorrectDates = checkDateFormat(study.domains.dm.getCol(dateCol), study.domains.dm.rowCount);
    if (incorrectDates.length) {
      let subjArray = incorrectDates.map(it => study.domains.dm.get(SUBJECT_ID, it));
      let formatErrorMessage = ui.info(`Subjects #${subjArray.join(',')} have incorrect format of ${SUBJ_REF_STDT}`);
      updateDivInnerHTML(this.enrollmentChart, formatErrorMessage);

    } else {
      let subjsPerDay = cumulativeEnrollemntByDay(study.domains.dm, dateCol, SUBJECT_ID, cumulativeCol);
      let refStartCol = subjsPerDay.col(dateCol);
      let lc = DG.Viewer.lineChart(subjsPerDay);
      lc.setOptions({ x: `${dateCol}`, yColumnNames: [cumulativeCol] });
      if (refStartCol.type != DG.TYPE.DATE_TIME) {
        await subjsPerDay.columns.addNewCalculated(`~${dateCol}`, `Date(\${${dateCol}}, 1, 1)`);
        lc.setOptions({ x: `~${dateCol}`, yColumnNames: [cumulativeCol] });
      } else {
        lc.setOptions({ x: `${dateCol}`, yColumnNames: [cumulativeCol] });
      }
      updateDivInnerHTML(this.enrollmentChart, lc.root);
      lc.root.prepend(ui.divText('Enrollment by day', viewerTitle));
    }
  }

  private createErrorsMap() {
    const errorsMap = {};
    const errorsMapWithCount = {};
    if (study.validationResults.rowCount) {
      const validationSummary = study.validationResults.groupBy(['Domain']).count().aggregate();
      for (let i = 0; i < validationSummary.rowCount; ++i) {
        const domain = validationSummary.get('Domain', i);
        const errorsCount = validationSummary.get('count', i);
        const link = ui.link(errorsCount, {}, '', { id: domain });
        link.addEventListener('click', (event) => {
          grok.shell.v = this.validationView;
          event.stopPropagation();
        });
        errorsMap[domain] = link;
        errorsMapWithCount[domain] = errorsCount;
      }
      return { withCount: errorsMapWithCount, withLinks: errorsMap };
    }
    return null;
  }

  override async propertyPanel() {
    const httpService = new HttpService();
    //let clinTrialsGovInfo = await httpService.getStudyData('R01NS050536', Object.keys(CLINICAL_TRIAL_GOV_FIELDS));
    let clinTrialsGovInfo = await httpService.getStudyData(this.studyId, Object.keys(CLINICAL_TRIAL_GOV_FIELDS));

    let acc = this.createAccWithTitle(this.studyId)

    if (clinTrialsGovInfo) {
      const summaryDict = {};
      Object.keys(clinTrialsGovInfo).forEach(key => {
        summaryDict[CLINICAL_TRIAL_GOV_FIELDS[key]] = clinTrialsGovInfo[key];
      })
      let studyLink = `${CLIN_TRIAL_GOV_SEARCH}${summaryDict['NCT ID']}`
      summaryDict[`Study link`] = ui.link('Go to study page', () => { window.open(studyLink, '_blank').focus(); })

      let acctable = ui.tableFromMap(summaryDict);
      acc.addPane('General', () => {
        $(acctable).find('tr').css('vertical-align', 'top');
        $(acctable).find('td').css('padding-bottom', '10px');
        $(acctable).find('.d4-entity-list>span').css('margin', '0px');
        return acctable
      }, true)
    } else {
      acc.addPane('General', () => {
        return ui.divText('Study not found on clinicaltrials.gov')
      }, true)
    }
    return acc.root;
  }
}
