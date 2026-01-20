import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {cumulativeEnrollemntByDay} from '../data-preparation/data-preparation';
import {CLINICAL_TRIAL_GOV_FIELDS} from '../constants/constants';
import {CLIN_TRIAL_GOV_SEARCH, HttpService} from '../services/http.service';
import {_package} from '../package';
import {AGE, BW_RES_N, RACE, SEX, SPECIES, SUBJECT_ID, SUBJ_REF_STDT,
  VISIT_DAY} from '../constants/columns-constants';
import {ClinicalCaseViewBase} from '../model/ClinicalCaseViewBase';
import $ from 'cash-dom';
import {checkDateFormat} from '../data-preparation/utils';
import {removeExtension, studyConfigToMap, updateDivInnerHTML} from '../utils/utils';
import {TRT_ARM_FIELD} from '../views-config';
import {checkColumnsAndCreateViewer} from '../utils/views-validation-utils';
import {CDISC_STANDARD} from '../utils/types';
import {addDomainAsTableView, studies} from '../utils/app-utils';


export class StudySummaryView extends ClinicalCaseViewBase {
  validationView: DG.ViewBase;
  errorsByDomain: any;
  errorsByDomainWithLinks: any;
  studyId: string;
  httpService = new HttpService();
  centralChart = ui.box();
  armChart = ui.box();
  sexChart = ui.box();
  raceChart = ui.box();
  ageOrWeightChart = ui.box();
  viewerTitle = {
    style: {
      'color': 'var(--grey-6)',
      'margin': '12px 0px 6px 12px',
      'font-size': '16px',
    },
  };
  test: any;
  validationErrorLinkHandler;
  isSend = false;

  constructor(name: string, studyId: string, errorLinkHandler?: () => void) {
    super(name, studyId);
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/summary.md`;
    this.path = '/summary';
    if (errorLinkHandler)
      this.validationErrorLinkHandler = errorLinkHandler;
  }

  createView() {
    this.isSend = studies[this.studyId].config.standard === CDISC_STANDARD.SEND;
    this.errorsByDomain = studies[this.studyId].errorsByDomain;
    this.errorsByDomainWithLinks = this.createErrorsMapWithLinks();
    this.buildView();
  }

  async buildView() {
    const summaryObj = {
      'subjects': studies[this.studyId].subjectsCount,
    };
    if (!this.isSend)
      summaryObj['sites'] = studies[this.studyId].sitesCount;

    const summary = ui.tableFromMap(summaryObj);

    const errorsSummary = this.errorsByDomainWithLinks ?
      ui.tableFromMap(this.errorsByDomainWithLinks) :
      ui.divText('No errors found', {style: {marginTop: '8px', marginLeft: '5px'}});

    const summaryStyle = {
      style: {
        'color': 'var(--grey-6)',
        'margin-top': '8px',
        'margin-left': '5px',
        'font-size': '16px',
      },
    };

    const addDomainToWorkspace = ui.icons.add(() => {
      const menu = DG.Menu.popup();
      for (const domain of studies[this.studyId].domains.all()) {
        const domainName = removeExtension(domain.name);
        menu.item(domainName, () => {
          addDomainAsTableView(this.studyId, domain, false);
        });
      }
      menu.show();
    }, 'Add domain to workspace');
    addDomainToWorkspace.classList.add('clinical-case-add-domain-to-workspace-icon');


    if (!this.isSend) { //in case of SDTM - show cumulative enrollment
      checkColumnsAndCreateViewer(
        studies[this.studyId].domains.dm,
        [SUBJ_REF_STDT],
        this.centralChart,
        async () => {await this.createCumulativeEnrollmentChart(this.viewerTitle);},
        'Cumulative enrollment');
    } else { //if SEND - show study summary
      const configMap = studyConfigToMap(studies[this.studyId].config);
      const summaryDiv = ui.divV([]);
      summaryDiv.append(ui.h1('Study summary',
        {style: {fontSize: '15px', marginLeft: '12px', position: 'absolute', width: '98%', background: 'white'}}));
      const summaryTableDiv = ui.div(ui.tableFromMap(configMap), {style: {marginTop: '15px'}});
      summaryDiv.append(summaryTableDiv);
      updateDivInnerHTML(this.centralChart, summaryDiv);
    }

    checkColumnsAndCreateViewer(
      studies[this.studyId].domains.dm,
      [studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD]],
      this.armChart, () => {
        const arm = DG.Viewer.barChart(studies[this.studyId].domains.dm,
          {split: studies[this.studyId].viewsConfig.config[this.name][TRT_ARM_FIELD], //@ts-ignore
            style: 'dashboard', barColor: DG.Color.lightBlue});
        arm.root.prepend(ui.divText('Treatment arm', this.viewerTitle));
        updateDivInnerHTML(this.armChart, arm.root);
      },
      'Treatment arm');

    checkColumnsAndCreateViewer(
      studies[this.studyId].domains.dm,
      [SEX],
      this.sexChart, () => {//@ts-ignore
        const sex = DG.Viewer.barChart(studies[this.studyId].domains.dm, {split: SEX, style: 'dashboard'});
        sex.root.prepend(ui.divText('Sex', this.viewerTitle));
        updateDivInnerHTML(this.sexChart, sex.root);
      },
      'Sex');

    checkColumnsAndCreateViewer(
      studies[this.studyId].domains.dm,
      [this.isSend ? SPECIES : RACE],
      this.raceChart, () => {//@ts-ignore
        const race = DG.Viewer.barChart(studies[this.studyId].domains.dm, //@ts-ignore
          {split: this.isSend ? SPECIES : RACE, style: 'dashboard'});
        race.root.prepend(ui.divText(this.isSend ? 'Species' : 'Race', this.viewerTitle));
        updateDivInnerHTML(this.raceChart, race.root);
      },
      this.isSend ? 'Species' : 'Race');

    //in case of SDTM - show age distribution
    if (!this.isSend) {
      checkColumnsAndCreateViewer(
        studies[this.studyId].domains.dm,
        [AGE],
        this.ageOrWeightChart, () => { //@ts-ignore
          const age = DG.Viewer.histogram(studies[this.studyId].domains.dm, {value: AGE, style: 'dashboard'});
          age.root.prepend(ui.divText('Age', this.viewerTitle));
          updateDivInnerHTML(this.ageOrWeightChart, age.root);
        },
        'Age');
    } else { //otherwise show initila body weight distribution
      //in some cases VISITDY column is missing - look for domain specific day column
      const visitDayCol = studies[this.studyId].domains.bw?.col(VISIT_DAY) ??
        studies[this.studyId].domains.bw?.col(`BWDY`);
      const visitDayName = visitDayCol ? visitDayCol.name : VISIT_DAY;
      checkColumnsAndCreateViewer(
        studies[this.studyId].domains.bw,
        [visitDayName, SUBJECT_ID, BW_RES_N],
        this.ageOrWeightChart, () => { //@ts-ignore
          const firstVisit = studies[this.studyId].domains.bw.col(visitDayName).stats.min;
          //extract bw at first visit
          const df = studies[this.studyId].domains.bw.groupBy([SUBJECT_ID, visitDayName, BW_RES_N])
            .where(`${visitDayName} = ${firstVisit}`).aggregate();
          const age = DG.Viewer.histogram(df, {value: BW_RES_N});
          age.root.prepend(ui.divText('Initial BW', this.viewerTitle));
          updateDivInnerHTML(this.ageOrWeightChart, age.root);
        },
        'Initial BW');
    }

    this.root.className = 'grok-view ui-box';
    const bottomCharts = ui.splitH([
      this.armChart,
      this.sexChart,
      this.raceChart,
      this.ageOrWeightChart,
    ]);

    const summaryDiv = ui.splitV([]);

    if (!this.isSend) {
      summaryDiv.append(ui.splitH([
        ui.panel([
          ui.divText(`${this.studyId} summary`, summaryStyle),
          summary,
        ]),
        ui.panel([
          ui.divText('Errors', summaryStyle),
          errorsSummary,
        ]),
        addDomainToWorkspace,
      ], {style: {maxHeight: '105px'}}));
    } else {
      summaryDiv.append(ui.div(addDomainToWorkspace,
        {style: {maxHeight: '15px', display: 'flex', alignSelf: 'flex-end'}}));
    }

    summaryDiv.append(this.centralChart);
    summaryDiv.append(bottomCharts);
    this.root.append(summaryDiv);
  }

  private async createCumulativeEnrollmentChart(viewerTitle: any) {
    const dateCol = SUBJ_REF_STDT;
    const cumulativeCol = 'CUMULATIVE_ENROLLMENT';
    const incorrectDates = checkDateFormat(studies[this.studyId].domains.dm.getCol(dateCol),
      studies[this.studyId].domains.dm.rowCount);
    if (incorrectDates.length) {
      const subjArray = incorrectDates.map((it) => studies[this.studyId].domains.dm.get(SUBJECT_ID, it));
      const formatErrorMessage = ui.info(`Subjects #${subjArray.join(',')} have incorrect format of ${SUBJ_REF_STDT}`);
      updateDivInnerHTML(this.centralChart, formatErrorMessage);
    } else {
      const subjsPerDay =
        cumulativeEnrollemntByDay(studies[this.studyId].domains.dm
          .clone(studies[this.studyId].domains.dm.filter), dateCol, SUBJECT_ID, cumulativeCol);
      const refStartCol = subjsPerDay.col(dateCol);
      const lc = DG.Viewer.lineChart(subjsPerDay);
      lc.setOptions({x: `${dateCol}`, yColumnNames: [cumulativeCol]});
      if (refStartCol.type != DG.TYPE.DATE_TIME) {
        await subjsPerDay.columns.addNewCalculated(`~${dateCol}`, `Date(\${${dateCol}}, 1, 1)`);
        lc.setOptions({x: `~${dateCol}`, yColumnNames: [cumulativeCol]});
      } else
        lc.setOptions({x: `${dateCol}`, yColumnNames: [cumulativeCol]});

      updateDivInnerHTML(this.centralChart, lc.root);
      lc.root.prepend(ui.divText('Enrollment by day', viewerTitle));
    }
  }

  private createErrorsMapWithLinks(): {[key: string]: HTMLAnchorElement} {
    const errorsByDomain = studies[this.studyId].errorsByDomain;
    if (errorsByDomain) {
      const errorsByDomainWithLinks: {[key: string]: HTMLAnchorElement} = {};
      Object.keys(errorsByDomain).forEach((domain) => {
        const link = ui.link(errorsByDomain[domain].toString(), {}, '', {id: domain});
        link.addEventListener('click', (event) => {
          if (this.validationErrorLinkHandler)
            this.validationErrorLinkHandler();
          else
            grok.shell.v = this.validationView;
          event.stopPropagation();
        });
        errorsByDomainWithLinks[domain] = link;
      });
      return errorsByDomainWithLinks;
    }
    return null;
  }

  override async propertyPanel() {
    if (this.isSend)
      return;
    const acc = this.createAccWithTitle(this.studyId);
    acc.addPane('General', () => {
      const configMap = studyConfigToMap(studies[this.studyId].config);
      return ui.tableFromMap(configMap);
    });
    this.getDataFromClinicalTrialsGov(acc);

    return acc.root;
  }

  getDataFromClinicalTrialsGov(acc: DG.Accordion): void {
    const httpService = new HttpService();
    // const clinTrialsGovInfo: {[key: string]: string | HTMLAnchorElement} =
    //   await httpService.getStudyData('R01NS050536', Object.keys(CLINICAL_TRIAL_GOV_FIELDS));

    httpService.getStudyData(this.studyId, Object.keys(CLINICAL_TRIAL_GOV_FIELDS))
      .then((clinTrialsGovInfo: {[key: string]: string | HTMLAnchorElement}) => {
        if (clinTrialsGovInfo) {
          const studyLink = `${CLIN_TRIAL_GOV_SEARCH}${clinTrialsGovInfo[CLINICAL_TRIAL_GOV_FIELDS.NCTId]}`;
          // eslint-disable-next-line max-len
          clinTrialsGovInfo[`Study link`] = ui.link('Go to study page', () => {window.open(studyLink, '_blank').focus();});
          const acctable = ui.tableFromMap(clinTrialsGovInfo);
          acc.addPane('ClinicalTrials.gov', () => {
            $(acctable).find('tr').css('vertical-align', 'top');
            $(acctable).find('td').css('padding-bottom', '10px');
            $(acctable).find('.d4-entity-list>span').css('margin', '0px');
            return acctable;
          }, true);
        } else {
          acc.addPane('ClinicalTrials.gov', () => {
            return ui.divText('Study not found on clinicaltrials.gov');
          }, true);
        }
      });
  }

  updateGlobalFilter(): void {
    checkColumnsAndCreateViewer(
      studies[this.studyId].domains.dm,
      [SUBJ_REF_STDT],
      this.centralChart,
      async () => {await this.createCumulativeEnrollmentChart(this.viewerTitle);},
      'Cumulative enrollment');
  }
}
