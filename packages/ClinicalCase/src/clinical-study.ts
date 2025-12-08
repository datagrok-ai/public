/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {vaidateAEDomain, vaidateDMDomain} from './sdtm-validation/services/validation-service';
import {createValidationDataFrame} from './sdtm-validation/validation-utils';
import {SITE_ID, STUDY_ID} from './constants/columns-constants';
import {addVisitDayFromTvDomain, createEventStartEndDaysCol} from './data-preparation/data-preparation';
import {createFilters} from './utils/utils';
import {createErrorsByDomainMap} from './utils/views-validation-utils';
import {CDISC_STANDARD, ClinStudyConfig} from './utils/types';
import {ClinicalCaseViewsConfig} from './views-config';
import {SUMMARY_VIEW_NAME} from './constants/view-names-constants';
import {funcs} from './package-api';
import {Subject} from 'rxjs';
import {ValidationResult} from './types/validation-result';

export class ClinicalDomains {
  ae: DG.DataFrame = null;
  ce: DG.DataFrame = null;
  ds: DG.DataFrame = null;
  dv: DG.DataFrame = null;
  ho: DG.DataFrame = null;
  mh: DG.DataFrame = null;
  cv: DG.DataFrame = null;
  da: DG.DataFrame = null;
  dd: DG.DataFrame = null;
  eg: DG.DataFrame = null;
  fa: DG.DataFrame = null;
  ft: DG.DataFrame = null;
  ie: DG.DataFrame = null;
  is: DG.DataFrame = null;
  lb: DG.DataFrame = null;
  mb: DG.DataFrame = null;
  mi: DG.DataFrame = null;
  mk: DG.DataFrame = null;
  mo: DG.DataFrame = null;
  ms: DG.DataFrame = null;
  nv: DG.DataFrame = null;
  oe: DG.DataFrame = null;
  pc: DG.DataFrame = null;
  pe: DG.DataFrame = null;
  pp: DG.DataFrame = null;
  qs: DG.DataFrame = null;
  re: DG.DataFrame = null;
  rp: DG.DataFrame = null;
  rs: DG.DataFrame = null;
  sc: DG.DataFrame = null;
  sr: DG.DataFrame = null;
  ss: DG.DataFrame = null;
  tr: DG.DataFrame = null;
  tu: DG.DataFrame = null;
  ur: DG.DataFrame = null;
  vs: DG.DataFrame = null;
  ag: DG.DataFrame = null;
  cm: DG.DataFrame = null;
  ec: DG.DataFrame = null;
  ex: DG.DataFrame = null;
  ml: DG.DataFrame = null;
  pr: DG.DataFrame = null;
  su: DG.DataFrame = null;
  co: DG.DataFrame = null;
  dm: DG.DataFrame = null;
  se: DG.DataFrame = null;
  sm: DG.DataFrame = null;
  sv: DG.DataFrame = null;
  ta: DG.DataFrame = null;
  td: DG.DataFrame = null;
  te: DG.DataFrame = null;
  ti: DG.DataFrame = null;
  tm: DG.DataFrame = null;
  ts: DG.DataFrame = null;
  tv: DG.DataFrame = null;
  bw: DG.DataFrame = null;
  bg: DG.DataFrame = null;
  cl: DG.DataFrame = null;

  all(): DG.DataFrame[] {
    return Object.keys(this).map((k) => this[k]).filter((v) => v != null);
  }

  static allClinicalDomainsNames(): string[] {
    return Object.keys(this);
  }
}

export class ClinicalStudy {
  description: string;
  domains: ClinicalDomains = new ClinicalDomains();
  subjectsCount: number;
  sitesCount: number;
  validationResults: ValidationResult;
  errorsByDomain: {[key: string]: number};
  dmFilters: DG.Viewer;
  studyId: string;
  validated = false;
  subjSitesCountsProcessed = false;
  initCompleted = false;
  config: ClinStudyConfig;
  viewsConfig = new ClinicalCaseViewsConfig();
  views: {[key: string]: DG.ViewBase} = {};
  loadingStudyData: boolean | null = null;
  currentViewName: string = SUMMARY_VIEW_NAME;
  validationCompleted = new Subject<boolean>();

  constructor(config: ClinStudyConfig) {
    this.config = config;
    this.studyId = config.name;
  }

  initFromWorkspace(): void {
    for (const t of grok.shell.tables) {
      const view = grok.shell.tableView(t.name);
      if (view != null)
        view.syncCurrentObject = false;

      if (t.name.toLowerCase() in this.domains)
        this.domains[t.name.toLowerCase()] = t;
    }

    this.process();

    this.validate();
  }

  init() {
    this.process();
    if (!this.validated)
      this.validate().then(() => this.initCompleted = true);
    else
      this.initCompleted = true;
  }


  private process(): void {
    if (!this.subjSitesCountsProcessed)
      this.processSitesAndSubjectCount();
    createEventStartEndDaysCol(this.studyId);
    addVisitDayFromTvDomain(this.studyId);
    if (this.domains.dm) {
      this.dmFilters = createFilters(this.domains.dm);
      this.domains.dm.onFilterChanged.subscribe(() => {
      });
      grok.shell.topMenu
        .group('Cohort')
        .item('Filter', () => {
          const dialog = ui.dialog({title: ''})
            .add(ui.div(this.dmFilters.root))
          //@ts-ignore
            .onOK(() => {})
            .show();
          dialog.root.addEventListener('mouseenter', (event) => {
            this.dmFilters.root.removeAttribute('data-widget');
          });
        });
    }
  }


  processSitesAndSubjectCount() {
    if (this.domains.dm != null) {
      this.subjectsCount = this.domains.dm.rowCount;
      if (this.domains.dm.col(SITE_ID))
        this.sitesCount = this.domains.dm.col(SITE_ID).stats.uniqueCount;

      if (!this.studyId && this.domains.dm.col(STUDY_ID))
        this.studyId = this.domains.dm.col(STUDY_ID).get(0);
    }
    this.subjSitesCountsProcessed = true;
  }

  async validate(): Promise<void> {
    let validationResStr = '';
    if (await grok.dapi.files
      .exists(`System:AppData/ClinicalCase/${this.config.standard!}/studies/${this.studyId}/validation_results.json`)) {
      validationResStr = await grok.dapi.files
        // eslint-disable-next-line max-len
        .readAsText(`System:AppData/ClinicalCase/${this.config.standard!}/studies/${this.studyId}/validation_results.json`);
    } else {
      validationResStr = await funcs.runCoreValidate(
        this.config.standard === CDISC_STANDARD.SEND ? 'sendig' : 'sdtmig',
        `ClinicalCase/${this.config.standard!}/studies/${this.studyId}`, '3.1', 'json', undefined);
      grok.dapi.files
        // eslint-disable-next-line max-len
        .writeAsText(`System:AppData/ClinicalCase/${this.config.standard!}/studies/${this.studyId}/validation_results.json`,
          validationResStr);
    }
    this.validationResults = JSON.parse(validationResStr) as ValidationResult;
    this.validated = true;
    this.validationCompleted.next(true);
  }
}

export class ClinRow {
  row: DG.Row;

  constructor(row: DG.Row) {
    this.row = row;
  }
}
