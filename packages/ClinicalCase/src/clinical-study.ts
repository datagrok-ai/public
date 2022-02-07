/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as meta from './sdtm-meta';
import { vaidateAEDomain, vaidateDMDomain } from './validation/services/validation-service';
import { createValidationDataFrame } from './validation/validation-utils';
import { AE_START_DAY, SITE_ID, STUDY_ID, VISIT_DAY, VISIT_NAME } from './columns-constants';
import { StudyVisit } from './model/study-visit';
import { getVisitNamesAndDays } from './data-preparation/utils';

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

  all(): DG.DataFrame[] {
    return Object.keys(this).map((k) => this[ k ]).filter((v) => v != null)
  }
}

export class ClinicalStudy {
  name: string;
  description: string;
  domains: ClinicalDomains = new ClinicalDomains();
  subjectsCount: number;
  sitesCount: number;
  validationResults: DG.DataFrame;
  labDataForCorelationMatrix: DG.DataFrame;
  visits: StudyVisit[] = [];

  initFromWorkspace(): void {
    for (let t of grok.shell.tables) {
      let view = grok.shell.tableView(t.name);
      if (view != null)
        view.syncCurrentObject = false;

      if (t.name.toLowerCase() in this.domains)
        this.domains[ t.name.toLowerCase() ] = t;
    }

    if (this.domains.dm != null) {
      this.subjectsCount = this.domains.dm.rowCount;
      if (this.domains.dm.col(SITE_ID)) {
        this.sitesCount = this.domains.dm.col(SITE_ID).stats.uniqueCount;
      }
      if (this.domains.dm.col(STUDY_ID)) {
        this.name = this.domains.dm.col(STUDY_ID).get(0);
      }
    }

    this.process();

    this.validate();

    this.createStudyConfig();

  }

  private process(): void {
    if (this.domains.ae != null) {
      if (this.domains.ae.col(AE_START_DAY)) {
        this.domains.ae.columns.addNewCalculated('week', `floor(\${${AE_START_DAY}} / 7)`);
      }
    }
  }

  private validate(): void {
    this.validationResults = createValidationDataFrame();
    if (this.domains.ae != null) {
      vaidateAEDomain(study.domains.ae, this.validationResults);
    }
    if (this.domains.dm != null) {
      vaidateDMDomain(study.domains.dm, this.validationResults);
    }
  }

  private createStudyConfig() {
    if (this.domains.sv) {
      let visits = getVisitNamesAndDays(this.domains.sv, true);
      visits.forEach(vis => this.visits.push(new StudyVisit(vis.name, vis.day)))
      this.visits.forEach((it, index) => it.num = index + 1);
    }
  }
}

export class ClinRow {
  row: DG.Row;

  constructor(row: DG.Row) {
    this.row = row;
  }
}

export let study: ClinicalStudy = new ClinicalStudy();
