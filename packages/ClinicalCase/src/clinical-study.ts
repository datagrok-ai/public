/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {COL_HAS_ERRORS_POSTFIX, DOMAIN, ERRORS_POSTFIX, HAS_VALIDATION_ERRORS_COL,
  SITE_ID, STUDY_ID, SUBJECT_ID, VIOLATED_RULES_COL} from './constants/columns-constants';
import {addVisitDayFromTvDomain, createEventStartEndDaysCol,
  calculateLBBaselineColumns} from './data-preparation/data-preparation';
import {createFilters, removeExtension} from './utils/utils';
import {CDISC_STANDARD, ClinStudyConfig} from './utils/types';
import {ClinicalCaseViewsConfig} from './views-config';
import {SUMMARY_VIEW_NAME} from './constants/view-names-constants';
import {funcs} from './package-api';
import {Subject} from 'rxjs';
import {ValidationResult, IssueDetail} from './types/validation-result';
import {COLUMN_FROM_DM_TAG} from './constants/constants';

export class ClinicalDomains {
  // Domains listed in alphabetical order
  ae: DG.DataFrame = null; // Adverse Events (SDTM-IG v3.4, Clinical)
  ag: DG.DataFrame = null; // Additional Genetic Tests (Custom/Extension, Genomics)
  apce: DG.DataFrame = null; // Associated Persons Clinical Events (SDTM-IG v3.4, Associated Persons)
  apcm: DG.DataFrame = null; // Associated Persons Concomitant Medications (SDTM-IG v3.4, Associated Persons)
  apeg: DG.DataFrame = null; // Associated Persons Electrocardiograms (SDTM-IG v3.4, Associated Persons)
  aplb: DG.DataFrame = null; // Associated Persons Laboratory Tests (SDTM-IG v3.4, Associated Persons)
  apmh: DG.DataFrame = null; // Associated Persons Medical History (SDTM-IG v3.4, Associated Persons)
  apvs: DG.DataFrame = null; // Associated Persons Vital Signs (SDTM-IG v3.4, Associated Persons)
  bg: DG.DataFrame = null; // Background Genetics (SEND-IG v3.2, SEND)
  bw: DG.DataFrame = null; // Body Weights (SEND-IG v3.2, SEND)
  ce: DG.DataFrame = null; // Clinical Events (SDTM-IG v3.4, Clinical)
  cg: DG.DataFrame = null; // Cardiovascular System Findings (SEND-IG v3.2, SEND)
  cl: DG.DataFrame = null; // Clinical Observations (SEND-IG v3.2, SEND)
  cm: DG.DataFrame = null; // Concomitant Medications (SDTM-IG v3.4, Clinical)
  co: DG.DataFrame = null; // Comments (SDTM-IG v3.4, Special Purpose)
  cv: DG.DataFrame = null; // Cardiovascular System Findings (SEND-IG v3.2, SEND)
  da: DG.DataFrame = null; // Drug Accountability (SDTM-IG v3.4, Clinical)
  dd: DG.DataFrame = null; // Death Details (SDTM-IG v3.4, Clinical)
  dm: DG.DataFrame = null; // Demographics (SDTM-IG v3.4, Special Purpose)
  ds: DG.DataFrame = null; // Disposition (SDTM-IG v3.4, Clinical)
  dv: DG.DataFrame = null; // Protocol Deviations (SDTM-IG v3.4, Clinical)
  ec: DG.DataFrame = null; // Exposure as Collected (SDTM-IG v3.4, Clinical)
  eg: DG.DataFrame = null; // Electrocardiograms (SDTM-IG v3.4, Clinical)
  ex: DG.DataFrame = null; // Exposure (SDTM-IG v3.4, Clinical)
  ey: DG.DataFrame = null; // Eye Examinations (SEND-IG v3.2, SEND)
  fa: DG.DataFrame = null; // Findings About (SDTM-IG v3.4, Clinical)
  ft: DG.DataFrame = null; // Functional Tests (SDTM-IG v3.4, Clinical)
  fw: DG.DataFrame = null; // Food and Water Consumption (SEND-IG v3.2, SEND)
  gf: DG.DataFrame = null; // Genetic Features (SDTM v2.0, Genomics)
  ho: DG.DataFrame = null; // Healthcare Encounters (SDTM-IG v3.4, Clinical)
  hp: DG.DataFrame = null; // Hematology Findings (SEND-IG v3.2, SEND)
  ie: DG.DataFrame = null; // Inclusion/Exclusion Criteria Not Met (SDTM-IG v3.4, Clinical)
  ig: DG.DataFrame = null; // In-Life Gestational Findings (SEND-IG v3.2, SEND)
  is: DG.DataFrame = null; // Immunogenicity Specimen Assessments (SDTM-IG v3.4, Clinical)
  la: DG.DataFrame = null; // Laboratory Findings (SEND-IG v3.2, SEND)
  lb: DG.DataFrame = null; // Laboratory Test Results (SDTM-IG v3.4, Clinical)
  ma: DG.DataFrame = null; // Macropathology Findings (SEND-IG v3.2, SEND)
  mb: DG.DataFrame = null; // Microbiology Findings (SDTM-IG v3.4, Clinical)
  mh: DG.DataFrame = null; // Medical History (SDTM-IG v3.4, Clinical)
  mi: DG.DataFrame = null; // Microbiology Findings (SEND-IG v3.2, SEND)
  mk: DG.DataFrame = null; // Microscopic Findings (SEND-IG v3.2, SEND)
  ml: DG.DataFrame = null; // Morphology/Localization (Custom/Extension, Genomics)
  mo: DG.DataFrame = null; // Morphology Findings (SEND-IG v3.2, SEND)
  ms: DG.DataFrame = null; // Microbiology Susceptibility (SDTM-IG v3.4, Clinical)
  nc: DG.DataFrame = null; // Neurobehavioral/Cognitive Findings (SEND-IG v3.2, SEND)
  nv: DG.DataFrame = null; // Nervous System Findings (SEND-IG v3.2, SEND)
  oe: DG.DataFrame = null; // Objective Evidence (SDTM-IG v3.4, Clinical)
  om: DG.DataFrame = null; // Organ Measurements (SEND-IG v3.2, SEND)
  pc: DG.DataFrame = null; // Pharmacogenetics/Genomics (SDTM-IG v3.4, Clinical)
  pm: DG.DataFrame = null; // Pharmacometrics
  pe: DG.DataFrame = null; // Physical Examination (SDTM-IG v3.4, Clinical)
  po: DG.DataFrame = null; // Physical Observations (SEND-IG v3.2, SEND)
  pp: DG.DataFrame = null; // Pharmacokinetic Parameters (SDTM-IG v3.4, Clinical)
  pr: DG.DataFrame = null; // Procedures (SDTM-IG v3.4, Clinical)
  qs: DG.DataFrame = null; // Questionnaires (SDTM-IG v3.4, Clinical)
  re: DG.DataFrame = null; // Reproductive System Findings (SEND-IG v3.2, SEND)
  relrec: DG.DataFrame = null; // Related Records (SDTM-IG v3.4, Special Purpose)
  rp: DG.DataFrame = null; // Reproductive System Findings (SDTM-IG v3.4, Clinical)
  rs: DG.DataFrame = null; // Disease Response/Clinical Classification (SDTM-IG v3.4, Clinical)
  sc: DG.DataFrame = null; // Subject Characteristics (SDTM-IG v3.4, Clinical)
  se: DG.DataFrame = null; // Subject Elements (SDTM-IG v3.4, Special Purpose)
  sm: DG.DataFrame = null; // Subject Milestones (SDTM-IG v3.4, Special Purpose)
  sr: DG.DataFrame = null; // Short-Term Routine Observations (SEND-IG v3.2, SEND)
  su: DG.DataFrame = null; // Substance Use (SDTM-IG v3.4, Clinical)
  sv: DG.DataFrame = null; // Subject Visits (SDTM-IG v3.4, Special Purpose)
  ta: DG.DataFrame = null; // Trial Arms (SDTM-IG v3.4, Trial Design)
  td: DG.DataFrame = null; // Trial Disease Assessments (SDTM-IG v3.4, Trial Design)
  te: DG.DataFrame = null; // Trial Elements (SDTM-IG v3.4, Trial Design)
  ti: DG.DataFrame = null; // Trial Inclusion/Exclusion Criteria (SDTM-IG v3.4, Trial Design)
  tm: DG.DataFrame = null; // Trial Disease Milestones (SDTM-IG v3.4, Trial Design)
  tr: DG.DataFrame = null; // Tumor/Lesion Results (SDTM-IG v3.4, Clinical)
  ts: DG.DataFrame = null; // Trial Summary (SDTM-IG v3.4, Trial Design)
  tu: DG.DataFrame = null; // Tumor/Lesion Identification (SDTM-IG v3.4, Clinical)
  tv: DG.DataFrame = null; // Trial Visits (SDTM-IG v3.4, Trial Design)
  tx: DG.DataFrame = null; // SEND Trial Summary (SEND-IG v3.2, SEND)
  ur: DG.DataFrame = null; // Urinalysis Findings (Custom/Extension, Clinical)
  vs: DG.DataFrame = null; // Vital Signs (SDTM-IG v3.4, Clinical)
  supp: DG.DataFrame[] = [];

  all(): DG.DataFrame[] {
    const dfs = Object.keys(this).filter((it) => it !== 'supp').map((k) => this[k]).filter((v) => v != null);
    return dfs.concat(this.supp);
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
  private validationColumnsSubscriptionSet = false;
  changeViewToSummary = true;
  treatmentAndControlConfig: {treatment: string, control: string}[] = [];

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
    if (this.config.standard === CDISC_STANDARD.SDTM) {
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
    this.domains.all().forEach((it) => {
      if (it.name !== 'dm' && it.columns.names().includes(SUBJECT_ID)) {
        const savedName = it.name;
        const columnsFromDm = this.domains.dm.columns.names().filter((it) => it !== DOMAIN);
        for (const colName of it.columns.names()) {
          const colIdx = columnsFromDm.findIndex((it) => it === colName);
          if (colIdx !== -1)
            columnsFromDm.splice(colIdx, 1);
        }
        grok.data.joinTables(it, this.domains.dm, [SUBJECT_ID], [SUBJECT_ID], null,
          columnsFromDm, DG.JOIN_TYPE.LEFT, true);
        it.name = savedName;
        for (const col of columnsFromDm)
          it.col(col)!.setTag(COLUMN_FROM_DM_TAG, 'true');
      }
      for (const col of it.columns.names()) {
        if (this.config.fieldsDefinitions && this.config.fieldsDefinitions[col])
          it.col(col)!.setTag('Description', this.config.fieldsDefinitions[col]);
      }
    });

    // Calculate LB baseline and post-baseline columns if LB domain exists
    if (this.domains.lb)
      calculateLBBaselineColumns(this.domains.lb);
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
      .exists(`System:AppData/ClinicalCase/${this.config.standard!}/${this.studyId}/validation_results.json`)) {
      validationResStr = await grok.dapi.files
        // eslint-disable-next-line max-len
        .readAsText(`System:AppData/ClinicalCase/${this.config.standard!}/${this.studyId}/validation_results.json`);
    } else {
      validationResStr = await funcs.runCoreValidate(
        this.config.standard === CDISC_STANDARD.SEND ? 'sendig' : 'sdtmig',
        `ClinicalCase/${this.config.standard!}/${this.studyId}`, '3.1', 'json', undefined);
      grok.dapi.files
        // eslint-disable-next-line max-len
        .writeAsText(`System:AppData/ClinicalCase/${this.config.standard!}/${this.studyId}/validation_results.json`,
          validationResStr);
    }
    this.validationResults = JSON.parse(validationResStr) as ValidationResult;
    this.validated = true;
    this.validationCompleted.next(true);
  }

  /**
   * Ensures validation columns are added to domains.
   * If validation is already completed, adds columns immediately.
   * If validation is not completed yet, subscribes to validationCompleted to add columns when ready.
   */
  ensureValidationColumnsAdded(): void {
    if (this.validated && this.validationResults) {
      // Validation already completed, add columns immediately
      this.addValidationColumnsToDomains();
    } else if (!this.validationColumnsSubscriptionSet) {
      // Validation not completed yet, subscribe to add columns when validation completes
      // Only subscribe once to avoid multiple subscriptions
      this.validationCompleted.subscribe(() => {
        this.addValidationColumnsToDomains();
      });
      this.validationColumnsSubscriptionSet = true;
    }
  }

  /**
   * Adds validation rule violation columns to all existing domains:
   * - A boolean column indicating whether any rule is violated for the row
   * - A string column containing JSON array of IssueDetail objects for violated rules
   */
  private addValidationColumnsToDomains(): void {
    if (!this.validationResults || !this.validationResults.Issue_Details)
      return;

    // Group issues by domain (dataset)
    const issuesByDomain: {[domain: string]: {[row: number]: IssueDetail[]}} = {};

    for (const issue of this.validationResults.Issue_Details) {
      if (!issue.dataset)
        continue;
      const domain = removeExtension(issue.dataset.toLowerCase());

      // Parse row number - it can be string or number
      let rowNum: number;
      if (typeof issue.row === 'string') {
        rowNum = parseInt(issue.row, 10);
        if (isNaN(rowNum))
          continue;
      } else
        rowNum = issue.row;

      // Initialize domain entry if needed
      if (!issuesByDomain[domain])
        issuesByDomain[domain] = {};

      // Initialize row entry if needed
      if (!issuesByDomain[domain][rowNum])
        issuesByDomain[domain][rowNum] = [];

      // Add IssueDetail object to the row
      issuesByDomain[domain][rowNum].push(issue);
    }

    const addValidationColumns = (domainName: string, domain: DG.DataFrame) => {
      const domainLower = domainName.toLowerCase();
      const domainIssues = issuesByDomain[domainLower];

      if (!domainIssues) {
        // No issues for this domain, but still add columns with false/empty values
        this.addValidationColumnsToDomain(domain, {});
        return;
      }

      this.addValidationColumnsToDomain(domain, domainIssues);
    };

    // Process each domain and add columns
    for (const domain of this.domains.all())
      addValidationColumns(domain.name, domain);
  }

  /**
   * Adds validation columns to a specific domain dataframe
   * @param {DG.DataFrame} domain - The domain dataframe to add columns to
   * @param {Object} issuesByRow - Map of row index to array of IssueDetail objects
   */
  private addValidationColumnsToDomain(
    domain: DG.DataFrame,
    issuesByRow: {[row: number]: IssueDetail[]},
  ): void {
    // Create boolean column for validation errors
    const hasErrorsCol = domain.columns.addNewBool(HAS_VALIDATION_ERRORS_COL);

    // Create string column for violated rules (stored as JSON array of IssueDetail objects)
    const violatedRulesCol = domain.columns.addNewString(VIOLATED_RULES_COL);
    violatedRulesCol.semType = 'sdisc-rule-violation';


    // Maps to track variable columns (created lazily as we encounter variables)
    const variableErrorColumns: {[variable: string]: DG.Column} = {};
    const variableHasErrorsColumns: {[variable: string]: DG.Column} = {};

    // Helper function to create columns for a variable if they don't exist
    const ensureVariableColumns = (variable: string): void => {
      if (!variableErrorColumns[variable] && domain.col(variable)) {
        const hasErrorsColName = `${variable}${COL_HAS_ERRORS_POSTFIX}`;
        const hasErrorsCol = domain.columns.addNewBool(hasErrorsColName);
        variableHasErrorsColumns[variable] = hasErrorsCol;

        const errorColName = `${variable}${ERRORS_POSTFIX}`;
        const errorCol = domain.columns.addNewString(errorColName);
        variableErrorColumns[variable] = errorCol;
      }
    };

    // Iterate through rows that have issues
    for (const [rowIndexStr, rowIssues] of Object.entries(issuesByRow)) {
      const rowIndex = parseInt(rowIndexStr, 10) - 1; //rows in issues are numbered starting from 1
      if (isNaN(rowIndex) || !rowIssues || rowIssues.length === 0)
        continue;

      hasErrorsCol.set(rowIndex, true);
      // Serialize IssueDetail objects as JSON array string
      violatedRulesCol.set(rowIndex, JSON.stringify(rowIssues));

      // Process issues, collect errors by variable, create columns, and set values in one loop
      const errorsByVariable: {[variable: string]: Array<{ruleID: string, message: string, value: string}>} = {};

      for (const issue of rowIssues) {
        if (issue.variables && issue.values) {
          // Iterate through variables and values arrays
          for (let varIndex = 0; varIndex < issue.variables.length; varIndex++) {
            const variable = issue.variables[varIndex];

            if (domain.col(variable)?.getTag(COLUMN_FROM_DM_TAG))
              continue;

            const value = varIndex < issue.values.length ? issue.values[varIndex] : '';

            // Create columns for this variable if they don't exist
            ensureVariableColumns(variable);

            // Only process if this variable has an error column (exists in domain)
            if (variableErrorColumns[variable]) {
              if (!errorsByVariable[variable])
                errorsByVariable[variable] = [];
              errorsByVariable[variable].push({
                ruleID: issue.core_id,
                message: issue.message,
                value: value,
              });

              // Set values immediately - update JSON string and set boolean flag
              const errorCol = variableErrorColumns[variable];
              const varHasErrorsCol = variableHasErrorsColumns[variable];
              errorCol.set(rowIndex, JSON.stringify(errorsByVariable[variable]));
              varHasErrorsCol.set(rowIndex, true);
            }
          }
        }
      }
    }
  }
}


export class ClinRow {
  row: DG.Row;

  constructor(row: DG.Row) {
    this.row = row;
  }
}
