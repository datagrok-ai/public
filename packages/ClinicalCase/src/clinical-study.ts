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
  ae: DG.DataFrame = null; // AE – Adverse Events (SDTM-IG v3.4; limited use in SEND)
  ag: DG.DataFrame = null; // AG – Additional Genetic Tests (SEND-IG, Genomics)
  apce: DG.DataFrame = null; // APCE – Associated Persons Clinical Events (SDTM-IG v3.4)
  apcm: DG.DataFrame = null; // APCM – Associated Persons Concomitant Medications (SDTM-IG v3.4)
  apeg: DG.DataFrame = null; // APEG – Associated Persons ECG (SDTM-IG v3.4)
  aplb: DG.DataFrame = null; // APLB – Associated Persons Laboratory Tests (SDTM-IG v3.4)
  apmh: DG.DataFrame = null; // APMH – Associated Persons Medical History (SDTM-IG v3.4)
  apvs: DG.DataFrame = null; // APVS – Associated Persons Vital Signs (SDTM-IG v3.4)
  bg: DG.DataFrame = null; // BG – Background Genetics (SEND-IG, SEND-only)
  bw: DG.DataFrame = null; // BW – Body Weights (SEND-IG, SEND-only)
  ce: DG.DataFrame = null; // CE – Clinical Events (SDTM-IG v3.4)
  cl: DG.DataFrame = null; // CL – Clinical Observations (SEND-IG, SEND-only)
  cm: DG.DataFrame = null; // CM – Concomitant Medications (SDTM-IG v3.4)
  co: DG.DataFrame = null; // CO – Comments (SDTM & SEND, Special Purpose)
  da: DG.DataFrame = null; // DA – Drug Accountability (SDTM-IG v3.4; not used in SEND)
  dd: DG.DataFrame = null; // DD – Death Details (SDTM-IG v3.4)
  dm: DG.DataFrame = null; // DM – Demographics (SDTM & SEND)
  ds: DG.DataFrame = null; // DS – Disposition (SDTM-IG v3.4)
  dv: DG.DataFrame = null; // DV – Protocol Deviations (SDTM-IG v3.4)
  ec: DG.DataFrame = null; // EC – Exposure as Collected (SDTM & SEND)
  eg: DG.DataFrame = null; // EG – Electrocardiograms (SDTM & SEND)
  ex: DG.DataFrame = null; // EX – Exposure (SDTM & SEND)
  ey: DG.DataFrame = null; // EY – Eye Examinations (SEND-IG, SEND-only)
  fa: DG.DataFrame = null; // FA – Findings About (SDTM-IG v3.4; not used in SEND)
  ft: DG.DataFrame = null; // FT – Functional Tests (SDTM-IG v3.4)
  fw: DG.DataFrame = null; // FW – Food and Water Consumption (SEND-IG, SEND-only)
  gf: DG.DataFrame = null; // GF – Genetic Features (SDTM Genomics v2.0)
  ho: DG.DataFrame = null; // HO – Healthcare Encounters (SDTM-IG v3.4; not used in SEND)
  ie: DG.DataFrame = null; // IE – Inclusion/Exclusion Criteria Not Met (SDTM-IG v3.4)
  ig: DG.DataFrame = null; // IG – In-Life Gestational Findings (SEND-IG, SEND-only)
  is: DG.DataFrame = null; // IS – Immunogenicity Specimen Assessments (SDTM-IG v3.4)
  lb: DG.DataFrame = null; // LB – Laboratory Test Results (SDTM & SEND)
  ma: DG.DataFrame = null; // MA – Macroscopic Findings (SEND-IG, SEND-only)
  mb: DG.DataFrame = null; // MB – Microbiology Findings (SDTM-IG v3.4)
  mh: DG.DataFrame = null; // MH – Medical History (SDTM-IG v3.4; rarely sponsor-defined in SEND)
  mi: DG.DataFrame = null; // MI – SDTM: Microbiology Findings | SEND: Microscopic Findings (⚠ semantic conflict)
  ml: DG.DataFrame = null; // ML – Morphology / Localization (SEND Genomics extension)
  ms: DG.DataFrame = null; // MS – Microbiology Susceptibility (SDTM-IG v3.4)
  nc: DG.DataFrame = null; // NC – Neurobehavioral Findings (SEND-IG, SEND-only)
  oe: DG.DataFrame = null; // OE – Objective Evidence (SDTM-IG v3.4)
  om: DG.DataFrame = null; // OM – Organ Measurements (SEND-IG, SEND-only)
  pc: DG.DataFrame = null; // PC – Pharmacokinetic Concentrations (SDTM & SEND)
  pe: DG.DataFrame = null; // PE – Physical Examination (SDTM-IG v3.4; not used in SEND)
  pm: DG.DataFrame = null; // PM – SEND: Palpable Masses | SDTM: Pharmacometrics (⚠ semantic conflict)
  po: DG.DataFrame = null; // PO – Physical Observations (SEND-IG, SEND-only)
  pp: DG.DataFrame = null; // PP – Pharmacokinetic Parameters (SDTM & SEND)
  pr: DG.DataFrame = null; // PR – Procedures (SDTM-IG v3.4)
  qs: DG.DataFrame = null; // QS – Questionnaires (SDTM-IG v3.4)
  relrec: DG.DataFrame = null; // RELREC – Related Records (SDTM & SEND)
  rp: DG.DataFrame = null; // RP – Reproductive Findings (SDTM-IG v3.4)
  rs: DG.DataFrame = null; // RS – Disease Response / Clinical Classification (SDTM Oncology; not used in SEND)
  sc: DG.DataFrame = null; // SC – Subject Characteristics (SDTM-IG v3.4)
  se: DG.DataFrame = null; // SE – Subject Elements (SDTM-IG v3.4)
  sm: DG.DataFrame = null; // SM – Subject Milestones (SDTM-IG v3.4)
  sr: DG.DataFrame = null; // SR – Short-Term Routine Observations (SEND-IG, SEND-only)
  su: DG.DataFrame = null; // SU – Substance Use (SDTM-IG v3.4)
  sv: DG.DataFrame = null; // SV – Subject Visits (SDTM-IG v3.4)
  ta: DG.DataFrame = null; // TA – Trial Arms (SDTM-IG v3.4; not used in SEND)
  td: DG.DataFrame = null; // TD – Trial Disease Assessments (SDTM-IG v3.4)
  te: DG.DataFrame = null; // TE – Trial Elements (SDTM-IG v3.4)
  tf: DG.DataFrame = null; // TF – Tumor Findings (SEND-IG, SEND-only)
  ti: DG.DataFrame = null; // TI – Trial Inclusion/Exclusion Criteria (SDTM-IG v3.4)
  tm: DG.DataFrame = null; // TM – Trial Disease Milestones (SDTM-IG v3.4)
  tr: DG.DataFrame = null; // TR – Tumor/Lesion Results (SDTM Oncology; not used in SEND)
  ts: DG.DataFrame = null; // TS – Trial Summary (SDTM-IG v3.4)
  tu: DG.DataFrame = null; // TU – Tumor/Lesion Identification (SDTM Oncology; not used in SEND)
  tv: DG.DataFrame = null; // TV – Trial Visits (SDTM-IG v3.4)
  tx: DG.DataFrame = null; // TX – SEND Trial Summary (SEND-IG, SEND-only)
  vs: DG.DataFrame = null; // VS – Vital Signs (SDTM & SEND)
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
