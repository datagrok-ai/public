/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {COL_HAS_ERRORS_POSTFIX, DOMAIN, ERRORS_POSTFIX, HAS_VALIDATION_ERRORS_COL,
  SITE_ID, STUDY_ID, SUBJECT_ID, VIOLATED_RULES_COL} from './constants/columns-constants';
import {addVisitDayFromTvDomain, createEventStartEndDaysCol,
  calculateLBBaselineColumns} from './data-preparation/data-preparation';
import {createFilters, removeExtension} from './utils/utils';
import {ClinStudyConfig} from './utils/types';
import {ClinicalCaseViewsConfig} from './views-config';
import {SUMMARY_VIEW_NAME} from './constants/view-names-constants';
import {funcs} from './package-api';
import {Subject} from 'rxjs';
import {ValidationResult, IssueDetail} from './types/validation-result';
import {COLUMN_FROM_DM_TAG} from './constants/constants';

export class ClinicalDomains {
  // Domains listed in alphabetical order
  ae: DG.DataFrame = null; // AE – Adverse Events (SDTM-IG v3.4)
  apce: DG.DataFrame = null; // APCE – Associated Persons Clinical Events (SDTM-IG v3.4)
  apcm: DG.DataFrame = null; // APCM – Associated Persons Concomitant Medications (SDTM-IG v3.4)
  apeg: DG.DataFrame = null; // APEG – Associated Persons ECG (SDTM-IG v3.4)
  aplb: DG.DataFrame = null; // APLB – Associated Persons Laboratory Tests (SDTM-IG v3.4)
  apmh: DG.DataFrame = null; // APMH – Associated Persons Medical History (SDTM-IG v3.4)
  apvs: DG.DataFrame = null; // APVS – Associated Persons Vital Signs (SDTM-IG v3.4)
  ce: DG.DataFrame = null; // CE – Clinical Events (SDTM-IG v3.4)
  cm: DG.DataFrame = null; // CM – Concomitant Medications (SDTM-IG v3.4)
  co: DG.DataFrame = null; // CO – Comments (SDTM, Special Purpose)
  da: DG.DataFrame = null; // DA – Drug Accountability (SDTM-IG v3.4)
  dd: DG.DataFrame = null; // DD – Death Details (SDTM-IG v3.4)
  dm: DG.DataFrame = null; // DM – Demographics (SDTM)
  ds: DG.DataFrame = null; // DS – Disposition (SDTM-IG v3.4)
  dv: DG.DataFrame = null; // DV – Protocol Deviations (SDTM-IG v3.4)
  ec: DG.DataFrame = null; // EC – Exposure as Collected (SDTM)
  eg: DG.DataFrame = null; // EG – Electrocardiograms (SDTM)
  ex: DG.DataFrame = null; // EX – Exposure (SDTM)
  fa: DG.DataFrame = null; // FA – Findings About (SDTM-IG v3.4)
  ft: DG.DataFrame = null; // FT – Functional Tests (SDTM-IG v3.4)
  gf: DG.DataFrame = null; // GF – Genetic Features (SDTM Genomics v2.0)
  ho: DG.DataFrame = null; // HO – Healthcare Encounters (SDTM-IG v3.4)
  ie: DG.DataFrame = null; // IE – Inclusion/Exclusion Criteria Not Met (SDTM-IG v3.4)
  is: DG.DataFrame = null; // IS – Immunogenicity Specimen Assessments (SDTM-IG v3.4)
  lb: DG.DataFrame = null; // LB – Laboratory Test Results (SDTM)
  mb: DG.DataFrame = null; // MB – Microbiology Findings (SDTM-IG v3.4)
  mh: DG.DataFrame = null; // MH – Medical History (SDTM-IG v3.4)
  mi: DG.DataFrame = null; // MI – Microbiology Findings (SDTM)
  ms: DG.DataFrame = null; // MS – Microbiology Susceptibility (SDTM-IG v3.4)
  oe: DG.DataFrame = null; // OE – Objective Evidence (SDTM-IG v3.4)
  pc: DG.DataFrame = null; // PC – Pharmacogenetics / Genomics (SDTM)
  pe: DG.DataFrame = null; // PE – Physical Examination (SDTM-IG v3.4)
  pm: DG.DataFrame = null; // PM – Pharmacometrics (SDTM)
  pp: DG.DataFrame = null; // PP – Pharmacokinetic Parameters (SDTM)
  pr: DG.DataFrame = null; // PR – Procedures (SDTM-IG v3.4)
  qs: DG.DataFrame = null; // QS – Questionnaires (SDTM-IG v3.4)
  relrec: DG.DataFrame = null; // RELREC – Related Records (SDTM)
  rp: DG.DataFrame = null; // RP – Reproductive Findings (SDTM-IG v3.4)
  rs: DG.DataFrame = null; // RS – Disease Response / Clinical Classification (SDTM Oncology)
  sc: DG.DataFrame = null; // SC – Subject Characteristics (SDTM-IG v3.4)
  se: DG.DataFrame = null; // SE – Subject Elements (SDTM-IG v3.4)
  sm: DG.DataFrame = null; // SM – Subject Milestones (SDTM-IG v3.4)
  su: DG.DataFrame = null; // SU – Substance Use (SDTM-IG v3.4)
  sv: DG.DataFrame = null; // SV – Subject Visits (SDTM-IG v3.4)
  ta: DG.DataFrame = null; // TA – Trial Arms (SDTM-IG v3.4)
  td: DG.DataFrame = null; // TD – Trial Disease Assessments (SDTM-IG v3.4)
  te: DG.DataFrame = null; // TE – Trial Elements (SDTM-IG v3.4)
  ti: DG.DataFrame = null; // TI – Trial Inclusion/Exclusion Criteria (SDTM-IG v3.4)
  tm: DG.DataFrame = null; // TM – Trial Disease Milestones (SDTM-IG v3.4)
  tr: DG.DataFrame = null; // TR – Tumor/Lesion Results (SDTM Oncology)
  ts: DG.DataFrame = null; // TS – Trial Summary (SDTM-IG v3.4)
  tu: DG.DataFrame = null; // TU – Tumor/Lesion Identification (SDTM Oncology)
  tv: DG.DataFrame = null; // TV – Trial Visits (SDTM-IG v3.4)
  vs: DG.DataFrame = null; // VS – Vital Signs (SDTM)
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
  // Dataframes for the two big arrays in the report. Built once when validate()
  // loads the cached results; the validation view consumes them directly so we
  // never run DG.DataFrame.fromObjects on 50k+ rows on the main thread.
  issueSummaryDf: DG.DataFrame | null = null;
  issueDetailsDf: DG.DataFrame | null = null;
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
    const base = `System:AppData/ClinicalCase/${this.config.standard!}/${this.studyId}`;
    const d42Path = `${base}/validation_results.d42`;
    const jsonPath = `${base}/validation_results.json`;
    const summaryCsvPath = `${base}/issue_summary.csv`;
    const detailsCsvPath = `${base}/issue_details.csv`;

    // Fast path: per-study binary cache produced from a previous load. d42
    // deserialization is much faster than CSV/JSON parsing for large reports.
    if (await grok.dapi.files.exists(d42Path)) {
      const dfs = await grok.dapi.files.readBinaryDataFrames(d42Path);
      this.issueSummaryDf = dfs.find((df) => df.name === 'Issue_Summary') ?? dfs[0] ?? null;
      this.issueDetailsDf = dfs.find((df) => df.name === 'Issue_Details') ?? dfs[1] ?? null;
      // Best-effort metadata: tiny file, parse cost negligible.
      if (await grok.dapi.files.exists(jsonPath)) {
        const meta = await grok.dapi.files.readAsText(jsonPath);
        this.validationResults = JSON.parse(meta) as ValidationResult;
      }
      this.validated = true;
      this.validationCompleted.next(true);
      return;
    }

    // Container output is missing — run validation in the container. After this,
    // AppData should contain the slim JSON + the two CSV files (see app.py).
    if (!await grok.dapi.files.exists(jsonPath)) {
      await funcs.runCoreValidate(
        'sdtmig',
        `ClinicalCase/${this.config.standard!}/${this.studyId}`,
        resolveSdtmigVersion(this.config.standardVersion), 'json', undefined);
    }

    // Preferred: CSV-per-table format. Skips DG.DataFrame.fromObjects entirely.
    const haveCsvs = await Promise.all([
      grok.dapi.files.exists(summaryCsvPath),
      grok.dapi.files.exists(detailsCsvPath),
    ]);
    if (haveCsvs.every(Boolean)) {
      const [summaryCsv, detailsCsv, metaText] = await Promise.all([
        grok.dapi.files.readAsText(summaryCsvPath),
        grok.dapi.files.readAsText(detailsCsvPath),
        grok.dapi.files.readAsText(jsonPath),
      ]);
      this.issueSummaryDf = DG.DataFrame.fromCsv(summaryCsv);
      this.issueSummaryDf.name = 'Issue_Summary';
      this.issueDetailsDf = DG.DataFrame.fromCsv(detailsCsv);
      this.issueDetailsDf.name = 'Issue_Details';
      this.validationResults = JSON.parse(metaText) as ValidationResult;
      // Write the d42 cache so the next open is fast. Best-effort: failing to
      // write should not break the current view.
      grok.dapi.files
        .writeBinaryDataFrames(d42Path, [this.issueSummaryDf, this.issueDetailsDf])
        .catch((e) => grok.shell.warning(`Failed to write d42 cache: ${e?.message ?? e}`));
      this.validated = true;
      this.validationCompleted.next(true);
      return;
    }

    // Fallback: old-format report — single fat JSON with Issue_Summary / Issue_Details
    // embedded. Slowest path; only happens when the container hasn't been rebuilt yet.
    const validationResStr = await grok.dapi.files.readAsText(jsonPath);
    this.validationResults = JSON.parse(validationResStr) as ValidationResult;
    if (this.validationResults.Issue_Summary?.length)
      this.issueSummaryDf = DG.DataFrame.fromObjects(this.validationResults.Issue_Summary) ?? null;
    if (this.validationResults.Issue_Details?.length)
      this.issueDetailsDf = DG.DataFrame.fromObjects(this.validationResults.Issue_Details) ?? null;
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
    // Issue details now live on a DataFrame (CSV/d42-backed). Read columns directly
    // from the dataframe — no per-row IssueDetail object reconstruction.
    const df = this.issueDetailsDf;
    if (!df || df.rowCount === 0)
      return;

    const datasetCol = df.col('dataset');
    const rowCol = df.col('row');
    if (!datasetCol || !rowCol)
      return;

    // Group dataframe row indices by domain → row-in-domain. We store indices
    // (not records) and look up issue fields directly when applying them.
    const issuesByDomain: {[domain: string]: {[row: number]: number[]}} = {};
    for (let i = 0; i < df.rowCount; i++) {
      const dataset = datasetCol.get(i);
      if (!dataset)
        continue;
      const domain = removeExtension(String(dataset).toLowerCase());

      // The `row` column can be int, NaN, empty string, or a numeric string.
      const rawRow = rowCol.get(i);
      let rowNum = NaN;
      if (typeof rawRow === 'number')
        rowNum = rawRow;
      else if (typeof rawRow === 'string' && rawRow.length)
        rowNum = parseInt(rawRow, 10);
      if (Number.isNaN(rowNum))
        continue;

      if (!issuesByDomain[domain])
        issuesByDomain[domain] = {};
      if (!issuesByDomain[domain][rowNum])
        issuesByDomain[domain][rowNum] = [];
      issuesByDomain[domain][rowNum].push(i);
    }

    for (const domain of this.domains.all())
      this.addValidationColumnsToDomain(domain, issuesByDomain[domain.name.toLowerCase()] ?? {});
  }

  /** Reconstructs an IssueDetail object from a row in issueDetailsDf.
   * Issue_Details no longer lives on validationResults — it's only on the
   * CSV/d42-backed dataframe — so consumers go through this helper.
   * @param {number} row Row index in issueDetailsDf.
   * @return {IssueDetail} Reconstructed record (or an empty object if no df). */
  issueDetailFromRow(row: number): IssueDetail {
    const df = this.issueDetailsDf;
    if (!df)
      return {} as IssueDetail;
    const detail: any = {};
    for (const col of df.columns.names()) {
      const v = df.col(col)?.get(row);
      if (col === 'variables' || col === 'values')
        detail[col] = parseArrayCell(v);
      else
        detail[col] = v;
    }
    return detail as IssueDetail;
  }

  /**
   * Adds validation columns to a specific domain dataframe. Reads issue fields
   * directly from `issueDetailsDf` using the row indices supplied in `issuesByRow`
   * — no IssueDetail object allocation per row.
   * @param {DG.DataFrame} domain - The domain dataframe to add columns to
   * @param {Object} issuesByRow - Map of "row in domain" → array of issueDetailsDf row indices
   */
  private addValidationColumnsToDomain(
    domain: DG.DataFrame,
    issuesByRow: {[row: number]: number[]},
  ): void {
    const detailsDf = this.issueDetailsDf;
    if (!detailsDf)
      return;
    const coreIdCol = detailsDf.col('core_id');
    const messageCol = detailsDf.col('message');
    const variablesCol = detailsDf.col('variables');
    const valuesCol = detailsDf.col('values');

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
    for (const [rowIndexStr, issueIndices] of Object.entries(issuesByRow)) {
      const rowIndex = parseInt(rowIndexStr, 10) - 1; //rows in issues are numbered starting from 1
      if (isNaN(rowIndex) || !issueIndices || issueIndices.length === 0)
        continue;
      // CORE occasionally references row numbers beyond what the domain actually
      // has (stale validation results, mismatched ingest). Skip rather than blow up.
      if (rowIndex < 0 || rowIndex >= domain.rowCount)
        continue;

      hasErrorsCol.set(rowIndex, true);

      // Build the JSON payload for the violated-rules column by reading directly
      // from the issue-details dataframe rather than materialising IssueDetail objects.
      const rowIssuesJson = issueIndices.map((idx) => ({
        core_id: coreIdCol?.get(idx) ?? '',
        message: messageCol?.get(idx) ?? '',
        variables: parseArrayCell(variablesCol?.get(idx)),
        values: parseArrayCell(valuesCol?.get(idx)),
      }));
      violatedRulesCol.set(rowIndex, JSON.stringify(rowIssuesJson));

      // Process issues, collect errors by variable, create columns, and set values in one loop
      const errorsByVariable: {[variable: string]: Array<{ruleID: string, message: string, value: string}>} = {};

      for (const issueIdx of issueIndices) {
        const variables = parseArrayCell(variablesCol?.get(issueIdx));
        const values = parseArrayCell(valuesCol?.get(issueIdx));
        const coreId = String(coreIdCol?.get(issueIdx) ?? '');
        const message = String(messageCol?.get(issueIdx) ?? '');
        if (!variables.length)
          continue;

        for (let varIndex = 0; varIndex < variables.length; varIndex++) {
          const variable = variables[varIndex];
          if (domain.col(variable)?.getTag(COLUMN_FROM_DM_TAG))
            continue;

          const value = varIndex < values.length ? values[varIndex] : '';

          // Create columns for this variable if they don't exist
          ensureVariableColumns(variable);

          // Only process if this variable has an error column (exists in domain)
          if (variableErrorColumns[variable]) {
            if (!errorsByVariable[variable])
              errorsByVariable[variable] = [];
            errorsByVariable[variable].push({
              ruleID: coreId,
              message: message,
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


export class ClinRow {
  row: DG.Row;

  constructor(row: DG.Row) {
    this.row = row;
  }
}

// CDISC CORE only ships rule sets for SDTMIG 3.2, 3.3, 3.4 (and uses `list-rule-sets` to expose them).
// We default to 3.4 when define.xml/standardVersion is missing, and clamp anything < 3.2 up to 3.2
// so that pre-3.2 studies still get validated against the oldest available rule set rather than
// silently producing an empty report.
const MIN_SDTMIG_VERSION = '3.2';
const DEFAULT_SDTMIG_VERSION = '3.4';

function resolveSdtmigVersion(declared?: string): string {
  if (!declared)
    return DEFAULT_SDTMIG_VERSION;
  return compareVersions(declared, MIN_SDTMIG_VERSION) < 0 ? MIN_SDTMIG_VERSION : declared;
}

function compareVersions(a: string, b: string): number {
  const aParts = a.split('.').map((p) => parseInt(p, 10) || 0);
  const bParts = b.split('.').map((p) => parseInt(p, 10) || 0);
  const len = Math.max(aParts.length, bParts.length);
  for (let i = 0; i < len; i++) {
    const diff = (aParts[i] ?? 0) - (bParts[i] ?? 0);
    if (diff !== 0)
      return diff;
  }
  return 0;
}

/** Reads a cell that should hold an array. Returns it as-is when already an
 * array (old fat-JSON fallback path), JSON-parses string cells (CSV / d42 path),
 * and returns [] for empty / unparseable cells.
 * @param {any} v Raw cell value.
 * @return {string[]} The decoded array, or [] on empty/unparseable input. */
function parseArrayCell(v: any): string[] {
  if (Array.isArray(v))
    return v;
  if (typeof v !== 'string' || v.length === 0)
    return [];
  try {
    const parsed = JSON.parse(v);
    return Array.isArray(parsed) ? parsed : [];
  } catch {
    return [];
  }
}
