import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ARMCD, COL_HAS_ERRORS_POSTFIX, DOMAIN, ERRORS_POSTFIX, HAS_VALIDATION_ERRORS_COL,
  PLANNED_TRT_ARM,
  STUDY_ID, SUBJECT_ID, VIOLATED_RULES_COL} from './constants/columns-constants';
import {calculateLBBaselineColumns} from './data-preparation/data-preparation';
import {removeExtension} from './utils/utils';
import {StudyConfig} from './types/types';
import {SUMMARY_VIEW_NAME} from './constants/view-names-constants';
import {funcs} from './package-api';
import {Subject} from 'rxjs';
import {ValidationResult, IssueDetail} from './types/validation-result';
import {COLUMN_FROM_DM_TAG} from './constants/constants';

export class PreclinicalDomains {
  ag: DG.DataFrame | null = null;
  bg: DG.DataFrame | null = null;
  bw: DG.DataFrame | null = null;
  cl: DG.DataFrame | null = null;
  co: DG.DataFrame | null = null;
  dd: DG.DataFrame | null = null;
  dm: DG.DataFrame | null = null;
  ds: DG.DataFrame | null = null;
  ec: DG.DataFrame | null = null;
  eg: DG.DataFrame | null = null;
  ex: DG.DataFrame | null = null;
  ey: DG.DataFrame | null = null;
  fw: DG.DataFrame | null = null;
  ig: DG.DataFrame | null = null;
  lb: DG.DataFrame | null = null;
  ma: DG.DataFrame | null = null;
  mi: DG.DataFrame | null = null;
  ml: DG.DataFrame | null = null;
  nc: DG.DataFrame | null = null;
  om: DG.DataFrame | null = null;
  pc: DG.DataFrame | null = null;
  pm: DG.DataFrame | null = null;
  po: DG.DataFrame | null = null;
  pp: DG.DataFrame | null = null;
  relrec: DG.DataFrame | null = null;
  sc: DG.DataFrame | null = null;
  se: DG.DataFrame | null = null;
  sm: DG.DataFrame | null = null;
  sr: DG.DataFrame | null = null;
  sv: DG.DataFrame | null = null;
  ta: DG.DataFrame | null = null;
  te: DG.DataFrame | null = null;
  tf: DG.DataFrame | null = null;
  ts: DG.DataFrame | null = null;
  tx: DG.DataFrame | null = null;
  vs: DG.DataFrame | null = null;
  supp: DG.DataFrame[] = [];

  all(): DG.DataFrame[] {
    const dfs = Object.keys(this).filter((it) => it !== 'supp').map((k) => (this as any)[k]).filter((v) => v != null);
    return dfs.concat(this.supp);
  }
}

export class PreclinicalStudy {
  domains: PreclinicalDomains = new PreclinicalDomains();
  subjectsCount: number = 0;
  validationResults: ValidationResult | null = null;
  errorsByDomain: {[key: string]: number} = {};
  studyId: string;
  validated = false;
  initCompleted = false;
  config: StudyConfig;
  views: {[key: string]: DG.ViewBase} = {};
  loadingStudyData: boolean | null = null;
  currentViewName: string = SUMMARY_VIEW_NAME;
  validationCompleted = new Subject<boolean>();
  private validationColumnsSubscriptionSet = false;
  changeViewToSummary = true;
  treatmentAndControlConfig: {treatment: string, control: string}[] = [];
  armCodeToName: {[armcd: string]: string} = {};

  constructor(config: StudyConfig) {
    this.config = config;
    this.studyId = config.name!;
  }

  init() {
    this.process();
    if (!this.validated)
      this.validate().then(() => this.initCompleted = true);
    else
      this.initCompleted = true;
  }

  private process(): void {
    if (this.domains.dm != null) {
      this.subjectsCount = this.domains.dm.rowCount;
      if (!this.studyId && this.domains.dm.col(STUDY_ID))
        this.studyId = this.domains.dm!.col(STUDY_ID)!.get(0);
    }

    this.domains.all().forEach((it) => {
      if (it.name !== 'dm' && it.columns.names().includes(SUBJECT_ID)) {
        const savedName = it.name;
        const columnsFromDm = this.domains.dm!.columns.names().filter((it) => it !== DOMAIN);
        for (const colName of it.columns.names()) {
          const colIdx = columnsFromDm.findIndex((it) => it === colName);
          if (colIdx !== -1)
            columnsFromDm.splice(colIdx, 1);
        }
        grok.data.joinTables(it, this.domains.dm!, [SUBJECT_ID], [SUBJECT_ID], null,
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

    if (this.domains.ta) {
      const armcdCol = this.domains.ta.col(ARMCD);
      const armCol = this.domains.ta.col(PLANNED_TRT_ARM);
      if (armcdCol && armCol) {
        for (let i = 0; i < this.domains.ta.rowCount; i++) {
          const armcd = armcdCol.get(i);
          const arm = armCol.get(i);
          if (armcd && arm && !this.armCodeToName[armcd])
            this.armCodeToName[armcd] = arm;
        }
      }
    }

    if (this.domains.lb)
      calculateLBBaselineColumns(this.domains.lb);
  }

  async validate(): Promise<void> {
    let validationResStr = '';
    if (await grok.dapi.files
      .exists(`System:AppData/Preclinicalcase/SEND/${this.studyId}/validation_results.json`)) {
      validationResStr = await grok.dapi.files
        .readAsText(`System:AppData/Preclinicalcase/SEND/${this.studyId}/validation_results.json`);
    } else {
      validationResStr = await funcs.runCoreValidate(
        'sendig', `Preclinicalcase/SEND/${this.studyId}`, '3.1', 'json', undefined);
      grok.dapi.files
        .writeAsText(`System:AppData/Preclinicalcase/SEND/${this.studyId}/validation_results.json`,
          validationResStr);
    }
    this.validationResults = JSON.parse(validationResStr) as ValidationResult;
    this.validated = true;
    this.validationCompleted.next(true);
  }

  ensureValidationColumnsAdded(): void {
    if (this.validated && this.validationResults)
      this.addValidationColumnsToDomains();
    else if (!this.validationColumnsSubscriptionSet) {
      this.validationCompleted.subscribe(() => {
        this.addValidationColumnsToDomains();
      });
      this.validationColumnsSubscriptionSet = true;
    }
  }

  private addValidationColumnsToDomains(): void {
    if (!this.validationResults || !this.validationResults.Issue_Details)
      return;

    const issuesByDomain: {[domain: string]: {[row: number]: IssueDetail[]}} = {};

    for (const issue of this.validationResults.Issue_Details) {
      if (!issue.dataset)
        continue;
      const domain = removeExtension(issue.dataset.toLowerCase());

      let rowNum: number;
      if (typeof issue.row === 'string') {
        rowNum = parseInt(issue.row, 10);
        if (isNaN(rowNum))
          continue;
      } else
        rowNum = issue.row;

      if (!issuesByDomain[domain])
        issuesByDomain[domain] = {};
      if (!issuesByDomain[domain][rowNum])
        issuesByDomain[domain][rowNum] = [];
      issuesByDomain[domain][rowNum].push(issue);
    }

    const addValidationColumns = (domainName: string, domain: DG.DataFrame) => {
      const domainLower = domainName.toLowerCase();
      const domainIssues = issuesByDomain[domainLower];
      if (!domainIssues) {
        this.addValidationColumnsToDomain(domain, {});
        return;
      }
      this.addValidationColumnsToDomain(domain, domainIssues);
    };

    for (const domain of this.domains.all())
      addValidationColumns(domain.name, domain);
  }

  private addValidationColumnsToDomain(
    domain: DG.DataFrame,
    issuesByRow: {[row: number]: IssueDetail[]},
  ): void {
    const hasErrorsCol = domain.columns.addNewBool(HAS_VALIDATION_ERRORS_COL);
    const violatedRulesCol = domain.columns.addNewString(VIOLATED_RULES_COL);
    violatedRulesCol.semType = 'sdisc-rule-violation';

    const variableErrorColumns: {[variable: string]: DG.Column} = {};
    const variableHasErrorsColumns: {[variable: string]: DG.Column} = {};

    const ensureVariableColumns = (variable: string): void => {
      if (!variableErrorColumns[variable] && domain.col(variable)) {
        const hasErrorsColName = `${variable}${COL_HAS_ERRORS_POSTFIX}`;
        const hasErrCol = domain.columns.addNewBool(hasErrorsColName);
        variableHasErrorsColumns[variable] = hasErrCol;

        const errorColName = `${variable}${ERRORS_POSTFIX}`;
        const errorCol = domain.columns.addNewString(errorColName);
        variableErrorColumns[variable] = errorCol;
      }
    };

    for (const [rowIndexStr, rowIssues] of Object.entries(issuesByRow)) {
      const rowIndex = parseInt(rowIndexStr, 10) - 1;
      if (isNaN(rowIndex) || !rowIssues || rowIssues.length === 0)
        continue;

      hasErrorsCol.set(rowIndex, true);
      violatedRulesCol.set(rowIndex, JSON.stringify(rowIssues));

      const errorsByVariable: {[variable: string]: Array<{ruleID: string, message: string, value: string}>} = {};

      for (const issue of rowIssues) {
        if (issue.variables && issue.values) {
          for (let varIndex = 0; varIndex < issue.variables.length; varIndex++) {
            const variable = issue.variables[varIndex];

            if (domain.col(variable)?.getTag(COLUMN_FROM_DM_TAG))
              continue;

            const value = varIndex < issue.values.length ? issue.values[varIndex] : '';

            ensureVariableColumns(variable);

            if (variableErrorColumns[variable]) {
              if (!errorsByVariable[variable])
                errorsByVariable[variable] = [];
              errorsByVariable[variable].push({
                ruleID: issue.core_id,
                message: issue.message,
                value: value,
              });

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
