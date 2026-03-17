export interface ConformanceDetails {
  CORE_Engine_Version: string;
  Report_Generation: string;
  Total_Runtime: string;
  Standard: string;
  Version: string;
  Substandard: string | null;
  CT_Version: string;
  Define_XML_Version: string | null;
  Issue_Limit_Per_Rule: number;
  Issue_Limit_Per_Dataset: boolean;
  UNII_Version: string | null;
  'Med-RT_Version': string | null;
  Meddra_Version: string | null;
  WHODRUG_Version: string | null;
  LOINC_Version: string | null;
  SNOMED_Version: string | null;
}

export interface DatasetDetail {
  filename: string;
  label: string;
  path: string;
  modification_date: string;
  size_kb: number;
  length: number;
}

export interface IssueSummary {
  dataset: string;
  core_id: string;
  message: string;
  issues: number;
}

export interface IssueDetail {
  core_id: string;
  message: string;
  executability: string;
  dataset: string;
  USUBJID: string;
  row: string | number;
  SEQ: string | number;
  variables: string[];
  values: string[];
}

export interface RuleReport {
  core_id: string;
  version: string;
  cdisc_rule_id: string;
  fda_rule_id: string;
  message: string;
  status: 'SUCCESS' | 'SKIPPED' | 'FAILED';
}

export interface VariableError {
  ruleID: string;
  message: string;
  value: string;
  isContext: boolean;
}

export interface ValidationResult {
  Conformance_Details: ConformanceDetails;
  Dataset_Details: DatasetDetail[];
  Issue_Summary: IssueSummary[];
  Issue_Details: IssueDetail[];
  Rules_Report: RuleReport[];
}
