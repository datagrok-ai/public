export const SCHEMA = 'artifacts';

export const T_PROGRAM = `${SCHEMA}.program`;
export const T_COMPOUND = `${SCHEMA}.compound`;
export const T_PROGRAM_COMPOUND = `${SCHEMA}.program_compound`;
export const T_STUDY = `${SCHEMA}.study`;
export const T_TAG = `${SCHEMA}.tag`;
export const T_ALIGNMENT = `${SCHEMA}.alignment`;
export const T_ALIGNMENT_COMPOUND = `${SCHEMA}.alignment_compound`;
export const T_ALIGNMENT_TAG = `${SCHEMA}.alignment_tag`;
export const T_HISTORY = `${SCHEMA}.alignment_history`;

export type PublicationStatus = 'pending' | 'approved' | 'rejected';

export const ARTIFACT_TYPE_WORKFLOW_RUN = 'workflow-run';
export const ARTIFACT_TYPE_FUNCTION_RUN = 'function-run';
export type ArtifactType = typeof ARTIFACT_TYPE_WORKFLOW_RUN | typeof ARTIFACT_TYPE_FUNCTION_RUN;

// FuncCall option keys stamped on frozen copies
export const OPT_PUBLICATION_ID = 'artifactsPublicationId';
export const OPT_FROZEN = 'artifactsFrozen';
export const OPT_PARENT_CALL_ID = 'parentCallId';

// Columns whose writes are gated through per-column property schemas
export const APPROVAL_COLUMNS = ['status', 'approved_by', 'approved_on', 'reject_reason'] as const;
export const CURATION_COLUMNS = ['path'] as const;

export const UMBRELLA_APPROVERS = 'Artifacts: Approvers';
export const UMBRELLA_CONTRIBUTORS = 'Artifacts: Contributors';

export const programGroupNames = (code: string) => ({
  viewers: `Artifacts: ${code} Viewers`,
  contributors: `Artifacts: ${code} Contributors`,
  approvers: `Artifacts: ${code} Approvers`,
});

export interface AlignmentRow {
  id: string;
  version: number;
  created_on?: any;
  updated_on?: any;
  author_id?: string;
  publication_id: string;
  revision: number;
  name: string;
  artifact_id: string;
  source_artifact_id?: string;
  program_id: string;
  study_id?: string | null;
  workstream?: string | null;
  artifact_type?: string;
  artifact_author?: string | null;
  artifact_created_on?: any;
  description?: string | null;
  status: PublicationStatus;
  approved_by?: string | null;
  approved_on?: any;
  reject_reason?: string | null;
  path?: string | null;
  tags?: {id: string, name: string}[];
  compounds?: {id: string, name: string}[];
}

export interface HistoryRow {
  id?: string;
  publication_id: string;
  revision: number;
  name: string;
  artifact_id: string;
  source_artifact_id?: string | null;
  program_id: string;
  study_id?: string | null;
  workstream?: string | null;
  artifact_type?: string;
  artifact_author?: string | null;
  artifact_created_on?: any;
  description?: string | null;
  final_status: PublicationStatus;
  approved_by?: string | null;
  approved_on?: any;
  reject_reason?: string | null;
  path?: string | null;
  tags_snapshot?: string | null;
  compounds_snapshot?: string | null;
  superseded_on?: any;
  original_row_id?: string | null;
}
