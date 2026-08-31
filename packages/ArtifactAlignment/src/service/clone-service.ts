import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {deepCopy, getStartedOrNull} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {
  CONFIG_PATH, saveInstanceState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {
  ArtifactType, ARTIFACT_TYPE_FUNCTION_RUN, ARTIFACT_TYPE_WORKFLOW_RUN,
  OPT_FROZEN, OPT_PUBLICATION_ID,
} from '../domain/constants';

dayjs.extend(utc);

export interface CloneOptions {
  /** Group granted View on dataframes the freeze uploads; All users when omitted. */
  audience?: DG.Group;
  publicationId: string;
  title?: string;
  description?: string;
}

export interface CloneResult {
  /** The frozen top-level FuncCall: a new meta call for a workflow run, or the
   * re-saved run itself for a single function run. */
  metaCall: DG.FuncCall;
  artifactType: ArtifactType;
  nqName: string;
  version?: string;
  sourceStartedOn?: dayjs.Dayjs;
  sourceAuthorId?: string;
}

/** Options every frozen call carries; Compute2's history listings filter on OPT_FROZEN
 * (doc § Platform tasks, sub-entity handling). */
function frozenOptions(publicationId: string): Record<string, string> {
  return {
    [OPT_FROZEN]: 'true',
    [OPT_PUBLICATION_ID]: publicationId,
  };
}

/** Freezes a run into a catalog artifact via the personal-history save paths. A saved
 * run whose options carry a pipeline config is a workflow run (a stamped meta call
 * over the SAME config — steps shared, not copied); anything else is a function run
 * (deep copy, audience-granted dataframes) — docs/artifact-alignment.html § Publish &
 * republish flow and § Known design gaps. */
export async function cloneRun(source: string | DG.FuncCall, options: CloneOptions): Promise<CloneResult> {
  const sourceCall = typeof source === 'string' ? await historyUtils.loadRun(source) : source;
  return sourceCall.options[CONFIG_PATH] != null ?
    cloneWorkflowRun(sourceCall, options) :
    cloneFunctionRun(sourceCall, options);
}

async function cloneFunctionRun(sourceCall: DG.FuncCall, options: CloneOptions): Promise<CloneResult> {
  const fc = deepCopy(sourceCall);
  fc.newId();
  for (const [key, value] of Object.entries(frozenOptions(options.publicationId)))
    fc.options[key] = value;
  if (options.title != null)
    fc.options['title'] = options.title;
  if (options.description != null)
    fc.options['description'] = options.description;
  const saved = await historyUtils.saveRun(fc, {audience: options.audience});
  return {
    metaCall: saved,
    artifactType: ARTIFACT_TYPE_FUNCTION_RUN,
    nqName: sourceCall.func.nqName,
    sourceStartedOn: getStartedOrNull(sourceCall),
    sourceAuthorId: sourceCall.author?.id ?? grok.shell.user?.id,
  };
}

async function cloneWorkflowRun(sourceCall: DG.FuncCall, options: CloneOptions): Promise<CloneResult> {
  const state: {version?: string} = deserialize(sourceCall.options[CONFIG_PATH] ?? '{}');
  const savedMeta = await saveInstanceState(sourceCall.func.nqName, state,
    {title: options.title, description: options.description}, state.version,
    {callOptions: frozenOptions(options.publicationId)});
  return {
    metaCall: savedMeta,
    artifactType: ARTIFACT_TYPE_WORKFLOW_RUN,
    nqName: sourceCall.func.nqName,
    version: state.version,
    sourceStartedOn: getStartedOrNull(sourceCall),
    sourceAuthorId: sourceCall.author?.id,
  };
}
