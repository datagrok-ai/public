import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {deepCopy, getStartedOrNull} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {
  CONFIG_PATH, makeMetaCall, saveInstanceState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {
  ArtifactType, ARTIFACT_TYPE_FUNCTION_RUN, ARTIFACT_TYPE_WORKFLOW_RUN,
  OPT_FROZEN, OPT_PARENT_CALL_ID, OPT_PUBLICATION_ID,
} from '../domain/constants';

dayjs.extend(utc);

interface SerializedNode {
  type?: string;
  funcCallId?: string;
  isReadonly?: boolean;
  steps?: SerializedNode[];
  nqName?: string;
  version?: string;
}

async function walkFuncCallNodes(node: SerializedNode, fn: (leaf: SerializedNode) => Promise<void>): Promise<void> {
  if (node == null || typeof node !== 'object')
    return;
  if (node.type === 'funccall') {
    await fn(node);
    return;
  }
  node.isReadonly = true;
  for (const step of node.steps ?? [])
    await walkFuncCallNodes(step, fn);
}

export interface CloneOptions {
  /** Group granted View on the frozen copy's uploaded dataframes; All users when omitted. */
  audience?: DG.Group;
  publicationId: string;
  title?: string;
  description?: string;
}

export interface CloneResult {
  /** The frozen top-level FuncCall: the meta call of a workflow run, or the run itself
   * for a single function run. */
  metaCall: DG.FuncCall;
  artifactType: ArtifactType;
  nqName: string;
  version?: string;
  childIdMap: {[oldId: string]: string};
  sourceStartedOn?: dayjs.Dayjs;
  sourceAuthorId?: string;
}

/** The catalog options every frozen call carries. Compute2's history listings filter
 * out runs carrying OPT_FROZEN, so frozen copies stay reachable by id only (how the
 * catalog addresses them); proper hiding awaits core sub-entity handling. */
function frozenOptions(publicationId: string): Record<string, string> {
  return {
    [OPT_FROZEN]: 'true',
    [OPT_PUBLICATION_ID]: publicationId,
  };
}

function stampFrozen(fc: DG.FuncCall, publicationId: string): void {
  for (const [key, value] of Object.entries(frozenOptions(publicationId)))
    fc.options[key] = value;
}

/** Freezes a Compute2 run into a frozen copy and returns it together with the
 * detected artifact type. A saved run id whose options carry a serialized pipeline
 * config is a workflow run: every step FuncCall is re-saved under a new id with the
 * audience-scoped dataframe grants, and the config is rewritten to the new ids with
 * `isReadonly` forced on every node. Any other source — a saved run without a config,
 * or a live in-memory FuncCall (an RFV model run or a workflow step, never mutated) —
 * is a function run: a deep copy is saved the same way, with no config to rewrite.
 *
 * Consistency-restriction dataframe refs inside INPUT_RESTRICTIONS are kept as-is:
 * they point at the tables the source run uploaded (already granted at source-save
 * time), so the frozen copy stays loadable without re-uploading them.
 *
 * The frozen funccalls themselves have no server-side ACL — FuncCall is not a full
 * entity. When the platform adds funccall permissions, the grants slot in here. */
export async function cloneRun(source: string | DG.FuncCall, options: CloneOptions): Promise<CloneResult> {
  const sourceCall = typeof source === 'string' ? await historyUtils.loadRun(source) : source;
  return sourceCall.options[CONFIG_PATH] != null ?
    cloneWorkflowRun(sourceCall, options) :
    cloneFunctionRun(sourceCall, options);
}

async function cloneFunctionRun(sourceCall: DG.FuncCall, options: CloneOptions): Promise<CloneResult> {
  const fc = deepCopy(sourceCall);
  fc.newId();
  stampFrozen(fc, options.publicationId);
  if (options.title != null)
    fc.options['title'] = options.title;
  if (options.description != null)
    fc.options['description'] = options.description;
  const saved = await historyUtils.saveRun(fc, {audience: options.audience});
  return {
    metaCall: saved,
    artifactType: ARTIFACT_TYPE_FUNCTION_RUN,
    nqName: sourceCall.func.nqName,
    childIdMap: {},
    sourceStartedOn: getStartedOrNull(sourceCall),
    sourceAuthorId: sourceCall.author?.id ?? grok.shell.user?.id,
  };
}

async function cloneWorkflowRun(sourceCall: DG.FuncCall, options: CloneOptions): Promise<CloneResult> {
  const state: SerializedNode = deserialize(sourceCall.options[CONFIG_PATH] ?? '{}');
  // minted up front — the frozen children carry it as their parentCallId
  const metaCall = await makeMetaCall(sourceCall.func.nqName);
  metaCall.newId();
  const frozenMetaId = metaCall.id;

  const childIdMap: {[oldId: string]: string} = {};
  await walkFuncCallNodes(state, async (leaf) => {
    leaf.isReadonly = true;
    if (leaf.funcCallId == null)
      return;
    const oldId = leaf.funcCallId;
    if (childIdMap[oldId] == null) {
      const fc = await historyUtils.loadRun(oldId);
      fc.newId();
      stampFrozen(fc, options.publicationId);
      fc.options[OPT_PARENT_CALL_ID] = frozenMetaId;
      const saved = await historyUtils.saveRun(fc, {audience: options.audience});
      childIdMap[oldId] = saved.id;
    }
    leaf.funcCallId = childIdMap[oldId];
  });

  const savedMeta = await saveInstanceState(sourceCall.func.nqName, state,
    {title: options.title, description: options.description}, state.version,
    {metaCall, callOptions: frozenOptions(options.publicationId)});
  return {
    metaCall: savedMeta,
    artifactType: ARTIFACT_TYPE_WORKFLOW_RUN,
    nqName: sourceCall.func.nqName,
    version: state.version,
    childIdMap,
    sourceStartedOn: getStartedOrNull(sourceCall),
    sourceAuthorId: sourceCall.author?.id,
  };
}
