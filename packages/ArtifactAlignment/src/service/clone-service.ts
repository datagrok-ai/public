import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {getStartedOrNull} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {
  CONFIG_PATH, makeMetaCall,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {serialize, deserialize} from '@datagrok-libraries/utils/src/json-serialization';
import {
  ArtifactType, ARTIFACT_TYPE_FUNCTION_RUN, ARTIFACT_TYPE_WORKFLOW_RUN,
  OPT_FROZEN, OPT_PARENT_CALL_ID, OPT_PUBLICATION_ID,
} from '../domain/constants';

dayjs.extend(utc);

// The Compute2 per-step history marker; a frozen copy must not inherit it, or it
// would show up in the author's step-history panel.
const STEP_HISTORY_OPTION = 'compute2StepHistory';

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

/** Stamps the catalog options every frozen call carries. `isDeleted` reuses the
 * compute history soft-delete convention so frozen copies never appear in the
 * author's history panels (workflow, RFV, or step history) — they stay loadable
 * by id, which is the only way the catalog addresses them. */
function stampFrozen(fc: DG.FuncCall, publicationId: string): void {
  fc.options[OPT_FROZEN] = 'true';
  fc.options[OPT_PUBLICATION_ID] = publicationId;
  fc.options['isDeleted'] = 'true';
  // options is a MapProxy — assign rather than delete; the step-history panel
  // matches on the exact string 'true'.
  if (fc.options[STEP_HISTORY_OPTION] != null)
    fc.options[STEP_HISTORY_OPTION] = 'false';
}

/** Deep-clones a saved Compute2 run into a frozen copy and returns it together with
 * the detected artifact type. A run whose options carry a serialized pipeline config
 * is a workflow run: every step FuncCall is re-saved under a new id with the
 * audience-scoped dataframe grants, and the config is rewritten to the new ids with
 * `isReadonly` forced on every node. Any other saved run (an RFV model run or an
 * individually saved workflow step) is a function run: the call itself is re-saved
 * the same way, with no config to rewrite.
 *
 * Consistency-restriction dataframe refs inside INPUT_RESTRICTIONS are kept as-is:
 * they point at the tables the source run uploaded (already granted at source-save
 * time), so the frozen copy stays loadable without re-uploading them.
 *
 * The frozen funccalls themselves have no server-side ACL — FuncCall is not a full
 * entity. When the platform adds funccall permissions, the grants slot in here. */
export async function cloneRun(sourceCallId: string, options: CloneOptions): Promise<CloneResult> {
  const sourceCall = await historyUtils.loadRun(sourceCallId);
  return sourceCall.options[CONFIG_PATH] != null ?
    cloneWorkflowRun(sourceCall, options) :
    cloneFunctionRun(sourceCall, options);
}

async function cloneFunctionRun(sourceCall: DG.FuncCall, options: CloneOptions): Promise<CloneResult> {
  const sourceStartedOn = getStartedOrNull(sourceCall);
  const sourceAuthorId = sourceCall.author?.id;
  sourceCall.newId();
  stampFrozen(sourceCall, options.publicationId);
  if (options.title != null)
    sourceCall.options['title'] = options.title;
  if (options.description != null)
    sourceCall.options['description'] = options.description;
  const saved = await historyUtils.saveRun(sourceCall, {audience: options.audience});
  return {
    metaCall: saved,
    artifactType: ARTIFACT_TYPE_FUNCTION_RUN,
    nqName: sourceCall.func.nqName,
    childIdMap: {},
    sourceStartedOn,
    sourceAuthorId,
  };
}

async function cloneWorkflowRun(sourceCall: DG.FuncCall, options: CloneOptions): Promise<CloneResult> {
  const state: SerializedNode = deserialize(sourceCall.options[CONFIG_PATH] ?? '{}');
  const nqName = sourceCall.func.nqName;

  const metaCall = await makeMetaCall(nqName);
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

  metaCall.options[CONFIG_PATH] = serialize(state, {useJsonDF: false});
  if (options.title != null)
    metaCall.options['title'] = options.title;
  if (options.description != null)
    metaCall.options['description'] = options.description;
  if (state.version != null)
    metaCall.options['version'] = state.version;
  stampFrozen(metaCall, options.publicationId);

  // await grok.dapi.permissions.grant(metaCall, options.audience, false);
  //   ^ requires core work: FuncCall has no per-entity ACL; the frozen copy is
  //     unlisted-but-public-by-id until the platform enforces funccall permissions.
  // await catalogService.transferOwnership(metaCall, serviceUser);
  //   ^ requires core work: service-identity hand-off; the clone currently stays
  //     owned by the publishing user's session.

  await metaCall.call(undefined, undefined, {processed: true, report: false});
  metaCall.started = dayjs();
  const savedMeta = await historyUtils.saveRun(metaCall);
  return {
    metaCall: savedMeta,
    artifactType: ARTIFACT_TYPE_WORKFLOW_RUN,
    nqName,
    version: state.version,
    childIdMap,
    sourceStartedOn: getStartedOrNull(sourceCall),
    sourceAuthorId: sourceCall.author?.id,
  };
}
