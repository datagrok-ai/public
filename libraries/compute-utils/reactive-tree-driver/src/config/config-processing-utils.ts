import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AbstractPipelineActionConfiguration, AbstractPipelineDynamicConfiguration, AbstractPipelineStaticConfiguration, LoadedPipeline, DataActionConfiguraion, NestedItemContext, PipelineConfigurationInitial, PipelineConfigurationDynamicInitial, PipelineConfigurationStaticInitial, PipelineInitConfiguration, PipelineLinkConfigurationBase, PipelineMutationConfiguration, PipelineRefInitial, PipelineSelfRef, PipelineStepConfiguration, FuncCallActionConfiguration, PipelineReturnConfiguration, PipelineDynamicItem} from './PipelineConfiguration';
import {isDynamicType, ItemId, LinkSpecString, NqName} from '../data/common-types';
import {callHandler, indexFromEnd} from '../utils';
import {LinkIOParsed, LinkSelectorSegment, parseLinkIO} from './LinkSpec';
import {normalizeIdRef} from './PipelineInstance';
import wu from 'wu';
import {getViewersHook} from '../../../shared-utils/utils';
import {DriverLogger, reportError} from '../data/Logger';

//
// Internal config processing
//

export type FuncCallIODescription = {
  id: ItemId;
  type: string;
  nullable: boolean;
  direction: 'input' | 'output';
}

type PipelineStepConfigurationInitial = PipelineStepConfiguration<never>;
type ConfigInitialTraverseItem = PipelineConfigurationInitial | PipelineStepConfigurationInitial | AbstractPipelineActionConfiguration;

export type PipelineConfigurationStaticProcessed = AbstractPipelineStaticConfiguration<FuncCallIODescription[]>;
export type PipelineConfigurationDynamicProcessed = AbstractPipelineDynamicConfiguration<FuncCallIODescription[]>;
export type PipelineConfigurationProcessed = PipelineConfigurationStaticProcessed | PipelineConfigurationDynamicProcessed;

export type IOType = 'input' | 'output' | 'base' | 'actions' | 'not' | 'showWhen' | 'hideWhen';

function isPipelineStaticInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationStaticInitial {
  return !!((c as PipelineConfigurationStaticInitial).type === 'static');
}

function isPipelineDynamicInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationDynamicInitial {
  return isDynamicType((c as PipelineConfigurationDynamicInitial).type);
}

function isPipelineRefInitial(c: ConfigInitialTraverseItem): c is PipelineRefInitial {
  return !!((c as PipelineRefInitial).type === 'ref');
}

function isActionConfigInitial(c: ConfigInitialTraverseItem): c is AbstractPipelineActionConfiguration {
  return (c as AbstractPipelineActionConfiguration).type === 'action';
}

function isStepConfigInitial(c: ConfigInitialTraverseItem): c is PipelineStepConfigurationInitial {
  return !isPipelineStaticInitial(c) && !isPipelineDynamicInitial(c) && !isPipelineRefInitial(c) && !isActionConfigInitial(c);
}

function isPipelineConfigInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationInitial {
  return isPipelineStaticInitial(c) || isPipelineDynamicInitial(c);
}

export async function getProcessedConfig(conf: PipelineConfigurationInitial, logger?: DriverLogger): Promise<PipelineConfigurationProcessed> {
  const pconf = await configProcessing(conf, new Map(), logger);
  return pconf as PipelineConfigurationProcessed;
}

export type PipelineRefStore<T> = Map<string, Map<string | undefined, T>>;

export function addPipelineRef<T>(store: PipelineRefStore<T>, nqName: string, version: string | undefined, item: T) {
  if (!store.has(nqName))
    store.set(nqName, new Map());
  const versionsStore = store.get(nqName)!;
  versionsStore.set(version, item);
}

export function getPipelineRef<T>(store: PipelineRefStore<T>, nqName: string, version: string | undefined): T | undefined {
  if (store.has(nqName))
    return store.get(nqName)!.get(version);
}

export function containsPipelineRef<T>(store: PipelineRefStore<T>, nqName: string, version: string | undefined) {
  if (store.has(nqName))
    return store.get(nqName)!.has(version);
  return false;
}

async function configProcessing(
  conf: ConfigInitialTraverseItem,
  loadedPipelines: PipelineRefStore<null>,
  logger?: DriverLogger,
): Promise<PipelineConfigurationProcessed | PipelineStepConfiguration<FuncCallIODescription[]> | AbstractPipelineActionConfiguration | PipelineSelfRef> {
  if (isPipelineConfigInitial(conf) && !isPipelineRefInitial(conf) && conf.nqName)
    addPipelineRef(loadedPipelines, conf.nqName, conf.version, null);

  if (isActionConfigInitial(conf)) {
    return processActionConfig(conf);
  } else if (isStepConfigInitial(conf)) {
    const pconf = await processStepConfig(conf, logger);
    return pconf;
  } else if (isPipelineStaticInitial(conf)) {
    const pconf = processStaticConfig(conf, logger);
    const steps = await Promise.all(conf.steps.map(async (step) => {
      processUIFlags(step);
      const sconf = await configProcessing(step, loadedPipelines, logger);
      return sconf;
    }));
    checkUniqId(steps, logger);
    return {...pconf, steps};
  } else if (isPipelineDynamicInitial(conf)) {
    const pconf = processDynamicConfig(conf, logger);
    const stepTypes = await Promise.all(conf.stepTypes.map(async (item) => {
      processUIFlags(item);
      const nconf = await configProcessing(item, loadedPipelines, logger);
      return nconf;
    }));
    checkUniqId(stepTypes, logger);
    return {...pconf, stepTypes};
  } else if (isPipelineRefInitial(conf)) {
    const pconf = await callHandler<LoadedPipeline>(conf.provider, conf).toPromise();
    if (containsPipelineRef(loadedPipelines, pconf.nqName, pconf.version))
      return {id: pconf.id, nqName: pconf.nqName, version: pconf.version, type: 'selfRef'};
    addPipelineRef(loadedPipelines, pconf.nqName, pconf.version, null);
    return configProcessing(pconf, loadedPipelines, logger);
  }
  throw new Error(`Pipeline configuration node type matching failed: ${conf}`);
}

function processUIFlags(item: PipelineDynamicItem<never>) {
  if (item.disableUIControlls) {
    item.disableUIAdding = true;
    item.disableUIDragging = true;
    item.disableUIRemoving = true;
  }
}

function processStaticConfig(conf: PipelineConfigurationStaticInitial, logger?: DriverLogger) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? [], logger);
  const onInit = processInitHook(conf.onInit);
  const onReturn = processReturnHook(conf.onReturn);
  const states = conf.states?.map((s) => normalizeIdRef(s));
  return {...conf, links, actions, onInit, onReturn, states};
}

function processDynamicConfig(conf: PipelineConfigurationDynamicInitial, logger?: DriverLogger) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? [], logger);
  const onInit = processInitHook(conf.onInit);
  const onReturn = processReturnHook(conf.onReturn);
  const initialSteps = conf.initialSteps?.map((s) => normalizeIdRef(s));
  const states = conf.states?.map((s) => normalizeIdRef(s));
  return {...conf, actions, links, onInit, onReturn, initialSteps, states};
}

async function processStepConfig(conf: PipelineStepConfiguration<never>, logger?: DriverLogger) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processStepActions(conf.actions ?? [], logger);
  const io = getFuncCallIO(conf.nqName);
  const func = DG.Func.byName(conf.nqName);
  const viewersHookMakerName = getViewersHook(func);
  let viewersHook = conf.viewersHook;
  if (!viewersHook && viewersHookMakerName) {
    const hookMaker = DG.Func.byName(viewersHookMakerName);
    viewersHook = await hookMaker.apply();
  }
  const states = conf.states?.map((s) => normalizeIdRef(s));
  return {...conf, viewersHook, io, links, actions, states};
}

function processActionConfig(conf: AbstractPipelineActionConfiguration & NestedItemContext): PipelineConfigurationStaticProcessed {
  return {
    id: conf.id,
    type: 'static',
    nqName: conf.nqName,
    friendlyName: conf.friendlyName,
    description: conf.description,
    tags: conf.tags,
    steps: [],
    links: [],
    actions: [],
    disableHistory: true,
    isActionStep: true,
    disableUIAdding: conf.disableUIAdding,
    disableUIDragging: conf.disableUIDragging,
    disableUIRemoving: conf.disableUIRemoving,
  } as PipelineConfigurationStaticProcessed;
}

function getFuncCallIO(nqName: NqName): FuncCallIODescription[] {
  const func = DG.Func.byName(nqName);
  if (!func)
    throw new Error(`Function '${nqName}' not found`);
  const fc = func.prepare();
  const inputs = wu(fc.inputParams.values()).map((p) => (
    {id: p.property.name, type: p.property.propertyType as any, direction: 'input' as const, nullable: isOptional(p.property)}
  ));
  const outputs = wu(fc.outputParams.values()).map((p) => (
    {id: p.property.name, type: p.property.propertyType as any, direction: 'output' as const, nullable: false}
  ));
  return [...inputs, ...outputs];
}

function isOptional(prop: DG.Property) {
  return prop.options.optional === 'true';
}

function processPipelineActions(actionsInput: (DataActionConfiguraion<LinkSpecString> | PipelineMutationConfiguration<LinkSpecString> | FuncCallActionConfiguration<LinkSpecString>)[], logger?: DriverLogger) {
  checkUniqId(actionsInput, logger);
  const actions = actionsInput.map((action) => ({...processLinkData(action), ...processActionVisibility(action)}));
  return actions;
}

function processStepActions(actionsInput: (DataActionConfiguraion<LinkSpecString> | FuncCallActionConfiguration<LinkSpecString>)[], logger?: DriverLogger) {
  checkUniqId(actionsInput, logger);
  const actions = actionsInput.map((action) => ({...processLinkData(action), ...processActionVisibility(action)}));
  return actions;
}

function processActionVisibility(action: {showWhen?: LinkSpecString, hideWhen?: LinkSpecString}) {
  const showWhen = processLink(action.showWhen ?? [], 'showWhen');
  const hideWhen = processLink(action.hideWhen ?? [], 'hideWhen');
  return {showWhen, hideWhen};
}

function processReturnHook(hooksInput?: PipelineReturnConfiguration<LinkSpecString>) {
  const hook = (hooksInput ? processLinkData(hooksInput) : undefined) as PipelineReturnConfiguration<LinkIOParsed[]> | undefined;
  return hook;
}

function processInitHook(hooksInput?: PipelineInitConfiguration<LinkSpecString>) {
  const hook = (hooksInput ? processLinkData(hooksInput) : undefined) as PipelineInitConfiguration<LinkIOParsed[]> | undefined;
  return hook;
}

function processLinkData<L extends PipelineLinkConfigurationBase<LinkSpecString>>(link: L) {
  const from = expandDeferredIOs(processLink(link.from ?? [], 'input'), link.id);
  const to = expandDeferredIOs(processLink(link.to ?? [], 'output'), link.id);
  const base = processLink(link.base ?? [], 'base');
  if (base.length > 1)
    throw new Error(`Link ${link.id}: base accepts at most one entry, got ${base.length}.`);
  const not = processLink(link.not ?? [], 'not');
  const actions = processLink(link.actions ?? [], 'actions');
  const linkType = (link as any).type;
  if (linkType !== 'funccall') {
    for (const io of from) {
      if (io.flags?.includes('call') && !io.flags?.includes('optional'))
        throw new Error(`Link ${link.id}: input ${io.name} uses (call) without (optional); only funccall actions allow that.`);
    }
    for (const io of to) {
      if (io.flags?.includes('call'))
        throw new Error(`Link ${link.id}: output ${io.name} uses (call); only funccall actions allow (call) on outputs.`);
    }
  }
  return {...link, from, to, base, not, actions};
}

function processLink(io: LinkSpecString, ioType: IOType) {
  if (Array.isArray(io))
    return io.flatMap((item) => parseLinkIO(item, ioType));
  else if (io)
    return parseLinkIO(io, ioType);
  else
    return [];
}

function checkUniqId(items: {id: string}[], logger?: DriverLogger) {
  const ids = new Set<string>();
  for (const item of items) {
    if (ids.has(item.id))
      reportError('warning', 'configProcessing', `Id ${item.id} is not unique`, logger);
    ids.add(item.id);
  }
}

// ---------------------------------------------------------------------------
// Deferred IO selector expansion
// ---------------------------------------------------------------------------

function expandDeferredIOs(ioList: LinkIOParsed[], linkId: string): LinkIOParsed[] {
  const seenTemplateNames = new Set<string | number>();
  let anonIdx = 0;
  return ioList.flatMap((io) => {
    const lastSeg = indexFromEnd(io.segments);
    if (!lastSeg || lastSeg.type !== 'selector' || !lastSeg.ioExpand) return [io];
    const isAnonymous = io.name === '_';
    const templateName: string | number = isAnonymous ? anonIdx++ : io.name;
    if (!isAnonymous && seenTemplateNames.has(templateName))
      throw new Error(`Link ${linkId}: multiple (template) operators on the same side use the link-io name "${io.name}". Give them distinct prefixes (e.g. a_(template), b_(template)) so they can be addressed individually.`);
    seenTemplateNames.add(templateName);
    const direction: 'input' | 'output' = lastSeg.ioExpand === 'inputs' ? 'input' : 'output';
    const excludeSet = new Set(lastSeg.excludeIds ?? []);
    let targetIO: FuncCallIODescription[];
    try {
      targetIO = getFuncCallIO(lastSeg.nqName!);
    } catch (e) {
      throw new Error(`Link ${linkId}: ${(e as Error).message}`);
    }
    return targetIO
      .filter((d) => d.direction === direction && !excludeSet.has(d.id))
      .map((d) => {
        const nname = isAnonymous ? d.id : io.name + d.id;
        const nlastSegment: LinkSelectorSegment = {type: 'selector', selector: 'first', ids: [d.id], stopIds: []};
        return {name: nname, segments: [...io.segments.slice(0, -1), nlastSegment], flags: io.flags, templateName};
      });
  });
}
