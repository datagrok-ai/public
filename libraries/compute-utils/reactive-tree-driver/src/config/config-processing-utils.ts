import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AbstractPipelineActionConfiguration, AbstractPipelineDynamicConfiguration, AbstractPipelineStaticConfiguration, LoadedPipeline, DataActionConfiguraion, PipelineConfigurationInitial, PipelineConfigurationDynamicInitial, PipelineConfigurationStaticInitial, PipelineInitConfiguration, PipelineLinkConfigurationBase, PipelineMutationConfiguration, PipelineRefInitial, PipelineSelfRef, PipelineStepConfiguration, FuncCallActionConfiguration, PipelineReturnConfiguration, PipelineDynamicItem} from './PipelineConfiguration';
import {ItemId, LinkSpecString, NqName} from '../data/common-types';
import {callHandler} from '../utils';
import {LinkIOParsed, parseLinkIO} from './LinkSpec';
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

type PipelineStepConfigurationInitial = PipelineStepConfiguration<LinkSpecString, never>;
type ConfigInitialTraverseItem = PipelineConfigurationInitial | PipelineStepConfigurationInitial | AbstractPipelineActionConfiguration;

export type PipelineConfigurationStaticProcessed = AbstractPipelineStaticConfiguration<LinkIOParsed[], FuncCallIODescription[], PipelineSelfRef>;
export type PipelineConfigurationDynamicProcessed = AbstractPipelineDynamicConfiguration<LinkIOParsed[], FuncCallIODescription[], PipelineSelfRef>;
/** @deprecated Use PipelineConfigurationDynamicProcessed */
export type PipelineConfigurationParallelProcessed = PipelineConfigurationDynamicProcessed;
/** @deprecated Use PipelineConfigurationDynamicProcessed */
export type PipelineConfigurationSequentialProcessed = PipelineConfigurationDynamicProcessed;
export type PipelineConfigurationProcessed = PipelineConfigurationStaticProcessed | PipelineConfigurationDynamicProcessed;

export type IOType = 'input' | 'output' | 'base' | 'actions' | 'not';

function isPipelineStaticInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationStaticInitial {
  return !!((c as PipelineConfigurationStaticInitial).type === 'static');
}

function isPipelineDynamicInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationDynamicInitial {
  const type = (c as PipelineConfigurationDynamicInitial).type;
  return type === 'dynamic' || type === 'parallel' || type === 'sequential';
}

/** @deprecated Use isPipelineDynamicInitial */
const isPipelineParallelInitial = isPipelineDynamicInitial;
/** @deprecated Use isPipelineDynamicInitial */
const isPipelineSequentialInitial = isPipelineDynamicInitial;

function isPipelineRefInitial(c: ConfigInitialTraverseItem): c is PipelineRefInitial {
  return !!((c as PipelineRefInitial).type === 'ref');
}

function isActionConfigInitial(c: ConfigInitialTraverseItem): c is AbstractPipelineActionConfiguration {
  return (c as AbstractPipelineActionConfiguration).type === 'action';
}

function isStepConfigInitial(c: ConfigInitialTraverseItem): c is PipelineStepConfigurationInitial {
  return !isPipelineStaticInitial(c) && !isPipelineParallelInitial(c) && !isPipelineSequentialInitial(c) && !isPipelineRefInitial(c) && !isActionConfigInitial(c);
}

function isPipelineConfigInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationInitial {
  return isPipelineStaticInitial(c) || isPipelineParallelInitial(c) || isPipelineSequentialInitial(c);
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
): Promise<PipelineConfigurationProcessed | PipelineStepConfiguration<LinkIOParsed[], FuncCallIODescription[]> | AbstractPipelineActionConfiguration | PipelineSelfRef> {
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

function processUIFlags(item: PipelineDynamicItem<LinkSpecString, never, PipelineRefInitial>) {
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
  return {...conf, links, actions, onInit, onReturn};
}

function processDynamicConfig(conf: PipelineConfigurationDynamicInitial, logger?: DriverLogger) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? [], logger);
  const onInit = processInitHook(conf.onInit);
  const onReturn = processReturnHook(conf.onReturn);
  return {...conf, actions, links, onInit, onReturn};
}

async function processStepConfig(conf: PipelineStepConfiguration<LinkSpecString, never>, logger?: DriverLogger) {
  const actions = processStepActions(conf.actions ?? [], logger);
  const io = await getFuncCallIO(conf.nqName);
  const func = DG.Func.byName(conf.nqName);
  const viewersHookMakerName = getViewersHook(func);
  let viewersHook = conf.viewersHook;
  if (!viewersHook && viewersHookMakerName) {
    const hookMaker = DG.Func.byName(viewersHookMakerName);
    viewersHook = await hookMaker.apply();
  }
  return {...conf, viewersHook, io, actions};
}

function processActionConfig(conf: AbstractPipelineActionConfiguration): PipelineConfigurationStaticProcessed {
  return {
    id: conf.id,
    type: 'static',
    friendlyName: conf.friendlyName,
    tags: conf.tags,
    steps: [],
    links: [],
    actions: [],
    disableHistory: true,
    isActionStep: true,
  } as PipelineConfigurationStaticProcessed;
}

async function getFuncCallIO(nqName: NqName): Promise<FuncCallIODescription[]> {
  const func = DG.Func.byName(nqName);
  const fc = func.prepare();
  const inputs = wu(fc.inputParams.values()).map((input) => (
    {id: input.property.name, type: input.property.propertyType as any, direction: 'input' as const, nullable: isOptional(input.property)}
  ));
  const outputs = wu(fc.outputParams.values()).map((output) => (
    {id: output.property.name, type: output.property.propertyType as any, direction: 'output' as const, nullable: false}
  ));
  const io = [...inputs, ...outputs];
  return io;
}

function isOptional(prop: DG.Property) {
  return prop.options.optional === 'true';
}

function processPipelineActions(actionsInput: (DataActionConfiguraion<LinkSpecString> | PipelineMutationConfiguration<LinkSpecString> | FuncCallActionConfiguration<LinkSpecString>)[], logger?: DriverLogger) {
  checkUniqId(actionsInput, logger);
  const actions = actionsInput.map((action) => ({...processLinkData(action)}));
  return actions;
}

function processStepActions(actionsInput: (DataActionConfiguraion<LinkSpecString> | FuncCallActionConfiguration<LinkSpecString>)[], logger?: DriverLogger) {
  checkUniqId(actionsInput, logger);
  const actions = actionsInput.map((action) => ({...processLinkData(action)}));
  return actions;
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
  const from = processLink(link.from ?? [], 'input');
  const to = processLink(link.to ?? [], 'output');
  const base = processLink(link.base ?? [], 'base');
  const not = processLink(link.not ?? [], 'not');
  const actions = processLink(link.actions ?? [], 'actions');
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
