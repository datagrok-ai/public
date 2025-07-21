import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AbstractPipelineParallelConfiguration, AbstractPipelineSequentialConfiguration, AbstractPipelineStaticConfiguration, LoadedPipeline, DataActionConfiguraion, PipelineConfigurationInitial, PipelineConfigurationParallelInitial, PipelineConfigurationSequentialInitial, PipelineConfigurationStaticInitial, PipelineInitConfiguration, PipelineLinkConfigurationBase, PipelineMutationConfiguration, PipelineRefInitial, PipelineSelfRef, PipelineStepConfiguration, FuncCallActionConfiguration, PipelineReturnConfiguration} from './PipelineConfiguration';
import {ItemId, LinkSpecString, NqName} from '../data/common-types';
import {callHandler} from '../utils';
import {LinkIOParsed, parseLinkIO} from './LinkSpec';
import wu from 'wu';

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
type ConfigInitialTraverseItem = PipelineConfigurationInitial | PipelineStepConfigurationInitial;

export type PipelineConfigurationStaticProcessed = AbstractPipelineStaticConfiguration<LinkIOParsed[], FuncCallIODescription[], PipelineSelfRef>;
export type PipelineConfigurationParallelProcessed = AbstractPipelineParallelConfiguration<LinkIOParsed[], FuncCallIODescription[], PipelineSelfRef>;
export type PipelineConfigurationSequentialProcessed = AbstractPipelineSequentialConfiguration<LinkIOParsed[], FuncCallIODescription[], PipelineSelfRef>;
export type PipelineConfigurationProcessed = PipelineConfigurationStaticProcessed | PipelineConfigurationParallelProcessed | PipelineConfigurationSequentialProcessed;

export type IOType = 'input' | 'output' | 'base' | 'actions' | 'not';

function isPipelineStaticInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationStaticInitial {
  return !!((c as PipelineConfigurationStaticInitial).type === 'static');
}

function isPipelineParallelInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationParallelInitial {
  return !!((c as PipelineConfigurationParallelInitial).type === 'parallel');
}

function isPipelineSequentialInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationSequentialInitial {
  return !!((c as PipelineConfigurationSequentialInitial).type === 'sequential');
}

function isPipelineRefInitial(c: ConfigInitialTraverseItem): c is PipelineRefInitial {
  return !!((c as PipelineRefInitial).type === 'ref');
}

function isStepConfigInitial(c: ConfigInitialTraverseItem): c is PipelineStepConfigurationInitial {
  return !isPipelineStaticInitial(c) && !isPipelineParallelInitial(c) && !isPipelineSequentialInitial(c) && !isPipelineRefInitial(c);
}

function isPipelineConfigInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationInitial {
  return isPipelineStaticInitial(c) || isPipelineParallelInitial(c) || isPipelineSequentialInitial(c);
}

export async function getProcessedConfig(conf: PipelineConfigurationInitial): Promise<PipelineConfigurationProcessed> {
  const pconf = await configProcessing(conf, new Map());
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
): Promise<PipelineConfigurationProcessed | PipelineStepConfiguration<LinkIOParsed[], FuncCallIODescription[]> | PipelineSelfRef> {
  if (isPipelineConfigInitial(conf) && !isPipelineRefInitial(conf) && conf.nqName)
    addPipelineRef(loadedPipelines, conf.nqName, conf.version, null);

  if (isStepConfigInitial(conf)) {
    const pconf = await processStepConfig(conf);
    return pconf;
  } else if (isPipelineStaticInitial(conf)) {
    const pconf = processStaticConfig(conf);
    const steps = await Promise.all(conf.steps.map(async (step) => {
      const sconf = await configProcessing(step, loadedPipelines);
      return sconf;
    }));
    checkUniqId(steps);
    return {...pconf, steps};
  } else if (isPipelineParallelInitial(conf)) {
    const pconf = processParallelConfig(conf);
    const stepTypes = await Promise.all(conf.stepTypes.map(async (item) => {
      const nconf = await configProcessing(item, loadedPipelines);
      return nconf;
    }));
    checkUniqId(stepTypes);
    return {...pconf, stepTypes};
  } else if (isPipelineSequentialInitial(conf)) {
    const pconf = processSequentialConfig(conf);
    const stepTypes = await Promise.all(conf.stepTypes.map(async (item) => {
      const nconf = await configProcessing(item, loadedPipelines);
      return nconf;
    }));
    checkUniqId(stepTypes);
    return {...pconf, stepTypes};
  } else if (isPipelineRefInitial(conf)) {
    const pconf = await callHandler<LoadedPipeline>(conf.provider, conf).toPromise();
    if (containsPipelineRef(loadedPipelines, pconf.nqName, pconf.version))
      return {id: pconf.id, nqName: pconf.nqName, version: pconf.version, type: 'selfRef'};
    addPipelineRef(loadedPipelines, pconf.nqName, pconf.version, null);
    return configProcessing(pconf, loadedPipelines);
  }
  throw new Error(`Pipeline configuration node type matching failed: ${conf}`);
}

function processStaticConfig(conf: PipelineConfigurationStaticInitial) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? []);
  const onInit = processInitHook(conf.onInit);
  const onReturn = processReturnHook(conf.onReturn);
  return {...conf, links, actions, onInit, onReturn};
}

function processParallelConfig(conf: PipelineConfigurationParallelInitial) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? []);
  const onInit = processInitHook(conf.onInit);
  const onReturn = processReturnHook(conf.onReturn);
  return {...conf, actions, links, onInit, onReturn};
}

function processSequentialConfig(conf: PipelineConfigurationSequentialInitial) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? []);
  const onInit = processInitHook(conf.onInit);
  const onReturn = processReturnHook(conf.onReturn);
  return {...conf, links, actions, onInit, onReturn};
}

async function processStepConfig(conf: PipelineStepConfiguration<LinkSpecString, never>) {
  const actions = processStepActions(conf.actions ?? []);
  const io = await getFuncCallIO(conf.nqName);
  return {...conf, io, actions};
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

function processPipelineActions(actionsInput: (DataActionConfiguraion<LinkSpecString> | PipelineMutationConfiguration<LinkSpecString> | FuncCallActionConfiguration<LinkSpecString>)[]) {
  checkUniqId(actionsInput);
  const actions = actionsInput.map((action) => ({...processLinkData(action)}));
  return actions;
}

function processStepActions(actionsInput: (DataActionConfiguraion<LinkSpecString> | FuncCallActionConfiguration<LinkSpecString>)[]) {
  checkUniqId(actionsInput);
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

function checkUniqId(items: {id: string}[]) {
  const ids = new Set<string>();
  for (const item of items) {
    if (ids.has(item.id)) {
      const msg = `Id ${item.id} is not unique`;
      console.error(msg);
      grok.shell.error(msg);
    }
    ids.add(item.id);
  }
}
