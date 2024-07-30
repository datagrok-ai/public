import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AbstractPipelineParallelConfiguration, AbstractPipelineSequentialConfiguration, AbstractPipelineStaticConfiguration, PipelineActionConfiguraion, PipelineConfigurationInitial, PipelineConfigurationParallelInitial, PipelineConfigurationSequentialInitial, PipelineConfigurationStaticInitial, PipelineGlobalId, PipelineHooks, PipelineLinkConfigurationBase, PipelineRefInitial, PipelineSelfRef, PipelineStepConfiguration} from './PipelineConfiguration';
import {ItemId, ItemPath, ItemPathArray, NqName} from '../data/common-types';
import {callHandler} from '../utils';

//
// Internal config processing
//

function parsePath(path: ItemPath): ItemPathArray {
  return path.split('/').filter((x) => x);
}

function normalizePaths(paths: ItemPath | ItemPath[]): ItemPathArray[] {
  if (Array.isArray(paths))
    return paths.map(parsePath);
  else
    return [parsePath(paths)];
}

export type FuncallStateItem = {
  id: ItemId;
  type: string;
  direction: 'input' | 'output';
}


type ConfigInitialTraverseItem = PipelineConfigurationInitial | PipelineStepConfiguration<never>;

export type PipelineConfigurationStaticProcessed = AbstractPipelineStaticConfiguration<ItemPathArray[], FuncallStateItem[], PipelineSelfRef, string | undefined>;
export type PipelineConfigurationParallelProcessed = AbstractPipelineParallelConfiguration<ItemPathArray[], FuncallStateItem[], PipelineSelfRef, string | undefined>;
export type PipelineConfigurationSequentialProcessed = AbstractPipelineSequentialConfiguration<ItemPathArray[], FuncallStateItem[], PipelineSelfRef, string | undefined>;
export type PipelineConfigurationProcessed = PipelineConfigurationStaticProcessed | PipelineConfigurationParallelProcessed | PipelineConfigurationSequentialProcessed;


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

function isStepConfigInitial(c: ConfigInitialTraverseItem): c is PipelineStepConfiguration<never> {
  return !isPipelineStaticInitial(c) && !isPipelineParallelInitial(c) && isPipelineSequentialInitial(c) && !isPipelineRefInitial(c);
}


export async function getProcessedConfig(conf: PipelineConfigurationInitial): Promise<PipelineConfigurationProcessed> {
  const pconf = await configProcessing(conf, new Set());
  return pconf as PipelineConfigurationProcessed;
}

async function configProcessing(conf: ConfigInitialTraverseItem, loadedPipelines: Set<string>): Promise<PipelineConfigurationProcessed | PipelineStepConfiguration<FuncallStateItem[]> | PipelineSelfRef> {
  if (isStepConfigInitial(conf)) {
    const pconf = await processStepConfig(conf);
    return pconf;
  } else if (isPipelineStaticInitial(conf)) {
    const pconf = processStaticConfig(conf);
    const steps = await Promise.all(conf.steps.map(async (step) => {
      const sconf = await configProcessing(step, loadedPipelines);
      return sconf as any;
    }));
    return {...pconf, steps};
  } else if (isPipelineParallelInitial(conf)) {
    const pconf = processParallelConfig(conf);
    const stepTypes = await Promise.all(conf.stepTypes.map(async (item) => {
      const nconf = await configProcessing(item, loadedPipelines);
      return {...item, ...nconf} as any;
    }));
    return {...pconf, stepTypes};
  } else if (isPipelineSequentialInitial(conf)) {
    const pconf = processSequentialConfig(conf);
    const stepTypes = await Promise.all(conf.stepTypes.map(async (item) => {
      const nconf = await configProcessing(item, loadedPipelines);
      return {...item, ...nconf} as any;
    }));
    return {...pconf, stepTypes};
  } else if (isPipelineRefInitial(conf)) {
    const pconf = await callHandler<PipelineConfigurationInitial & PipelineGlobalId>(conf.provider, conf).toPromise();
    if (loadedPipelines.has(pconf.globalId))
      return {selfRef: pconf.globalId, type: 'selfRef'};
    loadedPipelines.add(pconf.globalId);
    return configProcessing(pconf, loadedPipelines);
  }
  throw new Error(`configProcessing type matching failed: ${conf}`);
}

function processStaticConfig(conf: PipelineConfigurationStaticInitial) {
  const links = conf.links?.map((link) => processLinkBase(link));
  const actions = processActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, links, actions, hooks};
}

function processParallelConfig(conf: PipelineConfigurationParallelInitial) {
  const actions = processActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, actions, hooks, stepTypes: []};
}

function processSequentialConfig(conf: PipelineConfigurationSequentialInitial) {
  const links = conf.links?.map((link) => processLinkBase(link));
  const actions = processActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, links, actions, hooks};
}

async function processStepConfig(conf: PipelineStepConfiguration<never>) {
  if (conf.io == null) {
    const io = await getFuncCallIO(conf.nqName);
    return {...conf, io};
  }
  return conf;
}

async function getFuncCallIO(nqName: NqName): Promise<FuncallStateItem[]> {
  const func: DG.Func = await grok.functions.eval(nqName);
  const inputs = func.inputs.map((input) => ({id: input.name, type: input.propertyType as any, direction: 'input'} as FuncallStateItem));
  const outputs = func.outputs.map((output) => ({id: output.name, type: output.propertyType as any, direction: 'output'} as FuncallStateItem));
  const io = [...inputs, ...outputs];
  return io;
}

function processActions<P extends ItemPath | ItemPath[]>(actionsInput: PipelineActionConfiguraion<P>[]) {
  const actions = actionsInput.map((action) => ({...processLinkBase(action), path: normalizePaths(action.path)}));
  return actions;
}

function processHooks<P extends ItemPath | ItemPath[]>(hooksInput: PipelineHooks<P>) {
  const hooks = Object.fromEntries(Object.entries(hooksInput ?? {})?.map(([name, hooks]) => [name, hooks.map((link) => processLinkBase(link))] as const));
  return hooks;
}

function processLinkBase<L extends Partial<PipelineLinkConfigurationBase<ItemPath | ItemPath[]>>>(link: L) {
  const from = normalizePaths(link.from ?? []);
  const to = normalizePaths(link.to ?? []);
  const dataFrameMutations = Array.isArray(link.dataFrameMutations) ? normalizePaths(link.dataFrameMutations) : !!link.dataFrameMutations;
  const allTypeNodes = Array.isArray(link.allTypeNodes) ? normalizePaths(link.allTypeNodes) : !!link.allTypeNodes;
  return {...link, from, to, dataFrameMutations, allTypeNodes};
}
