import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AbstractPipelineConfiguration, AbstractPipelineParallelConfiguration, AbstractPipelineSequentialConfiguration, AbstractPipelineStaticConfiguration, PipelineActionConfiguraion, PipelineConfiguration, PipelineHooks, PipelineLinkConfigurationBase, PipelineRef, PipelineStepConfiguration} from './PipelineConfiguration';
import {ItemId, ItemPath, NqName} from '../data/common-types';
import {callHandler} from '../utils';

//
// Internal config processing
//

export type ItemPathArray = string[];

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

export type PipelineConfigurationProcessed = AbstractPipelineConfiguration<ItemPathArray[], FuncallStateItem[], never>;
type ConfigTypes = PipelineConfiguration | PipelineConfigurationProcessed;

type CArg1<C extends ConfigTypes> = C extends AbstractPipelineConfiguration<infer A, any, any> ? A : never;
type CArg2<C extends ConfigTypes> = C extends AbstractPipelineConfiguration<any, infer A, any> ? A : never;
type CArg3<C extends ConfigTypes> = C extends AbstractPipelineConfiguration<any, any, infer A> ? A : never;

export type StepConfigType<C extends ConfigTypes> = PipelineStepConfiguration<CArg2<C>>;
export type TraverseItem<C extends ConfigTypes> = C | StepConfigType<C> | PipelineRef;

export type PipelineStaticConfiguration<C extends ConfigTypes> = AbstractPipelineStaticConfiguration<CArg1<C>, CArg2<C>, CArg3<C>>;
export type PipelineParallelConfiguration<C extends ConfigTypes> = AbstractPipelineParallelConfiguration<CArg1<C>, CArg2<C>, CArg3<C>>;
export type PipelineSequentialConfiguration<C extends ConfigTypes> = AbstractPipelineSequentialConfiguration<CArg1<C>, CArg2<C>, CArg3<C>>;

export function isPipelineStaticConfig<C extends ConfigTypes>(c: TraverseItem<C>): c is PipelineStaticConfiguration<C> {
  return !!((c as PipelineStaticConfiguration<C>).steps);
}

export function isPipelineParallelConfig<C extends ConfigTypes>(c: TraverseItem<C>): c is PipelineParallelConfiguration<C> {
  return !!((c as PipelineParallelConfiguration<C>).dynamic === 'parallel');
}

export function isPipelineSequentialConfig<C extends ConfigTypes>(c: TraverseItem<C>): c is PipelineSequentialConfiguration<C> {
  return !!((c as PipelineSequentialConfiguration<C>).dynamic === 'sequential');
}

export function isPipelineConfig<C extends ConfigTypes>(c: TraverseItem<C>): c is C {
  return isPipelineStaticConfig(c) || isPipelineParallelConfig(c) || isPipelineSequentialConfig(c);
}

export function isPipelineRef<C extends ConfigTypes>(c: TraverseItem<C>): c is PipelineRef {
  return !!((c as PipelineRef).provider);
}

export function isStepConfig<C extends ConfigTypes>(c: TraverseItem<C>): c is StepConfigType<C> {
  return !isPipelineConfig(c) && !isPipelineRef(c);
}

export async function getProcessedConfig(conf: PipelineConfiguration): Promise<PipelineConfigurationProcessed> {
  const pconf = await configProcessing(conf);
  return pconf as PipelineConfigurationProcessed;
}

async function configProcessing<C extends TraverseItem<PipelineConfiguration>>(conf: C): Promise<TraverseItem<PipelineConfigurationProcessed>> {
  if (isStepConfig(conf)) {
    const pconf = await processStepConfig(conf);
    return pconf;
  } else if (isPipelineStaticConfig(conf)) {
    const pconf = await processStaticConfig(conf);
    pconf.steps = await Promise.all(conf.steps.map(async (step) => {
      const sconf = await configProcessing(step);
      return sconf as any;
    }));
    return pconf;
  } else if (isPipelineParallelConfig(conf)) {
    const pconf = await processParallelConfig(conf);
    pconf.stepType = await Promise.all(conf.stepType.map(async (item) => {
      const nconf = await configProcessing(item.config);
      return {...item, ...nconf} as any;
    }));
    return pconf;
  } else if (isPipelineSequentialConfig(conf)) {
    const pconf = await processSequentialConfig(conf);
    pconf.stepTypes = await Promise.all(conf.stepTypes.map(async (item) => {
      const nconf = await configProcessing(item.config);
      return {...item, ...nconf} as any;
    }));
    return pconf;
  } else if (isPipelineRef(conf)) {
    const pconf = await callHandler(conf.provider, conf).toPromise();
    return configProcessing(pconf);
  }
  throw new Error(`configProcessing type matching failed: ${conf}`);
}

async function processStaticConfig<C extends PipelineConfiguration>(
  conf: PipelineStaticConfiguration<C>,
): Promise<PipelineStaticConfiguration<PipelineConfigurationProcessed>> {
  const links = conf.links?.map((link) => processLinkBase(link));
  const actions = processActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, links, actions, hooks, steps: []};
}

async function processParallelConfig<C extends PipelineConfiguration>(
  conf: PipelineParallelConfiguration<C>,
): Promise<PipelineParallelConfiguration<PipelineConfigurationProcessed>> {
  const actions = processActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, actions, hooks, stepType: []};
}

async function processSequentialConfig<C extends PipelineConfiguration>(
  conf: PipelineSequentialConfiguration<C>,
): Promise<PipelineSequentialConfiguration<PipelineConfigurationProcessed>> {
  const links = conf.links?.map((link) => processLinkBase(link));
  const actions = processActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, links, actions, hooks, stepTypes: []};
}

async function processStepConfig(conf: StepConfigType<ConfigTypes>): Promise<StepConfigType<PipelineConfigurationProcessed>> {
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
