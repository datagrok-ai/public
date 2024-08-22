import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AbstractPipelineParallelConfiguration, AbstractPipelineSequentialConfiguration, AbstractPipelineStaticConfiguration, LoadedPipeline, PipelineActionConfiguraion, PipelineConfigurationInitial, PipelineConfigurationParallelInitial, PipelineConfigurationSequentialInitial, PipelineConfigurationStaticInitial, PipelineHooks, PipelineLinkConfigurationBase, PipelineRefInitial, PipelineSelfRef, PipelineStepConfiguration, StepActionConfiguraion} from './PipelineConfiguration';
import { ItemId, LinkSpecString, NqName} from '../data/common-types';
import {callHandler} from '../utils';
import {LinkParsed, parseLinkIO} from './LinkSpec';

//
// Internal config processing
//

export type FuncallStateItem = {
  id: ItemId;
  type: string;
  direction: 'input' | 'output';
}

type PipelineStepConfigurationInitial = PipelineStepConfiguration<LinkSpecString, never>;
type ConfigInitialTraverseItem = PipelineConfigurationInitial | PipelineStepConfigurationInitial;

export type PipelineConfigurationStaticProcessed = AbstractPipelineStaticConfiguration<LinkParsed[], FuncallStateItem[], PipelineSelfRef>;
export type PipelineConfigurationParallelProcessed = AbstractPipelineParallelConfiguration<LinkParsed[], FuncallStateItem[], PipelineSelfRef>;
export type PipelineConfigurationSequentialProcessed = AbstractPipelineSequentialConfiguration<LinkParsed[], FuncallStateItem[], PipelineSelfRef>;
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

function isStepConfigInitial(c: ConfigInitialTraverseItem): c is PipelineStepConfigurationInitial {
  return !isPipelineStaticInitial(c) && !isPipelineParallelInitial(c) && !isPipelineSequentialInitial(c) && !isPipelineRefInitial(c);
}

function isPipelineConfigInitial(c: ConfigInitialTraverseItem): c is PipelineConfigurationInitial {
  return isPipelineStaticInitial(c) || isPipelineParallelInitial(c) || isPipelineSequentialInitial(c);
}

export async function getProcessedConfig(conf: PipelineConfigurationInitial): Promise<PipelineConfigurationProcessed> {
  const pconf = await configProcessing(conf, new Set());
  return pconf as PipelineConfigurationProcessed;
}

async function configProcessing(conf: ConfigInitialTraverseItem, loadedPipelines: Set<string>): Promise<PipelineConfigurationProcessed | PipelineStepConfiguration<LinkParsed[], FuncallStateItem[]> | PipelineSelfRef> {
  if (isPipelineConfigInitial(conf) && !isPipelineRefInitial(conf) && conf.nqName)
    loadedPipelines.add(conf.nqName);

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
    const pconf = await callHandler<LoadedPipeline>(conf.provider, conf).toPromise();
    if (loadedPipelines.has(pconf.nqName))
      return {id: pconf.id, selfRef: pconf.nqName, type: 'selfRef'};
    loadedPipelines.add(pconf.nqName);
    return configProcessing(pconf, loadedPipelines);
  }
  throw new Error(`Pipeline configuration node type matching failed: ${conf}`);
}

function processStaticConfig(conf: PipelineConfigurationStaticInitial) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, links, actions, hooks};
}

function processParallelConfig(conf: PipelineConfigurationParallelInitial) {
  const actions = processPipelineActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, actions, hooks, stepTypes: []};
}

function processSequentialConfig(conf: PipelineConfigurationSequentialInitial) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? []);
  const hooks = processHooks(conf.hooks ?? {});
  return {...conf, links, actions, hooks};
}

async function processStepConfig(conf: PipelineStepConfiguration<LinkSpecString, never>) {
  const actions = processStepActions(conf.actions ?? []);
  const io = await getFuncCallIO(conf.nqName);
  return {...conf, io, actions};
}

async function getFuncCallIO(nqName: NqName): Promise<FuncallStateItem[]> {
  const func: DG.Func = await grok.functions.eval(nqName);
  const inputs = func.inputs.map((input) => ({id: input.name, type: input.propertyType as any, direction: 'input'} as FuncallStateItem));
  const outputs = func.outputs.map((output) => ({id: output.name, type: output.propertyType as any, direction: 'output'} as FuncallStateItem));
  const io = [...inputs, ...outputs];
  return io;
}

function processPipelineActions(actionsInput: PipelineActionConfiguraion<LinkSpecString>[]) {
  const actions = actionsInput.map((action) => ({...processLinkData(action)}));
  return actions;
}

function processStepActions(actionsInput: StepActionConfiguraion<LinkSpecString>[]) {
  const actions = actionsInput.map((action) => ({...processLinkData(action)}));
  return actions;
}

function processHooks(hooksInput: PipelineHooks<LinkSpecString>) {
  const hooks = Object.fromEntries(
    Object.entries(hooksInput ?? {})?.map(
      ([name, hooks]) => [
        name,
        hooks.map((link) => ({link, ...processLinkData(link)})),
      ] as const,
    ),
  );
  return hooks;
}

function processLinkData<L extends Partial<PipelineLinkConfigurationBase<LinkSpecString>>>(link: L) {
  const from = processLink(link.from ?? [], false, true);
  const to = processLink(link.to ?? [], false, false);
  const base = processLink(link.base ?? [], true, false);
  return {...link, from, to, base};
}

function processLink(io: LinkSpecString, isBase: boolean, isInput: boolean) {
  if (Array.isArray(io))
    return io.map(item => parseLinkIO(item, isBase, isInput));
  else if (io)
    return [parseLinkIO(io, isBase, isInput)];
  else
    return [];
}
