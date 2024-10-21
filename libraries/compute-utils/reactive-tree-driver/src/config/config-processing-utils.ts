import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {AbstractPipelineParallelConfiguration, AbstractPipelineSequentialConfiguration, AbstractPipelineStaticConfiguration, LoadedPipeline, PipelineActionConfiguraion, PipelineConfigurationInitial, PipelineConfigurationParallelInitial, PipelineConfigurationSequentialInitial, PipelineConfigurationStaticInitial, PipelineHookConfiguration, PipelineLinkConfigurationBase, PipelineMutationConfiguration, PipelineRefInitial, PipelineSelfRef, PipelineStepConfiguration, StepActionConfiguraion} from './PipelineConfiguration';
import {ItemId, LinkSpecString, NqName} from '../data/common-types';
import {callHandler} from '../utils';
import {LinkIOParsed, parseLinkIO} from './LinkSpec';

//
// Internal config processing
//

export type FuncallStateItem = {
  id: ItemId;
  type: string;
  nullable: boolean;
  direction: 'input' | 'output';
}

type PipelineStepConfigurationInitial = PipelineStepConfiguration<LinkSpecString, never>;
type ConfigInitialTraverseItem = PipelineConfigurationInitial | PipelineStepConfigurationInitial;

export type PipelineConfigurationStaticProcessed = AbstractPipelineStaticConfiguration<LinkIOParsed[], FuncallStateItem[], PipelineSelfRef>;
export type PipelineConfigurationParallelProcessed = AbstractPipelineParallelConfiguration<LinkIOParsed[], FuncallStateItem[], PipelineSelfRef>;
export type PipelineConfigurationSequentialProcessed = AbstractPipelineSequentialConfiguration<LinkIOParsed[], FuncallStateItem[], PipelineSelfRef>;
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

async function configProcessing(
  conf: ConfigInitialTraverseItem,
  loadedPipelines: Set<string>,
): Promise<PipelineConfigurationProcessed | PipelineStepConfiguration<LinkIOParsed[], FuncallStateItem[]> | PipelineSelfRef> {
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
  const onInit = processHook(conf.onInit);
  return {...conf, links, actions, onInit};
}

function processParallelConfig(conf: PipelineConfigurationParallelInitial) {
  const actions = processPipelineActions(conf.actions ?? []);
  const onInit = processHook(conf.onInit);
  return {...conf, actions, onInit, stepTypes: []};
}

function processSequentialConfig(conf: PipelineConfigurationSequentialInitial) {
  const links = conf.links?.map((link) => processLinkData(link));
  const actions = processPipelineActions(conf.actions ?? []);
  const onInit = processHook(conf.onInit);
  return {...conf, links, actions, onInit};
}

async function processStepConfig(conf: PipelineStepConfiguration<LinkSpecString, never>) {
  const actions = processStepActions(conf.actions ?? []);
  const io = await getFuncCallIO(conf.nqName);
  return {...conf, io, actions};
}

async function getFuncCallIO(nqName: NqName): Promise<FuncallStateItem[]> {
  const func: DG.Func = await grok.functions.eval(nqName);
  // force all io to be non-null
  const inputs = func.inputs.map((input) => (
    {id: input.name, type: input.propertyType as any, direction: 'input' as const, nullable: false}
  ));
  const outputs = func.outputs.map((output) => (
    {id: output.name, type: output.propertyType as any, direction: 'output' as const, nullable: false}
  ));
  const io = [...inputs, ...outputs];
  return io;
}

function processPipelineActions(actionsInput: (PipelineActionConfiguraion<LinkSpecString> | PipelineMutationConfiguration<LinkSpecString>)[]) {
  const actions = actionsInput.map((action) => ({...processLinkData(action)}));
  return actions;
}

function processStepActions(actionsInput: StepActionConfiguraion<LinkSpecString>[]) {
  const actions = actionsInput.map((action) => ({...processLinkData(action)}));
  return actions;
}

function processHook(hooksInput?: PipelineHookConfiguration<LinkSpecString>) {
  return hooksInput ? processLinkData(hooksInput) : undefined;
}

function processLinkData<L extends Partial<PipelineLinkConfigurationBase<LinkSpecString>>>(link: L) {
  const from = processLink(link.from ?? [], 'input');
  const to = processLink(link.to ?? [], 'output');
  const base = processLink(link.base ?? [], 'base');
  const actions = processLink(link.actions ?? [], 'actions');
  return {...link, from, to, base, actions};
}

function processLink(io: LinkSpecString, ioType: 'input' | 'output' | 'base' | 'actions') {
  if (Array.isArray(io))
    return io.map((item) => parseLinkIO(item, ioType));
  else if (io)
    return [parseLinkIO(io, ioType)];
  else
    return [];
}
