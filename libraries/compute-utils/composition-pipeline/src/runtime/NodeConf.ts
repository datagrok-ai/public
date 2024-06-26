import {FixedPipelineConfiguration, PipelineStepConfiguration, PipelineActionConfiguraion} from '../config/PipelineConfiguration';

export class Aborted extends Error { }

export type SubNodeConfTypes = 'action' | 'step';
export type SubNodeConf = PipelineActionConfiguraion | PipelineStepConfiguration;
export type NodeConfTypes = SubNodeConfTypes | 'pipeline';
export type NodeConf = SubNodeConf | FixedPipelineConfiguration;
