import {PipelineConfiguration, PipelineStepConfiguration, PipelinePopupConfiguration, PipelineActionConfiguraion} from '../PipelineConfiguration';

export class Aborted extends Error { }

export type SubNodeConfTypes = 'action' | 'popup' | 'step';
export type SubNodeConf = PipelineActionConfiguraion | PipelinePopupConfiguration | PipelineStepConfiguration;
export type NodeConfTypes = SubNodeConfTypes | 'pipeline';
export type NodeConf = SubNodeConf | PipelineConfiguration;
