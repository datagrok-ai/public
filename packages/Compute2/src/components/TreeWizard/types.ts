import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {Stat} from '@he-tree/vue/dist/v3/components/TreeProcessorVue';

type RestrictedStat = Stat<PipelineState>

export type AugmentedStat = RestrictedStat & {
  parent: AugmentedStat | null;
  children: AugmentedStat[];
};


export type Status = 'next' | 'next warn' | 'next error' | 'pending' | 'pending executed' | 'running' | 'succeeded' | 'failed' | 'succeeded warn' | 'succeeded inconsistent';
