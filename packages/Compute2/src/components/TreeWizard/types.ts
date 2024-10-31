import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {Stat} from '@he-tree/vue/dist/v3/components/TreeProcessorVue';


export type RuntimeData = {}

export type Data = RuntimeData & PipelineState

type NoStringIndex<T> = { [K in keyof T as string extends K ? never : K]: T[K] };

type RestrictedStat = NoStringIndex<Stat<Data>>

export type AugmentedStat = RestrictedStat & {
  parent: AugmentedStat | null;
  children: AugmentedStat[];
};


export type Status = 'locked' | 'next' | 'pending' | 'running' | 'succeeded' | 'failed' | 'partially succeeded';
