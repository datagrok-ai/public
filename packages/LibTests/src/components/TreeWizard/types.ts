import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Stat} from '@he-tree/vue/types/src/components/TreeProcessorVue';
import {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';


export type RuntimeData = {
  isHovered: boolean,
}

export type Data = RuntimeData & PipelineState

type NoStringIndex<T> = { [K in keyof T as string extends K ? never : K]: T[K] };

type RestrictedStat = NoStringIndex<Stat<Data>>

export type AugmentedStat = RestrictedStat & {
  parent: AugmentedStat | null;
  children: AugmentedStat[];
};


export type Status = 'locked' | `didn't run` | 'running' | 'succeeded' | 'failed' | 'partially succeeded';

export type HueTree = {
  add: Function,
  remove: Function,
  isDraggable: Function,
  isDroppable: Function,
  statsFlat: AugmentedStat[],
  $el: HTMLElement,
};
