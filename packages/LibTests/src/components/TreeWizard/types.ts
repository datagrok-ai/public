import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Stat} from '@he-tree/vue/types/src/components/TreeProcessorVue';

export type StepConfig = {
  text: string,
  children: StepConfig[]
}

type NoStringIndex<T> = { [K in keyof T as string extends K ? never : K]: T[K] };

type RestrictedStat = NoStringIndex<Stat<Data>>

export type Data = {
  status: Status,
  isHovered: boolean,
  funcCall?: string,
  text: string,
  order?: Order,
}

type Order = 'sequental' | 'parallel';

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
};
