import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type TreeNodeType = {
  text: string,
  children: TreeNodeType[]
}

export type Data = {
  funcCall?: string,
  text: string,
  order?: Order,
}


type Order = 'sequental' | 'parallel';

export type AugmentedStat = {
  isHovered: boolean,
  status: Status,
  data: Data;
  open: boolean;
  parent: AugmentedStat | null;
  children: AugmentedStat[];
  level: number;
};


export type Status = 'locked' | `didn't run` | 'running' | 'succeeded' | 'failed' | 'partially succeeded';

export type HueTree = {
  add: Function,
  remove: Function,
  isDraggable: Function,
  isDroppable: Function,
};
