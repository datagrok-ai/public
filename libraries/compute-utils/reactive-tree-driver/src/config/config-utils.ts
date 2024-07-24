import {ItemPathArray} from '../data/common-types';
import {pathJoin} from '../utils';

export type StateTraverseNode<C> = {
  id: string;
  steps: StateTraverseNode<C>[];
} & C;

export type StateTraverseItem<C> = {
  path: ItemPathArray;
} & StateTraverseNode<C>;

export function traverseConfig<C, T>(
  start: StateTraverseNode<C>,
  handler: (
    acc: T,
    conf: StateTraverseItem<C>,
    path: ItemPathArray,
    stop: () => void,
  ) => T,
  acc: T,
  path: ItemPathArray = [],
): T {
  const q: StateTraverseItem<C>[] = [{...start, path}];
  let stop = false;
  const signal = () => stop = true;
  while (q.length) {
    const item = q.shift()!;
    acc = handler(acc, item, item.path, signal);
    if (stop)
      return acc;
    const items = getNext(item, path);
    q.push(...items);
  }
  return acc;
}

function getNext<C>(config: StateTraverseNode<C>, path: ItemPathArray): StateTraverseItem<C>[] {
  const steps = config.steps ?? [];
  return steps.map((step) => ({...step, path: pathJoin(path, [step.id])}));
}
