import {TraverseHandler} from './common-types';

export function buildTraverseB<T, P, S>(startPath: P, getNext: (item: T, path: P, state?: S) => (readonly [T, P, S?])[], startState?: S) {
  return function traverse<A>(start: T, handler: TraverseHandler<A, T, P>, acc: A) {
    let stop = false;
    const signal = () => stop = true;
    const q: (readonly [T, P, S?])[] = [[start, startPath, startState]] as const;
    while (q.length) {
      const [item, path, state] = q.shift()!;
      acc = handler(acc, item, path, signal);
      if (stop)
        return acc;
      const next = getNext(item, path, state);
      if (next?.length)
        q.push(...next);
    }
    return acc;
  };
}

export function buildTraverseD<T, P, S>(startPath: P, getNext: (item: T, path: P, state?: S) => (readonly [T, P, S?])[], startState?: S) {
  return function traverse<A>(start: T, handler: TraverseHandler<A, T, P>, acc: A) {
    let stop = false;
    const signal = () => stop = true;
    const stk: (readonly [T, P, S?])[] = [[start, startPath, startState]] as const;
    while (stk.length) {
      const [item, path, state] = stk.pop()!;
      acc = handler(acc, item, path, signal);
      if (stop)
        return acc;
      const next = getNext(item, path, state);
      if (next?.length) {
        const rnext = [...next].reverse();
        stk.push(...rnext);
      }
    }
    return acc;
  };
}
