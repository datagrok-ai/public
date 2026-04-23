import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Observable, OperatorFunction, defer, from, merge, of} from 'rxjs';
import {concatMap, distinctUntilChanged, filter, map, reduce, share, withLatestFrom, windowToggle} from 'rxjs/operators';
import dayjs from 'dayjs';
import {deepEqual, createCustomEqual, TypeEqualityComparator} from 'fast-equals';
import {HandlerBase} from './config/PipelineConfiguration';
import {ValidationResult} from './data/common-types';
import {NodeAddressSegment, NodePathSegment, TreeNode} from './data/BaseTree';
import {StateTreeNode} from './runtime/StateTreeNodes';

/**
 * Buffers [key, value] pairs during a lock period, deduplicating by key (latest wins).
 * On unlock, all buffered entries are emitted. While unlocked, values pass through.
 */
export function bufferKeysDuringLock<K, V>(
  lock$: Observable<boolean>,
): OperatorFunction<readonly [K, V], readonly [K, V]> {
  return (source$) => {
    const shared$ = source$.pipe(share());
    const lockDistinct$ = lock$.pipe(distinctUntilChanged());
    const lockOn$ = lockDistinct$.pipe(filter((l) => l));
    const lockOff$ = lockDistinct$.pipe(filter((l) => !l));

    return merge(
      shared$.pipe(
        withLatestFrom(lockDistinct$),
        filter(([, locked]) => !locked),
        map(([item]) => item),
      ),
      shared$.pipe(
        windowToggle(lockOn$, () => lockOff$),
        concatMap((window$) => window$.pipe(
          reduce((acc, [k, v]) => acc.set(k, v), new Map<K, V>()),
          concatMap((buffered) => from(
            [...buffered].map(([k, v]) => [k, v] as const),
          )),
        )),
      ),
    );
  };
}

export function callHandler<R, P = any>(handler: HandlerBase<P, R>, params: P): Observable<R> {
  if (typeof handler === 'string') {
    return defer(async () => {
      const f = DG.Func.byName(handler);
      const call = f.prepare({params});
      await call.call(undefined, undefined, {processed: true, report: false});
      const res = call.getOutputParamValue() as R;
      return res;
    });
  } else {
    return defer(() => {
      const res = handler(params);
      if (res instanceof Observable || res instanceof Promise)
        return res;
      else
        return of(res);
    });
  }
}

export function mergeValidationResults(...results: (ValidationResult | undefined)[]): ValidationResult {
  const errors = [];
  const warnings = [];
  const notifications = [];

  for (const result of results) {
    if (result) {
      errors.push(...(result.errors ?? []));
      warnings.push(...(result.warnings ?? []));
      notifications.push(...(result.notifications ?? []));
    }
  }
  return {errors, warnings, notifications};
}

export function pathToUUID(
  rnode: TreeNode<StateTreeNode>,
  path: readonly NodeAddressSegment[] | readonly NodePathSegment[],
): string[] {
  let acc = {node: rnode, uuids: [rnode.getItem().uuid]};

  path.forEach(({idx}) => {
    const nnode = acc.node.getChild({idx});
    const {uuid} = nnode.getItem();

    acc = {node: nnode, uuids: [...acc.uuids, uuid]};
  });
  return acc.uuids;
}

export function indexFromEnd<T>(arr: Readonly<T[]>, offset = 0): T | undefined {
  return arr[arr.length - offset - 1];
}

const areObjectsEqual: TypeEqualityComparator<Record<any, any>> = (a, b) => {
  if (a instanceof DG.DataFrame && b instanceof DG.DataFrame) {
    if (a.rowCount !== b.rowCount || [...a.columns].length !== [...b.columns].length)
      return false;
    for (const columnA of a.columns) {
      const columnB = b.columns.byName(columnA.name);

      if (columnA.type !== columnB.type || columnA.name !== columnB.name)
        return false;

      for (let i = 0; i < a.rowCount; i++) {
        const valueA = columnA.get(i);
        const valueB = columnB.get(i);
        if (!customDeepEqual(valueA, valueB))
          return false;
      }
    }
    return true;
  }

  if (dayjs.isDayjs(a) && dayjs.isDayjs(b))
    return a.isSame(b);


  if (!deepEqual(a, b))
    return false;

  return true;
};

const FLOAT_TOLERANCE = 0.0001;

const areNumbersEqual: TypeEqualityComparator<number> = (a, b) => {
  if (isNaN(a) && isNaN(b))
    return true;
  return Math.abs(a - b) < FLOAT_TOLERANCE;
};

const customFastEqual = createCustomEqual({
  createCustomConfig: () => ({areNumbersEqual, areObjectsEqual}),
});

// TODO: try using lib apis
export const customDeepEqual = (a: any, b: any) => (a == null && b == null) || customFastEqual(a, b);

/** Format a single path segment: `id[0]`, `id[2]`. */
export function formatPathSegment(seg: {idx: number; id: string}): string {
  return `${seg.id}[${seg.idx}]`;
}

/** Format a node path as `step[0]/child[1]::ioName`.
 *  Segments joined by `/`, IO name appended after `::`. */
export function formatNodePath(path: readonly {idx: number; id: string}[], ioName?: string): string {
  const str = path.map(formatPathSegment).join('/');
  return ioName ? `${str}::${ioName}` : str;
}

/** Format a mutation path: `parent/+id[2]` (add), `parent/-id[2]` (remove), `parent/id[2→3]` (move). */
export function formatMutationPath(
  path: readonly {idx: number; id: string}[] | undefined,
  addIdx?: number, removeIdx?: number, id?: string,
): string {
  let str = path?.length ? path.map(formatPathSegment).join('/') : '';
  const append = (suffix: string) => str += str ? `/${suffix}` : suffix;
  if (addIdx != null && removeIdx != null)
    append(`${id}[${removeIdx}\u2192${addIdx}]`);
  else if (removeIdx != null)
    append(`-${id}[${removeIdx}]`);
  else if (addIdx != null)
    append(`+${id}[${addIdx}]`);
  else if (id)
    append(id);
  return str;
}
