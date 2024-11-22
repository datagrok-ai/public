import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {deserialize, serialize} from '@datagrok-libraries/utils/src/json-serialization';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {TestScheduler} from 'rxjs/testing';
import {Observable, of} from 'rxjs';
import {StateTree} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTree';
import {map, mapTo, startWith, switchMap} from 'rxjs/operators';
import {isFuncCallNode, StateTreeNode} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {NodePath} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/BaseTree';
import {PipelineConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';

declare global {
  var SNAPSHOTS_UPDATE_MODE: boolean;
}

const snapshotsPath = 'System:AppData/LibTests/snapshots/';

export async function snapshotCompare(actual: any, snapshotName: string) {
  if (globalThis.SNAPSHOTS_UPDATE_MODE) {
    const data = serialize(actual, {useJsonDF: true, space: 2});
    const name = snapshotName + '.json';
    const blob = new Blob([data]);
    DG.Utils.download(name, blob);
  } else {
    const data = await grok.dapi.files.readAsText(snapshotsPath + snapshotName + '.json');
    const expected = deserialize(data);
    expectDeepEqual(actual, expected);
  }
}

export type ExpectObservableNamed = (name: string, observable: Observable<any>, subscriptionMarbles?: string) => void;

export async function runRXTreeSnapshotTest(testName: string, fn: (expectObservable: ExpectObservableNamed, cold: typeof TestScheduler.prototype.createColdObservable) => void) {
  let current: {data: any, expectName?: string}[] = [];
  let snapshotName = testName;

  let testScheduler = new TestScheduler((data, fakeExp) => {
    const expectName = fakeExp[0].notification.value;
    current.push({expectName, data});
  });

  testScheduler.run((helpers) => {
    function expectObservable(expectName: string, observable: Observable<any>, subscriptionMarbles?: string) {
      return helpers.expectObservable(observable, subscriptionMarbles).toBe('a', {a: expectName});
    }
    fn(expectObservable, helpers.cold);
  });

  if (globalThis.SNAPSHOTS_UPDATE_MODE) {
    const data = serialize(current, {useJsonDF: true, space: 2});
    const name = snapshotName + '.json';
    const blob = new Blob([data]);
    DG.Utils.download(name, blob);
  } else {
    const data = await grok.dapi.files.readAsText(snapshotsPath + snapshotName + '.json');
    const allExpected = deserialize(data) as {data: any, expectName?: string}[];
    for (let idx = 0; idx < allExpected.length; idx++) {
      const expected = allExpected[idx];
      const actual = current[idx];
      expectDeepEqual(actual.expectName, expected.expectName, { prefix: `${idx} test name`});
      expectDeepEqual(actual.data, expected.data, {prefix: `${idx} ${expected.expectName}`});
    }
  }
}

export function getTreeStates(config: PipelineConfigurationProcessed): [StateTree, Observable<StateTree>] {
  const initialTree = StateTree.fromPipelineConfig({config, mockMode: true});
  initialTree.init().subscribe();
  return [
    initialTree,
    initialTree.makeStateRequests$.pipe(
      startWith(initialTree),
      mapTo(initialTree),
    )
  ] as const;
}

function expectTreeStates(tree$: Observable<StateTree>, expectObservable: ExpectObservableNamed, prefix: string, lister: (node: StateTreeNode) => string[], getter: (node: StateTreeNode, name: string) => Observable<any>) {
  const addedUUIDs = new Set<string>();
  tree$.pipe(
    map((tree, idx) => {
      const items = tree.nodeTree.traverse(tree.nodeTree.root, (acc, node, path) => {
        const item = node.getItem();
        if (addedUUIDs.has(item.uuid))
          return acc;
        addedUUIDs.add(item.uuid);
        return [...acc, {item, path}];
      }, [] as {item: StateTreeNode, path: NodePath}[]);
      for (const item of items) {
        for (const name of lister(item.item))
          expectObservable(`type: ${prefix}, gen: ${idx}, path: ${JSON.stringify(item.path)}, state: ${name}`, getter(item.item, name));
      }
    }),
  ).subscribe();
}

export function expectTreeData(tree$: Observable<StateTree>, expectObservable: ExpectObservableNamed) {
  expectTreeStates(
    tree$,
    expectObservable,
    'data',
    (node) => node.getStateStore().getStateNames(),
    (node, name) => node.getStateStore().getStateChanges(name)
  );
}

export function expectTreeMeta(tree$: Observable<StateTree>, expectObservable: ExpectObservableNamed) {
  expectTreeStates(
    tree$,
    expectObservable,
    'meta',
    (node) => node.getStateStore().getStateNames(),
    (node, name) => {
      if (isFuncCallNode(node)) {
        return node.metaInfo$.pipe(
          switchMap((meta) => {
            return meta[name] ?? of(undefined);
          })
        );
      }
      return of(undefined);
    }
  );
}

export function expectTreeValidations(tree$: Observable<StateTree>, expectObservable: ExpectObservableNamed) {
  expectTreeStates(
    tree$,
    expectObservable,
    'validations',
    (node) => node.getStateStore().getStateNames(),
    (node, name) => {
      if (isFuncCallNode(node)) {
        return node.validationInfo$.pipe(
          map((val) => {
            return val[name];
          })
        );
      }
      return of(undefined);
    }
  );
}


export function expectTreeConsistency(tree$: Observable<StateTree>, expectObservable: ExpectObservableNamed) {
  expectTreeStates(
    tree$,
    expectObservable,
    'consistency',
    (node) => node.getStateStore().getStateNames(),
    (node, name) => {
      if (isFuncCallNode(node)) {
        return node.consistencyInfo$.pipe(
          map((val) => {
            return val[name];
          })
        );
      }
      return of(undefined);
    }
  );
}
