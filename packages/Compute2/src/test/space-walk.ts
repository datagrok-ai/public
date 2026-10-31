import * as DG from 'datagrok-api/dg';
import {defer, of, throwError} from 'rxjs';
import {TestScheduler} from 'rxjs/testing';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {ListSpace, SpaceListing, makeIndex, reloadSpaces, walkSpaces} from '../components/RunStorage/space-index';

const space = (id: string) => ({id, name: id, friendlyName: id, isSpace: true}) as unknown as DG.Project;
const listingOf = (tree: Record<string, string[]>, id: string): SpaceListing =>
  ({subspaces: (tree[id] ?? []).map(space), runs: []});
const json = (value: unknown) => JSON.stringify(value);

/** Walks the tree from `r` in virtual time, each listing answering after its marble, and tells when each listing
 * started and each node came, as `frame:id`. */
function walkTree(tree: Record<string, string[]>,
  {marbles = {}, unsubscribeAt}: {marbles?: Record<string, string>, unsubscribeAt?: number} = {}) {
  const scheduler = new TestScheduler((actual, expected) => expect(json(actual), json(expected)));
  const started: string[] = [];
  const nodes: string[] = [];
  let completedAt: number | undefined;
  scheduler.run(({cold}) => {
    const list: ListSpace = (listed) => defer(() => {
      started.push(`${scheduler.now()}:${listed.id}`);
      return cold(marbles[listed.id] ?? '--(a|)', {a: listingOf(tree, listed.id)}, new Error(`${listed.id} is gone`));
    });
    const walk = walkSpaces([space('r')], list).subscribe({
      next: (node) => nodes.push(`${scheduler.now()}:${node.space.id}${node.error ? ` (${node.error})` : ''}`),
      complete: () => completedAt = scheduler.now(),
    });
    if (unsubscribeAt !== undefined)
      scheduler.schedule(() => walk.unsubscribe(), unsubscribeAt);
  });
  return {started, nodes, completedAt};
}

const sixSubspaces = {r: ['c1', 'c2', 'c3', 'c4', 'c5', 'c6'], c1: [], c2: [], c3: [], c4: [], c5: [], c6: []};

category('Run storage: space walk', () => {
  test('The walk keeps four listings going at every level', async () => {
    const {started, nodes, completedAt} = walkTree(sixSubspaces);
    expect(json(started), json(['0:r', '2:c1', '2:c2', '2:c3', '2:c4', '4:c5', '4:c6']));
    expect(json(nodes), json(['2:r', '4:c1', '4:c2', '4:c3', '4:c4', '6:c5', '6:c6']));
    expect(completedAt, 6);
  });

  test('A space under several parents is listed once', async () => {
    const {started, nodes} = walkTree({r: ['a', 'b'], a: ['shared'], b: ['shared'], shared: []});
    expect(json(started), json(['0:r', '2:a', '2:b', '4:shared']));
    expect(nodes.length, 4);
  });

  test('A space that can\'t be listed keeps its error, and the others load', async () => {
    const {started, nodes, completedAt} = walkTree({r: ['a', 'b'], a: ['below'], b: [], below: []},
      {marbles: {a: '---#'}});
    expect(json(started), json(['0:r', '2:a', '2:b']));
    expect(json(nodes), json(['2:r', '4:b', '5:a (a is gone)']));
    expect(completedAt, 5);
  });

  test('Unsubscribing stops the walk', async () => {
    const {started, nodes, completedAt} = walkTree(sixSubspaces, {unsubscribeAt: 3});
    expect(json(started), json(['0:r', '2:c1', '2:c2', '2:c3', '2:c4']));
    expect(json(nodes), json(['2:r']));
    expect(completedAt === undefined, true);
  });

  test('A reload lists the spaces again with their new subspaces, and a failing one keeps what it had', async () => {
    const node = (id: string, children: string[]) => ({space: space(id), children, runs: []});
    const index = makeIndex(['r', 'gone'], new Map([['r', node('r', ['a'])], ['a', node('a', [])],
      ['gone', node('gone', ['g1'])], ['g1', node('g1', [])]]));
    const now: Record<string, string[]> = {r: ['a', 'n'], n: []};
    const listed: string[] = [];
    const list: ListSpace = (relisted) => defer(() => {
      listed.push(relisted.id);
      return relisted.id === 'gone' ? throwError(new Error('gone is gone')) : of(listingOf(now, relisted.id));
    });
    let reloaded: ReturnType<typeof makeIndex> | undefined;
    reloadSpaces(index, ['r', 'gone'], list).subscribe((next) => reloaded = next);

    expect(json(listed.sort()), json(['gone', 'n', 'r']));
    expect(reloaded!.nodes.get('a') === index.nodes.get('a'), true);
    expect(json(reloaded!.parents.get('n')), json(['r']));
    const gone = reloaded!.nodes.get('gone')!;
    expect(json([gone.children, gone.error]), json([['g1'], 'gone is gone']));
  });
});
