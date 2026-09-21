import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {bySortKey, composeMenuGroups, hoverText, panelsEqual, reconcileByKey, resolveClick}
  from '@datagrok-libraries/webcomponents-vue/src/ViewService/ribbon-core';
import {RibbonItemBase} from '@datagrok-libraries/webcomponents-vue/src/ViewService/ribbon-types';
import {menuItem} from './ribbon-test-utils';

interface FakeState {
  key: string;
  item: {key: string};
  created: boolean;
}

function runReconcile(prev: FakeState[], nextKeys: string[]) {
  return reconcileByKey<FakeState, {key: string}>(
    prev, nextKeys.map((key) => ({key})),
    (item) => item.key, (state) => state.key,
    (item, key) => ({key, item, created: true}),
    (state, item) => {
      state.item = item;
    });
}

const states = (keys: string[]): FakeState[] =>
  keys.map((key) => ({key, item: {key}, created: false}));

category('WebComponentsVue: Ribbon core', () => {
  test('sorts by ascending priority with unset last', async () => {
    const input = [
      {priority: undefined, seq: 0},
      {priority: 10, seq: 1},
      {priority: 0, seq: 2},
    ];
    const sorted = [...input].sort(bySortKey);
    expectDeepEqual(sorted.map((e) => e.seq), [2, 1, 0]);
  });

  test('breaks priority ties and unset by registration order', async () => {
    const input = [
      {priority: undefined, seq: 3},
      {priority: 5, seq: 2},
      {priority: undefined, seq: 1},
      {priority: 5, seq: 0},
    ];
    const sorted = [...input].sort(bySortKey);
    expectDeepEqual(sorted.map((e) => e.seq), [0, 2, 1, 3]);
  });

  test('merges same name groups in registration order', async () => {
    const groups = composeMenuGroups([
      {seq: 0, name: 'G', items: [menuItem({text: 'a'})]},
      {seq: 1, name: 'H', items: [menuItem({text: 'h'})]},
      {seq: 2, name: 'G', items: [menuItem({text: 'b'})]},
    ]);
    expectDeepEqual(groups.map((g) => g.name), ['G', 'H']);
    expectDeepEqual(groups[0].items.map((i) => i.text), ['a', 'b']);
  });

  test('ranks a merged group by its first seen contributor', async () => {
    const groups = composeMenuGroups([
      {seq: 0, name: 'G', items: [menuItem({text: 'a'})]},
      {seq: 1, priority: 0, name: 'H', items: [menuItem({text: 'h'})]},
      {seq: 2, priority: -5, name: 'G', items: [menuItem({text: 'b'})]},
    ]);
    expectDeepEqual(groups.map((g) => g.name), ['H', 'G']);
    expectDeepEqual(groups[1].items.map((i) => i.text), ['a', 'b']);
  });

  test('drops empty groups', async () => {
    const groups = composeMenuGroups([
      {seq: 0, name: 'Empty', items: []},
      {seq: 1, name: 'Full', items: [menuItem()]},
    ]);
    expectDeepEqual(groups.map((g) => g.name), ['Full']);
  });

  test('keeps lower states when an upper item is removed', async () => {
    const prev = states(['a', 'b', 'c']);
    const {states: next, changed} = runReconcile(prev, ['b', 'c']);
    expect(changed);
    expect(next[0] === prev[1]);
    expect(next[1] === prev[2]);
  });

  test('creates only the inserted item', async () => {
    const prev = states(['b', 'c']);
    const {states: next, changed} = runReconcile(prev, ['a', 'b', 'c']);
    expect(changed);
    expectDeepEqual(next.map((s) => s.created), [true, false, false]);
    expect(next[1] === prev[0]);
  });

  test('reorders without recreating', async () => {
    const prev = states(['a', 'b']);
    const {states: next, changed} = runReconcile(prev, ['b', 'a']);
    expect(changed);
    expect(next[0] === prev[1]);
    expect(next[1] === prev[0]);
  });

  test('reports no change for an identical list', async () => {
    const prev = states(['a', 'b']);
    const {states: next, changed} = runReconcile(prev, ['a', 'b']);
    expect(changed, false);
    expect(next[0] === prev[0] && next[1] === prev[1]);
  });

  test('pools duplicate keys in order', async () => {
    const prev = states(['a', 'a']);
    const {states: next, changed} = runReconcile(prev, ['a', 'a', 'a']);
    expect(changed);
    expect(next[0] === prev[0] && next[1] === prev[1]);
    expect(next[2].created);
  });

  test('resolves clicks by disabled state and reason mode', async () => {
    const enabled: RibbonItemBase = {onClick: () => {}};
    expectDeepEqual(resolveClick(enabled), {kind: 'run'});
    const silent: RibbonItemBase = {onClick: () => {}, disabled: true, disabledReason: 'busy'};
    expectDeepEqual(resolveClick(silent), {kind: 'ignore'});
    const popup: RibbonItemBase =
      {onClick: () => {}, disabled: true, disabledReason: 'busy', disabledReasonMode: 'popup'};
    expectDeepEqual(resolveClick(popup), {kind: 'warn', reason: 'busy'});
    const both: RibbonItemBase =
      {onClick: () => {}, disabled: true, disabledReason: () => 'thunked', disabledReasonMode: 'both'};
    expectDeepEqual(resolveClick(both), {kind: 'warn', reason: 'thunked'});
    const noReason: RibbonItemBase = {onClick: () => {}, disabled: true, disabledReasonMode: 'popup'};
    expectDeepEqual(resolveClick(noReason), {kind: 'ignore'});
  });

  test('resolves hover text by disabled state and reason mode', async () => {
    expectDeepEqual(hoverText({onClick: () => {}, tooltip: 'tip'}), 'tip');
    expectDeepEqual(hoverText(
      {onClick: () => {}, tooltip: 'tip', disabled: true, disabledReason: 'busy'}), 'busy');
    expectDeepEqual(hoverText(
      {onClick: () => {}, tooltip: 'tip', disabled: true}), 'tip');
    expectDeepEqual(hoverText(
      {onClick: () => {}, tooltip: 'tip', disabled: true, disabledReason: 'busy', disabledReasonMode: 'popup'}),
    'tip');
    expectDeepEqual(hoverText({onClick: () => {}}), null);
  });

  test('compares panel element matrices', async () => {
    const a = document.createElement('div');
    const b = document.createElement('div');
    expect(panelsEqual([[a], [b]], [[a], [b]]));
    expect(panelsEqual([[a]], [[b]]), false);
    expect(panelsEqual([[a]], [[a], [b]]), false);
    expect(panelsEqual([[a, b]], [[a]]), false);
  });
});
