import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {
  STEP_HISTORY_OPTION, listMyModelRuns, loadRunAuthors, notifyRunsChanged, placeRunFilter,
} from '../components/RunStorage/run-storage';
import {useSpaceIndex} from '../components/RunStorage/use-space-index';
import {
  SpaceIndex, loadSpaceIndex, makeIndex, parentsOf, reloadSpaces, subtreeCounts, subtreeRuns,
} from '../components/RunStorage/space-index';
import {SPACES, locationKey, useSpaceNavigation} from '../components/RunStorage/use-space-navigation';
import {NO_FILTERS, RunFilters, useRunFilters} from '../components/RunStorage/use-run-filters';
import {RunFiltersBar} from '../components/RunStorage/RunFiltersBar';
import {useRunDetails} from '../components/RunStorage/use-run-details';
import {useRunSelection} from '../components/RunStorage/use-run-selection';
import {useMyRuns} from '../components/RunStorage/use-my-runs';
import {useRunAuthors} from '../components/RunStorage/use-run-authors';
import {RunList, spaceFolders} from '../components/RunStorage/RunList';
import {ModelSelector} from '../components/RunStorage/ModelSelector';
import {PlacesTree} from '../components/RunStorage/PlacesTree';
import {SaveRunDialog} from '../components/RunStorage/SaveRunDialog';
import {RunBrowser} from '../components/RunStorage/RunBrowser';
import {DAY_CALCULATION, uniqueName, withRunStorage} from './run-storage';

const dayFilter = () => placeRunFilter(DG.Func.byName(DAY_CALCULATION));
const runIds = (index: SpaceIndex, spaceIds: string[]) =>
  subtreeRuns(index, spaceIds, dayFilter()).map(({run}) => run.id).sort();

// a plain index for the parts that don't need the server
const space = (id: string) => ({id, name: id, friendlyName: id, isSpace: true}) as unknown as DG.Project;
const plainIndex = (tree: Record<string, string[]>, top: string[]): SpaceIndex =>
  makeIndex(top, new Map(Object.entries(tree).map(([id, children]) => [id, {space: space(id), children, runs: []}])));

// the id of a run that was never saved
const missingId = () => {
  const call = DG.Func.byName(DAY_CALCULATION).prepare({});
  call.newId();
  return call.id;
};

function mount(component: Vue.Component, props: Record<string, any>) {
  const root = ui.div();
  document.body.append(root);
  const app = Vue.createApp(component, props);
  app.mount(root);
  return {root, unmount: () => {
    app.unmount();
    root.remove();
  }};
}

category('Run storage: space index', () => {
  test('The index lists linked spaces once, counts each run once by the filter, and reloads one space', () =>
    withRunStorage(async (f) => {
      // a > b, d links b, and the run in a is linked into b too
      const a = await f.space();
      const b = await f.subspace(a);
      const d = await f.space();
      await grok.dapi.spaces.id(d.id).addSubspace(b, true);
      const inA = await f.saveDayRun(uniqueName('In a'), a.id);
      const inB = await f.saveDayRun(uniqueName('In b'), b.id);
      await grok.dapi.spaces.id(b.id).addEntity(inA.id, true);
      const step = await f.dayRun();
      step.options[STEP_HISTORY_OPTION] = 'true';
      await f.save(step, {metadata: {name: uniqueName('Step')}, target: {spaceId: a.id}});

      const index = await loadSpaceIndex().toPromise();
      expect(index.top.includes(a.id) && index.top.includes(d.id), true);
      expect(JSON.stringify(index.nodes.get(a.id)!.children), JSON.stringify([b.id]));
      expect(index.nodes.get(d.id)!.children.includes(b.id), true);
      expect(JSON.stringify(parentsOf(index, b.id).map((parent) => parent.id).sort()),
        JSON.stringify([a.id, d.id].sort()));
      const both = JSON.stringify([inA.id, inB.id].sort());
      expect(JSON.stringify(runIds(index, [d.id])), both);
      const day = DG.Func.byName(DAY_CALCULATION);
      const counts = subtreeCounts(index, placeRunFilter(day));
      expect(JSON.stringify([counts.get(a.id), counts.get(b.id), counts.get(d.id)]), JSON.stringify([3, 2, 2]));
      expect(subtreeCounts(index, placeRunFilter(day, true)).get(a.id), 1);

      const added = await f.saveDayRun(uniqueName('Added'), b.id);
      const reloaded = await reloadSpaces(index, [b.id]).toPromise();
      expect(reloaded.nodes.get(b.id)!.runs.some((run) => run.id === added.id), true);
      expect(reloaded.nodes.get(a.id) === index.nodes.get(a.id), true);
    }));

  test('Two quick changes to different spaces both reach the index', () =>
    withRunStorage(async (f) => {
      const a = await f.space();
      const b = await f.space();
      let spaces!: ReturnType<typeof useSpaceIndex>;
      const {unmount} = mount({setup: () => {
        spaces = useSpaceIndex();
        return () => null;
      }}, {});
      try {
        await awaitCheck(() => !!spaces.index.value, 'the index should load', 15000);
        const inA = await f.saveDayRun(uniqueName('In a'), a.id);
        const inB = await f.saveDayRun(uniqueName('In b'), b.id);
        notifyRunsChanged({spaceIds: [a.id]});
        notifyRunsChanged({spaceIds: [b.id]});
        const listed = (spaceId: string, runId: string) =>
          !!spaces.index.value?.nodes.get(spaceId)?.runs.some((run) => run.id === runId);
        await awaitCheck(() => listed(b.id, inB.id), 'the second change should reach the index', 15000);
        expect(listed(a.id, inA.id), true);
      } finally {
        unmount();
      }
    }));

  test('Navigation goes into spaces, back up the way it came and to a breadcrumb', async () => {
    const index = Vue.shallowRef(plainIndex({a: ['b'], b: ['c'], c: [], d: ['b']}, ['a', 'd']));
    const scope = Vue.effectScope();
    const nav = scope.run(() => useSpaceNavigation(index, SPACES))!;
    try {
      expect(JSON.stringify(nav.folders.value), JSON.stringify(['a', 'd']));
      nav.enter('d');
      nav.enter('b');
      expect(nav.crumbs.value.map((crumb) => crumb.label).join('>'), 'Spaces>d>b');
      nav.up();
      expect(nav.spaceId.value, 'd');
      nav.enter('b');
      nav.enter('c');
      nav.goTo(1);
      expect(locationKey(nav.location.value), 'spaces:d');
      nav.enter('b');
      // b stops being linked into d
      index.value = plainIndex({a: ['b'], b: ['c'], c: [], d: []}, ['a', 'd']);
      await Vue.nextTick();
      expect(locationKey(nav.location.value), 'spaces:d');
    } finally {
      scope.stop();
    }
  });

  test('Filters combine search, tags and authors, and count their options', async () => {
    const make = (name: string, tags: string[]) => {
      const call = DG.Func.byName(DAY_CALCULATION).prepare({});
      call.newId();
      call.options['title'] = name;
      call.options['tags'] = tags;
      return call;
    };
    const runs = [make('Early feed', ['feed']), make('Late feed', ['feed', 'late']), make('Baseline', [])];
    const authors = Vue.shallowRef(new Map([
      [runs[0].id, {id: 'u1', name: 'Maria'}], [runs[1].id, {id: 'u2', name: 'Tom'}]]));
    const scope = Vue.effectScope();
    const {filters, isActive, tagOptions, matches} = scope.run(() => useRunFilters(Vue.shallowRef(runs), authors))!;
    try {
      expect(isActive.value, false);
      expect(tagOptions.value.map((option) => `${option.value}:${option.detail}`).join(), 'feed:2,late:1');
      filters.value = {...NO_FILTERS, search: 'feed'};
      expect(runs.filter(matches).length, 2);
      filters.value = {...filters.value, tags: ['late']};
      expect(runs.filter(matches).map(getRunName).join(), 'Late feed');
      filters.value = {...filters.value, tags: [], authors: ['u1']};
      expect(runs.filter(matches).map(getRunName).join(), 'Early feed');
      expect(isActive.value, true);
      filters.value = NO_FILTERS;
      expect(isActive.value, false);
    } finally {
      scope.stop();
    }
  });

  test('A click selects one run, a ctrl-click adds or removes it, and a double click opens unless it toggled',
    async () => {
      const [a, b, c] = ['a', 'b', 'c'].map((name) => {
        const call = DG.Func.byName(DAY_CALCULATION).prepare({});
        call.newId();
        call.options['title'] = name;
        return call;
      });
      const listed = Vue.shallowRef([a, b, c]);
      const scope = Vue.effectScope();
      const selection = scope.run(() => useRunSelection(() => listed.value))!;
      const selected = () => selection.selection.value.map(getRunName).join();
      const plain = (detail: number) => ({ctrlKey: false, metaKey: false, detail});
      const ctrl = (detail: number) => ({ctrlKey: true, metaKey: false, detail});
      try {
        expect(selection.click(a, plain(1)), false);
        expect(selection.click(b, ctrl(1)), false);
        expect(selected(), 'a,b');
        expect(selection.actionRuns(a).map(getRunName).join(), 'a,b');
        expect(selection.actionRuns(c).map(getRunName).join(), 'c');
        selection.click(a, ctrl(1));
        expect(selected(), 'b');

        expect(selection.click(c, plain(1)) || selection.click(c, plain(2)), true);
        expect(selected(), 'c');
        expect(selection.click(a, ctrl(1)) || selection.click(a, plain(2)), false);
        expect(selection.click(b, plain(1)) || selection.click(b, ctrl(2)), false);

        selection.click(a, plain(1));
        selection.click(c, ctrl(1));
        listed.value = [a, b];
        await Vue.nextTick();
        expect(selected(), 'a');

        selection.keepOnly([b]);
        expect(selected(), 'b');
        selection.clear();
        expect(selected(), '');
      } finally {
        scope.stop();
      }
    });

  test('The filters bar reports each change as a new value and keeps the one it got', async () => {
    const value: RunFilters = {search: '', tags: ['feed', 'late'], authors: ['u1']};
    const changes: RunFilters[] = [];
    const {root, unmount} = mount(RunFiltersBar, {
      'modelValue': value, 'onUpdate:modelValue': (v: RunFilters) => changes.push(v),
      'tagOptions': [{value: 'feed', label: 'feed'}, {value: 'late', label: 'late'}],
      'authorOptions': [{value: 'u1', label: 'Maria'}],
    });
    try {
      const search = root.querySelector<HTMLInputElement>('.c2-runs-search')!;
      search.value = 'a';
      search.dispatchEvent(new Event('input'));
      const chips = () => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-tag'));
      expect(chips().map((chip) => chip.textContent?.trim()).join(), 'feed ×,late ×,Maria ×');
      chips()[0].click();
      chips()[2].click();
      expect(JSON.stringify(changes), JSON.stringify([{search: 'a', tags: ['feed', 'late'], authors: ['u1']},
        {search: '', tags: ['late'], authors: ['u1']}, {search: '', tags: ['feed', 'late'], authors: []}]));
      expect(JSON.stringify(value), JSON.stringify({search: '', tags: ['feed', 'late'], authors: ['u1']}));
    } finally {
      unmount();
    }
  });

  test('The list shows folders first, opens them and goes up, and shows results with their location', async () => {
    const run = DG.Func.byName(DAY_CALCULATION).prepare({});
    run.newId();
    run.options['title'] = 'Listed run';
    const events: string[] = [];
    const props = Vue.shallowReactive({
      crumbs: [{label: 'Spaces', depth: 0}, {label: 'Team', depth: 1}], canGoUp: true,
      folders: [{id: 'f1', name: 'Study', count: 2}], runs: [{run, spaceId: 's2'}],
      resultsLabel: undefined as string | undefined,
      isSelected: () => false, favorites: new Set<string>(), locationOf: () => 'Study',
      onEnter: (id: string) => events.push(`enter:${id}`), onUp: () => events.push('up'),
      onClearFilters: () => events.push('clear'),
    });
    const {root, unmount} = mount({render: () => Vue.h(RunList, props)}, {});
    try {
      const names = () => Array.from(root.querySelectorAll('.c2-runs-list .c2-runs-row-name'))
        .map((el) => el.textContent?.trim());
      expect(names().join('|'), '..|Study (2)|Listed run');
      root.querySelector('.c2-runs-folder')!.dispatchEvent(new MouseEvent('dblclick', {bubbles: true}));
      root.querySelector('.c2-runs-list')!
        .dispatchEvent(new KeyboardEvent('keydown', {key: 'Backspace', bubbles: true}));
      expect(events.join(), 'enter:f1,up');
      props.resultsLabel = 'Results in Team and its subspaces';
      await Vue.nextTick();
      expect(names().join('|'), 'Listed run');
      expect(root.querySelector('.c2-runs-row-location')?.textContent, 'Study');
      root.querySelector('.c2-runs-list')!
        .dispatchEvent(new KeyboardEvent('keydown', {key: 'Backspace', bubbles: true}));
      root.querySelector<HTMLElement>('.c2-runs-results .c2-runs-link')!.click();
      expect(events.join(), 'enter:f1,up,clear');
    } finally {
      unmount();
    }
  });

  test('A space that couldn\'t be listed shows a warning in the tree and as a folder', async () => {
    const index = makeIndex(['a'], new Map([['a', {space: space('a'), children: ['b'], runs: []}],
      ['b', {space: space('b'), children: [], runs: [], error: 'Permission denied'}]]));
    expect(spaceFolders(index, ['b'], 'a', undefined)[0].error, 'Permission denied');
    const {root, unmount} = mount(PlacesTree, {index, location: {myRuns: false, path: ['a', 'b']}});
    try {
      const warned = Array.from(root.querySelectorAll('.c2-runs-tree-row'))
        .filter((row) => row.querySelector('.c2-runs-error-icon'))
        .map((row) => row.querySelector('.c2-runs-tree-name')?.textContent);
      expect(JSON.stringify(warned), JSON.stringify(['b']));
    } finally {
      unmount();
    }
  });

  test('The details pane shows an edit once it is saved', () => withRunStorage(async (f) => {
    const name = uniqueName('Before edit');
    await f.saveDayRun(name);
    const {root, unmount} = mount(RunBrowser, {model: DG.Func.byName(DAY_CALCULATION)});
    try {
      const button = (label: string) => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-details button'))
        .find((item) => item.textContent?.trim() === label);
      const row = (text: string) => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-list .c2-runs-row'))
        .find((el) => el.textContent?.includes(text));
      await awaitCheck(() => !!row(name), 'the run should be listed', 15000);
      row(name)!.dispatchEvent(new MouseEvent('click', {bubbles: true, detail: 1}));
      await awaitCheck(() => !!button('Edit details'), 'the details should load', 15000);
      button('Edit details')!.click();
      await Vue.nextTick();

      const renamed = uniqueName('After edit');
      const input = root.querySelector<HTMLInputElement>('.c2-runs-form input')!;
      input.value = renamed;
      input.dispatchEvent(new Event('input'));
      button('Cancel')!.click();
      await Vue.nextTick();
      expect(root.querySelector('.c2-runs-details')!.textContent!.includes(name), true);

      button('Edit details')!.click();
      await Vue.nextTick();
      const again = root.querySelector<HTMLInputElement>('.c2-runs-form input')!;
      again.value = renamed;
      again.dispatchEvent(new Event('input'));
      button('Save details')!.click();
      await awaitCheck(() => !!button('Edit details') &&
        root.querySelector('.c2-runs-details')!.textContent!.includes(renamed), 'the pane should show the edit', 15000);
      await awaitCheck(() => !!row(renamed), 'the list should show the edit', 15000);
    } finally {
      unmount();
    }
  }));

  test('A search unselects the runs it hides, so the actions see only listed runs', () =>
    withRunStorage(async (f) => {
      const name = uniqueName('Searched');
      const kept = await f.saveDayRun(`${name} kept`);
      await f.saveDayRun(`${name} hidden`);
      const selections: string[][] = [];
      const {root, unmount} = mount(RunBrowser, {model: DG.Func.byName(DAY_CALCULATION),
        onSelectionChanged: (runs: DG.FuncCall[]) => selections.push(runs.map((run) => run.id))});
      try {
        const row = (suffix: string) => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-list .c2-runs-row'))
          .find((el) => el.textContent?.includes(`${name} ${suffix}`));
        await awaitCheck(() => !!row('kept') && !!row('hidden'), 'both runs should be listed', 15000);
        for (const suffix of ['kept', 'hidden'])
          row(suffix)!.dispatchEvent(new MouseEvent('click', {bubbles: true, ctrlKey: true, detail: 1}));
        await Vue.nextTick();
        expect(selections.at(-1)?.length, 2);
        const search = root.querySelector<HTMLInputElement>('.c2-runs-search')!;
        search.value = `${name} kept`;
        search.dispatchEvent(new Event('input'));
        await awaitCheck(() => !row('hidden'), 'the search should hide the other run', 5000);
        await Vue.nextTick();
        expect(JSON.stringify(selections.at(-1)), JSON.stringify([kept.id]));
      } finally {
        unmount();
      }
    }));

  test('Choosing runs quickly shows the last one, and choosing none clears it', () => withRunStorage(async (f) => {
    const first = await f.saveDayRun(uniqueName('Details first'));
    const second = await f.saveDayRun(uniqueName('Details second'));
    const scope = Vue.effectScope();
    const run = scope.run(() => useRunDetails(() => true, () => null))!;
    try {
      run.choose(first);
      run.choose(second);
      expect(run.isLoadingDetails.value, true);
      await awaitCheck(() => !run.isLoadingDetails.value, 'the details should load', 15000);
      expect(run.details.value?.run.id, second.id);

      run.choose(first);
      run.choose(undefined);
      expect(run.isLoadingDetails.value, false);
      expect(run.details.value == null, true);
      await historyUtils.shallowLoadRun(first.id);
      expect(run.details.value == null, true);

      run.choose({id: missingId()} as DG.FuncCall);
      await awaitCheck(() => !run.isLoadingDetails.value, 'the load should fail', 15000);
      expect(run.actionError.value?.startsWith('Couldn\'t load the run'), true);
      run.startEdit();
      expect(run.mode.value, 'edit');
      expect(run.actionError.value == null, true);
      run.stopEdit();
      expect(run.mode.value, 'view');
      run.choose(second);
      await awaitCheck(() => !!run.details.value, 'a later run should still load', 15000);
      expect(run.actionError.value == null, true);
    } finally {
      scope.stop();
    }
  }));

  test('My runs follow the models and run changes, and keep only the last load', () => withRunStorage(async (f) => {
    const day = DG.Func.byName(DAY_CALCULATION);
    const first = await f.saveDayRun(uniqueName('My runs first'));
    const models = Vue.shallowRef<DG.Func[]>([]);
    const scope = Vue.effectScope();
    const my = scope.run(() => useMyRuns(models, false))!;
    const isListed = (id: string) => my.runs.value.some((run) => run.id === id);
    try {
      expect(my.isLoading.value, false);
      models.value = [day];
      await awaitCheck(() => isListed(first.id), 'a chosen model should load its runs', 15000);
      expect(my.count.value, my.runs.value.length);

      const second = await f.saveDayRun(uniqueName('My runs second'));
      notifyRunsChanged({});
      await awaitCheck(() => isListed(second.id), 'a run change should load the runs again', 15000);

      notifyRunsChanged({});
      models.value = [];
      await Vue.nextTick();
      expect(my.isLoading.value, false);
      expect(my.count.value == null, true);
      expect(my.runs.value.length, 0);
      await listMyModelRuns(day);
      expect(my.runs.value.length, 0);
      expect(my.error.value == null, true);
    } finally {
      scope.stop();
    }
  }));

  test('Authors follow the runs and keep only the last load', () => withRunStorage(async (f) => {
    const saved = await f.saveDayRun(uniqueName('Authored'));
    const runs = Vue.shallowRef<DG.FuncCall[]>([saved]);
    const scope = Vue.effectScope();
    const {authors, error} = scope.run(() => useRunAuthors(runs))!;
    try {
      await awaitCheck(() => !!authors.value, 'the authors should load', 15000);
      expect(authors.value?.get(saved.id)?.id, DG.User.current().id);

      runs.value = [];
      await Vue.nextTick();
      expect(authors.value == null, true);

      runs.value = [saved];
      await Vue.nextTick();
      runs.value = [];
      await Vue.nextTick();
      await loadRunAuthors([saved]);
      expect(authors.value == null, true);
      expect(error.value == null, true);
    } finally {
      scope.stop();
    }
  }));

  test('The Save dialog lists the place chosen last', () =>
    withRunStorage(async (f) => {
      const [first, second] = [await f.space(), await f.space()];
      const firstRun = uniqueName('First listed');
      await f.saveDayRun(firstRun, first.id);
      const {root, unmount} = mount(SaveRunDialog, {func: DG.Func.byName(DAY_CALCULATION), name: 'New'});
      try {
        const treeRow = (place: DG.Project) => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-tree-row'))
          .find((row) => row.textContent?.includes(place.friendlyName || place.name));
        await awaitCheck(() => !!treeRow(first) && !!treeRow(second), 'the spaces should show in the tree', 15000);
        const listed = () => root.querySelector('.c2-runs-list')!.textContent ?? '';

        treeRow(first)!.click();
        await Vue.nextTick();
        treeRow(second)!.click();
        await awaitCheck(() => listed().includes('Nothing here yet'), 'the second space should be listed', 15000);
        await historyUtils.listSpace(first.id);
        expect(listed().includes(firstRun), false);
        const saveButton = () => root.querySelector<HTMLElement>('.ui-btn-ok')!;
        await awaitCheck(() => !saveButton().classList.contains('d4-disabled'),
          'the access check should allow the save', 5000);
      } finally {
        unmount();
      }
    }));

  test('The open model dropdown doesn\'t scroll sideways', async () => {
    const {root, unmount} = mount(ModelSelector, {models: [DG.Func.byName(DAY_CALCULATION)]});
    try {
      root.style.width = '400px';
      root.querySelector<HTMLElement>('.c2-runs-model-current')!.click();
      await Vue.nextTick();
      const dropdown = root.querySelector<HTMLElement>('.c2-runs-model-dropdown')!;
      expect(dropdown.scrollWidth, dropdown.clientWidth);
    } finally {
      unmount();
    }
  });

  test('Clicking a run of the space in the Save dialog offers its name to overwrite it', () =>
    withRunStorage(async (f) => {
      const space = await f.space();
      const taken = uniqueName('Taken');
      await f.saveDayRun(taken, space.id);
      const {root, unmount} = mount(SaveRunDialog, {func: DG.Func.byName(DAY_CALCULATION), name: 'My new run'});
      try {
        const spaceName = space.friendlyName || space.name;
        const treeRow = () => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-tree-row'))
          .find((row) => row.textContent?.includes(spaceName));
        await awaitCheck(() => !!treeRow(), 'the space should show in the tree', 15000);
      treeRow()!.click();
      const runRow = () => Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-list .c2-runs-row'))
        .find((row) => row.textContent?.includes(taken));
      await awaitCheck(() => !!runRow(), 'the space\'s run should be listed', 15000);
      runRow()!.click();
      await awaitCheck(() => !!root.querySelector('.c2-runs-warning'), 'clicking the run should warn', 3000);
      expect(root.querySelector('.c2-runs-warning')!.textContent!.includes(`Overwrite "${taken}"`), true);
      Array.from(root.querySelectorAll<HTMLElement>('.c2-runs-warning button'))
        .find((button) => button.textContent === 'Use this name')!.click();
      await Vue.nextTick();
      expect(root.querySelector<HTMLInputElement>('.c2-runs-form input')!.value, taken);
      expect(root.querySelector('.c2-runs-warning'), null);
      } finally {
        unmount();
      }
    }));
});
