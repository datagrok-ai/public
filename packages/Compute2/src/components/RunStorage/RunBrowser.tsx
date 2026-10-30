import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {ResizeHandle} from '@datagrok-libraries/webcomponents-vue';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {getModelFilter} from '@datagrok-libraries/compute-utils/model-catalog/src/model-handler';
import {PlacesTree} from './PlacesTree';
import {ModelSelector} from './ModelSelector';
import {addRunToSpace, addSeveralToSpace} from './add-to-space';
import {RunDetailsForm} from './RunDetailsForm';
import {RunFiltersBar} from './RunFiltersBar';
import {ListedRun, RunList, spaceFolders} from './RunList';
import {
  Place, bySavedDesc, deleteOrRemoveRuns, formatSaved, modelLabel, placeName, placeRunFilter, runVersion,
  notifyRunsChanged, spaceListingError,
} from './run-storage';
import {subtreeRuns} from './space-index';
import {useSpaceIndex} from './use-space-index';
import {locationKey, useSpaceNavigation} from './use-space-navigation';
import {NO_FILTERS, RunFilters, useRunFilters} from './use-run-filters';
import {useRunDetails} from './use-run-details';
import {useRunSelection} from './use-run-selection';
import {useMyRuns} from './use-my-runs';
import {useRunAuthors} from './use-run-authors';
import {openRunComparison} from '../RunComparison/comparison-views';
import './RunStorage.css';

/** Opens and manages saved runs of one model, in My runs and in spaces, browsed like folders. */
export const RunBrowser = Vue.defineComponent({
  name: 'RunBrowser',
  props: {
    // the model whose runs are listed; emits update:model when the user picks another one
    model: {type: DG.Func, required: false},
    // shows the model filter, with the catalog models and extraModels
    selectModel: {type: Boolean, default: false},
    // with no model chosen, lists the runs of every model in the selector
    allowAllModels: {type: Boolean, default: false},
    extraModels: {type: Array as Vue.PropType<DG.Func[]>, default: () => []},
    roleOnlyFilter: {type: Boolean, default: false},
    layout: {type: String as Vue.PropType<'dialog' | 'docked'>, default: 'docked'},
    stepHistory: {type: Boolean, default: false},
    allowCompare: {type: Boolean, default: true},
    allowOpen: {type: Boolean, default: true},
    // hosts that open each run in a new view open all the selected ones; the rest load one run into their view
    openSeveral: {type: Boolean, default: false},
    showDetails: {type: Boolean, default: true},
    openLabel: {type: String, default: 'Open'},
  },
  emits: {
    'open': (_run: DG.FuncCall) => true,
    'openSeveral': (_runs: DG.FuncCall[]) => true,
    'update:model': (_model: DG.Func | undefined) => true,
    'selectionChanged': (_runs: DG.FuncCall[]) => true,
    'cancel': () => true,
  },
  setup(props, {emit}) {
    const catalogModels = Vue.shallowRef<DG.Func[]>([]);
    const isLoadingModels = Vue.ref(false);
    const model = Vue.shallowRef<DG.Func | undefined>(props.model);
    Vue.watch(() => props.model, (func) => model.value = func);
    const modelList = Vue.computed(() => [
      ...catalogModels.value,
      ...props.extraModels.filter((extra) => !catalogModels.value.some((func) => func.nqName === extra.nqName)),
    ]);
    const chooseModel = (func: DG.Func | undefined) => {
      model.value = func;
      emit('update:model', func);
    };
    const listedModels = Vue.computed(() =>
      model.value ? [model.value] : (props.allowAllModels ? modelList.value : []));
    // space listings come without the package, so a run is matched to a model by its function name
    const modelOf = (run: DG.FuncCall) =>
      model.value ?? modelList.value.find((func) => historyUtils.isRunOf(func, run));
    const runFilter = Vue.computed(() => listedModels.value.length > 0 ?
      placeRunFilter(listedModels.value, props.stepHistory) : undefined);

    const spaces = useSpaceIndex();
    const nav = useSpaceNavigation(spaces.index);
    const isMyRuns = Vue.computed(() => nav.location.value.myRuns);
    const nodeOf = (spaceId?: string) => spaceId ? spaces.index.value?.nodes.get(spaceId) : undefined;
    const folders = Vue.computed(() => spaces.index.value ?
      spaceFolders(spaces.index.value, nav.folders.value, nav.spaceId.value, runFilter.value) : []);
    const spaceError = () => {
      const error = spaces.error.value ?? nodeOf(nav.spaceId.value)?.error;
      return error ? spaceListingError(error) : undefined;
    };

    const {runs: myRuns, count: myRunsCount, isLoading, error: listError} =
      useMyRuns(() => listedModels.value, () => props.stepHistory);
    const favorites = Vue.ref(new Set<string>());
    const treeHeight = Vue.ref(160);

    // what a filter searches: My runs, or the location's runs and those of everything below it
    const scopeRuns = Vue.computed<ListedRun[]>(() => {
      if (isMyRuns.value)
        return myRuns.value.map((run) => ({run}));
      const index = spaces.index.value;
      return index && runFilter.value ? subtreeRuns(index, nav.scope.value, runFilter.value)
        .sort((a, b) => bySavedDesc(a.run, b.run)) : [];
    });
    const folderRuns = Vue.computed<ListedRun[]>(() => {
      if (isMyRuns.value)
        return scopeRuns.value;
      const id = nav.spaceId.value;
      return scopeRuns.value.filter((item) => item.spaceId === id);
    });
    const {current, details, isLoadingDetails, mode, actionError, choose, startEdit, stopEdit} = useRunDetails(
      () => props.showDetails, (run) => spaceOf(run));
    const {authors, error: authorsError} =
      useRunAuthors(() => isMyRuns.value ? [] : scopeRuns.value.map((item) => item.run));
    const {filters, tagOptions, authorOptions, isActive, matches} =
      useRunFilters(() => scopeRuns.value.map((item) => item.run), authors);
    const listed = Vue.computed(() => isActive.value ?
      scopeRuns.value.filter(({run}) => matches(run)) : folderRuns.value);
    const {selection, isSelected, click, actionRuns, clear, keepOnly} =
      useRunSelection(() => listed.value.map(({run}) => run));
    const spaceOf = (run: DG.FuncCall): Place => isMyRuns.value ? null :
      nodeOf(scopeRuns.value.find((item) => item.run.id === run.id)?.spaceId ?? nav.spaceId.value)?.space ?? null;

    Vue.onMounted(async () => {
      grok.dapi.entities.getFavorites()
        .then((items) => favorites.value = new Set(items.filter((item) => item.entityType === 'FuncCall')
          .map((item) => item.id)))
        .catch(() => {});
      if (!props.selectModel)
        return;
      isLoadingModels.value = true;
      try {
        catalogModels.value = (await grok.dapi.functions.filter(getModelFilter(props.roleOnlyFilter)).list())
          .map((func) => Vue.markRaw(func));
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        isLoadingModels.value = false;
      }
    });

    Vue.watch(() => locationKey(nav.location.value), () => {
      choose(undefined);
      clear();
      filters.value = {...filters.value, authors: []};
    });
    // a run the listing hides stops being shown too
    Vue.watch(listed, (items) => {
      if (current.value && !items.some((item) => item.run.id === current.value!.id))
        choose(undefined);
    });
    Vue.watch(selection, (selected) => emit('selectionChanged', selected));

    const clickRun = (run: DG.FuncCall, e: MouseEvent) => {
      const opens = click(run, e);
      choose(run);
      if (opens)
        openRun(run);
    };

    const toggleFavorite = async (run: DG.FuncCall) => {
      const next = new Set(favorites.value);
      try {
        if (next.has(run.id)) {
          await DG.Favorites.remove(run);
          next.delete(run.id);
        } else {
          await DG.Favorites.add(run);
          next.add(run.id);
        }
        favorites.value = next;
      } catch (e: any) {
        grok.shell.error(`Couldn't change the favorite: ${e?.message ?? e}`);
      }
    };

    const runsToDelete = () => actionRuns(details.value?.run ?? current.value!);
    // runs found below the location are removed from the space each is in, so they must share one
    const removalSpace = (runs: DG.FuncCall[]): Place | undefined => {
      if (isMyRuns.value)
        return null;
      const ids = new Set(runs.map((run) => spaceOf(run)?.id));
      return ids.size === 1 ? spaceOf(runs[0]) ?? undefined : undefined;
    };

    const deleteRuns = async (runs: DG.FuncCall[], space: Place) => {
      const {failed} = await deleteOrRemoveRuns(runs, space);
      if (failed.length > 0) {
        const names = failed.map(({run}) => getRunName(run)).join(', ');
        grok.shell.error(`Couldn't ${space ? `remove ${names} from ${placeName(space)}` : `delete ${names}`}: ` +
          failed[0].error);
      }
      // the runs that failed stay selected, so they can be tried again
      keepOnly(failed.map(({run}) => run));
      choose(failed.length === 1 ? failed[0].run : undefined);
      // a run deleted from My runs stays in its spaces
      notifyRunsChanged({spaceIds: space ? [space.id] : []});
    };

    const askDelete = () => {
      const runs = runsToDelete();
      const space = removalSpace(runs);
      if (space === undefined)
        return;
      const several = runs.length > 1;
      const what = several ? `${runs.length} runs` : `"${getRunName(runs[0])}"`;
      const title = space ? `Remove from ${placeName(space)}` : `Delete ${several ? `${runs.length} runs` : 'run'}`;
      const text = space ?
        `Remove ${what} from ${placeName(space)}? ${several ? 'The runs themselves don\'t change and stay' :
          'The run itself doesn\'t change and stays'} in My runs of ${several ? 'their authors' : 'its author'}.` :
        `Delete ${what} from My runs? A deleted run can't be restored to My runs. ` +
          `Spaces ${several ? 'they were' : 'it was'} added to keep ${several ? 'them' : 'it'}.`;
      ui.dialog(title)
        .add(ui.divText(text))
        .onOK(() => deleteRuns(runs, space))
        .show({center: true, modal: true, width: 460});
    };

    const addToSpace = async () => {
      const runs = runsToDelete();
      const first = runs[0];
      const added = runs.length > 1 ?
        await addSeveralToSpace(modelOf(first) ?? first.func, runs) :
        await addRunToSpace(modelOf(first) ?? first.func, first);
      if (added)
        notifyRunsChanged({spaceIds: [added.id]});
    };

    // a comparison has one model to start from, so runs of several models aren't compared together
    const compareModel = Vue.computed(() => {
      const names = new Set(selection.value.map((run) => run.func?.name));
      return names.size === 1 ? modelOf(selection.value[0]) : undefined;
    });
    const compareSelected = () => {
      openRunComparison({runIds: selection.value.map((run) => run.id), model: compareModel.value});
      emit('cancel');
    };
    // the run as listed, without its tables: most hosts open runs by id, the rest load what they need
    const openRun = (run: DG.FuncCall) => {
      if (props.allowOpen)
        emit('open', run);
    };
    const canOpen = Vue.computed(() => !!current.value && (props.openSeveral || selection.value.length <= 1));
    const openSelected = () => {
      if (!canOpen.value)
        return;
      if (props.openSeveral && selection.value.length > 1)
        emit('openSeveral', selection.value);
      else
        openRun(current.value!);
    };
    const footerNote = () => {
      const count = selection.value.length;
      if (props.allowCompare && count > 0 && !compareModel.value)
        return 'Select runs of one model to compare them';
      const openOne = props.allowOpen && !props.openSeveral;
      if (count > 1)
        return openOne ? `${count} runs selected. Select one run to open it` : `${count} runs selected`;
      return props.allowCompare ? 'Ctrl+click to select several runs to compare' : '';
    };

    // what a search covers, as the list and the search box name it
    const scopeLabel = () => {
      const location = nav.location.value;
      if (location.myRuns)
        return 'My runs';
      if (location.path.length === 0)
        return 'all spaces';
      return `${nav.crumbs.value[nav.crumbs.value.length - 1].label} and its subspaces`;
    };
    const listStatus = () => {
      if (listedModels.value.length === 0)
        return 'Select a model';
      if (isMyRuns.value)
        return isLoading.value && myRuns.value.length === 0 ? 'Loading…' : undefined;
      return !spaces.index.value && !spaces.error.value ? 'Loading spaces…' : undefined;
    };

    const renderTree = () => (
      <div class='c2-runs-pane c2-runs-tree-pane'
        style={props.layout === 'docked' ? {height: `${treeHeight.value}px`} : {}}>
        <PlacesTree index={spaces.index.value} error={spaces.error.value} location={nav.location.value}
          runFilter={runFilter.value} myRunsCount={myRunsCount.value}
          onSelect={nav.go} onRetry={spaces.reload}/>
      </div>
    );

    const renderList = () => {
      const location = nav.location.value;
      return <RunList crumbs={nav.crumbs.value}
        canGoUp={!location.myRuns && location.path.length > 0}
        folders={folders.value}
        runs={listed.value}
        resultsLabel={isActive.value ? `Results in ${scopeLabel()}` : undefined}
        status={listStatus()} error={isMyRuns.value ? listError.value : spaceError()}
        isSelected={isSelected} favorites={favorites.value}
        modelOf={model.value ? undefined : (run: DG.FuncCall) => {
          const func = modelOf(run);
          return func ? modelLabel(func) : (run.func?.friendlyName || run.func?.name);
        }}
        authors={authors.value} showAuthors={!isMyRuns.value} authorsError={authorsError.value}
        locationOf={(spaceId?: string) =>
          spaceId === nav.spaceId.value ? 'here' : (nodeOf(spaceId) ? placeName(nodeOf(spaceId)!.space) : '')}
        onGoTo={nav.goTo} onEnter={nav.enter} onUp={nav.up}
        onClickRun={clickRun} onToggleFavorite={toggleFavorite} onClearFilters={() => filters.value = NO_FILTERS}>
        {{filters: () => <RunFiltersBar modelValue={filters.value}
          onUpdate:modelValue={(value: RunFilters) => filters.value = value}
          tagOptions={tagOptions.value} authorOptions={authorOptions.value}
          searchPlaceholder={`Search in ${scopeLabel()}`}/>}}
      </RunList>;
    };

    const renderDetails = () => {
      if (!current.value) {
        return <div class='c2-runs-pane c2-runs-details c2-runs-muted c2-runs-pad'>
          Select a run to see its details</div>;
      }
      const info = details.value;
      const run = info?.run ?? current.value;
      const {description, tags} = historyUtils.getRunMetadata(run);
      const version = runVersion(run);
      const runs = runsToDelete();
      const count = runs.length;
      const space = removalSpace(runs);
      const deleteButton = info && space !== undefined && (!space || info.canUnlink) &&
        <button class='ui-btn' onClick={askDelete}
          title={space ? `Remove from ${placeName(space)}` : 'Delete from My runs'}>
          {`${space ? 'Remove from space' : 'Delete'}${count > 1 ? ` (${count})` : ''}`}</button>;
      // several runs share only the actions that take them all; details are edited one run at a time
      if (count > 1) {
        return <div class='c2-runs-pane c2-runs-details'>
          <div class='c2-runs-details-name'>{`${count} runs selected`}</div>
          <div class='c2-runs-muted'>Select one run to see or edit its details.</div>
          { space === undefined && <div class='c2-runs-muted'>The runs are in different spaces.</div> }
          { actionError.value && <div class='c2-runs-error'>{actionError.value}</div> }
          <div class='c2-runs-actions'>
            { isMyRuns.value && <button class='ui-btn' onClick={addToSpace}>{`Add to space… (${count})`}</button> }
            { deleteButton }
          </div>
        </div>;
      }
      const runSpace = spaceOf(run);
      return <div class='c2-runs-pane c2-runs-details'>
        <div class='c2-runs-details-name'>{getRunName(run)}</div>
        { description && <div>{description}</div> }
        <div class='c2-runs-props'>
          <span>Saved</span><span>{formatSaved(run)}</span>
          <span>By</span><span>{isLoadingDetails.value ? '…' : (info?.author ?? '')}</span>
          { info?.modifiedBy && <><span>Modified by</span><span>{info.modifiedBy}</span></> }
          { runSpace && <><span>Space</span><span>{placeName(runSpace)}</span></> }
          { tags.length > 0 && <><span>Tags</span><span>{tags.join(', ')}</span></> }
          { version && <><span>Version</span><span>{version}</span></> }
        </div>
        { actionError.value && <div class='c2-runs-error'>{actionError.value}</div> }
        { info && mode.value === 'edit' &&
          <RunDetailsForm key={run.id} run={run} place={runSpace}
            onDone={(saved?: DG.FuncCall) => saved ? choose(saved) : stopEdit()}/> }
        { info && mode.value === 'view' &&
          <div class='c2-runs-actions'>
            { info.canEdit && <button class='ui-btn' onClick={startEdit}>Edit details</button> }
            { isMyRuns.value && <button class='ui-btn' onClick={addToSpace}>Add to space…</button> }
            { deleteButton }
          </div> }
      </div>;
    };

    const renderFooter = () => (
      <div class='c2-runs-buttons'>
        <span class='c2-runs-muted c2-runs-footer-note'>{footerNote()}</span>
        <div class='c2-runs-button-group'>
          { props.allowCompare &&
            <button class={['ui-btn', compareModel.value ? '' : 'd4-disabled']}
              onClick={() => compareModel.value && compareSelected()}>Compare</button> }
          { props.layout === 'dialog' && <button class='ui-btn' onClick={() => emit('cancel')}>Cancel</button> }
          { props.allowOpen &&
            <button class={['ui-btn', 'ui-btn-ok', canOpen.value ? '' : 'd4-disabled']}
              title={canOpen.value ? undefined : (current.value ? 'Select one run to open' : 'Select a run')}
              onClick={openSelected}>
              {props.openSeveral && selection.value.length > 1 ?
                `${props.openLabel} (${selection.value.length})` : props.openLabel}</button> }
        </div>
      </div>
    );

    return () => (
      <div class={['c2-runs-browser', `c2-runs-${props.layout}`]}>
        { props.selectModel &&
          <ModelSelector models={modelList.value} modelValue={model.value} isLoading={isLoadingModels.value}
            allLabel={props.allowAllModels ? 'All models' : undefined} onUpdate:modelValue={chooseModel}/> }
        { props.layout === 'dialog' ?
          <div class='c2-runs-split c2-runs-split-3'>
            {renderTree()}
            {renderList()}
            {renderDetails()}
          </div> :
          <>
            {renderTree()}
            <ResizeHandle axis='y' size={treeHeight.value} min={60} max={600}
              onUpdate:size={(size: number) => treeHeight.value = size}/>
            {renderList()}
            {props.showDetails && renderDetails()}
          </> }
        {renderFooter()}
      </div>
    );
  },
});
