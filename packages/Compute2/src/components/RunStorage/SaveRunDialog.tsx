import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {EMPTY, from} from 'rxjs';
import {watchExtractedObservable} from '@vueuse/rxjs';
import {
  historyUtils, RunMetadata, SaveTarget,
} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {PlacesTree} from './PlacesTree';
import {Breadcrumbs, FolderRows, spaceFolders} from './RunList';
import {useSpaceIndex} from './use-space-index';
import {MY_RUNS, SPACES, locationKey, useSpaceNavigation} from './use-space-navigation';
import {
  Place, SpaceAddOutcome, canEditSpace, entityKind, formatSaved, listMyModelRuns, modelLabel, parseTags, placeName,
  placeRunFilter, planSpaceAdd, spaceListingError,
} from './run-storage';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {showVueDialog} from '../../utils';
import './RunStorage.css';

export type SaveRunResult = {metadata: RunMetadata, target: SaveTarget, place: Place, outcomes?: SpaceAddOutcome[]};

type Listing = {runs: DG.FuncCall[], others: DG.Entity[]};

export const SaveRunDialog = Vue.defineComponent({
  name: 'SaveRunDialog',
  props: {
    func: {type: DG.Func, required: true},
    name: {type: String, default: ''},
    description: {type: String, default: ''},
    tags: {type: Array as Vue.PropType<string[]>, default: () => []},
    // Add to space: the run is linked as it is, so only a space is chosen
    linkedRun: {type: DG.FuncCall, required: false},
    // Add to space for several runs: each keeps its name
    linkedRuns: {type: Array as Vue.PropType<DG.FuncCall[]>, default: () => []},
  },
  emits: {
    done: (_result: SaveRunResult | null) => true,
  },
  setup(props, {emit}) {
    const isLinkingSeveral = props.linkedRuns.length > 0;
    const isLinking = !!props.linkedRun || isLinkingSeveral;
    const spaces = useSpaceIndex();
    // runs are added to spaces only, so linking starts at the spaces
    const nav = useSpaceNavigation(spaces.index, isLinking ? SPACES : MY_RUNS);
    const place = Vue.computed<Place>(() => nav.location.value.myRuns ? null : nav.space.value ?? null);
    const hasPlace = Vue.computed(() => nav.location.value.myRuns ? !isLinking : !!nav.space.value);
    const name = Vue.ref(props.name);
    const description = Vue.ref(props.description);
    const tagsText = Vue.ref(props.tags.join(', '));

    const listing = Vue.shallowRef<Listing>({runs: [], others: []});
    const listingError = Vue.ref<string>();
    // several runs may come from different models; the tree counts and the listing follow all of them
    const modelNames = new Set(isLinkingSeveral ?
      props.linkedRuns.map((run) => run.func?.name) : [props.func.name]);
    const isListedModel = (run: DG.FuncCall) => modelNames.has(run.func?.name);
    const runFilter = isLinkingSeveral ? isListedModel : placeRunFilter(props.func);
    const folders = Vue.computed(() => spaces.index.value ?
      spaceFolders(spaces.index.value, nav.folders.value, nav.spaceId.value, runFilter) : []);
    const myRunsCount = Vue.ref<number>();
    const isListing = Vue.ref(false);
    // undefined while the selected place's permission check runs
    const canWrite = Vue.ref<boolean | undefined>(true);
    const isChecking = Vue.ref(false);
    const error = Vue.ref<string>();
    // a run of the space picked in the listing, whose name the user may take to replace it
    const overwrite = Vue.shallowRef<DG.FuncCall>();

    const placeKey = Vue.computed(() => `${locationKey(nav.location.value)}/${place.value?.id}`);

    watchExtractedObservable(placeKey, () => {
      const selected = place.value;
      canWrite.value = selected ? undefined : true;
      return selected && hasPlace.value ? from(canEditSpace(selected).catch(() => false)) : EMPTY;
    }, (writable) => canWrite.value = writable, undefined, {immediate: true});

    const loadListing = async (selected: Place): Promise<Listing> => {
      const {runs, items} = selected ? await historyUtils.listSpace(selected.id) :
        {runs: await listMyModelRuns(props.func), items: [] as DG.Entity[]};
      return {
        runs: [...runs.filter(isListedModel), ...runs.filter((run) => !isListedModel(run))]
          .map((run) => Vue.markRaw(run)),
        // subspaces are the folders above the listing
        others: items.filter((item) => historyUtils.isNamedSpaceItem(item) &&
          !(item instanceof DG.Project && item.isSpace)).map((item) => Vue.markRaw(item)),
      };
    };

    watchExtractedObservable(placeKey, () => {
      error.value = undefined;
      overwrite.value = undefined;
      listingError.value = undefined;
      listing.value = {runs: [], others: []};
      isListing.value = hasPlace.value;
      return hasPlace.value ? from(loadListing(place.value)) : EMPTY;
    }, (loaded) => {
      listing.value = loaded;
      isListing.value = false;
      if (!place.value)
        myRunsCount.value = loaded.runs.length;
    }, {
      onError: (e: any) => {
        listingError.value = spaceListingError(String(e?.message ?? e));
        isListing.value = false;
      },
    }, {immediate: true});

    const alreadyLinked = Vue.computed(() =>
      !!props.linkedRun && listing.value.runs.some((run) => run.id === props.linkedRun!.id));

    const outcomes = Vue.computed(() => isLinkingSeveral && hasPlace.value && !isListing.value && !listingError.value ?
      planSpaceAdd(props.linkedRuns, {runs: listing.value.runs, items: listing.value.others}) : []);
    const addedCount = Vue.computed(() =>
      outcomes.value.filter((outcome) => outcome.kind === 'add' || outcome.kind === 'replace').length);

    const blocked = Vue.computed(() => {
      if (!hasPlace.value)
        return 'Choose a space';
      if (canWrite.value === undefined)
        return 'Checking access to this space…';
      if (!canWrite.value)
        return 'You can view this space but not save to it';
      if (alreadyLinked.value)
        return `This run is already in ${placeName(place.value)}`;
      if (isLinkingSeveral) {
        if (isListing.value || listingError.value)
          return isListing.value ? 'Loading the space…' : listingError.value;
        if (outcomes.value.some((outcome) => outcome.kind === 'blocked'))
          return 'Deselect the blocked runs, or choose another space';
        if (addedCount.value === 0)
          return `These runs are already in ${placeName(place.value)}`;
        return undefined;
      }
      if (!name.value.trim())
        return 'Enter a name';
      return undefined;
    });

    const save = async () => {
      if (blocked.value || isChecking.value)
        return;
      if (isLinkingSeveral) {
        const result: SaveRunResult = {metadata: {name: ''}, target: {spaceId: place.value!.id}, place: place.value,
          outcomes: outcomes.value};
        const replacing = outcomes.value.filter((outcome) => outcome.kind === 'replace');
        if (replacing.length > 0) {
          confirmThen(`Replace ${replacing.length} ${replacing.length > 1 ? 'runs' : 'run'}?`,
            `${placeName(place.value)} already has runs named ` +
            `${replacing.map((outcome) => `"${getRunName(outcome.run)}"`).join(', ')}. The selected runs take ` +
            'their places there. The old runs aren\'t deleted: they stay in their authors\' My runs and in any ' +
            'other space they were added to.', result);
        } else
          emit('done', result);
        return;
      }
      isChecking.value = true;
      error.value = undefined;
      try {
        const trimmed = name.value.trim();
        const metadata = {name: trimmed, description: description.value.trim(), tags: parseTags(tagsText.value)};
        // names repeat freely in My runs: only a space keeps them unique
        if (!place.value) {
          emit('done', {metadata, target: {}, place: null});
          return;
        }
        const conflict = await historyUtils.findNameConflict(trimmed, place.value.id);
        const where = placeName(place.value);
        if (conflict.other) {
          error.value = `A ${entityKind(conflict.other)} named ${trimmed} is already in ${where}. Choose another name.`;
          return;
        }
        const replace = conflict.runs.filter((run) => run.id !== props.linkedRun?.id);
        if (isLinking && conflict.runs.some((run) => run.id === props.linkedRun!.id) && replace.length === 0) {
          error.value = `This run is already in ${where}.`;
          return;
        }
        const result: SaveRunResult = {metadata, target: {spaceId: place.value.id, replace}, place: place.value};
        if (replace.length > 0)
          confirmReplace(replace, result);
        else
          emit('done', result);
      } catch (e: any) {
        error.value = String(e?.message ?? e);
      } finally {
        isChecking.value = false;
      }
    };

    // names repeat freely in My runs, so only a space's runs can be overwritten
    const pickRun = (run: DG.FuncCall) => {
      if (!isLinkingSeveral && place.value && isListedModel(run))
        overwrite.value = run;
    };
    const takeName = (run: DG.FuncCall) => {
      name.value = getRunName(run);
      overwrite.value = undefined;
      error.value = undefined;
    };

    const renderOverwrite = (run: DG.FuncCall) => (
      <div class='c2-runs-warning'>
        <span>{`Overwrite "${getRunName(run)}"? Saving under its name replaces it in ${placeName(place.value)}.`}</span>
        <div class='c2-runs-actions'>
          <button class='ui-btn' onClick={() => takeName(run)}>Use this name</button>
          <button class='ui-btn' onClick={() => overwrite.value = undefined}>Keep my name</button>
        </div>
      </div>
    );

    const renderListing = () => (
      <div class='c2-runs-list-pane'>
        <Breadcrumbs crumbs={nav.crumbs.value} onGoTo={nav.goTo}/>
        <div class='c2-runs-pane c2-runs-list'>
          <FolderRows canGoUp={!nav.location.value.myRuns && nav.location.value.path.length > 0}
            folders={folders.value}
            onEnter={nav.enter} onUp={nav.up}/>
          { !hasPlace.value && <div class='c2-runs-muted c2-runs-pad'>Choose a space</div> }
          { isListing.value && <div class='c2-runs-muted c2-runs-pad'>Loading…</div> }
          { listingError.value && <div class='c2-runs-error c2-runs-pad'>{listingError.value}</div> }
          { hasPlace.value && !isListing.value && !listingError.value && nav.folders.value.length === 0 &&
            listing.value.runs.length + listing.value.others.length === 0 &&
            <div class='c2-runs-muted c2-runs-pad'>Nothing here yet</div> }
          { listing.value.runs.map((run) => {
            const ofModel = isListedModel(run);
            return <div key={run.id}
              class={['c2-runs-row', ofModel ? '' : 'c2-runs-dimmed',
                overwrite.value?.id === run.id ? 'c2-runs-selected' : '']}
              onClick={() => pickRun(run)}>
              <span class='c2-runs-row-name'>{getRunName(run)}</span>
              <span class='c2-runs-muted'>
                {ofModel ? formatSaved(run) : `${run.func?.friendlyName || run.func?.name} run`}
              </span>
            </div>;
          })}
          { listing.value.others.map((item) =>
            <div key={item.id} class='c2-runs-row c2-runs-dimmed'>
              <span class='c2-runs-row-name'>{item.friendlyName || item.name}</span>
              <span class='c2-runs-muted'>{entityKind(item)}</span>
            </div>) }
        </div>
        { overwrite.value && renderOverwrite(overwrite.value) }
      </div>
    );

    // a modal confirmation keeps the Save dialog's fields as they were when Save was pressed
    const confirmThen = (title: string, text: string, result: SaveRunResult) => {
      const confirmation = ui.dialog(title);
      confirmation
        .add(ui.divText(text))
        .addButton('Replace', () => {
          confirmation.close();
          // the Save dialog closes on the next tick: chained centered dialogs otherwise leak a click-blocking overlay
          setTimeout(() => emit('done', result));
        })
        .show({modal: true, center: true, width: 460});
    };

    const confirmReplace = (replace: DG.FuncCall[], result: SaveRunResult) => {
      const name = result.metadata.name;
      const saved = replace.map(formatSaved).filter((date) => !!date);
      confirmThen(`Replace "${name}"?`,
        `${placeName(result.place)} already has ${replace.length > 1 ? `${replace.length} runs` :
          'a run'} named "${name}"${saved.length ? `, saved on ${saved.join(', ')}` : ''}. ` +
        'Your run takes its place there. The old run isn\'t deleted: it stays in its author\'s My runs and in ' +
        'any other space it was added to.', result);
    };

    const outcomeText = (outcome: SpaceAddOutcome) => {
      switch (outcome.kind) {
      case 'add': return 'Will be added';
      case 'skip': return 'Already there, skipped';
      case 'replace': return `Replaces the run with this name${outcome.replaces.map(formatSaved)
        .filter((date) => !!date).map((date) => `, saved on ${date}`).join('')}`;
      case 'blocked': return outcome.reason;
      }
    };

    const renderOutcomes = () => (
      <div class='c2-runs-outcomes'>
        <div class='c2-runs-label'>{`${props.linkedRuns.length} runs, each under its own name`}</div>
        { props.linkedRuns.map((run) => {
          const outcome = outcomes.value.find((item) => item.run.id === run.id);
          return <div key={run.id} class='c2-runs-outcome'>
            <span class='c2-runs-row-name'>{getRunName(run)}</span>
            { outcome && <span class={outcome.kind === 'blocked' ? 'c2-runs-error' : 'c2-runs-muted'}>
              {outcomeText(outcome)}</span> }
          </div>;
        }) }
      </div>
    );

    return () => (
      <div class='c2-runs-dialog'>
        <div class='c2-runs-split'>
          <div class='c2-runs-pane c2-runs-tree-pane'>
            <div class='c2-runs-label'>{isLinking ? 'Add to' : 'Save to'}</div>
            <PlacesTree index={spaces.index.value} error={spaces.error.value} location={nav.location.value}
              showMyRuns={!isLinking} runFilter={runFilter} myRunsCount={myRunsCount.value}
              onSelect={nav.go} onRetry={spaces.reload}/>
          </div>
          { renderListing() }
        </div>
        <div class='c2-runs-form'>
          { isLinkingSeveral && renderOutcomes() }
          { !isLinkingSeveral && <label class='c2-runs-field'>
            <span>Name</span>
            <input class='ui-input-editor' value={name.value}
              onInput={(e: Event) => {name.value = (e.target as HTMLInputElement).value; error.value = undefined;}}
              onKeydown={(e: KeyboardEvent) => e.key === 'Enter' && save()}/>
          </label> }
          { !isLinking && <label class='c2-runs-field'>
            <span>Description</span>
            <textarea class='ui-input-editor' rows={5} value={description.value}
              onInput={(e: Event) => description.value = (e.target as HTMLTextAreaElement).value}/>
          </label> }
          { !isLinking && <label class='c2-runs-field'>
            <span>Tags</span>
            <input class='ui-input-editor' placeholder='Comma-separated' value={tagsText.value}
              onInput={(e: Event) => tagsText.value = (e.target as HTMLInputElement).value}/>
          </label> }
          { error.value && <div class='c2-runs-error'>{error.value}</div> }
          { !error.value && hasPlace.value && canWrite.value === false &&
            <div class='c2-runs-error'>You can view this space but not save to it</div> }
          { !error.value && canWrite.value && alreadyLinked.value &&
            <div class='c2-runs-muted'>This run is already in {placeName(place.value)}. Choose another space.</div> }
        </div>
        <div class='c2-runs-buttons'>
          <span class='c2-runs-muted c2-runs-footer-note'>
            {isLinkingSeveral ? '' : modelLabel(props.func)}</span>
          <div class='c2-runs-button-group'>
            <button class='ui-btn' onClick={() => emit('done', null)}>Cancel</button>
            <button class={['ui-btn', 'ui-btn-ok', blocked.value || isChecking.value ? 'd4-disabled' : '']}
              title={blocked.value} onClick={save}>
              {isLinkingSeveral ? `Add (${outcomes.value.length ? addedCount.value : props.linkedRuns.length})` :
                (isLinking ? 'Add' : 'Save')}</button>
          </div>
        </div>
      </div>
    );
  },
});

export function openSaveRunDialog(options: {func: DG.Func, name: string, description?: string, tags?: string[]}) {
  return showVueDialog<SaveRunResult>('Save run', SaveRunDialog, options, {width: 760, height: 620});
}
