import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {IconFA} from '@datagrok-libraries/webcomponents-vue';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {RunAuthor, formatSaved, placeName, spaceListingError} from './run-storage';
import {RunFilter, SpaceIndex, parentsOf, subtreeRuns} from './space-index';
import {Crumb} from './use-space-navigation';

export type Folder = {id: string, name: string, count?: number, alsoIn?: string, error?: string};
export type ListedRun = {run: DG.FuncCall, spaceId?: string};

/** The subspaces of a location as folders, with the runs runFilter keeps in each and below it. */
export function spaceFolders(index: SpaceIndex, ids: string[], parentId: string | undefined,
  runFilter: RunFilter | undefined): Folder[] {
  return ids.flatMap((id) => {
    const node = index.nodes.get(id);
    if (!node)
      return [];
    const others = parentsOf(index, id).filter((parent) => parent.id !== parentId);
    return [{
      id, name: placeName(node.space),
      count: runFilter ? subtreeRuns(index, [id], runFilter).length : undefined,
      alsoIn: parentId && others.length > 0 ? others.map(placeName).join(', ') : undefined,
      error: node.error,
    }];
  });
}

/** The way back up from the location: a link per level. */
export const Breadcrumbs = Vue.defineComponent({
  name: 'Breadcrumbs',
  props: {crumbs: {type: Array as Vue.PropType<Crumb[]>, required: true}},
  emits: {goTo: (_depth: number) => true},
  setup(props, {emit}) {
    return () => (
      <div class='c2-runs-crumbs'>
        { props.crumbs.map((crumb, i) => <>
          { i > 0 && <span class='c2-runs-muted'>›</span> }
          { i < props.crumbs.length - 1 ?
            <span class='c2-runs-link' onClick={() => emit('goTo', crumb.depth)}>{crumb.label}</span> :
            <span>{crumb.label}</span> }
        </>) }
      </div>
    );
  },
});

/** `..` and the subspaces of the location, first in a listing; a double click or Enter opens one. */
export const FolderRows = Vue.defineComponent({
  name: 'FolderRows',
  props: {
    canGoUp: {type: Boolean, default: false},
    folders: {type: Array as Vue.PropType<Folder[]>, required: true},
  },
  emits: {enter: (_id: string) => true, up: () => true},
  setup(props, {emit}) {
    return () => <>
      { props.canGoUp &&
        <div class='c2-runs-row c2-runs-up' tabindex={0} title='Up one level' onDblclick={() => emit('up')}
          onKeydown={(e: KeyboardEvent) => e.key === 'Enter' && emit('up')}>
          <IconFA name='level-up-alt' class='c2-runs-row-icon c2-runs-folder-icon'/>
          <span class='c2-runs-row-name'>..</span>
        </div> }
      { props.folders.map((folder) =>
        <div key={folder.id} class='c2-runs-row c2-runs-folder' tabindex={0} title='Double-click to open'
          onDblclick={() => emit('enter', folder.id)}
          onKeydown={(e: KeyboardEvent) => e.key === 'Enter' && emit('enter', folder.id)}>
          <IconFA name='folder' class='c2-runs-row-icon c2-runs-folder-icon'/>
          <span class='c2-runs-row-name'>{folder.name}
            { folder.alsoIn && <IconFA name='link' class='c2-runs-muted' tooltip={`Also in ${folder.alsoIn}`}/> }
            { folder.error && <IconFA name='exclamation-triangle' class='c2-runs-error-icon'
              tooltip={spaceListingError(folder.error)}/> }
            { folder.count !== undefined && <span class='c2-runs-muted'> ({folder.count})</span> }
          </span>
        </div>) }
    </>;
  },
});

/** The location's folders and runs, or, while a filter is on, the matching runs of the location and below it. */
export const RunList = Vue.defineComponent({
  name: 'RunList',
  props: {
    crumbs: {type: Array as Vue.PropType<Crumb[]>, required: true},
    canGoUp: {type: Boolean, default: false},
    folders: {type: Array as Vue.PropType<Folder[]>, default: () => []},
    runs: {type: Array as Vue.PropType<ListedRun[]>, required: true},
    // shown while a filter is on, with a link that clears the filters; the filters themselves are the filters slot
    resultsLabel: {type: String, required: false},
    status: {type: String, required: false},
    error: {type: String, required: false},
    isSelected: {type: Function as Vue.PropType<(run: DG.FuncCall) => boolean>, required: true},
    favorites: {type: Set as Vue.PropType<Set<string>>, required: true},
    modelOf: {type: Function as Vue.PropType<(run: DG.FuncCall) => string>, required: false},
    // space runs only: in My runs every run is the user's
    authors: {type: Map as Vue.PropType<Map<string, RunAuthor>>, required: false},
    showAuthors: {type: Boolean, default: false},
    authorsError: {type: String, required: false},
    locationOf: {type: Function as Vue.PropType<(spaceId?: string) => string>, required: false},
  },
  emits: {
    goTo: (_depth: number) => true,
    enter: (_id: string) => true,
    up: () => true,
    clickRun: (_run: DG.FuncCall, _e: MouseEvent) => true,
    toggleFavorite: (_run: DG.FuncCall) => true,
    clearFilters: () => true,
  },
  setup(props, {emit, slots}) {
    const upOnBackspace = (e: KeyboardEvent) => {
      if (e.key === 'Backspace' && props.canGoUp && !props.resultsLabel) {
        e.preventDefault();
        emit('up');
      }
    };

    const renderRun = ({run, spaceId}: ListedRun) => {
      const isFavorite = props.favorites.has(run.id);
      return <div key={run.id} class={['c2-runs-row', props.isSelected(run) ? 'c2-runs-selected' : '']}
        onClick={(e: MouseEvent) => emit('clickRun', run, e)}>
        <IconFA name='star' faStyle={isFavorite ? 'fas' : 'fal'}
          class={['c2-runs-row-icon', isFavorite ? 'c2-runs-star c2-runs-star-on' : 'c2-runs-star']}
          tooltip={isFavorite ? 'Remove from favorites' : 'Add to favorites'}
          onClick={(e: MouseEvent) => {e.stopPropagation(); emit('toggleFavorite', run);}}/>
        <span class='c2-runs-row-name' title={getRunName(run)}>{getRunName(run)}</span>
        { props.modelOf && <span class='c2-runs-muted c2-runs-row-model'>{props.modelOf(run)}</span> }
        { props.resultsLabel && props.locationOf &&
          <span class='c2-runs-row-location' title={props.locationOf(spaceId)}>{props.locationOf(spaceId)}</span> }
        { props.showAuthors &&
          <span class='c2-runs-muted c2-runs-row-author' title={props.authors?.get(run.id)?.name}>
            {props.authors ? (props.authors.get(run.id)?.name ?? '') : (props.authorsError ? '' : '…')}</span> }
        <span class='c2-runs-muted'>{formatSaved(run)}</span>
      </div>;
    };

    // the columns the runs show, as grid tracks; folders span all but the icon
    const columns = () => [
      {title: '', track: '16px'},
      {title: 'Name', track: 'minmax(80px, 1fr)'},
      ...props.modelOf ? [{title: 'Model', track: 'minmax(60px, 120px)'}] : [],
      ...props.resultsLabel && props.locationOf ? [{title: 'Location', track: 'minmax(60px, 120px)'}] : [],
      ...props.showAuthors ? [{title: 'Author', track: 'minmax(60px, 110px)'}] : [],
      {title: 'Saved', track: '110px'},
    ];

    return () => {
      const showFolders = !props.resultsLabel;
      const shown = columns();
      const empty = !props.status && !props.error && props.runs.length === 0 &&
        !(showFolders && (props.folders.length > 0 || props.canGoUp));
      return <div class='c2-runs-list-pane'>
        <Breadcrumbs crumbs={props.crumbs} onGoTo={(depth: number) => emit('goTo', depth)}/>
        {slots.filters?.()}
        { props.resultsLabel &&
          <div class='c2-runs-results c2-runs-muted'>{props.resultsLabel} ·
            <span class='c2-runs-link' onClick={() => emit('clearFilters')}> Clear filters</span></div> }
        <div class='c2-runs-pane c2-runs-list c2-runs-grid' tabindex={0} onKeydown={upOnBackspace}
          style={{'--c2-runs-columns': shown.map((column) => column.track).join(' ')}}>
          { (props.runs.length > 0 || (showFolders && props.folders.length > 0)) &&
            <div class='c2-runs-header'>{shown.map((column) => <span>{column.title}</span>)}</div> }
          { props.status && <div class='c2-runs-muted c2-runs-pad'>{props.status}</div> }
          { props.error && <div class='c2-runs-error c2-runs-pad'>{props.error}</div> }
          { props.authorsError && <div class='c2-runs-error c2-runs-pad'>{props.authorsError}</div> }
          { showFolders && !props.status &&
            <FolderRows canGoUp={props.canGoUp} folders={props.folders}
              onEnter={(id: string) => emit('enter', id)} onUp={() => emit('up')}/> }
          { empty &&
            <div class='c2-runs-muted c2-runs-pad'>
              {props.resultsLabel ? 'No runs match the filters' : 'No runs here'}</div> }
          { props.runs.map(renderRun) }
        </div>
      </div>;
    };
  },
});
