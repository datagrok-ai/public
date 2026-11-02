import * as Vue from 'vue';
import {IconFA} from '@datagrok-libraries/webcomponents-vue';
import {placeName, spaceListingError} from './run-storage';
import {RunFilter, SpaceIndex, favoriteSpaces, parentsOf, subtreeCounts} from './space-index';
import {Favorites} from './use-favorites';
import {RunsLocation, locationKey} from './use-space-navigation';

const pathKey = (path: string[]) => locationKey({myRuns: false, path});

/** My runs and the spaces of an index, with the runs that runFilter keeps in each space and below it. */
export const PlacesTree = Vue.defineComponent({
  name: 'PlacesTree',
  props: {
    index: {type: Object as Vue.PropType<SpaceIndex>, required: false},
    error: {type: String, required: false},
    location: {type: Object as Vue.PropType<RunsLocation>, required: true},
    showMyRuns: {type: Boolean, default: true},
    runFilter: {type: Function as Vue.PropType<RunFilter>, required: false},
    // the host lists My runs itself, so it passes the count
    myRunsCount: {type: Number, required: false},
    favorites: {type: Object as Vue.PropType<Favorites>, required: true},
  },
  emits: {
    select: (_location: RunsLocation) => true,
    retry: () => true,
  },
  setup(props, {emit}) {
    const expanded = Vue.ref(new Set<string>());
    // the selected space's parents open, so a space entered from the list shows in the tree
    Vue.watch(() => props.location, (location) => {
      if (location.myRuns || location.path.length < 2)
        return;
      const next = new Set(expanded.value);
      for (let depth = 1; depth < location.path.length; depth++)
        next.add(pathKey(location.path.slice(0, depth)));
      expanded.value = next;
    }, {immediate: true});
    const toggle = (key: string) => {
      const next = new Set(expanded.value);
      next.has(key) ? next.delete(key) : next.add(key);
      expanded.value = next;
    };
    const counts = Vue.computed(() => props.index && props.runFilter ?
      subtreeCounts(props.index, props.runFilter) : undefined);
    const selectedKey = Vue.computed(() => locationKey(props.location));
    const onlyFavorites = Vue.ref(false);
    const kept = Vue.computed(() => onlyFavorites.value && props.index ?
      favoriteSpaces(props.index, props.favorites.ids.value) : undefined);

    const renderSpace = (index: SpaceIndex, id: string, parentPath: string[]): Vue.VNode[] => {
      const node = index.nodes.get(id);
      if (!node)
        return [];
      const path = [...parentPath, id];
      const key = pathKey(path);
      const children = kept.value ? node.children.filter((child) => kept.value!.has(child)) : node.children;
      // the favorites-only tree stays open: a closed parent would hide the favorite it is there for
      const canToggle = !kept.value && children.length > 0;
      const isOpen = !!kept.value || expanded.value.has(key);
      const parents = parentsOf(index, id);
      const count = counts.value?.get(id);
      const isFavorite = props.favorites.ids.value.has(id);
      return [
        <div key={key} class={['c2-runs-tree-row', selectedKey.value === key ? 'c2-runs-selected' : '']}
          style={{paddingLeft: `${8 + (path.length - 1) * 14}px`}}
          onClick={() => emit('select', {myRuns: false, path})}>
          <span class='c2-runs-tree-arrow' onClick={(e: MouseEvent) => {
            e.stopPropagation();
            if (canToggle)
              toggle(key);
          }}>{canToggle ? (isOpen ? '▾' : '▸') : ''}</span>
          <IconFA name='folder' class='c2-runs-folder-icon'/>
          <span class='c2-runs-tree-name'>{placeName(node.space)}</span>
          { node.error && <IconFA name='exclamation-triangle' class='c2-runs-error-icon'
            tooltip={spaceListingError(node.error)}/> }
          { parents.length > 1 &&
            <IconFA name='link' class='c2-runs-muted' tooltip={`Also in ${parents
              .filter((parent) => parent.id !== parentPath[parentPath.length - 1]).map(placeName).join(', ')}`}/> }
          { count !== undefined &&
            <span class='c2-runs-muted' title='Runs in this space and its subspaces'>({count})</span> }
          <IconFA name='star' faStyle={isFavorite ? 'fas' : 'fal'}
            class={['c2-runs-row-icon', isFavorite ? 'c2-runs-star c2-runs-star-on' : 'c2-runs-star']}
            tooltip={isFavorite ? 'Remove from favorites' : 'Add to favorites'}
            onClick={(e: MouseEvent) => {e.stopPropagation(); props.favorites.toggle(node.space);}}/>
        </div>,
        ...(isOpen ? children.flatMap((child) => renderSpace(index, child, path)) : []),
      ];
    };

    const renderFavoritesNote = () => {
      if (props.favorites.error.value) {
        return <div class='c2-runs-error c2-runs-pad'>Couldn't load favorites: {props.favorites.error.value}
          <span class='c2-runs-link' onClick={() => props.favorites.load()}> Try again</span></div>;
      }
      return kept.value?.size === 0 ? <div class='c2-runs-muted c2-runs-pad'>No favorite spaces</div> : undefined;
    };

    return () => (
      <div class='c2-runs-tree'>
        { props.showMyRuns &&
          <div class={['c2-runs-tree-row', props.location.myRuns ? 'c2-runs-selected' : '']}
            style={{paddingLeft: '8px'}} onClick={() => emit('select', {myRuns: true})}>
            <span class='c2-runs-tree-arrow'/>
            <IconFA name='user' class='c2-runs-folder-icon'/>
            <span class='c2-runs-tree-name'>My runs</span>
            { props.myRunsCount !== undefined &&
              <span class='c2-runs-muted' title='Runs in My runs'>({props.myRunsCount})</span> }
          </div> }
        { props.error &&
          <div class='c2-runs-error c2-runs-pad'>Couldn't load spaces: {props.error}
            <span class='c2-runs-link' onClick={() => emit('retry')}> Try again</span></div> }
        { !props.index && !props.error && <div class='c2-runs-muted c2-runs-pad'>Loading spaces…</div> }
        <div class='c2-runs-label c2-runs-tree-root' onClick={() => emit('select', {myRuns: false, path: []})}>
          <span class='c2-runs-tree-name'>Spaces</span>
          <IconFA name='star' faStyle={onlyFavorites.value ? 'fas' : 'fal'}
            class={['c2-runs-row-icon c2-runs-star-filter', onlyFavorites.value ? 'c2-runs-star-on' : '']}
            tooltip={onlyFavorites.value ? 'Show all spaces' : 'Show only favorite spaces'}
            onClick={(e: MouseEvent) => {e.stopPropagation(); onlyFavorites.value = !onlyFavorites.value;}}/>
        </div>
        { props.index && props.index.top.filter((id) => !kept.value || kept.value.has(id))
          .flatMap((id) => renderSpace(props.index!, id, [])) }
        { onlyFavorites.value && renderFavoritesNote() }
        { !onlyFavorites.value && props.index && props.index.top.length === 0 &&
          <div class='c2-runs-muted c2-runs-pad'>No spaces available</div> }
      </div>
    );
  },
});

