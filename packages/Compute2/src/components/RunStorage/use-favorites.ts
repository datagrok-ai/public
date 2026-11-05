import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export type FavoritesApi = {
  list: () => Promise<DG.Entity[]>,
  add: (entity: DG.Entity) => Promise<void>,
  remove: (entity: DG.Entity) => Promise<void>,
};

const platformFavorites: FavoritesApi = {
  list: () => grok.dapi.entities.getFavorites(),
  add: (entity) => DG.Favorites.add(entity),
  remove: (entity) => DG.Favorites.remove(entity),
};

/** The ids of the user's favorite entities; a failed load keeps the ids it had, with the error. */
export function createFavorites(api: FavoritesApi = platformFavorites) {
  const ids = Vue.shallowRef(new Set<string>());
  const error = Vue.ref<string>();

  // one request at a time, so a load under way can't undo a toggle and quick toggles all count
  let queue = Promise.resolve();
  const enqueue = (task: () => Promise<void>) => queue = queue.then(task);

  // every star looks empty until the first load finishes, so a click before that can't say what it wants
  let isLoaded = false;

  const load = () => enqueue(async () => {
    try {
      ids.value = new Set((await api.list()).map((item) => item.id));
      error.value = undefined;
    } catch (e: any) {
      error.value = String(e?.message ?? e);
    }
    isLoaded = true;
  });

  const toggle = (entity: DG.Entity) => !isLoaded ? queue : enqueue(async () => {
    const isFavorite = ids.value.has(entity.id);
    try {
      await (isFavorite ? api.remove(entity) : api.add(entity));
      const next = new Set(ids.value);
      isFavorite ? next.delete(entity.id) : next.add(entity.id);
      ids.value = next;
    } catch (e: any) {
      grok.shell.error(`Couldn't change the favorite: ${e?.message ?? e}`);
    }
  });

  return {ids: ids as Readonly<typeof ids>, error: error as Readonly<typeof error>, load, toggle};
}

export type Favorites = ReturnType<typeof createFavorites>;

const sessionFavorites = createFavorites();

/** The favorites every run browser and Save dialog of the session share, loaded again when one mounts. */
export function useFavorites() {
  Vue.onMounted(sessionFavorites.load);
  return sessionFavorites;
}
