import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, before} from '@datagrok-libraries/test/src/test';
import {until, untilAsync} from '../helpers';

// Favorites surface: DG.Favorites.add/remove, grok.shell.favorites, grok.dapi.entities.getFavorites.
// Server round-trips — reuse an existing script entity, verify via polled getFavorites(), un-favorite in finally.
category('AI: App: Favorites JS API', () => {
  // Fetch a real, already-persisted entity to (un)favorite. Returns null if the server has none.
  async function anyEntity(): Promise<DG.Entity | null> {
    const scripts = await grok.dapi.scripts.list({pageSize: 1});
    if (scripts != null && scripts.length > 0)
      return scripts[0];
    const queries = await grok.dapi.queries.list({pageSize: 1});
    return queries != null && queries.length > 0 ? queries[0] : null;
  }

  const hasId = (e: DG.Entity, id: string): boolean =>
    e != null && e.id === id;

  const favoritesContain = (favs: DG.Entity[], id: string): boolean =>
    Array.isArray(favs) && favs.some((f) => f != null && f.id === id);

  let sharedEntity: DG.Entity | null = null;
  before(async () => {
    sharedEntity = await anyEntity();
  });

  test('add an existing entity to favorites and verify it appears in getFavorites', async () => {
    const e = sharedEntity;
    expect(e != null, true);
    const id = e!.id;
    try {
      await DG.Favorites.add(e!);
      const found = await untilAsync(async () => favoritesContain(await grok.dapi.entities.getFavorites(), id));
      expect(found, true);
    } finally {
      await DG.Favorites.remove(e!);
    }
  });

  test('remove a favorited entity and verify it is absent', async () => {
    const e = sharedEntity;
    expect(e != null, true);
    const id = e!.id;
    let added = false;
    try {
      await DG.Favorites.add(e!);
      added = true;
      expect(await untilAsync(async () => favoritesContain(await grok.dapi.entities.getFavorites(), id)), true);
      await DG.Favorites.remove(e!);
      added = false;
      const gone = await untilAsync(async () => !favoritesContain(await grok.dapi.entities.getFavorites(), id));
      expect(gone, true);
    } finally {
      if (added)
        await DG.Favorites.remove(e!);
    }
  });

  test('getFavorites returns Entity array with non-null ids', async () => {
    const e = sharedEntity;
    expect(e != null, true);
    try {
      await DG.Favorites.add(e!);
      expect(await untilAsync(async () => favoritesContain(await grok.dapi.entities.getFavorites(), e!.id)), true);
      const favs = await grok.dapi.entities.getFavorites();
      expect(Array.isArray(favs), true);
      expect(favs.length > 0, true);
      for (const f of favs) {
        expect(f instanceof DG.Entity, true);
        expect(f.id != null && f.id !== '', true);
      }
    } finally {
      await DG.Favorites.remove(e!);
    }
  });

  test('adding the same entity twice does not create a duplicate', async () => {
    const e = sharedEntity;
    expect(e != null, true);
    const id = e!.id;
    try {
      await DG.Favorites.add(e!);
      await DG.Favorites.add(e!);
      expect(await untilAsync(async () => favoritesContain(await grok.dapi.entities.getFavorites(), id)), true);
      const favs = await grok.dapi.entities.getFavorites();
      const matches = favs.filter((f) => hasId(f, id)).length;
      expect(matches, 1);
    } finally {
      await DG.Favorites.remove(e!);
    }
  });

  test('grok.shell.favorites reflects a just-added entity', async () => {
    const e = sharedEntity;
    expect(e != null, true);
    const id = e!.id;
    try {
      await DG.Favorites.add(e!);
      const snapshot = grok.shell.favorites;
      expect(Array.isArray(snapshot), true);
      // The shell snapshot may lag the server write; poll it synchronously.
      await until(() => favoritesContain(grok.shell.favorites, id));
      expect(favoritesContain(grok.shell.favorites, id), true);
    } finally {
      await DG.Favorites.remove(e!);
    }
  });

  // Dropped: grok.dapi.entities.getFavoritesForGroups relies on api.grok_EntitiesDataSource_GetFavoritesForGroups,
  // which is not registered on this client (predates the binding) — the call throws "is not a function" at runtime.
}, {owner: 'agolovko@datagrok.ai'});
