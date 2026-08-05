import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Exercises DG.Favorites.add/remove. Uses a Project entity for portability —
// the same API path is what scripts/dapi/favorites.js documents for FileInfo,
// which additionally needs FileInfo.save() to assign an entity id.
async function isBookmarked(entityId: string): Promise<boolean> {
  const list = await grok.dapi.entities.filter('starredBy = @current').list({pageSize: 200});
  return list.some((e) => e.id === entityId);
}

category('Dapi: favorites', () => {
  test('add and remove', async () => {
    const project = DG.Project.create();
    project.name = `apitests-fav-${DG.Utils.randomString(6)}`;
    let saved: _DG.Project | undefined;
    try {
      saved = await grok.dapi.projects.save(project);
      expect(saved.id != null, true, 'project has no id after save');

      await DG.Favorites.add(saved);
      expect(await isBookmarked(saved.id), true, 'entity not bookmarked after add');

      await DG.Favorites.remove(saved);
      expect(await isBookmarked(saved.id), false, 'entity still bookmarked after remove');
    } finally {
      try { if (saved) await DG.Favorites.remove(saved); } catch (_) {}
      try { if (saved) await grok.dapi.projects.delete(saved); } catch (_) {}
    }
  });
}, {owner: 'aparamonov@datagrok.ai'});
