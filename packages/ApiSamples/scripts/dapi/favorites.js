// Adds a FileInfo entity to the current user's favorites, then verifies and cleans up.

(async () => {
  const path = 'System:AppData/ApiSamples/test-favorites.txt';

  await grok.dapi.files.writeAsText(path, 'hello favorites');

  const files = await grok.dapi.files.list('System:AppData/ApiSamples', false, 'test-favorites');
  const file = files.find((f) => f.fileName === 'test-favorites.txt');
  if (!file)
    throw new Error('test file not found');

  // Files browsed from cloud storage have no entity id until saved as a FileInfo record.
  const saved = await file.save();
  await DG.Favorites.add(saved);

  const favs = await grok.dapi.entities.getFavorites();
  const found = favs.some((e) => e.id === saved.id);
  console.log(`favorited: ${found}`);

  await DG.Favorites.remove(saved);
  await grok.dapi.files.delete(path);
})();
