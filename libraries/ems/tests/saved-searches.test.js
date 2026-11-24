import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {SavedSearches} = await import('../src/dg/saved-searches.js');
const {userSettings} = await import('datagrok-api/grok');

/** What reached the server: grok.userSettings holds writes in memory until flush() or its 10 s sync. */
function server() {
  const saved = {stored: {}};
  const flush = userSettings.flush;
  userSettings.flush = async () => {
    saved.stored = structuredClone(userSettings.values);
  };
  return {saved, restore: () => userSettings.flush = flush};
}

test('saved searches reach the server at once: save, remove and clear are flushed, not left to the 10 s sync', async () => {
  userSettings.values = {};
  const {saved, restore} = server();
  try {
    const high = {q: 'priority = "high"', search: '', sort: ''};
    await SavedSearches.save('grit.issue', 'High', high);
    assert.equal(saved.stored['ems-searches|grit.issue']?.High, JSON.stringify(high),
      'a saved search was not flushed: a reload within 10 s loses it');

    await SavedSearches.remove('grit.issue', 'High');
    assert.equal(saved.stored['ems-searches|grit.issue']?.High, undefined,
      'a removed saved search was not flushed: a reload within 10 s brings it back');

    SavedSearches.record('grit.issue', high);
    await SavedSearches.save('grit.issue', 'Any', high);
    SavedSearches.clearRecent('grit.issue');
    assert.equal(saved.stored['ems-recent-searches']?.['grit.issue'], undefined,
      'cleared recent searches were not flushed: a reload within 10 s brings them back');
  } finally {
    restore();
  }
});

test('a search saved in this tab keeps the ones another tab saved since this one loaded', async () => {
  userSettings.values = {};
  const {saved, restore} = server();
  const {dapi} = await import('datagrok-api/grok');
  const get = dapi.userDataStorage.get;
  dapi.userDataStorage.get = async () => ({Alpha: '{"q":"","search":"alpha","sort":""}'});
  try {
    await SavedSearches.save('grit.issue', 'Beta', {q: '', search: 'beta', sort: ''});
    assert.deepEqual(Object.keys(saved.stored['ems-searches|grit.issue']), ['Alpha', 'Beta']);
    await SavedSearches.remove('grit.issue', 'Beta');
    assert.deepEqual(Object.keys(saved.stored['ems-searches|grit.issue']), ['Alpha']);
  } finally {
    dapi.userDataStorage.get = get;
    restore();
  }
});
