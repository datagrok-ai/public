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

test('saved searches reach the server at once: save, remove and clear are flushed, not left to the 10 s sync', () => {
  userSettings.values = {};
  const {saved, restore} = server();
  try {
    const high = {q: 'priority = "high"', search: '', sort: ''};
    SavedSearches.save('grit.issue', 'High', high);
    assert.equal(saved.stored['ems-searches|grit.issue']?.High, JSON.stringify(high),
      'a saved search was not flushed: a reload within 10 s loses it');

    SavedSearches.remove('grit.issue', 'High');
    assert.equal(saved.stored['ems-searches|grit.issue']?.High, undefined,
      'a removed saved search was not flushed: a reload within 10 s brings it back');

    SavedSearches.record('grit.issue', high);
    SavedSearches.save('grit.issue', 'Any', high);
    SavedSearches.clearRecent('grit.issue');
    assert.equal(saved.stored['ems-recent-searches']?.['grit.issue'], undefined,
      'cleared recent searches were not flushed: a reload within 10 s brings them back');
  } finally {
    restore();
  }
});
