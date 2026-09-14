/* `SharedSession` (2-4a) over two memory-backed sources: the aggregate state, one `saveAll` for a
   draft parent and a draft child referencing it, a refusal that leaves both pending, discard
   fanning out, the ambient session `runWith` hands every source built inside it, and the
   `confirmDiscard` gate through the dialog DOM. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {backends} from '../src/sources/backends.js';
import {DomainSource} from '../src/sources/domain-source.js';
import {SharedSession, confirmDiscard} from '../src/sources/session.js';
import {Rows} from '../src/sources/rows-like.js';
import {notify} from '../src/components/display/notify.js';
import {backend} from './domain-fixtures.mjs';

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const env = {designTime: false, subBinds: {}, resolve: () => null};

async function pair(session) {
  const projects = new DomainSource({table: 'grit.project', session}, env);
  const issues = new DomainSource({table: 'grit.issue', session}, env);
  projects.start();
  issues.start();
  await flush();
  return {projects, issues};
}

const balloon = () => document.body.querySelector('.u2-notify-info')?.textContent ?? '';
const buttonNamed = (text) => [...document.body.querySelectorAll('.u2-dialog button')].find((b) => b.textContent === text);

scoped('aggregate: dirty, count, validity and the summary across sources; a source leaves at disposal', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  assert.deepEqual(session.sources.value, [projects, issues]);
  assert.equal(session.isDirty.value, false);
  assert.equal(session.summary.value, '');
  projects.rows.byKey('p1').name = 'Grit 2';
  assert.equal(session.isDirty.value, true);
  assert.equal(session.changeCount.value, 1);
  assert.equal(session.summary.value, '1 unsaved change');
  issues.rows.byKey('i1').title = 'A';
  issues.rows.byKey('i2').title = 'B';
  assert.equal(session.changeCount.value, 3);
  assert.equal(session.summary.value, '3 unsaved changes in 2 tables');
  assert.equal(session.validity.value, null);
  issues.rows.byKey('i2').title = '';
  assert.equal(session.validity.value, 'Value can\'t be empty', 'the first blocking problem');
  issues.dispose();
  assert.deepEqual(session.sources.value, [projects]);
  assert.equal(session.changeCount.value, 1);
  projects.dispose();
  assert.deepEqual(session.sources.value, []);
});

scoped('one saveAll: a draft parent and a draft child referencing it land as one transaction', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  const heard = [];
  const sub = session.onSaved.subscribe(() => heard.push('saved'));
  const parent = projects.newRow({key: 'NEW', name: 'New project'});
  const child = issues.newRow({project_id: parent.id, title: 'First issue'});
  assert.equal(Rows.isDraft(child.project_id), true, 'the child holds the parent\'s draft id');
  issues.currentRow.value = child;
  const saving = session.save();
  assert.equal(session.isSaving.value, true);
  assert.equal(await saving, true);
  sub.unsubscribe();
  assert.deepEqual(heard, ['saved']);
  assert.equal(session.isDirty.value, false);
  assert.equal(projects.isDirty.value, false);
  assert.equal(issues.isDirty.value, false);
  const savedParent = projects.currentRow.value;
  assert.match(savedParent.id, /^[0-9a-f-]{36}$/, 'the parent has its real id');
  assert.equal(savedParent.name, 'New project');
  const savedChild = issues.currentRow.value;
  assert.match(savedChild.id, /^[0-9a-f-]{36}$/, 'the current draft is re-pointed to the row it became');
  assert.equal(savedChild.project_id, savedParent.id, 'the child\'s FK is the parent\'s real id');
  assert.equal(projects.total.value, 3, 'totals counted again');
  assert.equal(issues.total.value, 4);
  assert.match(balloon(), /2 changes saved in 2 tables/);
  const [p, i] = ['grit.project', 'grit.issue'].map((t) => backends.domain.tableSync(t).history.at(-1));
  assert.equal(p.tx_id, i.tx_id, 'one tx_id in the audit');
  notify.closeAll();
  issues.rows.byKey('i1').title = 'Aspirin 100';
  assert.equal(await session.save(), true);
  assert.match(balloon(), /Issue saved/, 'one source, one change: the row');
  assert.equal(await session.save(), false, 'nothing pending: nothing landed');
  projects.dispose();
  issues.dispose();
});

scoped('a pristine parent draft a child refers to is saved with it: one saveAll, both inserted, the FK real', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  let calls = 0;
  const real = backends.domain.saveAll.bind(backends.domain);
  backends.domain.saveAll = (edits) => {
    calls++;
    return real(edits);
  };
  const parent = projects.newRow({key: 'NEW', name: 'Pristine'}, {pristine: true});
  assert.equal(projects.isDirty.value, false, 'created with defaults only, never edited');
  const child = issues.newRow({project_id: parent.id, title: 'Child'});
  assert.equal(session.changeCount.value, 1, 'the pristine parent counts for nothing');
  assert.equal(session.summary.value, '1 unsaved change');
  assert.equal(await session.save(), true);
  assert.equal(calls, 1, 'one saveAll for both');
  const savedParent = projects.currentRow.value;
  assert.match(savedParent.id, /^[0-9a-f-]{36}$/, 'the parent was inserted');
  assert.equal(savedParent.name, 'Pristine');
  assert.equal(issues.rows.byKey(child.id), undefined, 'the child is re-keyed');
  const savedChild = backends.domain.tableSync('grit.issue').rows.find((r) => r.title === 'Child');
  assert.equal(savedChild.project_id, savedParent.id, 'the child\'s FK is the parent\'s real id');
  assert.equal(backends.domain.tableSync('grit.project').rows.length, 3);
  assert.equal(session.isDirty.value, false);
  assert.match(balloon(), /1 change saved in 2 tables/);
  projects.dispose();
  issues.dispose();
});

scoped('a parent draft the frame filter hides is still saved with the child referring to it', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  const parent = projects.newRow({key: 'NEW', name: 'Hidden'}, {pristine: true});
  const at = projects.rows.items.value.findIndex((r) => r.id === parent.id);
  const df = projects.df.value;
  df.filter = {get: (i) => i !== at};
  df.onFilterChanged.fire(undefined);
  assert.equal(projects.rows.byKey(parent.id), undefined, 'the draft is out of sight');
  issues.newRow({project_id: parent.id, title: 'Child of a hidden parent'});
  assert.equal(await session.save(), true, 'the reference is discovered over the pending rows, not the visible ones');
  const savedChild = backends.domain.tableSync('grit.issue').rows.find((r) => r.title === 'Child of a hidden parent');
  const savedParent = backends.domain.tableSync('grit.project').rows.find((r) => r.name === 'Hidden');
  assert.equal(savedChild.project_id, savedParent.id, 'the FK is the parent\'s real id');
  projects.dispose();
  issues.dispose();
});

scoped('confirmDiscard refuses through the write-back interval, when the changes already read clean', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  let release;
  const real = backends.domain.saveAll.bind(backends.domain);
  backends.domain.saveAll = async (edits) => {
    const landed = await real(edits);
    await new Promise((resolve) => release = resolve);
    return landed;
  };
  issues.rows.byKey('i1').title = 'Aspirin 100';
  const saving = session.save();
  await flush();
  assert.equal(session.isDirty.value, false, 'the batch landed: the rows read clean already');
  assert.equal(session.isSaving.value, true, 'and the save is still in flight');
  assert.equal(await confirmDiscard(session), false, 'no navigation through that window');
  assert.equal(document.body.querySelector('.u2-dialog'), null, 'and nothing to answer');
  assert.match(document.body.querySelector('.u2-notify-warning')?.textContent ?? '', /Wait for the batch/);
  release();
  assert.equal(await saving, true);
  assert.equal(await confirmDiscard(session), true, 'once it is done, a clean session passes');
  projects.dispose();
  issues.dispose();
});

scoped('a refusal in the child leaves both pending with the error on both; a guard names its field', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  projects.rows.byKey('p1').name = 'Grit 2';
  issues.rows.byKey('i1').title = 'Aspirin 100';
  backends.domain.tableSync('grit.issue').rows[0].version = 9;
  assert.equal(await session.save(), false);
  assert.equal(session.isDirty.value, true, 'nothing landed');
  assert.equal(projects.error.value.code, 'version-conflict', 'the backend\'s refusal reaches every dirty source');
  assert.equal(issues.error.value.code, 'version-conflict');
  assert.match(projects.summary.value, /expected 1/);
  assert.equal(backends.domain.tableSync('grit.project').rows[0].name, 'Grit', 'the parent did not land either');
  backends.domain.tableSync('grit.issue').rows[0].version = 1;
  issues.newRow({project_id: 'p1', title: 'x', priority: 'medium'});
  assert.equal(await session.save(), false, 'the writer\'s own validity refuses before anything is sent');
  assert.match(issues.error.value.message, /Must be one of: low, high/);
  assert.equal(projects.error.value, undefined, 'and names only its own source');
  issues.discard();
  assert.equal(issues.isDirty.value, false, 'a source\'s own discard is the session\'s');
  assert.equal(projects.isDirty.value, false);

  projects.rows.byKey('p1').name = 'Grit 3';
  issues.rows.byKey('i1').title = 'A';
  const unguard = issues.guard(() => 'Title is required');
  assert.equal(await session.save(), false);
  assert.equal(issues.error.value.message, 'Cannot save: Title is required');
  assert.equal(projects.error.value, undefined, 'a guard is its own source\'s problem');
  assert.equal(session.isDirty.value, true);
  unguard();
  assert.equal(await session.save(), true);
  projects.dispose();
  issues.dispose();
});

scoped('discard fans out to every source and announces once', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  const heard = [];
  const sub = session.onDiscarded.subscribe(() => heard.push('discarded'));
  projects.rows.byKey('p1').name = 'x';
  issues.newRow({title: 'draft'});
  session.discard();
  assert.deepEqual(heard, ['discarded']);
  assert.equal(projects.rows.byKey('p1').name, 'Grit');
  assert.equal(issues.rows.items.value.length, 3);
  assert.equal(session.isDirty.value, false);
  sub.unsubscribe();
  projects.dispose();
  issues.dispose();
});

scoped('ambient: every source built under runWith joins that session; outside, each has its own', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const inside = SharedSession.runWith(session, () => [
    new DomainSource({table: 'grit.project'}, env), new DomainSource({table: 'grit.issue'}, env)]);
  assert.equal(SharedSession.ambient, undefined, 'restored');
  assert.deepEqual(session.sources.value, inside);
  assert.equal(inside[0].session, session);
  const alone = new DomainSource({table: 'grit.issue'}, env);
  assert.notEqual(alone.session, session);
  assert.ok(alone.session instanceof SharedSession);
  const explicit = SharedSession.runWith(session, () => new DomainSource({table: 'grit.issue', session: alone.session}, env));
  assert.equal(explicit.session, alone.session, 'an explicit session wins over the ambient one');
  for (const s of [...inside, alone, explicit])
    s.dispose();
});

scoped('confirmDiscard: clean passes silently; SAVE, DISCARD and CANCEL through the dialog; saving refuses', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  assert.equal(await confirmDiscard(session), true);
  assert.equal(document.body.querySelector('.u2-dialog'), null, 'no dialog for a clean session');

  projects.rows.byKey('p1').name = 'Grit 2';
  issues.rows.byKey('i1').title = 'A';
  let answer = confirmDiscard(session, {action: 'leave this page', subject: 'the issue'});
  await flush();
  const dialog = document.body.querySelector('.u2-dialog');
  assert.equal(dialog.querySelector('.u2-dialog-title-text').textContent, 'Unsaved changes');
  assert.deepEqual([...dialog.querySelectorAll('p')].map((p) => p.textContent),
    ['2 unsaved changes in the issue.', 'Save them, discard them, or cancel and do not leave this page.']);
  assert.deepEqual([...dialog.querySelectorAll('.u2-dialog-buttons button')].map((b) => b.textContent),
    ['CANCEL', 'SAVE', 'DISCARD'], 'the u2 dialog keeps CANCEL first');
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await answer, false);
  assert.equal(document.body.querySelector('.u2-dialog'), null);
  assert.equal(session.isDirty.value, true, 'cancel keeps everything');

  answer = confirmDiscard(session);
  await flush();
  assert.match(document.body.querySelector('.u2-dialog p').textContent, /in this view\./);
  fire(buttonNamed('SAVE'), 'click');
  assert.equal(await answer, true);
  assert.equal(session.isDirty.value, false, 'saved');
  assert.equal(backends.domain.tableSync('grit.project').rows[0].name, 'Grit 2');

  issues.rows.byKey('i1').title = 'B';
  answer = confirmDiscard(session);
  await flush();
  fire(buttonNamed('DISCARD'), 'click');
  assert.equal(await answer, true);
  assert.equal(issues.rows.byKey('i1').title, 'A', 'discarded');

  issues.rows.byKey('i1').title = 'C';
  issues.guard(() => 'Nope');
  answer = confirmDiscard(session);
  await flush();
  fire(buttonNamed('SAVE'), 'click');
  assert.equal(await answer, false, 'a refused save does not let the caller proceed');
  assert.equal(session.isDirty.value, true);

  const stuck = {isDirty: {peek: () => true}, isSaving: {peek: () => true}, changeCount: {peek: () => 1}};
  assert.equal(await confirmDiscard(stuck), false);
  assert.equal(document.body.querySelector('.u2-dialog'), null, 'no dialog while a batch is in flight');
  assert.match(document.body.querySelector('.u2-notify-warning')?.textContent ?? '', /Wait for the batch/);
  projects.dispose();
  issues.dispose();
});

scoped('a backend that refuses without throwing still says so in the summary, and says it once', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  projects.rows.byKey('p1').name = 'Grit 2';
  issues.rows.byKey('i1').title = 'Aspirin 100';
  // the platform editor reports the server's error itself and answers false
  backends.domain.saveAll = async () => false;
  assert.equal(await session.save(), false);
  assert.equal(session.isDirty.value, true, 'nothing landed');
  assert.equal(projects.error.value.code, 'refused');
  assert.equal(projects.summary.value, 'Cannot save: the changes were refused',
    'the summary says the changes are stuck, not how many there are');
  assert.equal(issues.summary.value, 'Cannot save: the changes were refused');
  projects.rows.byKey('p1').name = 'Grit 3';
  await flush();
  assert.equal(projects.error.value, undefined, 'the next edit takes the refusal back');
  projects.dispose();
  issues.dispose();
});

scoped('confirmDiscard agrees in number: one change is "it", several are "them"', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const {projects, issues} = await pair(session);
  issues.rows.byKey('i1').title = 'A';
  const answer = confirmDiscard(session, {action: 'leave this page'});
  await flush();
  assert.deepEqual([...document.body.querySelectorAll('.u2-dialog p')].map((p) => p.textContent),
    ['1 unsaved change in this view.', 'Save it, discard it, or cancel and do not leave this page.']);
  fire(buttonNamed('CANCEL'), 'click');
  assert.equal(await answer, false);
  projects.dispose();
  issues.dispose();
});

scoped('a child collection under a draft parent: the save puts the assigned id in its query and defaults', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const projects = new DomainSource({table: 'grit.project', session}, env);
  projects.start();
  await flush();
  const parent = projects.newRow({key: 'NEW', name: 'New project'});
  const issues = new DomainSource({table: 'grit.issue', session, defaults: {project_id: parent.id},
    query: `project_id = "${parent.id}"`}, env);
  issues.start();
  await flush();
  assert.equal(issues.state.value, 'ready', 'a query naming a draft id is not sent');
  assert.equal(issues.rows.items.value.length, 0);
  const child = issues.newRow({title: 'First issue'});
  assert.equal(child.project_id, parent.id, 'the draft takes the parent\'s draft id from the defaults');
  assert.equal(await session.save(), true);
  const saved = projects.currentRow.value;
  assert.match(saved.id, /^[0-9a-f-]{36}$/);
  assert.equal(issues.query.value, `project_id = "${saved.id}"`, 'the query names the row the parent became');
  assert.equal(issues.defaults.project_id, saved.id, 'and so do the defaults');
  const savedChild = issues.currentRow.value;
  assert.notEqual(savedChild, null, 'the re-read keeps the row current: the query matches it now');
  assert.equal(savedChild.project_id, saved.id, 'the child\'s FK is the parent\'s real id');
  assert.equal(Rows.isDraft(savedChild), false);
  assert.equal(issues.rows.items.value.length, 1);
  assert.equal(issues.total.value, 1, 'counted under the query as rewritten');
  assert.equal(issues.newRow({title: 'Second issue'}).project_id, saved.id, 'a later draft takes the real FK');
  projects.dispose();
  issues.dispose();
});

scoped('a query naming a draft id is never sent: a backend refusing the literal is not asked for rows', async () => {
  const inner = backend();
  const frames = [];
  const counts = [];
  // the platform refuses a `~new:` literal on a uuid column ("… is not a valid id for column")
  const refusing = (filter) => {
    if (JSON.stringify(filter ?? '').includes(Rows.DRAFT_PREFIX))
      throw new Error(`"${filter}" is not a valid id for column "project_id"`);
  };
  backends.domain = {
    saveAll: (edits) => inner.saveAll(edits),
    table: async (address) => {
      const t = await inner.table(address);
      return Object.assign(Object.create(Object.getPrototypeOf(t)), t, {
        frame: (spec) => {
          refusing(spec.filter);
          frames.push(spec);
          return t.frame(spec);
        },
        count: (filter, search) => {
          refusing(filter);
          counts.push(filter);
          return t.count(filter, search);
        },
      });
    },
  };
  const session = new SharedSession();
  const projects = new DomainSource({table: 'grit.project', session}, env);
  projects.start();
  await flush();
  const parent = projects.newRow({key: 'NEW', name: 'New project'});
  const issues = new DomainSource({table: 'grit.issue', session, defaults: {project_id: parent.id},
    query: `project_id = "${parent.id}"`}, env);
  issues.start();
  await flush();
  assert.equal(issues.state.value, 'ready', 'no query, no error');
  assert.equal(issues.error.value, undefined);
  assert.deepEqual(frames.at(-1), {limit: 0, offset: 0, withAccess: true},
    'the child asked for no rows, under no filter');
  assert.equal(counts.length, 1, 'and its total was not counted either');
  issues.newRow({title: 'First issue'});
  assert.equal(await session.save(), true);
  assert.equal(issues.state.value, 'ready');
  assert.equal(issues.rows.items.value.length, 1, 'the rows arrive once the query names the real id');
  assert.equal(issues.currentRow.value.project_id, projects.currentRow.value.id);
  projects.dispose();
  issues.dispose();
});

scoped('a source outside the batch: its defaults follow the assigned id, and nothing it holds is re-read away', async () => {
  backends.domain = backend();
  const session = new SharedSession();
  const projects = new DomainSource({table: 'grit.project', session}, env);
  projects.start();
  await flush();
  const parent = projects.newRow({key: 'NEW', name: 'New project'});
  // the children pane's shape under a draft parent: no query, the parent's id as the default
  const pane = new DomainSource({table: 'grit.issue', session, empty: true, defaults: {project_id: parent.id}}, env);
  pane.start();
  await flush();
  const pristine = pane.newRow({title: 'Later'}, {pristine: true});
  assert.equal(await session.save(), true);
  const saved = projects.currentRow.value;
  assert.equal(pane.defaults.project_id, saved.id, 'the defaults name the row the parent became');
  assert.notEqual(pane.rows.byKey(pristine.id), undefined, 'a defaults-only rewrite re-reads nothing');
  assert.equal(pane.newRow({title: 'Next'}).project_id, saved.id);
  projects.dispose();
  pane.dispose();
});
