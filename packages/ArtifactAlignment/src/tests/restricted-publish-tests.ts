import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {T_ALIGNMENT, T_PROGRAM} from '../domain/constants';
import {updateCuration} from '../service/publication-service';
import {
  addToGroup, cleanupTestPrograms, makeSavedRun, makeTestProgram, publishRun, withRestrictedUser,
} from './fixtures';

const alignment = () => grok.dapi.domains.table(T_ALIGNMENT);

category('ArtifactAlignment: restricted access', () => {
  test('a contributor sees the catalog and publishes pending; outsiders stay blind', async () => {
    const program = await makeTestProgram();
    const foreign = await makeTestProgram();
    const run = await makeSavedRun();
    await withRestrictedUser('aares', async (user) => {
      const programsSeen = (id: string) => user.asUser(async () =>
        (await grok.dapi.domains.table(T_PROGRAM).query(
          {filter: DG.cond('id', '=', id), columns: ['code']})).length);
      expect(await programsSeen(program.id), 0,
        'a user outside every artifact group must not see programs');
      await addToGroup(user.group, program.contributors);
      // Visibility is umbrella-wide (a master-row-scoped predicate is a core
      // ask): a program contributor sees the registry, foreign programs
      // included; per-program write isolation is the service guard, below.
      expect(await programsSeen(program.id), 1, 'a contributor must see their program');
      expect(await programsSeen(foreign.id), 1,
        'audience visibility is registry-wide until per-row trimming lands');
      // The auto-approve gate runs as the acting user: without approval rights
      // the column security rejects the flip and the version stays pending.
      const result = await user.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'RestrictedPublish'}));
      expect(result.status, 'pending', 'a non-approver publish must stay pending');
      const row: any = await alignment().get(result.rowId);
      expect(row.status, 'pending');
    });
  });

  test('publish service refuses a cross-program contributor', async () => {
    const program = await makeTestProgram();
    const foreign = await makeTestProgram();
    const run = await makeSavedRun();
    await withRestrictedUser('aacross', async (user) => {
      await addToGroup(user.group, program.contributors);
      let error = '';
      try {
        await user.asUser(() => publishRun({
          sourceMetaCallId: run.metaCallId, programId: foreign.id, name: 'CrossProgram'}));
      } catch (e: any) {
        error = `${e?.message ?? e}`;
      }
      expect(error.includes('Only members of the program contributors or approvers groups'), true,
        `got: '${error}'`);
      expect((await alignment().query({filter: DG.cond('program_id', '=', foreign.id)})).length, 0,
        'no version row must land in the foreign program');
    });
  });

  test('curation service refuses a cross-program contributor', async () => {
    const program = await makeTestProgram();
    const foreign = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: foreign.id, name: 'ForeignRow',
      skipAutoApprove: true});
    await withRestrictedUser('aacur', async (user) => {
      await addToGroup(user.group, program.contributors);
      let error = '';
      try {
        await user.asUser(() => updateCuration(pub.rowId, {path: 'hijacked'}));
      } catch (e: any) {
        error = `${e?.message ?? e}`;
      }
      expect(error.includes('Only members of the program contributors or approvers groups'), true,
        `got: '${error}'`);
      const row: any = await alignment().get(pub.rowId);
      expect(row.path == null, true, 'the denied curation must not stick');
    });
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 240000});
