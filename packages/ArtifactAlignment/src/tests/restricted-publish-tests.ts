import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {T_ALIGNMENT, T_PROGRAM} from '../domain/constants';
import {approvePublication, rejectPublication, updateCuration} from '../service/publication-service';
import {
  addToTier, cleanupTestPrograms, makeSavedRun, makeTestProgram, publishRun, withRestrictedUser,
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
      await addToTier(user.group, program, 'contributor');
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
      await addToTier(user.group, program, 'contributor');
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
      await addToTier(user.group, program, 'contributor');
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

  test('viewer-tier cannot publish into their own program', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    await withRestrictedUser('aavpub', async (user) => {
      await addToTier(user.group, program, 'viewer');
      let error = '';
      try {
        await user.asUser(() => publishRun({
          sourceMetaCallId: run.metaCallId, programId: program.id, name: 'ViewerPublish'}));
      } catch (e: any) {
        error = `${e?.message ?? e}`;
      }
      expect(error.includes('Only members of the program contributors or approvers groups'), true,
        `got: '${error}'`);
      expect((await alignment().query({filter: DG.cond('program_id', '=', program.id)})).length, 0,
        'no version row must land');
    });
  });

  test('viewer-tier cannot curate', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'ViewerCurate'});
    await withRestrictedUser('aavcur', async (user) => {
      await addToTier(user.group, program, 'viewer');
      let error = '';
      try {
        await user.asUser(() => updateCuration(pub.rowId, {path: 'nope'}));
      } catch (e: any) {
        error = `${e?.message ?? e}`;
      }
      expect(error.includes('Only members of the program contributors or approvers groups'), true,
        `got: '${error}'`);
      const row: any = await alignment().get(pub.rowId);
      expect(row.path == null, true, 'the denied curation must not stick');
    });
  });

  test('contributor-tier cannot approve', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrApprove',
      skipAutoApprove: true});
    await withRestrictedUser('aacappr', async (user) => {
      await addToTier(user.group, program, 'contributor');
      let error = '';
      try {
        await user.asUser(() => approvePublication(pub.rowId));
      } catch (e: any) {
        error = `${e?.message ?? e}`;
      }
      expect(error.includes('Only members of the program approvers group'), true, `got: '${error}'`);
      const row: any = await alignment().get(pub.rowId);
      expect(row.status, 'pending', 'the denied approval must not stick');
    });
  });

  test('viewer- and contributor-tier cannot reject', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'TierReject',
      skipAutoApprove: true});
    for (const role of ['viewer', 'contributor'] as const) {
      await withRestrictedUser(`aarj${role[0]}`, async (user) => {
        await addToTier(user.group, program, role);
        let error = '';
        try {
          await user.asUser(() => rejectPublication(pub.rowId, 'not allowed'));
        } catch (e: any) {
          error = `${e?.message ?? e}`;
        }
        expect(error.includes('Only members of the program approvers group'), true,
          `${role}: got '${error}'`);
      });
    }
    const row: any = await alignment().get(pub.rowId);
    expect(row.status, 'pending');
  });

  test('tiered self-approve — the author-approver is denied, a second approver succeeds', async () => {
    const program = await makeTestProgram();
    await withRestrictedUser('aaself', async (author) => {
      await addToTier(author.group, program, 'approver');
      // the author must own the source run for the author != approver rule to bite
      const run = await author.asUser(() => makeSavedRun());
      const pub = await author.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'TierSelf',
        skipAutoApprove: true}));
      let error = '';
      try {
        await author.asUser(() => approvePublication(pub.rowId));
      } catch (e: any) {
        error = `${e?.message ?? e}`;
      }
      expect(error.includes('Authors cannot approve'), true, `got: '${error}'`);
      await withRestrictedUser('aasecond', async (approver) => {
        await addToTier(approver.group, program, 'approver');
        await approver.asUser(() => approvePublication(pub.rowId));
        const row: any = await alignment().get(pub.rowId);
        expect(row.status, 'approved');
        expect(row.approved_by, approver.id);
      });
    });
  });

  test('viewer-tier raw insert into alignment is denied', async () => {
    const program = await makeTestProgram();
    await withRestrictedUser('aavins', async (user) => {
      await addToTier(user.group, program, 'viewer');
      let error: any = null;
      try {
        await user.asUser(() => alignment().insert({
          publication_id: crypto.randomUUID(), revision: 1, name: 'ViewerInsert',
          artifact_id: crypto.randomUUID(), program_id: program.id}));
      } catch (e) {
        error = e;
      }
      expect(error instanceof DG.DomainError, true, `expected a typed denial, got: ${error}`);
      expect((await alignment().query({filter: DG.cond('program_id', '=', program.id)})).length, 0,
        'the denied insert must not stick');
    });
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 240000});
