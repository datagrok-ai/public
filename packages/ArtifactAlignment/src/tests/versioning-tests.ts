import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {T_ALIGNMENT, T_HISTORY} from '../domain/constants';
import {
  approvePublication, discover, getLiveVersions, getPublicationHistory,
  publishWorkflowRun, rejectPublication,
} from '../service/publication-service';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram, makeTestStudy} from './fixtures';

category('ArtifactAlignment: versioning', () => {
  const alignment = () => grok.dapi.domains.table(T_ALIGNMENT);
  const history = () => grok.dapi.domains.table(T_HISTORY);

  test('first publish auto-approves and is discoverable', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const result = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'First',
      workstream: 'modeling',
    });
    expect(result.revision, 1);
    expect(result.status, 'approved');
    expect(result.artifactId === run.metaCallId, false, 'the frozen copy must be a new meta call');
    const discovered = await discover(program.id);
    expect(discovered.length, 1);
    expect(discovered[0].artifact_id, result.artifactId);
    expect(discovered[0].revision, 1);
  });

  test('republish under the same key bumps the revision and archives the old version on approval', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Repub'});
    const v2 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Repub'});
    expect(v1.publicationId, v2.publicationId, 'same key must address the same publication');
    expect(v2.revision, 2);
    const discovered = await discover(program.id);
    expect(discovered.length, 1, 'discovery must hold exactly one row per publication');
    expect(discovered[0].revision, 2);
    const archived = await getPublicationHistory(v1.publicationId);
    expect(archived.length, 1);
    expect(archived[0].revision, 1);
    expect(archived[0].final_status, 'approved');
    expect(archived[0].original_row_id, v1.rowId);
  });

  test('latest-approved-stays-live — a pending update never affects the audience', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Stays'});
    const v2 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Stays', skipAutoApprove: true});
    expect(v2.status, 'pending');
    const live = await getLiveVersions(v1.publicationId);
    expect(live.length, 2, 'approved + in-review rows coexist');
    const discovered = await discover(program.id);
    expect(discovered.length, 1);
    expect(discovered[0].revision, 1, 'discovery still serves v1 while v2 is in review');
    await approvePublication(v2.rowId, {auto: true});
    const after2 = await discover(program.id);
    expect(after2.length, 1);
    expect(after2[0].revision, 2);
    expect((await getPublicationHistory(v1.publicationId)).some((h) => h.revision === 1), true);
  });

  test('rejection leaves the approved version untouched; resubmission archives the rejected draft', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Rej'});
    const v2 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Rej', skipAutoApprove: true});
    await rejectPublication(v2.rowId, 'wrong cohort', {auto: true});
    const discovered = await discover(program.id);
    expect(discovered[0].revision, 1, 'audience keeps seeing v1 after the rejection');
    const v3 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Rej', skipAutoApprove: true});
    expect(v3.revision, 3);
    const archived = await getPublicationHistory(v1.publicationId);
    const rejected = archived.find((h) => h.revision === 2);
    expect(rejected?.final_status, 'rejected');
    expect(rejected?.reject_reason, 'wrong cohort');
    expect((await getLiveVersions(v1.publicationId)).length, 2, 'approved v1 + pending v3');
  });

  test('widened key — the same name under a different study is a different publication', async () => {
    const program = await makeTestProgram();
    const study = await makeTestStudy(program.id);
    const run = await makeSavedRun();
    const noStudy = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Keyed'});
    const withStudy = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, studyId: study.id, name: 'Keyed'});
    expect(noStudy.publicationId === withStudy.publicationId, false);
    expect(withStudy.revision, 1);
    expect((await discover(program.id)).length, 2);
  });

  test('server-enforced singleton — a second pending row of the same publication is rejected', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Singleton', skipAutoApprove: true});
    let error: any = null;
    try {
      await alignment().insert({
        publication_id: v1.publicationId, revision: 99, name: 'Singleton',
        artifact_id: crypto.randomUUID(), program_id: program.id,
      }, {errorOnDuplicate: true});
    } catch (e) {
      error = e;
    }
    expect(error instanceof DG.DomainValidationError, true, `unexpected: ${error}`);
  });

  test('approve flip is blocked while another approved row is live (direct write)', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Flip'});
    const v2 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Flip', skipAutoApprove: true});
    let error: any = null;
    try {
      await alignment().update(v2.rowId, {status: 'approved'});
    } catch (e) {
      error = e;
    }
    expect(error instanceof DG.DomainValidationError, true,
      'the raw flip must violate the (publication_id, status) live unique key');
    await approvePublication(v2.rowId, {auto: true});
    expect((await discover(program.id))[0].revision, 2);
  });

  test('history archive insert is idempotent on (publication_id, revision)', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Idem'});
    await publishWorkflowRun({sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Idem'});
    const archived = (await getPublicationHistory(v1.publicationId))[0];
    const [dup] = await history().insert({
      publication_id: archived.publication_id, revision: archived.revision,
      name: archived.name, artifact_id: archived.artifact_id,
      program_id: program.id, final_status: archived.final_status,
    });
    expect(dup.created, false);
    expect(dup.status, 'duplicate');
  });

  test('publication with only archived versions resumes identity and numbering', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Resume'});
    // Simulate an interrupted flip: archive the only live row by hand.
    const row: any = (await getLiveVersions(v1.publicationId))[0];
    await history().insert({
      publication_id: row.publication_id, revision: row.revision, name: row.name,
      artifact_id: row.artifact_id, program_id: row.program_id,
      final_status: row.status, original_row_id: row.id,
    });
    await alignment().delete(row.id);
    const v2 = await publishWorkflowRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Resume'});
    expect(v2.publicationId, v1.publicationId, 'identity resumes from history');
    expect(v2.revision, 2);
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
