import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {T_ALIGNMENT, T_STUDY} from '../domain/constants';
import {approvePublication} from '../service/publication-service';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram, makeTestStudy, publishRun} from './fixtures';

// Pins the platform mechanics the design doc (docs/artifact-alignment.html) marks as
// verified, so a platform upgrade that breaks one fails loudly here.
category('ArtifactAlignment: platform mechanics', () => {
  const alignment = () => grok.dapi.domains.table(T_ALIGNMENT);

  test('insert without status lands on the server default pending', async () => {
    const program = await makeTestProgram();
    const [ins] = await alignment().insert({
      publication_id: crypto.randomUUID(), revision: 1, name: 'Default status',
      artifact_id: crypto.randomUUID(), program_id: program.id,
    });
    const row: any = await alignment().get(ins.id);
    expect(row.status, 'pending', 'new versions must start unapproved structurally');
  });

  test('choices are validated server-side with a typed error', async () => {
    const program = await makeTestProgram();
    let error: any = null;
    try {
      await alignment().insert({
        publication_id: crypto.randomUUID(), revision: 1, name: 'Bad choice',
        artifact_id: crypto.randomUUID(), program_id: program.id, workstream: 'bogus',
      });
    } catch (e) {
      error = e;
    }
    expect(error instanceof DG.DomainValidationError, true, `unexpected: ${error}`);
  });

  test('onDelete restrict rejects deleting a study referenced by a live version row', async () => {
    const program = await makeTestProgram();
    const study = await makeTestStudy(program.id);
    await alignment().insert({
      publication_id: crypto.randomUUID(), revision: 1, name: 'Restricts',
      artifact_id: crypto.randomUUID(), program_id: program.id, study_id: study.id,
    });
    let error: any = null;
    try {
      await grok.dapi.domains.table(T_STUDY).delete(study.id);
    } catch (e) {
      error = e;
    }
    expect(error instanceof DG.DomainRestrictError, true, `unexpected: ${error}`);
  });

  test('the approval flip leaves an audit row with before/after diffs and the actor', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const result = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Audited',
      skipAutoApprove: true});
    await approvePublication(result.rowId, {auto: true});
    const trail = await alignment().audit(result.rowId);
    const flip = trail.find((e) => e.op === 'update' &&
      e.before?.status === 'pending' && e.after?.status === 'approved');
    expect(flip != null, true, `no approval diff in: ${trail.map((e) => e.op).join(',')}`);
    const user = await grok.dapi.users.current();
    expect(flip!.actor_id, user.id, 'the approval audit row must carry the actor');
  });

  test('dotted reference paths filter through the program registry', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const result = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Dotted'});
    const rows = await alignment().query({filter: [
      {property: 'program_id.code', operator: '=', value: program.code},
      {property: 'status', operator: '=', value: 'approved'},
    ]});
    expect(rows.some((r: any) => r.id === result.rowId), true);
  });

  test('set-level molecule questions — compounds.id = null finds unlinked publications', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const linked = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Linked',
      compounds: [`GRK-SET-${Date.now()}`]});
    const unlinked = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Unlinked'});
    const rows = await alignment().query({filter: [
      {property: 'program_id', operator: '=', value: program.id},
      {property: 'compounds.id', operator: '=', value: null},
    ]});
    const ids = rows.map((r: any) => r.id);
    expect(ids.includes(unlinked.rowId), true, 'the compound-less publication must match');
    expect(ids.includes(linked.rowId), false, 'the linked publication must not match');
  });

  test('published artifact ids resolve through batch entity resolution', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const result = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Resolvable'});
    const [entity] = await grok.dapi.getEntities([result.artifactId]);
    expect(entity != null, true, 'the frozen meta call id must resolve to an entity');
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
