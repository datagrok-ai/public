import * as grok from 'datagrok-api/grok';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {T_ALIGNMENT, T_COMPOUND, T_TAG} from '../domain/constants';
import {
  approvePublication, discover, getLiveVersions, getPublicationHistory, rejectPublication, updateCuration,
} from '../service/publication-service';
import {
  addToTier, cleanupTestPrograms, makeSavedRun, makeTestProgram, publishRun, withRestrictedUser,
} from './fixtures';

const alignment = () => grok.dapi.domains.table(T_ALIGNMENT);

// The positive workflows under real role tiers (viewer ⊂ contributor ⊂ approver) —
// every state flip runs under the acting user's own session, not the admin's.
category('ArtifactAlignment: multiuser versioning', () => {
  test('contributor republishes over an approved v1', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrRepub', path: 'PK/original'});
    expect(v1.status, 'approved');
    await withRestrictedUser('aactr1', async (probe) => {
      await addToTier(probe.group, program, 'contributor');
      const v2 = await probe.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrRepub'}));
      expect(v2.publicationId, v1.publicationId, 'same key must address the same publication');
      expect(v2.revision, 2);
      expect(v2.status, 'pending', 'a contributor republish must stay pending');
      const discovered = await discover(program.id);
      expect(discovered[0].revision, 1, 'the audience keeps seeing v1');
      const row: any = await alignment().get(v2.rowId);
      expect(row.path, 'PK/original', 'curation must copy forward under the contributor session');
    });
  });

  test('contributor republishes over their own pending draft', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    await withRestrictedUser('aactr2', async (probe) => {
      await addToTier(probe.group, program, 'contributor');
      const v1 = await probe.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrDraft'}));
      expect(v1.status, 'pending');
      const v2 = await probe.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrDraft'}));
      expect(v2.revision, 2);
      expect((await getLiveVersions(v1.publicationId)).length, 1, 'the draft must be superseded');
      const archived = await getPublicationHistory(v1.publicationId);
      expect(archived.length, 1);
      expect(archived[0].revision, 1);
      expect(archived[0].final_status, 'pending');
    });
  });

  test('contributor resubmits after an approver rejection', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    await withRestrictedUser('aactr3', async (contributor) => {
      await addToTier(contributor.group, program, 'contributor');
      const v1 = await contributor.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrResub'}));
      await withRestrictedUser('aaapp3', async (approver) => {
        await addToTier(approver.group, program, 'approver');
        await approver.asUser(() => rejectPublication(v1.rowId, 'wrong cohort'));
      });
      const v2 = await contributor.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrResub'}));
      expect(v2.revision, 2);
      expect(v2.status, 'pending');
      const rejected = (await getPublicationHistory(v1.publicationId)).find((h) => h.revision === 1);
      expect(rejected?.final_status, 'rejected');
      expect(rejected?.reject_reason, 'wrong cohort');
    });
  });

  test('three-actor lifecycle: publish, approve, discover, republish, approve, discover', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    await withRestrictedUser('aalcc', async (contributor) => {
      await addToTier(contributor.group, program, 'contributor');
      await withRestrictedUser('aalca', async (approver) => {
        await addToTier(approver.group, program, 'approver');
        await withRestrictedUser('aalcv', async (viewer) => {
          await addToTier(viewer.group, program, 'viewer');
          const v1 = await contributor.asUser(() => publishRun({
            sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Lifecycle'}));
          expect(v1.status, 'pending');
          await approver.asUser(() => approvePublication(v1.rowId));
          let seen = await viewer.asUser(() => discover(program.id));
          expect(seen.length, 1);
          expect(seen[0].revision, 1);
          const v2 = await contributor.asUser(() => publishRun({
            sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Lifecycle'}));
          expect(v2.status, 'pending');
          await approver.asUser(() => approvePublication(v2.rowId));
          seen = await viewer.asUser(() => discover(program.id));
          expect(seen.length, 1);
          expect(seen[0].revision, 2);
          const history = await viewer.asUser(() => getPublicationHistory(v1.publicationId));
          expect(history.some((h) => h.revision === 1 && h.final_status === 'approved'), true,
            'the viewer must see v1 in the archive');
        });
      });
    });
  });

  test('approver-tier publish auto-approves under their own rights', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    await withRestrictedUser('aaapub', async (probe) => {
      await addToTier(probe.group, program, 'approver');
      const result = await probe.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'ApproverPublish'}));
      expect(result.status, 'approved');
      const row: any = await alignment().get(result.rowId);
      expect(row.approved_by, probe.id, 'the auto-approval must stamp the acting approver');
    });
  });

  test('contributor curates via the service under their own session', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const marker = `mu${Date.now()}`;
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrCurate'});
    await withRestrictedUser('aacur6', async (probe) => {
      await addToTier(probe.group, program, 'contributor');
      await probe.asUser(() => updateCuration(pub.rowId, {
        tags: [`${marker}-t1`, `${marker}-t2`], compounds: [`GRK-${marker}`],
        path: 'PK/contributed', description: 'curated by contributor'}));
    });
    const row: any = (await getLiveVersions(pub.publicationId))[0];
    expect(row.path, 'PK/contributed');
    expect(row.tags.map((t: any) => t.name).sort().join(','), `${marker}-t1,${marker}-t2`);
    expect(row.compounds.map((c: any) => c.name).join(','), `GRK-${marker}`);
    expect(row.description, 'curated by contributor');
  });

  test('contributor publish creates missing tags and compounds in the registry', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const marker = `mr${Date.now()}`;
    await withRestrictedUser('aareg7', async (probe) => {
      await addToTier(probe.group, program, 'contributor');
      const pub = await probe.asUser(() => publishRun({
        sourceMetaCallId: run.metaCallId, programId: program.id, name: 'CtrRegistry',
        tags: [`${marker}-tag`], compounds: [`GRK-${marker}`]}));
      const tag = await grok.dapi.domains.table(T_TAG).getByKey({name: `${marker}-tag`});
      expect(tag != null, true, 'the tag must be created in the registry');
      const compound = await grok.dapi.domains.table(T_COMPOUND)
        .getByKey({registration_code: `GRK-${marker}`});
      expect(compound != null, true, 'the compound must be created in the registry');
      const row: any = (await getLiveVersions(pub.publicationId))[0];
      expect(row.tags.map((t: any) => t.name).join(','), `${marker}-tag`);
    });
  });

  test('contributor-published function run keeps its dataframes viewer-readable', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun(4);
    await withRestrictedUser('aafr8', async (contributor) => {
      await addToTier(contributor.group, program, 'contributor');
      const pub = await contributor.asUser(() => publishRun({
        sourceMetaCallId: run.stepCallId, programId: program.id, name: 'CtrFrozen'}));
      const frozen = await grok.dapi.functions.calls.allPackageVersions()
        .include('inputs, outputs').find(pub.artifactId);
      const tableId = Object.values(frozen.outputs).find((v) => typeof v === 'string') as string;
      expect(tableId != null, true, 'the frozen function run must reference an uploaded table');
      await withRestrictedUser('aavw8', async (viewer) => {
        await addToTier(viewer.group, program, 'viewer');
        const table = await viewer.asUser(() => grok.dapi.tables.getTable(tableId));
        expect(table.rowCount, 4, 'a program viewer must read the contributor-frozen dataframe');
      });
    });
  });

  test('viewer reads the discovery slice and publication history through the services', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'ViewerRead'});
    await publishRun({sourceMetaCallId: run.metaCallId, programId: program.id, name: 'ViewerRead'});
    await withRestrictedUser('aavr9', async (probe) => {
      await addToTier(probe.group, program, 'viewer');
      const slice = await probe.asUser(() => discover(program.id));
      expect(slice.length, 1);
      expect(slice[0].revision, 2);
      const history = await probe.asUser(() => getPublicationHistory(v1.publicationId));
      expect(history.length, 1);
      expect(history[0].revision, 1);
    });
  });

  test('approver-tier curation succeeds — the tier includes the contributor column grants', async () => {
    // updateCuration accepts approvers service-side, but the curation column Edit is
    // granted to the contributors umbrella only — this passes because an approver is
    // also a contributor by the tier model, and pins exactly that assumption.
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'ApprCurate'});
    await withRestrictedUser('aaac10', async (probe) => {
      await addToTier(probe.group, program, 'approver');
      await probe.asUser(() => updateCuration(pub.rowId, {path: 'PK/approver-tier'}));
      const row: any = await alignment().get(pub.rowId);
      expect(row.path, 'PK/approver-tier');
    });
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 360000});
