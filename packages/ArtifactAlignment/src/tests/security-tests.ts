import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  APPROVAL_COLUMNS, CURATION_COLUMNS, T_ALIGNMENT, T_COMPOUND, T_HISTORY, T_PROGRAM, T_STUDY, T_TAG,
} from '../domain/constants';
import {currentUserIn, ensureProgram, findGroup, setupSchemaSecurity} from '../domain/security';
import {
  approvePublication, discover, getPublicationHistory, rejectPublication,
} from '../service/publication-service';
import {
  addToTier, cleanupTestPrograms, makeSavedRun, makeTestProgram, makeTestStudy, publishRun,
  withRestrictedUser,
} from './fixtures';

const alignment = () => grok.dapi.domains.table(T_ALIGNMENT);

category('ArtifactAlignment: security', () => {
  test('program provisioning is idempotent and wires groups into the row', async () => {
    const program = await makeTestProgram();
    const again = await ensureProgram({code: program.code, name: 'renamed'});
    expect(again.id, program.id);
    expect(again.viewers.id, program.viewers.id);
    const row: any = await grok.dapi.domains.table(T_PROGRAM).get(program.id);
    expect(row.viewers_group, program.viewers.id);
    expect(row.contributors_group, program.contributors.id);
    expect(row.approvers_group, program.approvers.id);
    expect(row.name, 'renamed');
  });

  test('schema security setup is idempotent and leaves approval/curation columns writable for admin', async () => {
    await setupSchemaSecurity();
    await setupSchemaSecurity();
    const caps = await alignment().capabilities();
    for (const column of [...APPROVAL_COLUMNS, ...CURATION_COLUMNS]) {
      expect(caps.writableColumns.includes(column), true,
        `'${column}' not writable for admin after column restriction`);
    }
  });

  test('per-program groups are members of the umbrella groups', async () => {
    const program = await makeTestProgram();
    const approversUmbrella = await findGroup('Artifacts: Approvers');
    const membershipIds = [
      ...(await grok.dapi.groups.find(program.approvers.id)).memberships,
      ...(await grok.dapi.groups.find(program.approvers.id)).adminMemberships,
    ].map((g) => g.id);
    expect(membershipIds.includes(approversUmbrella!.id), true,
      'program approvers must be included in the umbrella that holds column Edit');
  });

  test('transitive membership resolves through the program group into the umbrella', async () => {
    const program = await makeTestProgram();
    // The creator becomes an admin member of the groups they create, so the current
    // user is transitively in the fresh approvers group — and through its umbrella
    // membership, in the umbrella too.
    expect(await currentUserIn(program.approvers), true);
    const umbrella = await findGroup('Artifacts: Approvers');
    expect(await currentUserIn(umbrella), true, 'membership must be transitive through the program group');
  });

  test('service-side reviewer checks — author cannot approve or reject their own publication', async () => {
    // The current user is in the approvers group as its creator, so the membership
    // check passes and the author != approver rule is what rejects.
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const result = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'SelfApprove',
      skipAutoApprove: true});
    let approveError = '';
    try {
      await approvePublication(result.rowId);
    } catch (e: any) {
      approveError = `${e?.message ?? e}`;
    }
    expect(approveError.includes('Authors cannot approve'), true, `got: '${approveError}'`);
    let rejectError = '';
    try {
      await rejectPublication(result.rowId, 'no');
    } catch (e: any) {
      rejectError = `${e?.message ?? e}`;
    }
    expect(rejectError.includes('Authors cannot review'), true, `got: '${rejectError}'`);
  });

  test('non-member of the approvers group cannot approve (service check)', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const result = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'OutsiderApprove',
      skipAutoApprove: true});
    await withRestrictedUser('aaoutappr', async (probe) => {
      // viewer tier: can read the row, holds no approval authority
      await addToTier(probe.group, program, 'viewer');
      let error = '';
      try {
        await probe.asUser(() => approvePublication(result.rowId));
      } catch (e: any) {
        error = `${e?.message ?? e}`;
      }
      expect(error.includes('Only members of the program approvers group'), true, `got: '${error}'`);
      const row: any = await alignment().get(result.rowId);
      expect(row.status, 'pending', 'the denied approval must not stick');
    });
  });

  test('non-approver typed denial on a status write (column security)', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const result = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'ColumnGate',
      skipAutoApprove: true});
    await withRestrictedUser('aacolsec', async (probe) => {
      await addToTier(probe.group, program, 'contributor');
      // a curation column IS writable for the contributor — what follows is the
      // column gate on the approval columns, not a blanket write denial
      await probe.asUser(() => alignment().update(result.rowId, {path: 'probe/ok'} as any));
      let error: any = null;
      try {
        await probe.asUser(() => alignment().update(result.rowId, {status: 'approved'} as any));
      } catch (e) {
        error = e;
      }
      expect(error instanceof DG.DomainError, true, `expected a typed denial, got: ${error}`);
      const row: any = await alignment().get(result.rowId);
      expect(row.status, 'pending', `the denied status write must not stick ('${error?.message}')`);
      expect(row.path, 'probe/ok', 'the curation write must have gone through');
    });
  });

  test('viewers group member sees the program and its publications; outsiders do not', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Visible'});
    await withRestrictedUser('aavis', async (probe) => {
      const programsSeen = () => probe.asUser(async () =>
        (await grok.dapi.domains.table(T_PROGRAM).query({filter: DG.cond('id', '=', program.id)})).length);
      const versionsSeen = () => probe.asUser(async () =>
        (await alignment().query({filter: DG.cond('id', '=', pub.rowId)})).length);
      expect(await programsSeen(), 0, 'an outsider must not see the program');
      expect(await versionsSeen(), 0, 'an outsider must not see its publications');
      await addToTier(probe.group, program, 'viewer');
      expect(await programsSeen(), 1, 'a viewers member must see the program');
      expect(await versionsSeen(), 1, 'a viewers member must see its publications');
    });
  });

  test('non-audience user cannot read the frozen run dataframes', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun(4);
    // a function-run publish re-uploads the dataframes with the viewers-only grant
    const pub = await publishRun({
      sourceMetaCallId: run.stepCallId, programId: program.id, name: 'FrozenTables'});
    const frozen = await grok.dapi.functions.calls.allPackageVersions()
      .include('inputs, outputs').find(pub.artifactId);
    const tableId = Object.values(frozen.outputs).find((v) => typeof v === 'string') as string;
    expect(tableId != null, true, 'the frozen function run must reference an uploaded table');
    await withRestrictedUser('aatbl', async (probe) => {
      const readTable = () => probe.asUser(() => grok.dapi.tables.getTable(tableId));
      let outsiderRead: DG.DataFrame | null = null;
      try {
        outsiderRead = await readTable();
      } catch (_) {/* denied */}
      expect(outsiderRead == null, true, 'an outsider must not read the frozen dataframe');
      await addToTier(probe.group, program, 'viewer');
      expect((await readTable()).rowCount, 4, 'a viewers member must read the frozen dataframe');
    });
  });

  test('an approver who is not the author approves, and the approval stamps them', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'TwoActorApprove',
      skipAutoApprove: true});
    await withRestrictedUser('aaappr', async (probe) => {
      await addToTier(probe.group, program, 'approver');
      await probe.asUser(() => approvePublication(pub.rowId));
      const row: any = await alignment().get(pub.rowId);
      expect(row.status, 'approved');
      expect(row.approved_by, probe.id, 'the approval must stamp the acting approver');
    });
  });

  test('an approver who is not the author rejects with a reason', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'TwoActorReject',
      skipAutoApprove: true});
    await withRestrictedUser('aarej', async (probe) => {
      await addToTier(probe.group, program, 'approver');
      await probe.asUser(() => rejectPublication(pub.rowId, 'not aligned'));
      const row: any = await alignment().get(pub.rowId);
      expect(row.status, 'rejected');
      expect(row.reject_reason, 'not aligned');
      expect(row.approved_by, probe.id, 'the rejection must stamp the acting reviewer');
    });
  });

  test('outsider gets an empty discovery slice and empty history', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'OutDisc'});
    await publishRun({sourceMetaCallId: run.metaCallId, programId: program.id, name: 'OutDisc'});
    await withRestrictedUser('aaodisc', async (probe) => {
      expect((await probe.asUser(() => discover(program.id))).length, 0,
        'discovery must be empty for an outsider');
      expect((await probe.asUser(() => getPublicationHistory(v1.publicationId))).length, 0,
        'history must be empty for an outsider');
    });
  });

  test('outsider cannot read the study registry', async () => {
    const program = await makeTestProgram();
    const study = await makeTestStudy(program.id);
    await withRestrictedUser('aastdy', async (probe) => {
      const seen = await probe.asUser(() => grok.dapi.domains.table(T_STUDY)
        .query({filter: DG.cond('id', '=', study.id)}));
      expect(seen.length, 0, 'the study registry must be hidden from outsiders');
    });
  });

  test('outsider cannot read alignment_history', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const v1 = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'OutHist'});
    await publishRun({sourceMetaCallId: run.metaCallId, programId: program.id, name: 'OutHist'});
    expect((await getPublicationHistory(v1.publicationId)).length, 1, 'admin must see the archived row');
    await withRestrictedUser('aahist', async (probe) => {
      const seen = await probe.asUser(() => grok.dapi.domains.table(T_HISTORY)
        .query({filter: DG.cond('publication_id', '=', v1.publicationId)}));
      expect(seen.length, 0, 'the archive must be hidden from outsiders');
    });
  });

  test('compound and tag registries stay readable for everyone (by design)', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const marker = `pub${Date.now()}`;
    await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'OpenReg',
      tags: [`${marker}-tag`], compounds: [`GRK-${marker}`]});
    await withRestrictedUser('aaoreg', async (probe) => {
      const tag = await probe.asUser(() =>
        grok.dapi.domains.table(T_TAG).getByKey({name: `${marker}-tag`}));
      expect(tag != null, true, 'tags are all-users readable by design');
      const compound = await probe.asUser(() => grok.dapi.domains.table(T_COMPOUND)
        .getByKey({registration_code: `GRK-${marker}`}));
      expect(compound != null, true, 'compounds are all-users readable by design');
    });
  });

  test('cross-program raw table write is open — canary until core row-scoped predicates', async () => {
    // Table Edit is umbrella-wide: a contributor of program A can raw-insert a row
    // into program B; per-program isolation exists only in the services (a
    // master-row-scoped write predicate is a core ask). When this canary fails, the
    // predicate landed: turn it into a denial test.
    const program = await makeTestProgram();
    const foreign = await makeTestProgram();
    await withRestrictedUser('aaxw', async (probe) => {
      await addToTier(probe.group, program, 'contributor');
      const inserted = await probe.asUser(() => alignment().insert({
        publication_id: crypto.randomUUID(), revision: 1, name: 'CrossWrite',
        artifact_id: crypto.randomUUID(), program_id: foreign.id})).catch(() => null);
      expect(inserted != null, true, 'cross-program raw writes are now denied — ' +
        'core row-scoped predicate landed; write the real denial test');
    });
  });

  test('workflow-run step dataframes are open to outsiders — canary until scoped step grants', async () => {
    // Workflow freezes share the source run's step calls; their dataframes keep the
    // personal-history default grant. When this canary fails, step dataframes got
    // audience-scoped: write the real denial test mirroring the function-run one above.
    const program = await makeTestProgram();
    const run = await makeSavedRun(3);
    await publishRun({sourceMetaCallId: run.metaCallId, programId: program.id, name: 'StepTables'});
    const step = await grok.dapi.functions.calls.allPackageVersions()
      .include('inputs, outputs').find(run.stepCallId);
    const tableId = Object.values(step.outputs).find((v) => typeof v === 'string') as string;
    expect(tableId != null, true, 'the shared step must reference an uploaded table');
    await withRestrictedUser('aastep', async (probe) => {
      const table = await probe.asUser(() => grok.dapi.tables.getTable(tableId)).catch(() => null);
      expect(table != null, true, 'step dataframes are no longer world-readable — ' +
        'scoped grants landed; write the real denial test');
    });
  });

  test('frozen funccall shells are loadable by outsiders — canary until core per-call ACL', async () => {
    // FuncCall is not a full entity: no per-funccall ACL exists, so frozen calls are
    // unlisted but public-by-id. When this canary fails, core has landed the ACL:
    // write the real denial test and enable the commented grants in clone-service.ts.
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'AclCanary',
      skipAutoApprove: true});
    await withRestrictedUser('aafc', async (probe) => {
      const loaded = await probe.asUser(() =>
        grok.dapi.functions.calls.allPackageVersions().find(pub.artifactId).catch(() => null));
      expect(loaded != null, true,
        'outsiders can no longer load frozen calls by id — core FuncCall ACL landed; write the real test');
    });
  });

  test('frozen copy stays author-mutable — canary until core ownership transfer', async () => {
    // The publishing user's session owns the frozen call and can re-save it by id.
    // When this canary fails, core has landed service-identity ownership transfer:
    // write the real immutability test.
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const pub = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'OwnerCanary',
      skipAutoApprove: true});
    const frozen = await grok.dapi.functions.calls.allPackageVersions().find(pub.artifactId);
    frozen.options['aa-ownership-canary'] = 'touched';
    const saved = await grok.dapi.functions.calls.allPackageVersions().save(frozen).catch(() => null);
    expect(saved != null, true,
      'frozen calls are no longer author-mutable — core ownership transfer landed; write the real test');
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 240000});
