import * as grok from 'datagrok-api/grok';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {APPROVAL_COLUMNS, CURATION_COLUMNS, T_ALIGNMENT, T_PROGRAM} from '../domain/constants';
import {currentUserIn, ensureProgram, findGroup, setupSchemaSecurity} from '../domain/security';
import {approvePublication, publishWorkflowRun, rejectPublication} from '../service/publication-service';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram} from './fixtures';

const SECOND_USER = 'requires a second user session (restricted-account run)';
const CORE_FUNCCALL_ACL = 'requires core: FuncCall per-entity ACL is not enforced yet';
const CORE_OWNERSHIP = 'requires core: service-identity ownership transfer for funccalls';

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
    const caps = await grok.dapi.domains.table(T_ALIGNMENT).capabilities();
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
    const result = await publishWorkflowRun({
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
    // The group creator is always a transitive member, so exercising the
    // membership denial needs a session for a user outside the group.
  }, {skipReason: SECOND_USER});

  test('non-approver typed denial on a status write (column security)', async () => {
    // Verified manually on the original fixtures ('Column not writable'); asserting it
    // headlessly needs a session for a user outside the approvers umbrella.
  }, {skipReason: SECOND_USER});

  test('viewers group member sees the program and its publications; outsiders do not', async () => {
    // Row-level trimming through the program delegate chain; needs a restricted account.
  }, {skipReason: SECOND_USER});

  test('non-audience user cannot read the frozen run dataframes', async () => {
    // The publish clone grants View on the frozen tables to the viewers group only;
    // asserting the denial needs a session outside that group.
  }, {skipReason: SECOND_USER});

  test('frozen funccall shells are not loadable by outsiders', async () => {
    // FuncCall is not a full entity: no per-funccall ACL exists, so the frozen meta
    // call and step calls are unlisted but public-by-id. The commented grants in
    // clone-service.ts mark where enforcement slots in.
  }, {skipReason: CORE_FUNCCALL_ACL});

  test('frozen copy is immutable for the author', async () => {
    // Immutability-by-ownership needs the service identity to own the frozen calls;
    // today the publishing user's session owns them and could re-save them by id.
  }, {skipReason: CORE_OWNERSHIP});

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
