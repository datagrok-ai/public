import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {expect} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {saveInstanceState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {T_ALIGNMENT, T_HISTORY, T_PROGRAM, T_STUDY} from '../domain/constants';
import {ensureProgram, ProgramInfo} from '../domain/security';
import {publishWorkflowRun, PublishRequest, PublishResult} from '../service/publication-service';

export const STEP_NQ = 'ArtifactAlignment:AaTestStep';
export const PROVIDER_NQ = 'ArtifactAlignment:AaTestProvider';

export interface SavedRunFixture {
  metaCallId: string;
  stepCallId: string;
}

/** Every FuncCall a test creates registers here so cleanup deletes exactly this
 * run's artifacts by id — no server-wide listing. */
const trackedCallIds = new Set<string>();

export function trackCalls(...ids: string[]): void {
  for (const id of ids)
    trackedCallIds.add(id);
}

/** Publishes through the service and tracks the frozen artifact for cleanup.
 * Tests use this instead of calling publishWorkflowRun directly. */
export async function publishRun(req: PublishRequest): Promise<PublishResult> {
  const result = await publishWorkflowRun(req);
  trackCalls(result.artifactId);
  return result;
}

/** A genuine saved compute2-style workflow run, built headlessly: an executed step
 * plus a meta call carrying the pipeline state that references it. */
export async function makeSavedRun(a: number = 3): Promise<SavedRunFixture> {
  const stepFunc = DG.Func.byName(STEP_NQ);
  const fc = stepFunc.prepare({a});
  await fc.call(false, undefined, {processed: true, report: false});
  fc.newId();
  const savedStep = await historyUtils.saveRun(fc);
  const state = {
    type: 'static',
    configId: 'root',
    uuid: crypto.randomUUID(),
    isReadonly: false,
    nqName: PROVIDER_NQ,
    friendlyName: 'AA fixture workflow',
    steps: [{
      type: 'funccall',
      configId: 'step1',
      uuid: crypto.randomUUID(),
      isReadonly: false,
      nqName: STEP_NQ,
      funcCallId: savedStep.id,
      friendlyName: 'Step 1',
    }],
  };
  const metaCall = await saveInstanceState(PROVIDER_NQ, state, {title: 'AA fixture run'});
  trackCalls(savedStep.id, metaCall.id);
  return {metaCallId: metaCall.id, stepCallId: savedStep.id};
}

const createdPrograms: ProgramInfo[] = [];

export async function makeTestProgram(): Promise<ProgramInfo> {
  const code = `ZZT-${Date.now()}-${Math.floor(Math.random() * 1e5)}`;
  const info = await ensureProgram({code, name: `Test program ${code}`});
  createdPrograms.push(info);
  return info;
}

export async function makeTestStudy(programId: string): Promise<{id: string, protocolCode: string}> {
  const protocolCode = `ZZS-${Date.now()}-${Math.floor(Math.random() * 1e5)}`;
  const [ins] = await grok.dapi.domains.table(T_STUDY).insert(
    {program_id: programId, protocol_code: protocolCode, phase: 'ph2', status: 'ongoing'});
  return {id: ins.id, protocolCode};
}

/** Deletes one tracked run together with the dataframes it uploaded, all by id. */
async function deleteTrackedRun(callId: string): Promise<void> {
  const source = () => grok.dapi.functions.calls.allPackageVersions();
  const fc = await source().include('inputs, outputs').find(callId).catch(() => null);
  if (fc == null)
    return;
  // non-materialized dataframe params hold the uploaded table id
  const tableIds: string[] = [];
  for (const [params, values] of
    [[fc.inputParams, fc.inputs], [fc.outputParams, fc.outputs]] as const) {
    for (const p of params.values() as Iterable<DG.FuncCallParam>) {
      if (p.property.propertyType === DG.TYPE.DATA_FRAME && typeof values[p.name] === 'string')
        tableIds.push(values[p.name]);
    }
  }
  for (const id of tableIds) {
    try {
      const table = await grok.dapi.tables.find(id);
      if (table != null)
        await grok.dapi.tables.delete(table);
    } catch (_) {/* already gone */}
  }
  // the plain calls source silently no-ops the delete — allPackageVersions is required
  await source().delete(fc);
  const remains = await source().find(callId).catch(() => null);
  if (remains != null)
    throw new Error(`tracked run ${callId} survived its delete`);
}

/** Deletes every tracked run by id; frozen workflow clones share their source's
 * steps and dataframes, so sources + frozen ids cover everything. */
export async function cleanupTrackedRuns(): Promise<void> {
  const ids = [...trackedCallIds];
  trackedCallIds.clear();
  await Promise.all(ids.map((id) => deleteTrackedRun(id)));
}

/** Manual utility draining fixture runs that escaped tracking; not wired into any
 * category. */
export async function sweepFixtureBacklog(budgetMs: number = 45000): Promise<number> {
  const started = Date.now();
  let deleted = 0;
  for (const nqName of [PROVIDER_NQ, STEP_NQ]) {
    const funcName = nqName.split(':')[1];
    while (Date.now() - started < budgetMs) {
      const page = await grok.dapi.functions.calls.allPackageVersions()
        .filter(`func.name="${funcName}"`).list({pageSize: 50});
      if (page.length === 0)
        break;
      await Promise.all(page.map((fc) => deleteTrackedRun(fc.id)));
      deleted += page.length;
    }
  }
  return deleted;
}

/** A freshly signed-up user with no privileges, and the way to act as them. */
export interface RestrictedUser {
  login: string;
  id: string;
  /** The user's personal group — what a grant or membership is addressed to. */
  group: string;
  /** The user's session token — for raw authenticated fetches when a test
   * asserts on the HTTP surface itself; asUser covers everything else. */
  token: string;
  /** Runs [action] under THAT user's session, restoring the admin one after. */
  asUser<T>(action: () => Promise<T>): Promise<T>;
}

/** Adds the probe's personal group into [group] (admin session). Test groups are
 * deleted wholesale by cleanupTestPrograms, so memberships need no undo. */
export async function addToGroup(memberGroupId: string, group: DG.Group): Promise<void> {
  const g = await grok.dapi.groups.find(group.id);
  g.addMember(await grok.dapi.groups.find(memberGroupId));
  await grok.dapi.groups.saveRelations(g);
}

export type ProgramRole = 'viewer' | 'contributor' | 'approver';

/** Assigns [role] the way production does: the role's group plus every lower
 * tier (viewer ⊂ contributor ⊂ approver), so an approver probe is also a
 * contributor and a viewer. */
export async function addToTier(memberGroupId: string, program: ProgramInfo,
  role: ProgramRole): Promise<void> {
  const tiers: DG.Group[] = [program.viewers];
  if (role !== 'viewer')
    tiers.push(program.contributors);
  if (role === 'approver')
    tiers.push(program.approvers);
  for (const group of tiers)
    await addToGroup(memberGroupId, group);
}

/** Runs [body] with a throwaway restricted user (adapted from ApiTests'
 * domain-lifecycle helper — packages cannot share test code). Impersonation goes
 * through `grok.dapi.impersonationToken` (the session cookie is HttpOnly); the
 * finally BLOCKS the user (not API-deletable). Throws where self-signup is
 * unavailable — a multiuser test must fail loudly, not pass vacuously. */
export async function withRestrictedUser<T>(prefix: string,
  body: (user: RestrictedUser) => Promise<T>): Promise<T> {
  const restoreAdmin = () => grok.dapi.impersonationToken = null;
  const stamp = `${Date.now()}${Math.floor(Math.random() * 1e4)}`;
  const login = `${prefix}${stamp}`;
  let user: any = null;
  try {
    // credentials: 'omit' is load-bearing — same-origin fetch would otherwise
    // store the signup's Set-Cookie over the admin's HttpOnly auth cookie.
    const signup = await (await fetch(`${grok.dapi.root}/users/signup`, {
      method: 'POST', credentials: 'omit', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({login: login, email: `${login}@test.datagrok.ai`,
        password: btoa(`Pw-${stamp}`), firstName: 'AaTests', lastName: 'Probe'}),
    })).json();
    if (!signup?.token) {
      throw new Error(`restricted-user signup unavailable (${signup?.reason ?? 'no token'}) — ` +
        'multiuser tests need self-signup enabled on the stand');
    }
    user = await grok.dapi.users.filter(`login = "${login}"`).include('group').first();
    if (user == null)
      throw new Error('signed-up user not found');
    return await body({
      login: login,
      id: user.id,
      group: user.group.id,
      token: signup.token,
      asUser: async <R>(action: () => Promise<R>): Promise<R> => {
        grok.dapi.impersonationToken = signup.token;
        try {
          return await action();
        } finally {
          restoreAdmin();
        }
      },
    });
  } finally {
    restoreAdmin();
    if (user?.id == null) {
      try {
        user = await grok.dapi.users.filter(`login = "${login}"`).first();
      } catch (x) {
        console.error(`test user ${login} could not be looked up — block it manually: ${x}`);
      }
    }
    if (user?.id != null) {
      const blocked = await fetch(`${grok.dapi.root}/users/block`, {
        method: 'POST', credentials: 'include', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({'#type': 'User', 'id': user.id, 'login': login}),
      });
      expect(blocked.ok, true, `test user ${login} was not blocked: ${blocked.status}`);
    }
  }
}

/** Removes everything the fixtures created: tracked runs and their dataframes,
 * version rows, history, studies, program rows, and the per-program groups. */
export async function cleanupTestPrograms(): Promise<void> {
  await cleanupTrackedRuns();
  for (const program of createdPrograms.splice(0)) {
    const byProgram = DG.cond('program_id', '=', program.id);
    while ((await grok.dapi.domains.table(T_ALIGNMENT).deleteWhere(byProgram)).hasMore);
    while ((await grok.dapi.domains.table(T_HISTORY).deleteWhere(byProgram)).hasMore);
    while ((await grok.dapi.domains.table(T_STUDY).deleteWhere(byProgram)).hasMore);
    await grok.dapi.domains.table(T_PROGRAM).delete(program.id);
    for (const group of [program.viewers, program.contributors, program.approvers]) {
      try {
        await grok.dapi.groups.delete(group);
      } catch (_) {/* platform may refuse deleting groups still referenced */}
    }
  }
}
