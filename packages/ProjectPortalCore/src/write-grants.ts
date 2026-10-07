/** Server-side table grants: Edit on the workflow tables (program_role, parked_folder, approval,
 * sync_log) to the portfolios' groups the app lets write them, View on parked_folder to the
 * approvers, and the tables that leaked across portfolios closed to All users. The rules are in
 * ./model/write-grants.
 *
 * A schema redeploy never revokes a grant (the manifest's `grants` are additive), so the schema
 * declaring View alone is not enough on a server that already has All users' Edit: the revoke
 * runs here. Changing grants needs Share on the tables, which the publisher and admins hold, so
 * this runs when portfolios are loaded and from the securePortalTables function. */
import * as grok from 'datagrok-api/grok';
import {portalDb, PortfolioRow} from './generated/db';
import {describeWriteGrants, Grant, MANAGED_TABLES, planWriteGrants, PortfolioGroups, WriteGrantPlan} from './model/write-grants';
import {ROW_LIMIT} from './rows';

const table = (name: string) => grok.dapi.domains.table(`portal.${name}`);

/** The grants on each managed table now. Needs Share on the tables. */
export async function currentWriteGrants(): Promise<Record<string, Grant[]>> {
  const out: Record<string, Grant[]> = {};
  for (const name of MANAGED_TABLES)
    out[name] = (await table(name).grants()).map((g) => [g.group.id, String(g.permission)] as Grant);
  return out;
}

const groupsOf = (p: PortfolioRow): PortfolioGroups => ({code: p.code, approvers: p.approvers, contributors: p.contributors});

/** Gives the writer and reader groups of [portfolios] (every portfolio when omitted) their grants
 * on the managed tables and takes All users' away. Idempotent; with [dryRun] only reports. A
 * user without Share on the tables gets a report line instead of an error: the tables stay as
 * they were. */
export async function secureWriteTables(portfolios?: PortfolioGroups[], dryRun = false): Promise<string[]> {
  let plan: WriteGrantPlan;
  try {
    const all = portfolios ?? (await portalDb.portfolios.query({limit: ROW_LIMIT})).map(groupsOf);
    plan = planWriteGrants(all, await currentWriteGrants());
  } catch (e: any) {
    return [`write permissions not checked: ${e?.message ?? e}`];
  }
  const names = new Map<string, string>();
  for (const id of new Set(plan.grant.map((g) => g.group))) {
    const g = await grok.dapi.groups.find(id).catch(() => null);
    names.set(id, g ? g.friendlyName || g.name : id);
  }
  const lines = describeWriteGrants(plan, (id) => names.get(id) ?? id);
  if (dryRun || !lines.length)
    return lines.length ? lines.map((l) => `would set ${l}`) : [];
  try {
    // grants first: a failure part-way leaves the old, wider grants rather than nobody able to write
    for (const g of plan.grant)
      await table(g.table).grant(g.group, g.permission);
    for (const r of plan.revoke)
      await table(r.table).revoke(r.group, r.permission);
  } catch (e: any) {
    return [`write permissions not set: ${e?.message ?? e}`];
  }
  return lines;
}
