/** Parking orphan folders (FR-22, the action half of the orphan report). The folders the report
 * lists (root folders no program claims, folders of programs that ended or left their feed) are
 * moved, never deleted, into `_Parked/<date>` under the portfolio's root, each move recorded so it
 * can be undone. Deleting is out of scope: a parked folder stays until a person removes it by hand.
 *
 * Pure: the caller reads the report, the parked records and who the user is. */
import {isParkingFolder, isSpaceNamed, Orphan, PARKED_FOLDER, SpaceRef} from './workspace-plan';

export {isParkingFolder, PARKED_FOLDER};

/** The day folder a park on [when] goes into: `2026-10-02`, in UTC. */
export const parkDay = (when: Date): string => when.toISOString().slice(0, 10);

export interface ParkMove {
  /** The folder's space id. */
  id: string;
  name: string;
  kind: Orphan['kind'];
  reason: string;
  /** Where it is now and where it goes, as paths under the root. */
  from: string;
  to: string;
}

export interface ParkPlan {
  root: string;
  day: string;
  moves: ParkMove[];
  /** Orphans the plan leaves where they are, and why. */
  skipped: {path: string, reason: string}[];
}

/** What parking the orphans [orphans] of [root] on [day] does: each goes to `root/_Parked/day/`.
 * [only] narrows it to those paths or ids; a folder already parked is skipped. */
export function planParking(root: string, orphans: Orphan[], day: string, parkedIds: Iterable<string> = [],
  only?: string[]): ParkPlan {
  const parked = new Set(parkedIds);
  const plan: ParkPlan = {root, day, moves: [], skipped: []};
  for (const o of orphans) {
    if (only?.length && !only.includes(o.path) && !only.includes(o.id))
      continue;
    const name = o.path.startsWith(`${root}/`) ? o.path.slice(root.length + 1) : o.path;
    if (parked.has(o.id)) {
      plan.skipped.push({path: o.path, reason: 'already parked'});
      continue;
    }
    plan.moves.push({id: o.id, name, kind: o.kind, reason: o.reason, from: o.path, to: `${root}/${PARKED_FOLDER}/${day}/${name}`});
  }
  return plan;
}

/** Why the current user may not park or unpark, empty when they may: the portfolio's approvers. */
export function checkParking(c: {isApprover: boolean, dryRun: boolean}): string[] {
  return c.dryRun || c.isApprover ? [] : ['Only the portfolio\'s approvers park or unpark folders; preview with dryRun'];
}

/** A parked folder as recorded. */
export interface ParkedRef {
  id: string;
  name: string;
  status: 'parked' | 'restored';
  fromPath: string;
  parkedPath: string;
}

/** Why [p] cannot go back, empty when it can: it is parked, and no folder of the same name has
 * taken its place under the root since (the builder may have made a new one for its program). */
export function checkUnpark(p: ParkedRef, rootChildren: SpaceRef[]): string[] {
  const refused: string[] = [];
  if (p.status !== 'parked')
    refused.push(`${p.name} is not parked (${p.status})`);
  else if (rootChildren.some((c) => isSpaceNamed(c, p.name)))
    refused.push(`a folder named ${p.name} is back under the root; move or rename it first`);
  return refused;
}

/** The plan as text: a line per move, then what is left where it is. */
export function summarizeParking(code: string, plan: ParkPlan, dryRun: boolean, errors: string[] = []): string {
  const head = `${code} workspace "${plan.root}": ${plan.moves.length} folder(s) ${dryRun ? 'would be parked' : 'parked'} ` +
    `in ${PARKED_FOLDER}/${plan.day}; nothing deleted`;
  return [head,
    ...plan.moves.map((m) => `  ${dryRun ? '→' : '✓'} ${m.from} → ${m.to} (${m.kind === 'park' ? 'program ended: ' : ''}${m.reason})`),
    ...plan.skipped.map((s) => `  = ${s.path}: ${s.reason}`),
    ...errors.map((e) => `  ! ${e}`)].join('\n');
}
