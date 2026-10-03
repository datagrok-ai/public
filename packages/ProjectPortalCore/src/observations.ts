/** Observations (F-26): curated conclusions on a program, its work units and its compounds, stored
 * in `portal.observation`. The rules are pure (src/model/observations.ts); this reads and writes. */
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {ObservationRow, portalDb, PortfolioRow} from './generated/db';
import {compoundsOf, enabledWorkUnits, workUnitsOf} from './portal';
import {programAccess} from './program-roles';
import {checkObservationAccess, checkRetraction, evidenceOfObservation, ObservationInput, ObservationTargets,
  ObservationView, sortObservations, subjectLabel, validateObservation} from './model/observations';
import {filterValue, ROW_LIMIT} from './rows';

const iso = (d: any): string | null => d ? dayjs(d).toISOString() : null;

/** What the program holds that an observation may be about: its series, campaigns and studies
 * (the kinds its portfolio enables) and its compounds. */
export async function observationTargetsOf(portfolio: PortfolioRow, programId: string): Promise<ObservationTargets> {
  const out: ObservationTargets = {};
  for (const kind of enabledWorkUnits(portfolio)) {
    if (kind === 'series' || kind === 'campaign' || kind === 'study')
      out[kind] = (await workUnitsOf(kind, [programId]).catch(() => [])).map((u) => String(u.code));
  }
  out.compound = (await compoundsOf(programId).catch(() => [])).map((c) => String(c.compound.code));
  return out;
}

export function observationView(r: ObservationRow): ObservationView {
  return {id: r.id, subject: r.subject_kind, code: r.subject_code ?? null, title: r.title, finding: r.finding,
    evidence: evidenceOfObservation(r.evidence), status: r.status ?? 'current', by: r.observed_by ?? r.author_id ?? null,
    on: iso(r.observed_on ?? r.created_on), retractedBy: r.retracted_by ?? null, retractedOn: iso(r.retracted_on),
    retraction: r.retraction ?? null};
}

/** The program's observations: current first, each group newest first. */
export async function observationsOf(programId: string): Promise<ObservationView[]> {
  const rows = await portalDb.observations.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT});
  return sortObservations(rows.map(observationView));
}

/** Records an observation as the current user, after checking who they are and what it is about. */
export async function addObservation(programId: string, o: ObservationInput): Promise<ObservationView> {
  const program = await portalDb.programs.get(programId);
  const portfolio = await portalDb.portfolios.get(program.portfolio_id);
  const access = await programAccess(portfolio, programId);
  const errors = [...checkObservationAccess({isApprover: access.approver, isContributor: access.contributor}),
    ...validateObservation(o, await observationTargetsOf(portfolio, programId))];
  if (errors.length)
    throw new Error(errors.join('; '));
  const [row] = await portalDb.observations.insert({program_id: programId, subject_kind: o.subject,
    subject_code: o.subject === 'program' ? undefined : o.code!.trim(), title: o.title.trim(), finding: o.finding.trim(),
    // a json column holds a map: the list goes under items
    evidence: {items: (o.evidence ?? []).map((e) => ({kind: e.kind, target: e.target.trim(), name: e.name.trim()}))},
    status: 'current', observed_by: DG.User.current().id, observed_on: dayjs()});
  return observationView(await portalDb.observations.get(row.id));
}

/** Retracts an observation: kept, marked, with who, when and why. Its author or an approver. */
export async function retractObservation(id: string, reason: string): Promise<ObservationView> {
  const row = await portalDb.observations.get(id);
  const program = await portalDb.programs.get(row.program_id);
  const portfolio = await portalDb.portfolios.get(program.portfolio_id);
  const access = await programAccess(portfolio, program.id);
  const view = observationView(row);
  const errors = checkRetraction({status: view.status, by: view.by}, DG.User.current().id, access.approver, reason);
  if (errors.length)
    throw new Error(errors.join('; '));
  await portalDb.observations.update(id, {status: 'retracted', retracted_by: DG.User.current().id, retracted_on: dayjs(),
    retraction: reason.trim()});
  return observationView(await portalDb.observations.get(id));
}

/** The observations as a table. */
export function observationsFrame(list: ObservationView[], names: Map<string, string>, name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('About', list.map((o) => subjectLabel(o.subject, o.code))),
    DG.Column.fromStrings('Title', list.map((o) => o.title)),
    DG.Column.fromStrings('Finding', list.map((o) => o.finding)),
    DG.Column.fromStrings('Evidence', list.map((o) => o.evidence.map((e) => e.name).join('; '))),
    DG.Column.fromStrings('Status', list.map((o) => o.status)),
    DG.Column.fromStrings('By', list.map((o) => names.get(o.by ?? '') ?? '')),
    DG.Column.fromStrings('On', list.map((o) => o.on?.slice(0, 10) ?? '')),
    DG.Column.fromStrings('Id', list.map((o) => o.id)),
  ]);
  df.name = name;
  return df;
}
