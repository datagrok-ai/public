/** Design sets and selection sets (F-23): stored in `design_set`, `design_proposal` and
 * `selection_set`. The rules are pure (src/model/design-sets.ts); this reads and writes them. */
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {DesignProposalRow, DesignSetRow, portalDb, PortfolioRow, ProgramRow, SelectionSetRow} from './generated/db';
import {compoundsOf, portfolioByCode, workUnitsOf} from './portal';
import {programAccess} from './program-roles';
import {checkDesignAccess, DesignSetInput, designSetCode, rankProposals, SelectionSetInput, selectionCompounds,
  validateDesignSet, validateSelectionSet} from './model/design-sets';
import {filterValue, ROW_LIMIT} from './rows';
import {newestFirstBy} from './model/time';

export interface DesignSetItem {
  set: DesignSetRow;
  /** Best rank first; unranked last, by label. */
  proposals: DesignProposalRow[];
}

async function accessOf(program: ProgramRow): Promise<string[]> {
  const portfolio = await portalDb.portfolios.get(program.portfolio_id);
  const a = await programAccess(portfolio, program.id);
  return checkDesignAccess({isApprover: a.approver, isContributor: a.contributor, programStatus: program.status});
}

/** Why the current user may not propose or choose designs in [program]; empty when they may. */
export async function designRefusal(program: ProgramRow): Promise<string[]> {
  return accessOf(program).catch((e) => [e?.message ?? String(e)]);
}

const byRank = (a: DesignProposalRow, b: DesignProposalRow) =>
  (a.rank ?? Infinity) - (b.rank ?? Infinity) || a.label.localeCompare(b.label);

/** The program's design sets, newest first, each with its proposals. */
export async function designSetsOf(programId: string): Promise<DesignSetItem[]> {
  const by = `program_id = ${filterValue(programId)}`;
  const [sets, proposals] = await Promise.all([portalDb.designSets.query({filter: by, limit: ROW_LIMIT}),
    portalDb.designProposals.query({filter: by, limit: 5000})]);
  return sets.sort(newestFirstBy((s) => s.designed_on ?? s.created_on)).map((set) =>
    ({set, proposals: proposals.filter((p) => p.design_set_id === set.id).sort(byRank)}));
}

/** Writes a design set and its proposals, ranked by score, as the current user. [extra] carries
 * what only an import of examples sets: a fixed code and dates. */
export async function createDesignSet(program: ProgramRow, d: DesignSetInput,
  extra: {code?: string, designedOn?: string, chosen?: Map<string, string>, compounds?: Map<string, string>} = {}):
  Promise<DesignSetItem> {
  const refused = await accessOf(program);
  const series = (await workUnitsOf('series', [program.id]).catch(() => [])).map((u) => String(u.code));
  const errors = [...refused, ...validateDesignSet(d, series)];
  if (errors.length)
    throw new Error(errors.join('; '));
  const taken = (await portalDb.designSets.query({filter: `program_id = ${filterValue(program.id)}`, limit: ROW_LIMIT})).map((s) => s.code);
  const code = extra.code ?? designSetCode(d.series.trim(), new Date(), taken);
  if (taken.includes(code))
    throw new Error(`${program.code} has a design set ${code} already`);
  const me = DG.User.current().id;
  const [set] = await portalDb.designSets.insert({program_id: program.id, code, name: d.name.trim(), series_code: d.series.trim(),
    source: d.source, score_name: d.scoreName ?? undefined, higher_is_better: d.higherIsBetter !== false,
    description: d.description ?? undefined, designed_by: me, designed_on: extra.designedOn ? dayjs(extra.designedOn) : dayjs()});
  const ranks = rankProposals(d.proposals, d.higherIsBetter !== false);
  await portalDb.designProposals.insert(d.proposals.map((p) => {
    const chosenOn = extra.chosen?.get(p.label.trim());
    return {program_id: program.id, design_set_id: set.id, label: p.label.trim(), structure: p.structure.trim(),
      score: p.score ?? undefined, rank: ranks.get(p.label) ?? undefined, note: p.note ?? undefined,
      compound_code: extra.compounds?.get(p.label.trim()),
      ...(chosenOn ? {chosen: true, chosen_by: me, chosen_on: dayjs(chosenOn)} : {chosen: false})};
  }));
  return (await designSetsOf(program.id)).find((x) => x.set.id === set.id)!;
}

/** Marks proposals of one set chosen for synthesis, or not; who and when of the last change are recorded. */
export async function chooseForSynthesis(program: ProgramRow, setId: string, labels: string[], chosen: boolean): Promise<number> {
  const refused = await accessOf(program);
  if (refused.length)
    throw new Error(refused.join('; '));
  const set = await portalDb.designSets.get(setId);
  if (set.program_id !== program.id)
    throw new Error('the design set belongs to another program');
  const rows = (await portalDb.designProposals.query({filter: `design_set_id = ${filterValue(setId)}`, limit: 5000}));
  const unknown = labels.filter((l) => !rows.some((r) => r.label === l));
  if (unknown.length)
    throw new Error(`${set.code} has no proposal ${unknown.join(', ')}`);
  let changed = 0;
  for (const r of rows.filter((x) => labels.includes(x.label) && !!x.chosen !== chosen)) {
    // who and when of the last change, either way
    await portalDb.designProposals.update(r.id, {chosen, chosen_by: DG.User.current().id, chosen_on: dayjs()});
    changed++;
  }
  return changed;
}

/** The program's selection sets, newest first. */
export async function selectionSetsOf(programId: string): Promise<SelectionSetRow[]> {
  return (await portalDb.selectionSets.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT}))
    .sort(newestFirstBy((s) => s.selected_on ?? s.created_on));
}

export async function createSelectionSet(program: ProgramRow, s: SelectionSetInput,
  extra: {code?: string, selectedOn?: string} = {}): Promise<SelectionSetRow> {
  const refused = await accessOf(program);
  const campaigns = (await workUnitsOf('campaign', [program.id]).catch(() => [])).map((u) => String(u.code));
  const compounds = (await compoundsOf(program.id).catch(() => [])).map((c) => String(c.compound.code));
  const errors = [...refused, ...validateSelectionSet(s, campaigns, compounds)];
  if (errors.length)
    throw new Error(errors.join('; '));
  const taken = (await selectionSetsOf(program.id)).map((x) => x.code);
  const code = extra.code ?? designSetCode(s.campaign || 'HITS', new Date(), taken);
  if (taken.includes(code))
    throw new Error(`${program.code} has a selection set ${code} already`);
  const [row] = await portalDb.selectionSets.insert({program_id: program.id, code, name: s.name.trim(),
    campaign_code: s.campaign || undefined, rationale: s.rationale.trim(), compounds: {items: s.compounds},
    selected_by: DG.User.current().id, selected_on: extra.selectedOn ? dayjs(extra.selectedOn) : dayjs()});
  return await portalDb.selectionSets.get(row.id);
}

/** The proposals as a table: set, series, label, structure (a Molecule column), score, rank, chosen. */
export function designSetsFrame(items: DesignSetItem[], name: string): DG.DataFrame {
  const rows = items.flatMap((i) => i.proposals.map((p) => ({i, p})));
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Design set', rows.map((r) => r.i.set.name)),
    DG.Column.fromStrings('Series', rows.map((r) => r.i.set.series_code)),
    DG.Column.fromStrings('Label', rows.map((r) => r.p.label)),
    DG.Column.fromStrings('Structure', rows.map((r) => r.p.structure)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Score', rows.map((r) => r.p.score ?? null)),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Rank', rows.map((r) => r.p.rank ?? null)),
    DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'Chosen', rows.map((r) => !!r.p.chosen)),
    DG.Column.fromStrings('Compound', rows.map((r) => r.p.compound_code ?? '')),
  ]);
  df.col('Structure')!.semType = DG.SEMTYPE.MOLECULE;
  df.col('Structure')!.meta.units = 'smiles';
  df.name = name;
  return df;
}

export function selectionSetsFrame(rows: SelectionSetRow[], name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Selection set', rows.map((r) => r.name)),
    DG.Column.fromStrings('Campaign', rows.map((r) => r.campaign_code ?? '')),
    DG.Column.fromStrings('Compounds', rows.map((r) => selectionCompounds(r.compounds).join(', '))),
    DG.Column.fromStrings('Rationale', rows.map((r) => r.rationale)),
    DG.Column.fromStrings('On', rows.map((r) => r.selected_on ? dayjs(r.selected_on as any).format('YYYY-MM-DD') : '')),
  ]);
  df.name = name;
  return df;
}

interface ExampleDesignSet {
  portfolio: string; program: string; code: string; name: string; series: string; source: DesignSetInput['source'];
  scoreName?: string; designedOn?: string; description?: string;
  proposals: {label: string, structure: string, score?: number, chosenOn?: string, compound?: string}[];
}
interface ExampleSelectionSet {
  portfolio: string; program: string; code: string; name: string; campaign?: string; rationale: string; compounds: string[];
  selectedOn?: string;
}

/** Loads the synthetic design and selection sets the package ships, EX-* portfolios only; a set
 * already there (by code) is left as it is. */
export async function loadExampleDesignSets(pkg: DG.Package): Promise<string> {
  const ex: {designSets: ExampleDesignSet[], selectionSets: ExampleSelectionSet[]} =
    JSON.parse(await pkg.files.readAsText('design-sets/examples.json'));
  const lines: string[] = [];
  const programs = new Map<string, ProgramRow | null>();
  const programOf = async (portfolio: string, code: string) => {
    const key = `${portfolio}/${code}`;
    if (!programs.has(key)) {
      const pf: PortfolioRow = await portfolioByCode(portfolio);
      programs.set(key, (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(pf.id)} and code = ${filterValue(code)}`,
        limit: ROW_LIMIT}))[0] ?? null);
    }
    return programs.get(key)!;
  };
  for (const e of [...ex.designSets, ...ex.selectionSets]) {
    if (!e.portfolio.startsWith('EX-')) {
      lines.push(`${e.portfolio} ${e.program}: skipped, examples go to EX-* portfolios only`);
      continue;
    }
    try {
      const program = await programOf(e.portfolio, e.program);
      if (!program) {
        lines.push(`${e.portfolio} ${e.program}: no such program (sync the portfolio first)`);
        continue;
      }
      if ('proposals' in e) {
        if ((await designSetsOf(program.id)).some((x) => x.set.code === e.code)) {
          lines.push(`${e.program}: design set ${e.code} already there`);
          continue;
        }
        await createDesignSet(program, {name: e.name, series: e.series, source: e.source, scoreName: e.scoreName,
          description: e.description, proposals: e.proposals}, {code: e.code, designedOn: e.designedOn,
          chosen: new Map(e.proposals.filter((p) => p.chosenOn).map((p) => [p.label, p.chosenOn!])),
          compounds: new Map(e.proposals.filter((p) => p.compound).map((p) => [p.label, p.compound!]))});
        lines.push(`${e.program}: design set ${e.code} loaded (${e.proposals.length} proposals)`);
      } else {
        if ((await selectionSetsOf(program.id)).some((x) => x.code === e.code)) {
          lines.push(`${e.program}: selection set ${e.code} already there`);
          continue;
        }
        await createSelectionSet(program, {name: e.name, campaign: e.campaign, rationale: e.rationale, compounds: e.compounds},
          {code: e.code, selectedOn: e.selectedOn});
        lines.push(`${e.program}: selection set ${e.code} loaded (${e.compounds.length} compounds)`);
      }
    } catch (err: any) {
      lines.push(`${e.program}: ${e.code} not loaded: ${err?.message ?? err}`);
    }
  }
  return lines.join('\n');
}
