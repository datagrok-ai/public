/** Compounds in programs (F-02, FR-03/04). A compound is a projection of the customer's registry,
 * keyed by registry and registration code: the portal never mints one. A program holds compounds
 * many-to-many, each with a role, and every role change is logged. Compounds are never folders
 * (FR-26) and never promoted entities (07 §2.3). Pure: the caller reads the feed and the portal's
 * rows, and applies the plan. */

export type CompoundRole = 'explored' | 'hit' | 'lead' | 'backup' | 'candidate' | 'tool' | 'discontinued';

export const COMPOUND_ROLES: CompoundRole[] = ['explored', 'hit', 'lead', 'backup', 'candidate', 'tool', 'discontinued'];

export interface CompoundFeedConfig {
  /** The compound registry; a compound is the same compound wherever this system and code meet. */
  system: string;
  connection: string;
  /** One row per program × compound: `program` and `compound` (the registration code) are
   * required; `role`, `name`, `smiles`, `aliases`, `batches` are mapped when present ('|' between
   * aliases and batches); anything else lands in the compound's `details`. */
  sql: string;
}

export function validateCompoundFeed(f: CompoundFeedConfig | undefined): string[] {
  if (f === undefined)
    return [];
  const errors: string[] = [];
  for (const k of ['system', 'connection', 'sql'] as const) {
    if (!f[k])
      errors.push(`compoundFeed: ${k} is required`);
  }
  return errors;
}

export interface CompoundValues {
  system: string;
  code: string;
  name: string | null;
  smiles: string | null;
  aliases: string | null;
  batches: string | null;
  details: Record<string, unknown> | null;
}

export interface CompoundRef {
  id: string;
  code: string;
  values: Partial<CompoundValues>;
}

export interface LinkRef {
  id: string;
  programId: string;
  compoundId: string;
  role: CompoundRole | null;
  source?: string | null;
  missingFromFeed?: boolean | null;
}

export interface CompoundPlan {
  portfolio: string;
  system: string;
  compoundInserts: CompoundValues[];
  compoundUpdates: {id: string, code: string, changes: Partial<CompoundValues>}[];
  /** New program × compound links; the compound by code, since a new compound has no id yet. */
  linkInserts: {programId: string, programCode: string, compoundCode: string, role: CompoundRole}[];
  roleChanges: {linkId: string, programId: string, programCode: string, compoundCode: string,
    from: CompoundRole | null, to: CompoundRole}[];
  /** Links the feed no longer delivers: flagged, kept. A link made or changed by hand is never
   * missing, and its role is not overwritten. */
  linkMissing: {id: string, key: string}[];
  linkReturned: {id: string, key: string}[];
  unchanged: number;
  rejected: {row: number, reason: string}[];
  warnings: string[];
}

export interface CompoundPlanInput {
  portfolio: string;
  system: string;
  rows: Record<string, unknown>[];
  programs: Map<string, string>;
  /** Compounds of this system the portal already holds, by code. */
  compounds: Map<string, CompoundRef>;
  /** Links under the portfolio's programs. */
  links: LinkRef[];
}

const blank = (v: unknown) => v === null || v === undefined || (typeof v === 'string' && v.trim() === '');
const text = (v: unknown) => blank(v) ? null : String(v).trim();
const list = (v: unknown) => {
  const items = blank(v) ? [] : String(v).split('|').map((s) => s.trim()).filter((s) => s.length > 0);
  return items.length ? [...new Set(items)].join('|') : null;
};
const MAPPED = new Set(['program', 'compound', 'role', 'name', 'smiles', 'aliases', 'batches']);
const same = (a: unknown, b: unknown) => JSON.stringify(a ?? null) === JSON.stringify(b ?? null);

export function planCompounds(input: CompoundPlanInput): CompoundPlan {
  const plan: CompoundPlan = {portfolio: input.portfolio, system: input.system, compoundInserts: [], compoundUpdates: [],
    linkInserts: [], roleChanges: [], linkMissing: [], linkReturned: [], unchanged: 0, rejected: [], warnings: []};
  const codeOfProgram = new Map([...input.programs.entries()].map(([c, id]) => [id, c]));
  const idToCompound = new Map([...input.compounds.values()].map((c) => [c.id, c.code]));
  const links = new Map(input.links.map((l) => [`${l.programId}/${idToCompound.get(l.compoundId) ?? l.compoundId}`, l]));
  const compoundSeen = new Map<string, CompoundValues>();
  const linkSeen = new Set<string>();

  input.rows.forEach((row, i) => {
    const program = text(row.program);
    const code = text(row.compound);
    if (!program || !code) {
      plan.rejected.push({row: i + 1, reason: 'program and compound are required'});
      return;
    }
    const programId = input.programs.get(program);
    if (!programId) {
      plan.rejected.push({row: i + 1, reason: `${program}/${code}: ${program} is not a program of ${input.portfolio}`});
      return;
    }
    const key = `${programId}/${code}`;
    if (linkSeen.has(key)) {
      plan.rejected.push({row: i + 1, reason: `${program}/${code} is delivered twice`});
      return;
    }
    linkSeen.add(key);
    let role = (text(row.role)?.toLowerCase() ?? 'explored') as CompoundRole;
    if (!COMPOUND_ROLES.includes(role)) {
      plan.warnings.push(`${program}/${code}: role "${row.role}" is not one of ${COMPOUND_ROLES.join(', ')}; taken as explored`);
      role = 'explored';
    }

    // the compound itself: one row per compound, whichever program delivers it first
    const details: Record<string, unknown> = {};
    for (const [k, v] of Object.entries(row)) {
      if (!MAPPED.has(k) && !blank(v))
        details[k] = v;
    }
    const values: CompoundValues = {system: input.system, code, name: text(row.name), smiles: text(row.smiles),
      aliases: list(row.aliases), batches: list(row.batches), details: Object.keys(details).length ? details : null};
    const first = compoundSeen.get(code);
    if (first) {
      if (!same({...first, details: null}, {...values, details: null}))
        plan.warnings.push(`${code}: two programs deliver it with different registry values; the first is kept`);
    } else {
      compoundSeen.set(code, values);
      const was = input.compounds.get(code);
      if (!was)
        plan.compoundInserts.push(values);
      else {
        const changes: Partial<CompoundValues> = {};
        for (const k of ['name', 'smiles', 'aliases', 'batches', 'details'] as const) {
          if (!same(was.values[k], values[k]))
            (changes as any)[k] = values[k];
        }
        if (Object.keys(changes).length)
          plan.compoundUpdates.push({id: was.id, code, changes});
      }
    }

    // the link
    const link = links.get(`${programId}/${code}`);
    if (!link)
      plan.linkInserts.push({programId, programCode: program, compoundCode: code, role});
    else if (link.source === 'manual') {
      // a role set by hand overrides the registry, the way a manual program_flow row overrides the rules
      if (link.role !== role)
        plan.warnings.push(`${program}/${code}: the registry says ${role}; kept ${link.role}, set by hand`);
      plan.unchanged++;
    } else {
      if (link.missingFromFeed)
        plan.linkReturned.push({id: link.id, key: `${program}/${code}`});
      if (link.role !== role)
        plan.roleChanges.push({linkId: link.id, programId, programCode: program, compoundCode: code, from: link.role, to: role});
      else if (!link.missingFromFeed)
        plan.unchanged++;
    }
  });

  for (const l of input.links) {
    const code = idToCompound.get(l.compoundId);
    if (l.source === 'feed' && !l.missingFromFeed && code && !linkSeen.has(`${l.programId}/${code}`))
      plan.linkMissing.push({id: l.id, key: `${codeOfProgram.get(l.programId) ?? '?'}/${code}`});
  }
  return plan;
}

export function summarizeCompoundPlan(p: CompoundPlan): string {
  const lines = [`${p.portfolio} compounds (${p.system}): ${p.compoundInserts.length} new, ${p.compoundUpdates.length} updated; ` +
    `links ${p.linkInserts.length} new, ${p.roleChanges.length} role changes, ${p.unchanged} unchanged, ` +
    `${p.linkMissing.length} missing from feed`];
  for (const c of p.roleChanges)
    lines.push(`  ${c.programCode}/${c.compoundCode}: ${c.from ?? '—'} → ${c.to}`);
  if (p.linkReturned.length)
    lines.push(`  back in the feed: ${p.linkReturned.map((r) => r.key).join(', ')}`);
  for (const r of p.rejected)
    lines.push(`  rejected row ${r.row}: ${r.reason}`);
  for (const w of p.warnings)
    lines.push(`  ${w}`);
  return lines.join('\n');
}
