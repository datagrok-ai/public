/** Who sees a program row. The row is private to its author (`defaultRowVisibility: none`), and the
 * server refuses a grant on it (GROK-21041). What the server does accept is a link: a promoted
 * domain row linked into a space is visible to everyone who may view that space, and to nobody
 * else; the server clamps that path to View, so a link never lets anyone write.
 *
 * So each program row is linked into its own workspace folder. The folder already carries the
 * portfolio's groups (from the root) and the program's roles, so the row reaches exactly them.
 * A program without a folder stays private to its author.
 *
 * Pure: the caller reads the programs and what each folder links. */

export interface ProgramFolderRef {
  id: string;
  code: string;
  spaceId: string | null;
}

export interface RowLinkPlan {
  /** Programs to link into their folders. */
  link: ProgramFolderRef[];
  /** Programs linked already. */
  linked: ProgramFolderRef[];
  /** Programs that cannot be linked, with why: no folder, or the folder is gone. */
  skipped: {program: ProgramFolderRef, reason: string}[];
}

/** What to link. [links]: each existing folder's linked entity ids; a folder missing from it is
 * gone. */
export function planRowLinks(programs: ProgramFolderRef[], links: Map<string, Set<string>>): RowLinkPlan {
  const plan: RowLinkPlan = {link: [], linked: [], skipped: []};
  for (const p of programs) {
    if (!p.spaceId)
      plan.skipped.push({program: p, reason: 'no workspace folder: private to its author'});
    else if (!links.has(p.spaceId))
      plan.skipped.push({program: p, reason: 'its workspace folder is gone'});
    else if (links.get(p.spaceId)!.has(p.id))
      plan.linked.push(p);
    else
      plan.link.push(p);
  }
  return plan;
}

/** One line for the report. */
export function describeRowLinks(plan: RowLinkPlan, dryRun = false): string[] {
  const lines: string[] = [];
  if (plan.link.length)
    lines.push(`${dryRun ? 'would link' : 'linked'} ${plan.link.length} program row(s) into their folders: ` +
      plan.link.map((p) => p.code).join(', '));
  if (plan.linked.length)
    lines.push(`${plan.linked.length} already linked`);
  for (const s of plan.skipped)
    lines.push(`${s.program.code}: ${s.reason}`);
  return lines;
}
