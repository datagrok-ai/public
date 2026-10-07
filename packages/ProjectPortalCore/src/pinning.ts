/** Pinned programs against the platform: a program in favorites is a favorite on the program's own
 * record (`DG.Favorites`, the `portal.program` row as a platform entity), the user's own or a
 * group's, so the Workspace on Home lists the program under its name and the "My programs" widget
 * gives it a tile. Nothing is stored in the portal. A favorite on the program's workspace folder (a
 * pin made before 0.3.0) still counts; removing the program removes both. The rules live in
 * ./model/pinning. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {isApprover, portfolioByCode} from './portal';
import {resolvePrincipal} from './program-roles';
import {checkGroupPin, FavoriteRef, isProgramRecordName, PinHolders, pinOf, pinsOf, whyNotPinnable} from './model/pinning';
import {filterValue, ROW_LIMIT} from './rows';

const asPin = (e: DG.Entity): FavoriteRef => ({id: e.id,
  kind: e instanceof DG.Project ? 'folder' : isProgramRecordName(e.name) ? 'program' : null});

/** The user's groups whose favorites count as theirs: every group they belong to but their own. */
async function myGroups(): Promise<DG.Group[]> {
  const mine = await grok.dapi.groups.currentUserGroups().catch(() => [] as DG.Group[]);
  const own = DG.User.current().group?.id;
  return mine.filter((g) => g.id !== own && !g.personal);
}

/** Every pin the user sees in their Workspace, by entity id (a program's id for its record, a space
 * id for a folder), and who pinned it. One call for the user's own favorites and one for all their
 * groups'. */
export async function myPins(): Promise<Map<string, PinHolders>> {
  const groups = await myGroups();
  const [personal, byGroup] = await Promise.all([
    grok.dapi.entities.getFavorites(DG.User.current().group),
    groups.length ? grok.dapi.entities.getFavoritesForGroups(groups).catch(() => ({} as {[gid: string]: DG.Entity[]})) : {},
  ]);
  return pinsOf(personal.map(asPin), groups.map((g) => ({group: g.friendlyName || g.name,
    entities: ((byGroup as {[gid: string]: DG.Entity[]})[g.id] ?? []).map(asPin)})));
}

/** The program's record as a platform entity, null when the platform does not find it; and why it
 * cannot be a favorite, when not. */
export async function pinnableRecord(program: ProgramRow): Promise<{record: DG.Entity | null, reason: string | null}> {
  const record = await grok.dapi.entities.find(program.id).catch(() => null);
  return {record: record ?? null, reason: whyNotPinnable({recordFound: !!record})};
}

/** Whether the program is pinned: by the user, and for which of their groups. */
export async function pinStateOf(program: ProgramRow): Promise<PinHolders> {
  return pinOf(await myPins(), {programId: program.id, spaceId: program.space_id}) ?? {mine: false, groups: []};
}

/** The ids among [ids] that [group]'s favorites hold (the user's own without a group). */
async function heldBy(group: DG.Group | null, ids: string[]): Promise<Set<string>> {
  const favorites = await grok.dapi.entities.getFavorites(group ?? DG.User.current().group).catch(() => [] as DG.Entity[]);
  return new Set(favorites.map((e) => e.id).filter((id) => ids.includes(id)));
}

const idsOf = (program: ProgramRow): string[] => [program.id, ...(program.space_id ? [program.space_id] : [])];

/** Adds the program to (or removes it from) the favorites of the user or, by an approver of its
 * portfolio, of a group. Adding favorites its record; removing takes out its record and, from a pin
 * made before 0.3.0, its folder. Adding twice, or removing what is not there, changes nothing. */
export async function setPinned(portfolio: PortfolioRow, program: ProgramRow, pinned: boolean, group?: DG.Group | null):
  Promise<string> {
  if (group) {
    const refused = checkGroupPin({isApprover: await isApprover(portfolio)});
    if (refused.length)
      throw new Error(refused.join('; '));
  }
  const label = `${program.code}${group ? ` for ${group.friendlyName || group.name}` : ''}`;
  // a group's favorites are read from the group itself: an approver outside it does not see them as their own
  const held = await heldBy(group ?? null, idsOf(program));
  if (pinned) {
    if (held.size)
      return `${label}: already in favorites`;
    const {record, reason} = await pinnableRecord(program);
    if (!record)
      throw new Error(`${program.code} cannot be added to favorites: ${reason}`);
    await DG.Favorites.add(record, group ?? undefined);
    return `${label}: added to favorites; it is in the Workspace and has a tile on Home`;
  }
  if (!held.size)
    return `${label}: was not in favorites`;
  for (const id of held) {
    const entity = id === program.id ? await grok.dapi.entities.find(id).catch(() => null) :
      await grok.dapi.spaces.find(id).catch(() => null);
    if (entity)
      await DG.Favorites.remove(entity, group ?? undefined);
  }
  return `${label}: removed from favorites`;
}

/** Whether [group] has [program] among its favorites: its record, or its folder. */
export async function groupHasPin(group: DG.Group, program: ProgramRow): Promise<boolean> {
  return (await heldBy(group, idsOf(program))).size > 0;
}

async function programByCode(portfolioCode: string, programCode: string): Promise<{portfolio: PortfolioRow, program: ProgramRow}> {
  const portfolio = await portfolioByCode(portfolioCode);
  const program = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)} and code = ${filterValue(programCode)}`, limit: ROW_LIMIT}))[0];
  if (!program)
    throw new Error(`${portfolioCode} has no program ${programCode}, or it is not shared with you`);
  return {portfolio, program};
}

export async function pinByCode(portfolioCode: string, programCode: string, pinned: boolean, groupName?: string | null):
  Promise<string> {
  const {portfolio, program} = await programByCode(portfolioCode, programCode);
  const group = groupName ? await resolvePrincipal(groupName) : null;
  return setPinned(portfolio, program, pinned, group);
}
