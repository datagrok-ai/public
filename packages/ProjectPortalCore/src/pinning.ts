/** Pinned programs against the platform: a pin is a favorite on the program's workspace folder
 * (`DG.Favorites`), the user's own or a group's, so the Workspace on Home lists the folder and the
 * "My programs" widget lists the program. Nothing is stored in the portal. The rules live in
 * ./model/pinning. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {isApprover, portfolioByCode} from './portal';
import {existingSpaces} from './workspace';
import {resolvePrincipal} from './program-roles';
import {checkGroupPin, PinHolders, pinsOf, whyNotPinnable} from './model/pinning';
import {filterValue, ROW_LIMIT} from './rows';

const asPin = (e: DG.Entity) => ({id: e.id, isSpace: e instanceof DG.Project});

/** The user's groups whose favorites count as theirs: every group they belong to but their own. */
async function myGroups(): Promise<DG.Group[]> {
  const mine = await grok.dapi.groups.currentUserGroups().catch(() => [] as DG.Group[]);
  const own = DG.User.current().group?.id;
  return mine.filter((g) => g.id !== own && !g.personal);
}

/** Every pinned folder the user sees in their Workspace, and who pinned it. One call for the user's
 * own favorites and one for all their groups'. */
export async function myPins(): Promise<Map<string, PinHolders>> {
  const groups = await myGroups();
  const [personal, byGroup] = await Promise.all([
    grok.dapi.entities.getFavorites(DG.User.current().group),
    groups.length ? grok.dapi.entities.getFavoritesForGroups(groups).catch(() => ({} as {[gid: string]: DG.Entity[]})) : {},
  ]);
  return pinsOf(personal.map(asPin), groups.map((g) => ({group: g.friendlyName || g.name,
    entities: ((byGroup as {[gid: string]: DG.Entity[]})[g.id] ?? []).map(asPin)})));
}

/** The program's folder, null when there is none to pin; and why not, when not. */
export async function pinnableFolder(program: ProgramRow): Promise<{space: DG.Project | null, reason: string | null}> {
  const live = program.space_id ? await existingSpaces([program.space_id]).catch(() => new Set<string>()) : new Set<string>();
  const reason = whyNotPinnable({spaceId: program.space_id, spaceLive: !!program.space_id && live.has(program.space_id)});
  if (reason)
    return {space: null, reason};
  const space = await grok.dapi.spaces.find(program.space_id!).catch(() => null);
  return space ? {space, reason: null} : {space: null, reason: whyNotPinnable({spaceId: program.space_id, spaceLive: false})};
}

/** Whether the program is pinned: by the user, and for which of their groups. */
export async function pinStateOf(program: ProgramRow): Promise<PinHolders> {
  const held = program.space_id ? (await myPins()).get(program.space_id) : undefined;
  return held ?? {mine: false, groups: []};
}

/** Pins (or unpins) the program's folder, for the user or, by an approver of its portfolio, for a
 * group. Pinning twice, or unpinning what is not pinned, changes nothing. */
export async function setPinned(portfolio: PortfolioRow, program: ProgramRow, pinned: boolean, group?: DG.Group | null):
  Promise<string> {
  if (group) {
    const refused = checkGroupPin({isApprover: await isApprover(portfolio)});
    if (refused.length)
      throw new Error(refused.join('; '));
  }
  const {space, reason} = await pinnableFolder(program);
  const label = `${program.code}${group ? ` for ${group.friendlyName || group.name}` : ''}`;
  if (!space) {
    if (!pinned)
      return `${label}: not pinned (no folder)`;
    throw new Error(`${program.code} cannot be pinned: ${reason}`);
  }
  // a group's pins read from the group itself: an approver outside it does not see them as their own
  const was = group ? await groupHasPin(group, space.id) : (await pinStateOf(program)).mine;
  if (was === pinned)
    return `${label}: ${pinned ? 'already pinned' : 'was not pinned'}`;
  if (pinned)
    await DG.Favorites.add(space, group ?? undefined);
  else
    await DG.Favorites.remove(space, group ?? undefined);
  return `${label}: ${pinned ? 'pinned; its folder is in the Workspace and its tile on Home' : 'unpinned'}`;
}

/** Whether [group] has [spaceId] among its favorites. */
export async function groupHasPin(group: DG.Group, spaceId: string): Promise<boolean> {
  return (await grok.dapi.entities.getFavorites(group).catch(() => [] as DG.Entity[])).some((e) => e.id === spaceId);
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
