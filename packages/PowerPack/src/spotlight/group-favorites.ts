import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export interface GroupFavorites {
  group: DG.Group;
  entities: DG.Entity[];
  isAdmin: boolean;
}

let _userGroupPromise: Promise<DG.Group | null> | null = null;

/** Cached lookup of the current user's group; null when it cannot be retrieved.
 *  A failed lookup is not cached, so a transient outage does not last the whole session. */
export function getCurrentUserGroup(): Promise<DG.Group | null> {
  return _userGroupPromise ??= loadCurrentUserGroup();
}

async function loadCurrentUserGroup(): Promise<DG.Group | null> {
  try {
    const group = await grok.dapi.groups.find(DG.User.current().group.id);
    if (group != null)
      return group;
  }
  catch (e) {
    console.warn('Failed to load the current user group', e);
  }
  _userGroupPromise = null;
  return null;
}

/** Returns all pinned objects across the current user's groups. */
export async function getMyGroupFavorites(): Promise<GroupFavorites[]> {
  const userGroup = await getCurrentUserGroup();
  if (userGroup == null)
    return [];
  const adminIds = new Set(userGroup.adminMemberships.map((g) => g.id));
  const groups = [...userGroup.memberships, ...userGroup.adminMemberships]
    .filter((g) => !g.personal && g.friendlyName);
  const seen = new Set<string>();
  const uniqueGroups = groups.filter((g) => {
    if (seen.has(g.id))
      return false;
    seen.add(g.id);
    return true;
  });

  const results: GroupFavorites[] = [];
  try {
    const byGroup = await grok.dapi.entities.getFavoritesForGroups(uniqueGroups);
    for (const g of uniqueGroups) {
      const entities = byGroup[g.id] ?? [];
      if (entities.length > 0)
        results.push({group: g, entities, isAdmin: adminIds.has(g.id)});
    }
  }
  catch (e) {
    console.warn('Failed to load group favorites', e);
  }
  return results;
}

/** Returns groups the current user is an admin of (excluding personal groups). */
export async function getAdminGroups(): Promise<DG.Group[]> {
  const userGroup = await getCurrentUserGroup();
  return userGroup?.adminMemberships.filter((g) => !g.personal && g.friendlyName) ?? [];
}

/** Alphabetical sort by `friendlyName`; returns a new array. */
export function sortGroupsByFriendlyName<T extends {friendlyName: string}>(groups: T[]): T[] {
  return [...groups].sort((a, b) => a.friendlyName.localeCompare(b.friendlyName));
}

/** Returns entities the current user pinned to their personal group ("Myself only"). */
export async function getMyPersonalFavorites(): Promise<DG.Entity[]> {
  try {
    return await grok.dapi.entities.getFavorites(DG.User.current().group);
  }
  catch (e) {
    console.warn('Failed to load personal favorites', e);
    return [];
  }
}

/**
 * Pins an entity to a group's favorites.
 * FileInfo entities must be registered on the server before they can be favorited,
 * so we save them first to ensure they have a persistent server-side ID.
 */
export async function pinEntityToGroup(entity: DG.Entity, group: DG.Group): Promise<void> {
  if (entity instanceof DG.FileInfo)
    await entity.save();
  await DG.Favorites.add(entity, group);
}