import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export interface GroupFavorites {
  group: DG.Group;
  entities: DG.Entity[];
  isAdmin: boolean;
}

/** Returns all pinned objects across the current user's groups. */
export async function getMyGroupFavorites(): Promise<GroupFavorites[]> {
  const userGroup = await grok.dapi.groups.find(DG.User.current().group.id);
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
  await Promise.all(uniqueGroups.map(async (g) => {
    try {
      const entities = await grok.dapi.entities.getFavorites(g);
      if (entities.length > 0)
        results.push({group: g, entities, isAdmin: adminIds.has(g.id)});
    }
    catch (e) {
      console.warn(`Failed to load group favorites for "${g.friendlyName}"`, e);
    }
  }));
  return results;
}

/** Returns groups the current user is an admin of (excluding personal groups). */
export async function getAdminGroups(): Promise<DG.Group[]> {
  const userGroup = await grok.dapi.groups.find(DG.User.current().group.id);
  return userGroup.adminMemberships.filter((g) => !g.personal && g.friendlyName);
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