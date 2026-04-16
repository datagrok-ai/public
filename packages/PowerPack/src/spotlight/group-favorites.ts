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
