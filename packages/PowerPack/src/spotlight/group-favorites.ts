import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


const PROJECT_PREFIX = 'Home: ';

function getProjectName(group: DG.Group): string {
  return `${PROJECT_PREFIX}${group.friendlyName}`;
}

async function getOrCreatePinnedProject(group: DG.Group): Promise<DG.Project> {
  const name = getProjectName(group);
  let project = await grok.dapi.projects.filter(`name = "${name}"`).first();
  if (!project) {
    project = DG.Project.create();
    project.name = name;
    await grok.dapi.projects.save(project);
    await grok.dapi.permissions.grant(project, group, false);
  }
  return project;
}

/** Pins an entity to a group's home favorites. */
export async function pin(entity: DG.Entity, group: DG.Group): Promise<void> {
  const project = await getOrCreatePinnedProject(group);
  project.addLink(entity);
  await grok.dapi.projects.save(project);
}

/** Unpins an entity from a group's home favorites. */
export async function unpin(entity: DG.Entity, group: DG.Group): Promise<void> {
  const name = getProjectName(group);
  const project = await grok.dapi.projects.filter(`name = "${name}"`).first();
  if (!project)
    return;
  project.removeLink(entity);
  await grok.dapi.projects.save(project);
}

/** Returns all pinned objects across the current user's groups. */
export async function getMyGroupFavorites(): Promise<{group: DG.Group; entities: DG.Entity[]}[]> {
  const userGroup = await grok.dapi.groups.find(DG.User.current().group.id);
  const groups = [...userGroup.memberships, ...userGroup.adminMemberships]
    .filter((g) => !g.personal);
  const seen = new Set<string>();
  const uniqueGroups = groups.filter((g) => {
    if (seen.has(g.id))
      return false;
    seen.add(g.id);
    return true;
  });

  const results: {group: DG.Group; entities: DG.Entity[]}[] = [];
  await Promise.all(uniqueGroups.map(async (g) => {
    const name = getProjectName(g);
    try {
      const project = await grok.dapi.projects.filter(`name = "${name}"`).first();
      if (project && project.links.length > 0) {
        const ids = project.links.map((l) => l.id);
        const resolved = await grok.dapi.getEntities(ids);
        results.push({group: g, entities: resolved.filter((e) => e != null)});
      }
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
  return userGroup.adminMemberships.filter((g) => !g.personal);
}
