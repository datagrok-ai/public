import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/**
 * Loads the teams eligible as reviewers — excludes hidden and personal groups
 * and the built-in 'All users' group. Shared by the Share for Review dialog and
 * the package settings editor so both offer exactly the same set.
 */
export async function loadReviewerGroups(): Promise<DG.Group[]> {
  const allGroups = await grok.dapi.groups.list();
  const allUsersId = (DG.Group as any).defaultGroupsIds?.['All users'];
  return (allGroups ?? []).filter((g: any) => {
    if (g == null) return false;
    if (g.hidden) return false;
    if (g.personal) return false;
    if (allUsersId && g.id === allUsersId) return false;
    return true;
  });
}

/** Display name used as the choice key — friendlyName, falling back to name. */
export function reviewerGroupName(g: DG.Group): string {
  return (g as any).friendlyName ?? (g as any).name ?? '';
}
