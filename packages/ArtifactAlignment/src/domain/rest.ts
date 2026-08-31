import * as grok from 'datagrok-api/grok';

/** Shares an entity (incl. promoted domain rows) with a group by NAME via the public
 * share route. Raw fetch is deliberate: the typed permissions API refuses promoted
 * rows — doc § Platform tasks, "Grant API parity for promoted domain rows". */
export async function shareEntity(entityId: string, groupName: string,
  access: 'View' | 'Edit'): Promise<void> {
  const url = `${grok.dapi.root}/public/v1/entities/${entityId}/shares` +
    `?groups=${encodeURIComponent(groupName)}&access=${access}`;
  const resp = await fetch(url, {method: 'POST', headers: {Authorization: grok.dapi.token}});
  if (!resp.ok)
    throw new Error(`Sharing ${entityId} with "${groupName}" failed: ${await resp.text()}`);
}
