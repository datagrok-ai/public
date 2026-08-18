import * as grok from 'datagrok-api/grok';

/** Shares an entity (including promoted domain rows) with a group by NAME through
 * the platform's public share route — the same call `grok s shares add` makes.
 *
 * Raw fetch is used deliberately: the typed dapi has no surface for entity shares
 * (`grok.dapi.permissions.grant` refuses promoted domain rows), which is a known
 * platform gap — see the design doc's core task "Typed JS API for property-schema
 * grants" / classic grant API parity. */
export async function shareEntity(entityId: string, groupName: string,
  access: 'View' | 'Edit'): Promise<void> {
  const url = `${grok.dapi.root}/public/v1/entities/${entityId}/shares` +
    `?groups=${encodeURIComponent(groupName)}&access=${access}`;
  const resp = await fetch(url, {method: 'POST', headers: {Authorization: grok.dapi.token}});
  if (!resp.ok)
    throw new Error(`Sharing ${entityId} with "${groupName}" failed: ${await resp.text()}`);
}
