/** The platform's own notifications (the bell): the one route the portal tells people by. Nothing
 * external. There is no JS API to create a notification (PLAT-18), so this posts to the server's
 * endpoint, which lets a user notify only themselves unless they may edit users. */
import * as grok from 'datagrok-api/grok';

/** Posts [text] to [userId]'s notifications; a word on how it went, never a throw. */
export async function notifyUser(userId: string, text: string): Promise<{ok: boolean, detail: string}> {
  try {
    const response = await fetch(`${grok.dapi.root}/users/notifications`, {
      method: 'POST', credentials: 'include',
      headers: {'Content-Type': 'application/json', 'Authorization': grok.dapi.token},
      body: JSON.stringify({user: {id: userId}, text, data: '', sender: 'Project Portal'}),
    });
    return response.ok ? {ok: true, detail: 'notified'} :
      {ok: false, detail: `not notified (${response.status} ${(await response.text()).slice(0, 120)})`};
  } catch (e: any) {
    return {ok: false, detail: `not notified (${e?.message ?? e})`};
  }
}
