import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { RevvityUser } from "./revvity-api";

export let users: {[key: string]: RevvityUser} | undefined = undefined;

export async function getRevvityUsers() {
    if (!users) {
        const usersStr = await grok.functions.call('RevvitySignalsLink:getUsers');
        users = JSON.parse(usersStr);
    }
    return users;
}