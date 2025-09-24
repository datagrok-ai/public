import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { RevvityUser } from "./revvity-api";

//key is revvity user Id
export let users: RevvityUser[] | undefined = undefined;
export const revvityToDatagrokUsersMapping: {[key: string]: DG.User} = {};

export async function getRevvityUsers(): Promise<RevvityUser[] | undefined> {
    if (!users) {
        const usersStr = await grok.functions.call('RevvitySignalsLink:getUsers');
        users = JSON.parse(usersStr);
    }
    return users;
}


export async function getUserIdByUserString(userString: string): Promise<string> {
    //extract username
    let emailArr = userString.split('(');
    if (emailArr.length < 2)
        return '';
    const email = emailArr[1].split(')')[0];
    const users = await getRevvityUsers();
    const user = users!.filter((it) => it.email === email);
    if (!user.length || !user[0].userId)
        return '';
    return user[0].userId!;
}

export async function getUserStringIdById(id: string): Promise<string> {
    const users = await getRevvityUsers();
    const user = users!.filter((it) => it.userId === id);
    if (!user.length || !user[0].userId)
       return '';
    return `${user[0].firstName} ${user[0].lastName}(${user[0].email})`;
}

export async function getUsersSuggestions(text: string): Promise<string[]> {
    const users = await getRevvityUsers();
    const usersSuggestions = users!.filter((it) => it.firstName?.toLowerCase().includes(text.toLowerCase()) || it.lastName?.toLowerCase().includes(text.toLowerCase()));
    return usersSuggestions.map((user) => `${user.firstName} ${user.lastName}(${user.email})`);
}