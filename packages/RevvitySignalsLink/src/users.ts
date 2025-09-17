import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { RevvityUser } from "./revvity-api";

//key is revvity user Id
export let users: {[key: string]: {revvityUser: RevvityUser, datagrokUser?: DG.User}} | undefined = undefined;
export const revvityToDatagrokUsersMapping: {[key: string]: DG.User} = {};

export async function getRevvityUsersWithMapping(): Promise<{[key: string]: {revvityUser: RevvityUser, datagrokUser?: DG.User}} | undefined> {
    if (!users) {
        const usersObjWithmapping: {[key: string]: {revvityUser: RevvityUser, datagrokUser?: DG.User}} = {};
        const usersStr = await grok.functions.call('RevvitySignalsLink:getUsers');
        const revvityUsers = JSON.parse(usersStr);
        const datagrokUsers = await grok.dapi.users.list();
        Object.keys(revvityUsers).forEach((id: string) => {
            const revvityUser = revvityUsers[id];
            usersObjWithmapping[id] = {revvityUser: revvityUser};
            const datagrokUserIdx = datagrokUsers.findIndex((it) => (revvityUser.firstName && it.firstName === revvityUser.firstName) &&
                (revvityUser.lastName && it.lastName === revvityUser.lastName));
            if (datagrokUserIdx !== -1)
                usersObjWithmapping[id].datagrokUser = datagrokUsers[datagrokUserIdx];
        });
        users = usersObjWithmapping;

    }
    return users;
}

export async function createRevvityToDatagrokUsersMapping(revvityUsers: RevvityUser[]) {
    const datagrokUsers = await grok.dapi.users.list();
    revvityUsers.forEach((user) => {
        const datagrokUserIdx = datagrokUsers.findIndex((it) => (user.firstName && it.firstName === user.firstName) &&
            (user.lastName && it.lastName === user.lastName));
        if (datagrokUserIdx !== -1)
            revvityToDatagrokUsersMapping[`${user.firstName} ${user.lastName}`] = datagrokUsers[datagrokUserIdx];
    })
}

export async function getUserIdByUserString(userString: string): Promise<string> {
    //extract username
    let usernameArr = userString.split('(');
    if (usernameArr.length < 2)
        return '';
    const username = usernameArr[1].split(')')[0];
    const users = await getRevvityUsersWithMapping();
    const user = Object.values(users!).filter((it) => it.revvityUser.userName === username);
    if (!user.length || !user[0].revvityUser.userId)
        return '';
    return user[0].revvityUser.userId!;
}

export async function getUserStringIdById(id: string): Promise<string> {
    const users = await getRevvityUsersWithMapping();
    const user = Object.values(users!).filter((it) => it.revvityUser.userId === id);
    if (!user.length || !user[0].revvityUser.userId)
       return '';
    return `${user[0].revvityUser.firstName} ${user[0].revvityUser.lastName}(${user[0].revvityUser.userName})`;
}

export async function getUsersSuggestions(text: string): Promise<string[]> {
    const users = await getRevvityUsersWithMapping();
    const usersSuggestions = Object.values(users!).filter((it) => it.revvityUser.firstName?.toLowerCase().includes(text.toLowerCase()) || it.revvityUser.lastName?.toLowerCase().includes(text.toLowerCase()));
    return usersSuggestions.map((user) => `${user.revvityUser.firstName} ${user.revvityUser.lastName}(${user.revvityUser.userName})`);
}