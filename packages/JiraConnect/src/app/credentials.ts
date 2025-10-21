
import * as constants from './constants';
import { _package } from '../package';

interface JiraCredsData {
    authKey: string;
    userName: string;
    host: string;
}

let credentialsPromise: Promise<any> | null = null;
let credentials: any | null = null;

export async function getJiraCreds(): Promise<JiraCredsData | null> {
    if (!credentialsPromise)
        credentialsPromise = _package.getCredentials();
    credentials ??= await credentialsPromise;
    if (!credentials)
        return null;
    let credentialsParams = credentials.parameters;
    let authKey = credentialsParams[constants.AUTH_KEY];
    let userName = credentialsParams[constants.USERNAME];
    let host = credentialsParams[constants.HOST_NAME];
    if (!(authKey && userName && host))
        return null;
    return { authKey: authKey, userName: userName, host: host };
}