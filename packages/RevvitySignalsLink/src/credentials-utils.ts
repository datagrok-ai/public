import { API_KEY_PARAM_NAME, API_URL_PARAM_NAME } from "./constants";
import { _package } from "./package";

const cache = new Map<string, string>();

async function getCredential(paramName: string, label: string): Promise<string> {
    const cached = cache.get(paramName);
    if (cached !== undefined)
        return cached;
    const credentials = await _package.getCredentials();
    if (!credentials || !credentials.parameters[paramName])
        throw new Error(`${label} is not set in package credentials`);
    const value = credentials.parameters[paramName];
    cache.set(paramName, value);
    return value;
}

export function getApiKey(): Promise<string> {
    return getCredential(API_KEY_PARAM_NAME, 'API key');
}

export function getApiUrl(): Promise<string> {
    return getCredential(API_URL_PARAM_NAME, 'API URL');
}