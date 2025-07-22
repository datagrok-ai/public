import { _package } from "./package";

let apiKey = '';
const API_KEY_PARAM_NAME = 'apiKey';

let apiUrl = '';
const API_URL_PARAM_NAME = 'apiUrl';

export async function getApiKey(): Promise<string> {
    if (apiKey === '') {
        const credentials = await _package.getCredentials();
        if (!credentials)
            throw new Error('API key is not set in package credentials');
        if (!credentials.parameters[API_KEY_PARAM_NAME])
            throw new Error('API key is not set in package credentials');
        apiKey = credentials.parameters[API_KEY_PARAM_NAME];
    }
    return apiKey;
}

export async function getApiUrl(): Promise<string> {
    if (apiUrl === '') {
        const credentials = await _package.getCredentials();
        if (!credentials)
            throw new Error('API URL is not set in package credentials');
        if (!credentials.parameters[API_URL_PARAM_NAME])
            throw new Error('API URL is not set in package credentials');
        apiUrl = credentials.parameters[API_URL_PARAM_NAME];
    }
    return apiUrl;
}