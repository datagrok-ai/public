// ⚠️ UNVERIFIED AUTH MECHANISM — re-check once real Benchling endpoints are reachable.
//
// This file assumes a single static API key stored as the `apiKey` package credential, injected
// as `Authorization: Bearer <apiKey>` by the commented HTTP blocks in each *Api.ts. That matches
// Benchling's v1-style API-token auth, but the swagger `swaggers/benchling_openapi.yaml` (see the
// `/token` endpoint near line 17482) describes OAuth2 with `basicClientIdSecretAuth` — i.e. a
// client-id/client-secret pair exchanged for a short-lived access token via POST /token.
//
// When wiring live HTTP, verify which scheme the target tenant actually requires:
//   - tenant API-token (static bearer): the current shape below is sufficient; just fill in
//     `apiKey` on the deployed package credentials.
//   - OAuth2 client credentials: this file must be rewritten to (a) store `clientId` + `clientSecret`
//     instead of `apiKey`, (b) POST them to `/token` to mint an access token, (c) cache the token
//     with its TTL (swagger `TokenResponse` includes `expires_in`), and (d) refresh on 401.
//
// Also reconsider the module-level `apiKey` cache: it holds the key forever until page reload even
// if the credential rotates on the server, and concurrent first calls race on `_package.getCredentials()`.
// For OAuth2 the TTL makes a time-boxed cache mandatory.
//
// Do not treat this file as load-bearing before that investigation — the whole auth layer may need
// to be re-implemented.
import { _package } from "./package";

let apiKey = '';
const API_KEY_PARAM_NAME = 'apiKey';

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