import {BROKER_BASE, BrokerMode, ProviderInfo} from './provider-env';

interface StatusResponse {
  mode: BrokerMode;
  authRequired?: boolean;
  region?: string;
  foundryResource?: string;
  models?: {opus?: string; sonnet?: string; haiku?: string};
}

export type ProviderStatus = ProviderInfo & {authRequired: boolean};

async function fetchStatus(): Promise<StatusResponse | undefined> {
  try {
    const resp = await fetch(`${BROKER_BASE}/status`);
    return resp.ok ? await resp.json() as StatusResponse : undefined;
  } catch {
    return undefined;
  }
}

export async function getProviderInfo(): Promise<ProviderStatus> {
  const s = await fetchStatus();
  return {
    mode: s?.mode ?? 'none',
    authRequired: s?.authRequired ?? false,
    region: s?.region,
    foundryResource: s?.foundryResource,
    opusModel: s?.models?.opus,
    sonnetModel: s?.models?.sonnet,
    haikuModel: s?.models?.haiku,
  };
}

export async function registerMcpSession(
  targetUrl: string, apiKey?: string, apiUrl?: string,
): Promise<string | undefined> {
  try {
    const resp = await fetch(`${BROKER_BASE}/mcp-session`, {
      method: 'POST',
      headers: {'content-type': 'application/json'},
      body: JSON.stringify({targetUrl, apiKey, apiUrl}),
    });
    if (!resp.ok)
      return undefined;
    const {token} = await resp.json() as {token: string};
    return `${BROKER_BASE}/mcp/${token}`;
  } catch {
    return undefined;
  }
}

export async function authStart(): Promise<{url?: string; error?: string}> {
  const resp = await fetch(`${BROKER_BASE}/auth/start`, {method: 'POST'});
  return await resp.json() as {url?: string; error?: string};
}

export async function authCode(code: string): Promise<{status: string; message?: string}> {
  const resp = await fetch(`${BROKER_BASE}/auth/code`, {
    method: 'POST', headers: {'content-type': 'application/json'}, body: JSON.stringify({code}),
  });
  return await resp.json() as {status: string; message?: string};
}
