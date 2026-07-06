export const BOLTZ_API_BASE_URL = 'https://api.boltz.bio';
export const BOLTZ_API_CONFIG_PATH = 'System:AppData/Boltz1/boltz-api';
export const BOLTZ_API_KEY_PARAM = 'apiKey';

export const STRUCTURE_AND_BINDING_MODEL = 'boltz-2.1';
export const ADME_MODEL = 'adme-v1';

export const TERMINAL_STATUSES: ReadonlySet<string> = new Set(['succeeded', 'failed', 'stopped']);

export const BOLTZ_POLL_INTERVAL_MS = 10000;

export interface BoltzJob {
  id: string;
  status: 'pending' | 'running' | 'succeeded' | 'failed' | 'stopped';
  error: {code: string; message: string} | null;
}
