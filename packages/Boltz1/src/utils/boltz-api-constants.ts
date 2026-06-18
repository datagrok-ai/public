export const BOLTZ_API_BASE_URL = 'https://api.boltz.bio';
export const BOLTZ_API_KEY_SETTING = 'apiKey';

export const STRUCTURE_AND_BINDING_MODEL = 'boltz-2.1';
export const ADME_MODEL = 'adme-v1';

export const TERMINAL_STATUSES: ReadonlySet<string> = new Set(['succeeded', 'failed', 'stopped']);

export interface BoltzJob {
  id: string;
  status: 'pending' | 'running' | 'succeeded' | 'failed' | 'stopped';
  error: {code: string; message: string} | null;
}
