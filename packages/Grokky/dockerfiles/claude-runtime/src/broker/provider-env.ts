export const BROKER_HOST = '127.0.0.1';
export const BROKER_PORT = Number(process.env['BROKER_PORT']) || 8377;
export const BROKER_BASE = `http://${BROKER_HOST}:${BROKER_PORT}`;

export const PLACEHOLDER_KEY = 'sk-ant-broker-placeholder';

export type BrokerMode = 'anthropic' | 'subscription' | 'bedrock' | 'foundry' | 'none';

export interface ProviderInfo {
  mode: BrokerMode;
  region?: string;
  foundryResource?: string;
  opusModel?: string;
  sonnetModel?: string;
  haikuModel?: string;
}

export function buildCliEnv(info: ProviderInfo): Record<string, string> {
  const env: Record<string, string> = {
    PATH: process.env['PATH'] ?? '/usr/local/bin:/usr/bin:/bin',
    HOME: process.env['HOME'] ?? '/home/grok',
    TERM: 'dumb',
  };

  if (info.mode === 'bedrock') {
    env['CLAUDE_CODE_USE_BEDROCK'] = '1';
    env['CLAUDE_CODE_SKIP_BEDROCK_AUTH'] = '1';
    env['ANTHROPIC_BEDROCK_BASE_URL'] = `${BROKER_BASE}/bedrock`;
    if (info.region)
      env['AWS_REGION'] = info.region;
  } else if (info.mode === 'foundry') {
    env['CLAUDE_CODE_USE_FOUNDRY'] = '1';
    env['CLAUDE_CODE_SKIP_FOUNDRY_AUTH'] = '1';
    env['ANTHROPIC_FOUNDRY_BASE_URL'] = `${BROKER_BASE}/foundry`;
    if (info.foundryResource)
      env['ANTHROPIC_FOUNDRY_RESOURCE'] = info.foundryResource;
  } else {
    env['ANTHROPIC_BASE_URL'] = BROKER_BASE;
    env['ANTHROPIC_API_KEY'] = PLACEHOLDER_KEY;
  }

  if (info.opusModel)
    env['ANTHROPIC_DEFAULT_OPUS_MODEL'] = info.opusModel;
  if (info.sonnetModel)
    env['ANTHROPIC_DEFAULT_SONNET_MODEL'] = info.sonnetModel;
  if (info.haikuModel)
    env['ANTHROPIC_DEFAULT_HAIKU_MODEL'] = info.haikuModel;

  return env;
}
