export const WORKSPACE = process.env['CLAUDE_WORKSPACE'] || '/workspace';

export function rewriteForDocker(url: string): string {
  return url.replace('//localhost', '//host.docker.internal');
}

export function apiUrlFromMcpUrl(mcpUrl: string): string | undefined {
  const idx = mcpUrl.indexOf('/docker/containers/proxy/');
  return idx > 0 ? mcpUrl.substring(0, idx) : undefined;
}

export const PLUGIN_DIR = process.env['GROKKY_PLUGIN_DIR'] || '/app/plugin';
export const AGENT_DIR = process.env['GROKKY_AGENT_DIR'] || '/app/agent';
export const AGENT_PYTHON = process.env['GROKKY_AGENT_PYTHON'] || '/opt/hermes/.venv/bin/python';
// Not under /tmp, which an agent may clear.
export const HERMES_HOMES = process.env['GROKKY_HERMES_HOMES'] || '/var/lib/grokky/hermes';

/** The keys of Grokky's package credentials handed to the agent's provider (agent/grokky_agent/providers). */
export const PROVIDER_CREDENTIALS = ['provider', 'apiKey', 'baseUrl', 'region', 'apiFormat', 'model', 'service',
  'TEST_UPSTREAM'];
