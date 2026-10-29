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

/** Container credentials handed to the agent's provider (agent/grokky_agent/providers). */
export const PROVIDER_CREDENTIALS = ['provider', 'apiKey', 'baseUrl', 'region', 'apiFormat', 'model', 'service',
  'TEST_UPSTREAM'];

/** Removes the named variables from `env` and returns them; an unset or empty one is left out. */
export function takeEnv(names: string[], env: NodeJS.ProcessEnv = process.env): Record<string, string> {
  const taken = Object.fromEntries(names.filter((k) => env[k]).map((k) => [k, env[k]!]));
  for (const k of names)
    delete env[k];
  return taken;
}

// Out of the environment at startup: every child the server runs (git for the workspace sync) would inherit them,
// and a child is dumpable again after its exec.
export const PROVIDER_ENV = takeEnv(PROVIDER_CREDENTIALS);
