import * as fs from 'node:fs';
import {spawn} from 'node:child_process';
import type {ChildProcess} from 'node:child_process';
import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {createNodeWebSocket} from '@hono/node-ws';
import {syncUserFiles} from './sync/orchestrator';
import {ensureUserDir} from './user/user-dir';
import {startWorkspaceSync} from './sync/workspace';
import {rewriteForDocker, apiUrlFromMcpUrl} from './query-options';
import {emit, handleMessage, handleAbort, handleInputResponse, handleDisconnect} from './session';
import type {WsSender} from './session';

const PORT = 5355;

// ---------------------------------------------------------------------------
// Subscription auth — relays the interactive `claude auth login` flow to the browser
// ---------------------------------------------------------------------------

let authProc: ChildProcess | null = null;

function handleAuthStart(ws: WsSender): void {
  if (authProc) {
    emit(ws, {type: 'auth_error', message: 'Authentication already in progress'});
    return;
  }
  authProc = spawn('claude', ['auth', 'login'], {
    env: {...process.env, TERM: 'dumb', FORCE_COLOR: '0'},
    stdio: ['pipe', 'pipe', 'pipe'],
  });

  let urlSent = false;
  const onData = (data: Buffer) => {
    if (urlSent)
      return;
    const text = data.toString();
    const m = text.match(/https:\/\/claude\.com\/cai\/oauth[^\s]+/);
    if (m) {
      urlSent = true;
      emit(ws, {type: 'auth_url', url: m[0]});
    }
  };
  authProc.stdout?.on('data', onData);
  authProc.stderr?.on('data', onData);

  authProc.on('exit', (code) => {
    authProc = null;
    if (code === 0)
      emit(ws, {type: 'auth_done'});
    else
      emit(ws, {type: 'auth_error', message: `claude auth login exited with code ${code}`});
  });

  authProc.on('error', (err: Error) => {
    authProc = null;
    emit(ws, {type: 'auth_error', message: err.message});
  });
}

function handleAuthCode(ws: WsSender, code: string): void {
  if (!authProc?.stdin) {
    emit(ws, {type: 'auth_error', message: 'No authentication in progress'});
    return;
  }
  authProc.stdin.write(code + '\n');
}

// ---------------------------------------------------------------------------
// HTTP + WebSocket transport — routes incoming ws messages to their handlers
// ---------------------------------------------------------------------------

const app = new Hono();
const {injectWebSocket, upgradeWebSocket} = createNodeWebSocket({app});

app.get('/health', (c) => c.json({status: 'ok'}));

app.get('/ws', upgradeWebSocket(() => {
  const sessionIds = new Set<string>();
  return {
    onMessage(evt: any, ws: any) {
      const sender = ws as unknown as WsSender;
      let data: any;
      try {
        data = JSON.parse(String(evt.data));
      } catch {
        return emit(sender, {type: 'error', sessionId: '', message: 'Invalid JSON'});
      }

      if (data.apiKey)
        ensureUserDir(data.apiKey).catch((e: any) => console.warn('user-dir pre-hook:', e.message));

      if (data.type === 'abort') {
        handleAbort(sender, data);
        return;
      }

      if (data.type === 'input_response') {
        handleInputResponse(sender, data);
        return;
      }

      if (data.type === 'auth_start') {
        handleAuthStart(sender);
        return;
      }

      if (data.type === 'auth_code') {
        handleAuthCode(sender, data.code ?? '');
        return;
      }

      if (data.type === 'sync_user_files') {
        const mcpUrl = rewriteForDocker(data.mcpServerUrl || '');
        const apiUrl = apiUrlFromMcpUrl(mcpUrl);
        if (apiUrl && data.apiKey) {
          const scope = data.scope || 'all';
          const packageName = data.packageName;
          console.log(`sync_user_files: scope=${scope}, packageName=${packageName ?? '<none>'}`);
          (async () => {
            try {
              const result = await syncUserFiles(apiUrl, data.apiKey, scope, packageName);
              console.log(`sync_user_files: synced ${result.files.length} file(s) (scope=${scope})`);
              emit(sender, {type: 'sync_status', status: 'done', files: result.files});
            } catch (e: any) {
              console.warn('sync_user_files failed:', e.message);
              emit(sender, {type: 'sync_status', status: 'error', message: e.message});
            }
          })();
        }
        return;
      }

      if (data.type !== 'user_message') {
        return emit(sender, {
          type: 'error', sessionId: data.sessionId ?? '',
          message: `Unknown type: ${data.type}`,
        });
      }

      sessionIds.add(data.sessionId ?? '');
      handleMessage(sender, data).catch((e: any) =>
        emit(sender, {type: 'error', sessionId: data.sessionId ?? '', message: String(e?.message ?? e)}));
    },
    onClose() {
      handleDisconnect(sessionIds);
    },
    onError() {
      handleDisconnect(sessionIds);
    },
  };
}));

app.notFound((c) => c.json({error: 'Not found'}, 404));
app.onError((err, c) => c.json({error: String(err)}, 500));

// ---------------------------------------------------------------------------
// Provider configuration — translates injected credentials into SDK env vars
// ---------------------------------------------------------------------------

// Provider config arrives as container env, forwarded from the Grokky package credentials.
// Here we translate those into the env vars the Claude Agent SDK (which wraps Claude Code) reads
// at spawn. Field-name -> SDK-env mappings below mirror Claude Code's documented provider setup:
//   Bedrock  -> CLAUDE_CODE_USE_BEDROCK + AWS_REGION + (AWS_BEARER_TOKEN_BEDROCK | AWS_* IAM creds)
//              https://code.claude.com/docs/en/amazon-bedrock
//   Foundry  -> CLAUDE_CODE_USE_FOUNDRY + ANTHROPIC_FOUNDRY_RESOURCE + (ANTHROPIC_FOUNDRY_API_KEY | Entra ID)
//              https://code.claude.com/docs/en/microsoft-foundry
//   Anthropic-> ANTHROPIC_API_KEY
// The model aliases buildOptions() passes (sonnet/opus/haiku) resolve per provider via
// ANTHROPIC_DEFAULT_{OPUS,SONNET,HAIKU}_MODEL — Bedrock needs inference-profile ids, Foundry deployment names.
// Translate the injected credential fields into the SDK provider env, collecting any
// missing-required-credential problems so they surface in the container logs at startup.
function applyProviderConfig(): void {
  const e = process.env;
  const problems: string[] = [];

  const provider = e['provider'] || 'Anthropic';
  if (provider === 'Bedrock') {
    e['CLAUDE_CODE_USE_BEDROCK'] = '1';
    if (e['region'])
      e['AWS_REGION'] = e['region'];
    if (e['awsBearerToken'])
      e['AWS_BEARER_TOKEN_BEDROCK'] = e['awsBearerToken'];
    if (e['awsAccessKeyId'])
      e['AWS_ACCESS_KEY_ID'] = e['awsAccessKeyId'];
    if (e['awsSecretAccessKey'])
      e['AWS_SECRET_ACCESS_KEY'] = e['awsSecretAccessKey'];
    if (e['awsSessionToken'])
      e['AWS_SESSION_TOKEN'] = e['awsSessionToken'];
    if (!e['awsBearerToken'] && !(e['awsAccessKeyId'] && e['awsSecretAccessKey']))
      problems.push('Bedrock selected but no credentials — set awsBearerToken, or awsAccessKeyId + awsSecretAccessKey');
  }
  else if (provider === 'Microsoft Foundry') {
    e['CLAUDE_CODE_USE_FOUNDRY'] = '1';
    if (e['foundryResource'])
      e['ANTHROPIC_FOUNDRY_RESOURCE'] = e['foundryResource'];
    if (e['foundryApiKey'])
      e['ANTHROPIC_FOUNDRY_API_KEY'] = e['foundryApiKey'];
    if (!e['foundryResource'])
      problems.push('Microsoft Foundry selected but foundryResource is missing — required to reach the endpoint');
    if (!e['foundryApiKey'])
      problems.push('Microsoft Foundry selected without foundryApiKey — falls back to Entra ID, which is not configured in this container');
  }
  else {
    if (e['apiKey'])
      e['ANTHROPIC_API_KEY'] = e['apiKey'];
  }

  if (e['opusModel'])
    e['ANTHROPIC_DEFAULT_OPUS_MODEL'] = e['opusModel'];
  if (e['sonnetModel'])
    e['ANTHROPIC_DEFAULT_SONNET_MODEL'] = e['sonnetModel'];
  if (e['haikuModel'])
    e['ANTHROPIC_DEFAULT_HAIKU_MODEL'] = e['haikuModel'];

  for (var p of problems)
    console.warn(`[provider-config] ${p}`);
}

applyProviderConfig();

const usingBedrock = process.env['CLAUDE_CODE_USE_BEDROCK'] === '1';
const usingFoundry = process.env['CLAUDE_CODE_USE_FOUNDRY'] === '1';
const hasApiKey = !!process.env['ANTHROPIC_API_KEY'];
// Subscription auth requires the host's ~/.claude/.credentials.json to be mounted into the container at this path.
const hasSubscription = fs.existsSync('/home/grok/.claude/.credentials.json');
if (usingBedrock)
  console.log('Claude auth: using Amazon Bedrock');
else if (usingFoundry)
  console.log('Claude auth: using Microsoft Foundry');
else if (hasApiKey)
  console.log('Claude auth: using ANTHROPIC_API_KEY');
else if (hasSubscription)
  console.log('Claude auth: using subscription credentials at ~/.claude/.credentials.json');
else
  console.warn('Claude auth: no provider configured (no Bedrock/Foundry/ANTHROPIC_API_KEY and no ~/.claude/.credentials.json) — API calls will fail');

// ---------------------------------------------------------------------------
// Startup
// ---------------------------------------------------------------------------

// Survive stray SDK rejections (e.g. abort races) instead of letting Node kill the container.
process.on('unhandledRejection', (reason) => console.error('unhandledRejection (survived):', reason));
process.on('uncaughtException', (err) => console.error('uncaughtException (survived):', err));

const server = serve({fetch: app.fetch, port: PORT});
injectWebSocket(server);
startWorkspaceSync();

console.log(`claude-runtime listening on :${PORT}`);
