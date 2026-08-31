import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {createNodeWebSocket} from '@hono/node-ws';
import {syncUserFiles} from './sync/orchestrator';
import {userDirs} from './user/user-dir';
import {startWorkspaceSync} from './sync/workspace';
import {buildHelpIndex} from './help-index';
import {WORKSPACE, rewriteForDocker} from './constants';
import {apiUrlFromMcpUrl} from './query-options';
import {BROKER_BASE} from './broker/provider-env';
import {getProviderInfo, authStart, authCode} from './broker/broker-client';
import {emit, handleMessage, handleAbort, handleInputResponse, handleDisconnect} from './session';
import {claimTask, releaseTask} from './tasks';
import type {WsSender} from './session';

const PORT = 5355;

let authOwner: WsSender | null = null;

async function handleAuthStart(ws: WsSender): Promise<void> {
  authOwner = ws;
  try {
    const r = await authStart();
    if (r.url)
      emit(ws, {type: 'auth_url', url: r.url});
    else
      emit(ws, {type: 'auth_error', message: r.error ?? 'Authentication failed to start'});
  } catch (e: any) {
    emit(ws, {type: 'auth_error', message: e.message});
  }
}

async function handleAuthCode(ws: WsSender, code: string): Promise<void> {
  try {
    const r = await authCode(code);
    if (r.status === 'done')
      emit(ws, {type: 'auth_done'});
    else
      emit(ws, {type: 'auth_error', message: r.message ?? 'Authentication failed'});
  } catch (e: any) {
    emit(ws, {type: 'auth_error', message: e.message});
  }
}

// ---------------------------------------------------------------------------
// HTTP + WebSocket transport — routes incoming ws messages to their handlers
// ---------------------------------------------------------------------------

const app = new Hono();
const {injectWebSocket, upgradeWebSocket} = createNodeWebSocket({app});

app.get('/health', async (c) => {
  try {
    if (!(await fetch(`${BROKER_BASE}/health`)).ok)
      return c.json({status: 'degraded', broker: 'unhealthy'}, 503);
  } catch {
    return c.json({status: 'degraded', broker: 'unreachable'}, 503);
  }
  return c.json({status: 'ok'});
});

// Admission-task long-poll for the queued route (see tasks.ts). The celery worker
// re-polls while the answer is {status: 'running'}; any other status ends its task.
app.post('/task/claim', async (c) => {
  const {taskId, sessionId} = await c.req.json();
  if (!taskId || !sessionId)
    return c.json({error: 'taskId and sessionId are required'}, 400);
  return c.json(await claimTask(taskId, sessionId));
});

// Task revoked platform-side: free the slot and abort the turn it was admitting.
app.post('/task/release', async (c) => {
  const {taskId} = await c.req.json();
  const sessionId = releaseTask(taskId ?? '');
  if (sessionId)
    handleDisconnect([sessionId]);
  return c.json({status: 'released'});
});

app.get('/ws', upgradeWebSocket(() => {
  const sessionIds = new Set<string>();
  let conn: WsSender | null = null;
  return {
    onMessage(evt: any, ws: any) {
      const sender = ws as unknown as WsSender;
      conn = sender;
      let data: any;
      try {
        data = JSON.parse(String(evt.data));
      } catch {
        return emit(sender, {type: 'error', sessionId: '', message: 'Invalid JSON'});
      }

      if (data.apiKey && (data.type === 'user_message' || data.type === 'sync_user_files')) {
        userDirs.ensureDir(data.apiKey, apiUrlFromMcpUrl(rewriteForDocker(data.mcpServerUrl || ''))).catch((e: any) =>
          console.warn(`user-dir pre-hook (${data.type}, session ${data.sessionId ?? '?'}):`, e.message));
      }

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
      if (authOwner === conn)
        authOwner = null;
      handleDisconnect(sessionIds);
    },
    onError() {
      if (authOwner === conn)
        authOwner = null;
      handleDisconnect(sessionIds);
    },
  };
}));

app.notFound((c) => c.json({error: 'Not found'}, 404));
app.onError((err, c) => c.json({error: String(err)}, 500));

getProviderInfo().then((info) => {
  console.log(`Claude auth: broker mode = ${info.mode}`);
  if (info.mode === 'none')
    console.warn('Claude auth: broker reports no provider configured — API calls will fail');
}).catch((e: any) => console.warn('Claude auth: could not reach broker /status:', e.message));

// ---------------------------------------------------------------------------
// Startup
// ---------------------------------------------------------------------------

// Survive stray SDK rejections (e.g. abort races) instead of letting Node kill the container.
process.on('unhandledRejection', (reason) => console.error('unhandledRejection (survived):', reason));
process.on('uncaughtException', (err) => console.error('uncaughtException (survived):', err));

const server = serve({fetch: app.fetch, port: PORT});
injectWebSocket(server);
startWorkspaceSync();
try {
  buildHelpIndex(WORKSPACE);
} catch (e: any) {
  console.warn('help-index: build failed:', e.message);
}

if (!process.env.DATAGROK_API_URL) {
  console.warn('DATAGROK_API_URL is not set: ' +
    'identity verification will trust the first client-supplied URL (dev only)');
}

console.log(`claude-runtime listening on :${PORT}`);
