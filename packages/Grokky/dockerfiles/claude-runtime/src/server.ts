import './privileges';
import {Hono} from 'hono';
import {serve} from '@hono/node-server';
import {createNodeWebSocket} from '@hono/node-ws';
import type {WSContext} from 'hono/ws';
import {syncUserFiles} from './sync/orchestrator';
import {userDirs} from './user/user-dir';
import {startWorkspaceSync} from './sync/workspace';
import {buildHelpIndex} from './help-index';
import {WORKSPACE, apiUrlFromMcpUrl, rewriteForDocker} from './constants';
import {emit} from './stream-filter';
import type {WsSender} from './stream-filter';
import {handleMessage, handleAbort, handleInputResponse, handleDisconnect} from './session';
import {claimTask, releaseTask} from './tasks';
import {disposeAll, startPool} from './agent/pool';
import {CONTAINER_TOKEN_HEADER, userCredentials} from './credentials';

const PORT = 5355;
// The life of the container token a connection brings: its credentials are not re-read after that.
const CREDENTIALS_TTL_MS = 15 * 60_000;

const app = new Hono();
const {injectWebSocket, upgradeWebSocket} = createNodeWebSocket({app});

app.get('/health', (c) => c.json({status: 'ok'}));

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

app.get('/ws', upgradeWebSocket((c) => {
  // The token comes only with the handshake. A connection whose lookup failed, or whose credentials have outlived
  // the token, is closed once idle: the browser's next prompt reconnects and the credentials are read again.
  const credentials = userCredentials(c.req.header(CONTAINER_TOKEN_HEADER));
  let stale = false;
  let turns = 0;
  let socket: WSContext | undefined;
  const closeIfStale = () => {
    if (stale && !turns)
      socket?.close(1000, 'Credentials refresh');
  };
  credentials.catch((e: Error) => {
    console.warn(`credentials: ${e.message}`);
    stale = true;
  });
  const expiry = setTimeout(() => {
    stale = true;
    closeIfStale();
  }, CREDENTIALS_TTL_MS);
  const sessionIds = new Set<string>();
  const disconnect = (_evt: unknown, ws: WsSender) => {
    clearTimeout(expiry);
    handleDisconnect(sessionIds, ws);
  };
  return {
    onOpen(_evt: Event, ws: WSContext) {
      socket = ws;
    },
    onMessage(evt: any, ws: any) {
      const sender = ws as unknown as WsSender;
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
      turns++;
      handleMessage(sender, data, credentials).catch((e: any) =>
        emit(sender, {type: 'error', sessionId: data.sessionId ?? '', message: String(e?.message ?? e)}))
        .finally(() => {
          turns--;
          closeIfStale();
        });
    },
    onClose: disconnect,
    onError: disconnect,
  };
}));

app.notFound((c) => c.json({error: 'Not found'}, 404));
app.onError((err, c) => c.json({error: String(err)}, 500));

// Survive stray rejections (e.g. abort races) instead of letting Node kill the container.
process.on('unhandledRejection', (reason) => console.error('unhandledRejection (survived):', reason));
process.on('uncaughtException', (err) => console.error('uncaughtException (survived):', err));
// Agent processes run in their own process groups, so they do not die with this one on their own.
for (const signal of ['SIGTERM', 'SIGINT'] as const) {
  process.on(signal, () => {
    disposeAll(signal);
    process.exit(0);
  });
}

const server = serve({fetch: app.fetch, port: PORT}, () => startPool());
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
