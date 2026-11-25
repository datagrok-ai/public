import {spawnSync} from 'node:child_process';
import {randomBytes} from 'node:crypto';
import {existsSync, mkdirSync, readFileSync, rmSync, writeFileSync} from 'node:fs';
import net from 'node:net';
import {homedir} from 'node:os';
import {dirname, join, relative, resolve} from 'node:path';
import {fileURLToPath} from 'node:url';

const USAGE = `The SQL Server of identifiers-catalogs.feature: one Docker container, any machine with Docker and Node.

  node catalogs-db.mjs up [--reload] [--server <address>] [--port <port>]
      creates the container if it is absent (starts it if it is stopped), waits until SQL Server answers,
      creates or repairs the login catalogs_demo, loads catalogs-fixture.sql when its databases are missing
      (or with --reload: its three working catalogs are made again), gives them to the login, and writes the
      connection with the login's password to .catalogs-db/catalogs-connection.json. The passwords are kept
      in .catalogs-db/secrets.json; the folder is git-ignored.
      --server/--port: the address the stand's Grok Connect reaches the database by (kept for later runs;
      default host.docker.internal:14330, which is Docker Desktop's)
  node catalogs-db.mjs connection [--host <alias>] [--server <address>] [--port <port>]
      saves the CatalogsDemo connection on a stand (grok s connections save --save-credentials; the
      running account's own connection of that name is updated) and tests it there. The host is an
      alias of ~/.grok/config.yaml, by default its default one; never "dev"
  node catalogs-db.mjs status
      the container, the volume, the login, the fixture's databases, the port, the local files
  node catalogs-db.mjs down
      stops the container; it, its volume and every database on it are kept
  node catalogs-db.mjs down --purge --yes
      removes the container AND its volume: every database on the server goes, also the ones other
      logins keep there, not only this fixture's. Without --yes it says what would go and stops`;

const HERE = dirname(fileURLToPath(import.meta.url));
const CONTAINER = 'dg-catalogs-db';
const IMAGE = 'mcr.microsoft.com/mssql/server:2022-CU27-GDR1-ubuntu-22.04';
const VOLUME = 'dg-catalogs-db-data';
const LABEL_KEY = 'ai.datagrok.fixture';
const LABEL_VALUE = 'catalogs-db';
const HOST_PORT = 14330;
const LOGIN = 'catalogs_demo';
const CONNECTION = 'CatalogsDemo';
const SQLCMD = '/opt/mssql-tools18/bin/sqlcmd';
const FIXTURE = join(HERE, 'catalogs-fixture.sql');
const TEMPLATE = join(HERE, 'catalogs-connection.json');
const LOCAL = join(HERE, '.catalogs-db');
const SECRETS = join(LOCAL, 'secrets.json');
const LOCAL_CONNECTION = join(LOCAL, 'catalogs-connection.json');
const GROK_JS = resolve(HERE, '../../../../../tools/bin/grok.js');
// the fixture's databases, and the only ones this script ever drops or changes
const OWN = ['research', 'archive', 'sales',
  ...Array.from({length: 200}, (_, i) => `tenant_${String(i + 1).padStart(3, '0')}`)];
const OWN_LIST = OWN.map((n) => `N'${n}'`).join(', ');
// what the login sees: its own databases, master and tempdb
const VISIBLE = OWN.length + 2;
const MINUTE = 60000;
const shown = (path) => relative(process.cwd(), path) || path;

class Stop extends Error {}
const stop = (message) => {
  throw new Stop(message);
};
const delay = (ms) => new Promise((r) => setTimeout(r, ms));

function run(cmd, args, {input, timeoutMs = MINUTE, env, inherit = false} = {}) {
  const r = spawnSync(cmd, args, {input, timeout: timeoutMs, env: env ? {...process.env, ...env} : process.env,
    encoding: 'utf8', maxBuffer: 64 * 1024 * 1024, windowsHide: true,
    stdio: inherit ? ['ignore', 'inherit', 'inherit'] : 'pipe'});
  if (r.error?.code === 'ENOENT')
    stop(`${cmd} was not found on PATH`);
  const timedOut = r.error?.code === 'ETIMEDOUT';
  return {ok: !r.error && r.status === 0, status: r.status, timedOut, stdout: r.stdout ?? '',
    stderr: (r.stderr ?? '') + (timedOut ? `\n(${cmd} ${args[0]} gave no answer in ${timeoutMs / 1000} s and was stopped)` : '')};
}

const docker = (args, options) => run('docker', args, options);

function must(r, what) {
  if (!r.ok)
    stop(`${what} failed: ${(r.stderr || r.stdout).trim().slice(0, 1500)}`);
  return r.stdout;
}

const alnum = (n) => {
  const chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
  return [...randomBytes(n)].map((b) => chars[b % chars.length]).join('');
};

function readSecrets() {
  try {
    return JSON.parse(readFileSync(SECRETS, 'utf8'));
  }
  catch {
    return {};
  }
}

function writeSecrets(patch) {
  mkdirSync(LOCAL, {recursive: true});
  writeFileSync(SECRETS, JSON.stringify({...readSecrets(), ...patch}, null, 2) + '\n', {mode: 0o600});
}

// --- Docker -------------------------------------------------------------------------------------

function container() {
  const format = `{{.State.Status}}|{{index .Config.Labels "${LABEL_KEY}"}}|{{.Config.Image}}|` +
    '{{json .HostConfig.PortBindings}}|{{.HostConfig.RestartPolicy.Name}}|{{.HostConfig.Memory}}';
  const r = docker(['container', 'inspect', '--format', format, CONTAINER]);
  if (!r.ok) {
    if (/no such (container|object)/i.test(r.stderr))
      return null;
    stop(`docker container inspect failed: ${r.stderr.trim()}`);
  }
  const [status, label, image, ports, restart, memory] = r.stdout.trim().split('|');
  const bindings = JSON.parse(ports || '{}')?.['1433/tcp'] ?? [];
  const differences = [
    image !== IMAGE && `image ${image}, not ${IMAGE}`,
    !(bindings.length === 1 && bindings[0].HostIp === '127.0.0.1' && bindings[0].HostPort === `${HOST_PORT}`) &&
      `port ${ports}, not 127.0.0.1:${HOST_PORT} -> 1433`,
    restart !== 'unless-stopped' && `restart policy ${restart}, not unless-stopped`,
    memory !== `${8 * 1024 ** 3}` && `memory limit ${memory} bytes, not 8g`,
  ].filter(Boolean);
  return {status, ours: label === LABEL_VALUE, differences};
}

function ownContainer() {
  const c = container();
  if (c && !c.ours)
    stop(`a container named ${CONTAINER} exists without the label ${LABEL_KEY}=${LABEL_VALUE}: it is not this ` +
      'fixture\'s, and the script leaves it alone');
  return c;
}

const volumeExists = () => docker(['volume', 'inspect', VOLUME]).ok;

/** The sa password is in the container's environment: it is read, never printed or put on a command line. */
function saPassword() {
  const r = docker(['container', 'inspect', '--format', '{{range .Config.Env}}{{println .}}{{end}}', CONTAINER]);
  if (!r.ok)
    stop(`reading the container's environment failed: ${r.stderr.trim()}`);
  return r.stdout.split(/\r?\n/).find((l) => l.startsWith('MSSQL_SA_PASSWORD='))?.slice('MSSQL_SA_PASSWORD='.length);
}

function createContainer(password) {
  if (!docker(['image', 'inspect', IMAGE]).ok) {
    console.log(`pulling ${IMAGE}`);
    must(docker(['pull', IMAGE], {timeoutMs: 30 * MINUTE, inherit: true}), `docker pull ${IMAGE}`);
  }
  // `-e MSSQL_SA_PASSWORD` without a value: docker takes it from its own environment, off the command line
  must(docker(['run', '-d', '--name', CONTAINER, '--label', `${LABEL_KEY}=${LABEL_VALUE}`,
    '--restart', 'unless-stopped', '--memory', '8g', '-p', `127.0.0.1:${HOST_PORT}:1433`,
    '-v', `${VOLUME}:/var/opt/mssql`, '-e', 'ACCEPT_EULA=Y', '-e', 'MSSQL_PID=Developer',
    '-e', 'MSSQL_MEMORY_LIMIT_MB=6144', '-e', 'MSSQL_SA_PASSWORD', IMAGE],
  {timeoutMs: 2 * MINUTE, env: {MSSQL_SA_PASSWORD: password}}), 'docker run');
}

// --- SQL ----------------------------------------------------------------------------------------

/** Runs a script as sa inside the container; the script comes on stdin, the password from the container. */
function sqlAsSa(script, {timeoutMs = MINUTE, vars = {}} = {}) {
  const v = Object.entries(vars).map(([k, x]) => ` -v ${k}=${x}`).join('');
  return docker(['exec', '-i', CONTAINER, 'sh', '-c',
    `${SQLCMD} -S localhost -U sa -P "$MSSQL_SA_PASSWORD" -C -b -l 10 -W -h -1 -s "|"${v} -i /dev/stdin`],
  {input: script, timeoutMs});
}

/** Runs a script as the fixture's login; its password goes in through docker's environment (SQLCMDPASSWORD). */
function sqlAsLogin(script, password) {
  return docker(['exec', '-i', '-e', 'SQLCMDPASSWORD', CONTAINER, SQLCMD, '-S', 'localhost', '-U', LOGIN, '-C', '-b',
    '-l', '10', '-W', '-h', '-1', '-s', '|', '-i', '/dev/stdin'], {input: script, env: {SQLCMDPASSWORD: password}});
}

const rows = (out) => out.split(/\r?\n/).map((l) => l.trim()).filter(Boolean).map((l) => l.split('|'));

/** A read asked again when it gave no answer in time: the first one after two hundred new databases once did not. */
function read(ask) {
  let r = ask();
  for (let attempt = 2; r.timedOut && attempt <= 3; attempt++) {
    console.log(`(no answer in time; asking again, ${attempt} of 3)`);
    r = ask();
  }
  return r;
}

async function waitForSql(deadlineMs = 4 * MINUTE) {
  const until = Date.now() + deadlineMs;
  let refused = 0;
  let last = '';
  while (Date.now() < until) {
    const c = container();
    if (!c || c.status === 'exited' || c.status === 'dead') {
      const logs = docker(['logs', '--tail', '20', CONTAINER]);
      stop(`the container stopped while SQL Server was starting:\n${(logs.stdout + logs.stderr).trim()}`);
    }
    const r = sqlAsSa('SELECT 1;', {timeoutMs: 30000});
    if (r.ok)
      return;
    last = (r.stderr || r.stdout).trim();
    // a starting server refuses with a reason ("script upgrade mode"); a wrong password has none
    refused = /Login failed for user 'sa'\.(?![^\n]*Reason)/.test(last) ? refused + 1 : 0;
    if (refused >= 5)
      stop('SQL Server refuses sa with the password in the container\'s environment: the volume was made with ' +
        'another one. Recreate the container with that password, or see `down --purge`.');
    await delay(2000);
  }
  stop(`SQL Server did not answer in ${deadlineMs / 1000} s: ${last.slice(0, 500)}`);
}

/** After a start, the server answers before its databases have recovered: the fixture's are waited for. */
async function waitForOnline(deadlineMs = 4 * MINUTE) {
  const until = Date.now() + deadlineMs;
  let left = '';
  while (Date.now() < until) {
    const r = sqlAsSa(`SET NOCOUNT ON; SELECT COUNT(*) FROM sys.databases WHERE name IN (${OWN_LIST}) AND state <> 0;`);
    left = r.ok ? rows(r.stdout)[0]?.[0] : (r.stderr || r.stdout).trim().slice(0, 300);
    if (left === '0')
      return;
    await delay(2000);
  }
  console.warn(`warning: the fixture's databases are not all online after ${deadlineMs / 1000} s (${left})`);
}

/** name -> owner of the fixture's databases that exist */
function fixtureDatabases() {
  const out = must(read(() => sqlAsSa(`SET NOCOUNT ON; SELECT name, ISNULL(SUSER_SNAME(owner_sid), '?') ` +
    `FROM sys.databases WHERE name IN (${OWN_LIST});`)), 'listing the fixture\'s databases');
  return new Map(rows(out).map(([name, owner]) => [name, owner]));
}

/** The rows the feature and the tests read, counted: an interrupted load leaves the three working catalogs short. */
function fixtureComplete() {
  const r = read(() => sqlAsSa('SET NOCOUNT ON; SELECT (SELECT COUNT(*) FROM research.chem.compounds), ' +
    '(SELECT COUNT(*) FROM research.bio.assays), (SELECT COUNT(*) FROM archive.chem.compounds), ' +
    '(SELECT COUNT(*) FROM sales.dbo.orders);'));
  return r.ok && rows(r.stdout)[0]?.join(',') === '12,4,12,4';
}

/** Databases of the fixture's names that another login owns are someone else's: nothing goes on. */
function stopOnOthers(databases) {
  const others = [...databases].filter(([, owner]) => owner !== LOGIN && owner !== 'sa');
  if (others.length)
    stop(`${others.map(([n, o]) => `${n} (owner ${o})`).join(', ')}: of the fixture's names, but owned by another ` +
      'login. The script does not touch them; nothing was changed.');
}

/** Creates the login or sets its password; the password is written to the git-ignored file before the server has it. */
function ensureLogin() {
  const existed = rows(must(read(() => sqlAsSa(`SET NOCOUNT ON; ` +
    `SELECT CASE WHEN SUSER_ID(N'${LOGIN}') IS NULL THEN 0 ELSE 1 END;`)), 'looking up the login'))[0]?.[0] === '1';
  let password = readSecrets().login;
  const fresh = !password;
  if (fresh) {
    password = alnum(32);
    writeSecrets({login: password});
  }
  must(sqlAsSa(`SET NOCOUNT ON;
IF SUSER_ID(N'${LOGIN}') IS NULL
  CREATE LOGIN [${LOGIN}] WITH PASSWORD = N'${password}', CHECK_EXPIRATION = OFF, CHECK_POLICY = OFF;
ELSE
  ALTER LOGIN [${LOGIN}] WITH PASSWORD = N'${password}', CHECK_EXPIRATION = OFF, CHECK_POLICY = OFF;
ALTER LOGIN [${LOGIN}] ENABLE;
DENY VIEW ANY DATABASE TO [${LOGIN}];
`), `setting up the login ${LOGIN}`);
  if (!existed)
    console.log(`login ${LOGIN}: created`);
  else if (fresh)
    console.log(`login ${LOGIN}: its password was not in ${shown(SECRETS)}, so it was reset; save the connection on ` +
      'the stand again (`connection`)');
  else
    console.log(`login ${LOGIN}: ok`);
  return password;
}

/** What the login sees, as Grok Connect will: [all catalogs, the fixture's] */
function visibleToLogin(password) {
  const r = read(() => sqlAsLogin(`SET NOCOUNT ON; SELECT COUNT(*), ` +
    `SUM(CASE WHEN name IN (${OWN_LIST}) THEN 1 ELSE 0 END) FROM sys.databases;`, password));
  if (!r.ok)
    return null;
  return rows(r.stdout)[0]?.map(Number) ?? null;
}

// --- the host side ------------------------------------------------------------------------------

/** A TDS pre-login sent to the published port: a wedged Docker port forward accepts the connection and never answers. */
function portAnswers(timeoutMs = 5000) {
  const prelogin = Buffer.from([0x12, 0x01, 0x00, 0x1A, 0x00, 0x00, 0x01, 0x00,
    0x00, 0x00, 0x0B, 0x00, 0x06, 0x01, 0x00, 0x11, 0x00, 0x01, 0xFF,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02]);
  return new Promise((done) => {
    const socket = net.connect({host: '127.0.0.1', port: HOST_PORT});
    const end = (answer) => {
      socket.destroy();
      done(answer);
    };
    socket.setTimeout(timeoutMs, () => end('connected, but no answer'));
    socket.on('error', (e) => end(e.code ?? e.message));
    socket.on('connect', () => socket.write(prelogin));
    socket.on('data', (d) => end(d[0] === 0x04 ? 'ok' : 'an answer that is not SQL Server\'s'));
  });
}

function readJson(path) {
  try {
    return JSON.parse(readFileSync(path, 'utf8'));
  }
  catch {
    return null;
  }
}

/** The tracked connection with the login's password; the address is the one given, else the one kept, else the template's. */
function writeConnection(password, opts) {
  const conn = readJson(TEMPLATE);
  if (!conn)
    stop(`${TEMPLATE} is missing or not JSON`);
  const kept = readJson(LOCAL_CONNECTION)?.parameters ?? {};
  conn.parameters.server = opts.server ?? kept.server ?? conn.parameters.server;
  conn.parameters.port = Number(opts.port ?? kept.port ?? conn.parameters.port);
  conn.credentials.parameters = {login: LOGIN, password};
  mkdirSync(LOCAL, {recursive: true});
  writeFileSync(LOCAL_CONNECTION, JSON.stringify(conn, null, 2) + '\n', {mode: 0o600});
  return conn;
}

function grokConfig() {
  try {
    return readFileSync(join(homedir(), '.grok', 'config.yaml'), 'utf8');
  }
  catch {
    return '';
  }
}

function grokHost(requested) {
  const config = grokConfig();
  const host = requested ?? /^default:\s*['"]?([^'"\s#]+)/m.exec(config)?.[1];
  if (!host)
    stop('no --host given and no default server in ~/.grok/config.yaml');
  const escaped = host.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const url = new RegExp(`^ {2}${escaped}:\\s*\\r?\\n(?: {4}.*\\r?\\n)*? {4}url:\\s*(\\S+)`, 'm').exec(config)?.[1] ?? host;
  if (host.toLowerCase() === 'dev' || /(^|\/\/|\.)dev\.datagrok\.ai/i.test(url))
    stop(`refusing the server "${host}" (${url}): the fixture is for a local stand`);
  return host;
}

function grok(args, host) {
  if (!existsSync(GROK_JS))
    stop(`${GROK_JS} is missing: run the script from a checkout of public with \`grok setup\` done`);
  return run(process.execPath, [GROK_JS, 's', ...args, '--host', host], {timeoutMs: 2 * MINUTE});
}

const parseJson = (text, what) => {
  try {
    return JSON.parse(text);
  }
  catch {
    stop(`${what}: not JSON: ${text.trim().slice(0, 500)}`);
  }
};

// --- commands -----------------------------------------------------------------------------------

async function up(opts) {
  let c = ownContainer();
  if (!c) {
    let password;
    if (volumeExists()) {
      password = readSecrets().sa;
      if (!password)
        stop(`the volume ${VOLUME} exists but no container uses it, and no sa password is cached in ` +
          `${shown(SECRETS)}: the password SQL Server keeps in the volume is unknown. Put the one the volume was ` +
          'made with into that file as {"sa": "..."} and run `up` again, or, if nothing on the server is needed, ' +
          '`down --purge --yes` deletes the volume with every database on it.');
      console.log(`creating the container ${CONTAINER} on the existing volume ${VOLUME} (sa password from the cache)`);
    }
    else {
      password = `Aa1${alnum(29)}`;
      writeSecrets({sa: password});
      console.log(`creating the container ${CONTAINER} and its volume ${VOLUME}`);
    }
    createContainer(password);
    c = ownContainer();
  }
  else {
    for (const d of c.differences)
      console.warn(`warning: the container differs from the contract: ${d} (the script does not change it)`);
    if (c.status !== 'running') {
      console.log(`starting the container ${CONTAINER} (${c.status})`);
      must(docker(['start', CONTAINER], {timeoutMs: 2 * MINUTE}), 'docker start');
    }
  }
  await waitForSql();
  await waitForOnline();
  const sa = saPassword();
  if (sa && sa !== readSecrets().sa)
    writeSecrets({sa});
  console.log('SQL Server answers');

  const password = ensureLogin();
  let databases = fixtureDatabases();
  stopOnOthers(databases);
  const missing = OWN.filter((n) => !databases.has(n));
  const complete = missing.length === 0 && fixtureComplete();
  if (opts.reload || !complete) {
    console.log(opts.reload ? 'loading catalogs-fixture.sql (--reload)' : missing.length ?
      `loading catalogs-fixture.sql: ${missing.length} of its ${OWN.length} databases are missing` :
      'loading catalogs-fixture.sql: its working catalogs are incomplete');
    const r = sqlAsSa(readFileSync(FIXTURE, 'utf8'), {timeoutMs: 20 * MINUTE, vars: {owner: LOGIN}});
    must(r, 'catalogs-fixture.sql');
    databases = fixtureDatabases();
  }
  else
    console.log(`the fixture's ${OWN.length} databases are there`);
  // a database the fixture made as sa and did not give away (a load stopped half-way, or one run by hand)
  const unowned = [...databases].filter(([, owner]) => owner === 'sa').map(([n]) => n);
  if (unowned.length) {
    must(sqlAsSa(unowned.map((n) => `ALTER AUTHORIZATION ON DATABASE::[${n}] TO [${LOGIN}];`).join('\n')),
      'giving the databases to the login');
    console.log(`gave ${unowned.length} databases to ${LOGIN}`);
    databases = fixtureDatabases();
  }
  const owned = [...databases.values()].filter((o) => o === LOGIN).length;
  if (owned !== OWN.length || !fixtureComplete())
    stop(`the fixture is not whole after loading: ${owned} of ${OWN.length} databases owned by ${LOGIN}`);
  const seen = visibleToLogin(password);
  if (!seen)
    stop(`${LOGIN} cannot sign in with the password in ${shown(SECRETS)}`);
  console.log(`${LOGIN} sees ${seen[0]} catalogs: the fixture's ${seen[1]}, master and tempdb` +
    (seen[0] === VISIBLE ? '' : ` (expected ${VISIBLE}: does the login own other databases?)`));

  const conn = writeConnection(password, opts);
  console.log(`connection written to ${shown(LOCAL_CONNECTION)} (${conn.parameters.server}:${conn.parameters.port})`);
  const port = await portAnswers();
  if (port !== 'ok')
    console.warn(`warning: 127.0.0.1:${HOST_PORT} on this machine: ${port}. A wedged Docker port forward is fixed ` +
      `by \`docker restart ${CONTAINER}\``);
  console.log('next: node features/identifiers/catalogs-db.mjs connection [--host <alias>]');
}

async function connection(opts) {
  const host = grokHost(opts.host);
  if (opts.server || opts.port) {
    const password = readSecrets().login;
    if (!password)
      stop('no login password yet: run `up` first');
    writeConnection(password, opts);
  }
  const body = readJson(LOCAL_CONNECTION);
  if (!body)
    stop(`${shown(LOCAL_CONNECTION)} is missing: run \`up\` first`);
  const me = parseJson(must(grok(['raw', 'GET', '/users/current', '--output', 'json'], host), 'reading the account'),
    'the account');
  const found = parseJson(must(grok(['connections', 'list', '--filter', `shortName = "${CONNECTION}"`,
    '--output', 'json'], host), 'looking up the connection'), 'the connections');
  const own = (Array.isArray(found) ? found : []).find((x) => x.author?.id === me.id);
  if (own)
    body.id = own.id;
  const file = join(LOCAL, `save-${process.pid}.json`);
  writeFileSync(file, JSON.stringify(body), {mode: 0o600});
  let saved;
  try {
    saved = parseJson(must(grok(['connections', 'save', '--json', file, '--save-credentials', '--output', 'json'],
      host), 'saving the connection'), 'the saved connection');
  }
  finally {
    rmSync(file, {force: true});
  }
  console.log(`${own ? 'updated' : 'created'} ${saved.namespace ?? ''}${saved.name} on ${host} ` +
    `(${body.parameters.server}:${body.parameters.port}, login ${LOGIN})`);
  const test = grok(['connections', 'test', saved.id], host);
  const answer = (test.stdout + test.stderr).trim();
  if (!test.ok) {
    process.exitCode = 1;
    console.error(`the stand cannot reach it: ${answer.slice(0, 800)}\nIs the container up (\`status\`)? After the ` +
      'database was restarted, Grok Connect may hold dead pooled connections: restart its container. On Linux, see ' +
      'the README for the address.');
  }
  else
    console.log(`connection test: ${answer || 'ok'}`);
}

async function status() {
  let healthy = true;
  const problem = (text) => {
    healthy = false;
    console.log(`  ! ${text}`);
  };
  const c = container();
  const volume = volumeExists();
  console.log(`container ${CONTAINER}: ${c ? `${c.status}${c.ours ? '' : ' (NOT labelled as this fixture)'}` : 'absent'}`);
  for (const d of c?.differences ?? [])
    problem(`differs from the contract: ${d}`);
  console.log(`volume ${VOLUME}: ${volume ? 'present' : 'absent'}`);
  const secrets = readSecrets();
  console.log(`${shown(SECRETS)}: ${existsSync(SECRETS) ? `sa password ${secrets.sa ? 'cached' : 'not cached'}, login ` +
    `password ${secrets.login ? 'kept' : 'missing'}` : 'absent'}`);
  console.log(`${shown(LOCAL_CONNECTION)}: ${existsSync(LOCAL_CONNECTION) ? 'present' : 'absent'}`);
  if (c?.status !== 'running' || !c.ours) {
    problem('the server is not running: `up` starts it');
    return healthy;
  }
  const ping = sqlAsSa('SELECT 1;', {timeoutMs: 30000});
  if (!ping.ok) {
    problem(`SQL Server does not answer: ${(ping.stderr || ping.stdout).trim().slice(0, 300)}`);
    return healthy;
  }
  const databases = fixtureDatabases();
  const byOwner = {};
  for (const owner of databases.values())
    byOwner[owner] = (byOwner[owner] ?? 0) + 1;
  console.log(`the fixture's databases: ${databases.size} of ${OWN.length}` +
    (databases.size ? ` (owners: ${Object.entries(byOwner).map(([o, n]) => `${o} ${n}`).join(', ')})` : ''));
  if (databases.size !== OWN.length || byOwner[LOGIN] !== OWN.length)
    problem('`up` loads the missing ones and gives them to the login');
  else if (!fixtureComplete())
    problem('the working catalogs are incomplete: `up --reload`');
  if (secrets.login) {
    const seen = visibleToLogin(secrets.login);
    if (!seen)
      problem(`${LOGIN} cannot sign in with the kept password: \`up\` resets it`);
    else
      console.log(`${LOGIN} sees ${seen[0]} catalogs (expected ${VISIBLE})`);
  }
  else
    problem('no login password kept: `up` resets it');
  const port = await portAnswers();
  console.log(`127.0.0.1:${HOST_PORT}: ${port}`);
  if (port !== 'ok')
    problem(`a wedged Docker port forward is fixed by \`docker restart ${CONTAINER}\``);
  return healthy;
}

function down(opts) {
  const c = ownContainer();
  if (!opts.purge) {
    if (!c)
      return console.log(`no container ${CONTAINER}`);
    if (c.status === 'running')
      must(docker(['stop', '-t', '30', CONTAINER], {timeoutMs: 2 * MINUTE}), 'docker stop');
    console.log(`stopped ${CONTAINER}: every database on it is down until \`up\`; the container and the volume are kept`);
    return;
  }
  if (!opts.yes) {
    let others = '';
    if (c?.status === 'running') {
      const r = sqlAsSa(`SET NOCOUNT ON; SELECT COUNT(*) FROM sys.databases WHERE database_id > 4 AND name NOT IN (${OWN_LIST});`);
      if (r.ok)
        others = ` It holds ${rows(r.stdout)[0]?.[0]} databases that are not this fixture's.`;
    }
    console.error(`--purge removes the container ${CONTAINER} and its volume ${VOLUME}: EVERY database on the server ` +
      `is deleted, also the ones other logins keep there, not only this fixture's.${others} Add --yes to do it.`);
    process.exitCode = 1;
    return;
  }
  if (c)
    must(docker(['rm', '-f', CONTAINER], {timeoutMs: 2 * MINUTE}), 'docker rm');
  if (volumeExists())
    must(docker(['volume', 'rm', VOLUME], {timeoutMs: 2 * MINUTE}), 'docker volume rm');
  rmSync(LOCAL, {recursive: true, force: true});
  console.log(`removed ${CONTAINER}, the volume ${VOLUME} with every database on it, and ${shown(LOCAL)}. A stand's saved ` +
    'connection holds a password that no longer works: `up`, then `connection`.');
}

function options(args) {
  const opts = {};
  for (let i = 0; i < args.length; i++) {
    const a = args[i];
    if (a === '--reload' || a === '--purge' || a === '--yes')
      opts[a.slice(2)] = true;
    else if (a === '--host' || a === '--server' || a === '--port') {
      if (!args[i + 1] || args[i + 1].startsWith('--'))
        stop(`${a} needs a value`);
      opts[a.slice(2)] = args[++i];
    }
    else
      stop(`unknown option ${a}\n\n${USAGE}`);
  }
  if (opts.port && !/^\d+$/.test(opts.port))
    stop(`--port ${opts.port} is not a port`);
  return opts;
}

const [command, ...rest] = process.argv.slice(2);
try {
  const opts = options(rest);
  if (command === 'up')
    await up(opts);
  else if (command === 'connection')
    await connection(opts);
  else if (command === 'status')
    process.exitCode = await status() ? 0 : 1;
  else if (command === 'down')
    down(opts);
  else {
    console.log(USAGE);
    process.exitCode = command && !['help', '--help', '-h'].includes(command) ? 1 : 0;
  }
}
catch (e) {
  if (!(e instanceof Stop))
    throw e;
  console.error(`catalogs-db: ${e.message}`);
  process.exitCode = 1;
}
