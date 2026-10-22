/* The AI panel's chat as a user drives it, against a live claude-runtime container and a real model.
   A turn ends on the runtime's own terminal event (final, error or aborted) for the panel's session,
   read from the runtime WebSocket, and the panel leaving its busy state; the socket closing first loses
   the turn. The tap keeps each frame's type, socket, session and tool, a tool's one-line summary and an
   error's text (both written by the runtime from the model's calls, key-shaped text redacted), a prompt's
   mode and whether it carried the panel's saved transcript; not the prompt or the answer. */
import * as fs from 'node:fs';
import * as path from 'node:path';
import {randomUUID} from 'node:crypto';
import {execFileSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {type Page, type WebSocket, test} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, escapeRegExp, exactText, expect, pollMs, serverRequests} from '@datagrok-libraries/bdd/runtime';
import {AI_PANEL, AI_PROMPT_BOX} from './elements.js';

declare const grok: any;

const RUNTIME_WS = /\/docker\/containers\/proxy(-ws)?\/[^/]+\/ws\b/;
const TERMINAL = ['final', 'error', 'aborted'];
// A turn is a scenario's work, not a check: unlike pollMs, its budget is not narrowed in a @known-failure scenario.
const TURN_MS = Number(process.env.GROKKY_BDD_TURN_MS ?? 300000);
const SECRET_RE = /sk-(?:ant|proj|svcacct)-[\w-]{20,}|dapi[0-9a-f]{20}|Bearer\s+\S{20}/;
const SECRET_TEXT = /(sk-(?:ant|proj|svcacct)-|dapi[0-9a-f]{20}|Bearer\s+)\S*/g;
const RESTORED_TRANSCRIPT = '[Conversation restored from saved history';
const AGENT_PROCESS = /-m grokky_agent\b/;
const CANARY_DIRS = '/users/00000000-0000-4bdd-8000-';
const SHOTS = process.env.GROKKY_BDD_SHOTS ?? fileURLToPath(new URL('../test-results/screenshots', import.meta.url));

/* The terminal commands an agent ran for a turn, from its Hermes session store: the runtime's activity
   summary of a terminal call is "Run a command". Arguments: the turn's start in the container's epoch
   seconds, and the start of its prompt, which the store keeps in the user message. */
const TURN_COMMANDS = `
import glob, json, sqlite3, sys
since, prompt = float(sys.argv[1]), sys.argv[2]
for db in glob.glob('/var/lib/grokky/hermes/*/state.db') + glob.glob('/var/lib/grokky/hermes/parked/*/state.db'):
    try:
        c = sqlite3.connect('file:' + db + '?mode=ro', uri=True)
        sessions = [r[0] for r in c.execute("select distinct session_id from messages where role = 'user' and timestamp >= ? and instr(content, ?) > 0", (since, prompt))]
        calls = [r[0] for s in sessions for r in c.execute("select tool_calls from messages where session_id = ? and timestamp >= ? and tool_calls is not null", (s, since))]
    except sqlite3.Error:
        continue
    for call in (x for listed in calls for x in json.loads(listed)):
        f = call.get('function') or {}
        if f.get('name') == 'terminal':
            print(json.dumps(json.loads(f.get('arguments') or '{}').get('command', '')))
`;

interface Frame {
  dir: 'in' | 'out';
  type: string;
  socket: number;
  sessionId?: string;
  tool?: string;
  summary?: string;
  error?: string;
  mode?: string;
  restored?: boolean;
}

interface PanelState {
  busy: boolean;
  responses: number;
  lastText: string;
  systemMessages: string[];
  handledNatively: number;
  question: boolean;
  confirmCards: number;
}

interface Turn {
  prompt: string;
  startedAt: number;
  from: number;
  before: PanelState;
  sessionId?: string;
  socket?: number;
  native?: boolean;
  end?: string;
  error?: string;
  answer: string;
  question?: {options: string[]; picked: string};
  confirmations: string[];
}

interface Process {
  pid: number;
  cmd: string;
}

interface Canary {
  note: string;
  dir: string;
  accessed: string;
}

interface PageColumn {
  name: string;
  type: string;
  isNumerical: boolean;
  get(row: number): number | string;
}

interface PageTable {
  name: string;
  rowCount: number;
  dart: unknown;
  columns: {names(): string[]; toList(): PageColumn[]};
  get(column: string, row: number): number;
}

const taps = new WeakMap<Page, Frame[]>();
const turns = new WeakMap<Page, Turn[]>();
const containers = new WeakMap<Page, string>();
const canaries = new WeakMap<Page, Canary>();

function tapOf(page: Page): Frame[] {
  const frames = taps.get(page);
  if (!frames)
    throw new Error('the runtime WebSocket is not recorded: put "Given the AI runtime traffic is recorded" in the Background, before "user is logged in"');
  return frames;
}

function turnsOf(page: Page): Turn[] {
  const list = turns.get(page) ?? [];
  turns.set(page, list);
  return list;
}

function lastTurn(page: Page): Turn {
  const list = turnsOf(page);
  if (list.length === 0)
    throw new Error('no prompt was sent to the AI panel in this scenario');
  return list[list.length - 1];
}

const panel = (page: Page): Promise<PanelState> => page.evaluate((sel) => {
  const roots = Array.from(document.querySelectorAll(sel)) as HTMLElement[];
  const root = roots.find((p) => p.offsetParent !== null) ?? roots[0];
  if (!root)
    throw new Error('the AI panel is not on the page');
  const responses = Array.from(root.querySelectorAll('.d4-ai-panel-output-area .d4-ai-assistant-response-container')) as HTMLElement[];
  const submit = (Array.from(root.querySelectorAll('.d4-ai-panel-output-area button')) as HTMLButtonElement[])
    .find((b) => b.textContent!.trim().toLowerCase() === 'submit' && !b.disabled);
  return {
    busy: !!root.querySelector('.d4-ai-panel-run-controls .fa-stop'),
    responses: responses.length,
    lastText: responses.length ? responses[responses.length - 1].innerText.trim() : '',
    systemMessages: (Array.from(root.querySelectorAll('.grokky-system-message')) as HTMLElement[]).map((e) => e.innerText.trim()),
    handledNatively: root.querySelectorAll('.d4-ai-handled-natively-icon').length,
    question: !!submit,
    confirmCards: root.querySelectorAll('.d4-ai-panel-output-area .grokky-confirm-card').length,
  };
}, AI_PANEL);

function framesOf(page: Page, turn: Turn): Frame[] {
  const frames = tapOf(page).slice(turn.from);
  const sent = frames.find((f) => f.dir === 'out' && f.type === 'user_message');
  turn.sessionId ??= sent?.sessionId;
  turn.socket ??= sent?.socket;
  return turn.sessionId ? frames.filter((f) => f.sessionId === turn.sessionId) : [];
}

const toolsOf = (page: Page, turn: Turn): string[] => framesOf(page, turn)
  .filter((f) => f.dir === 'in' && (f.type === 'tool_activity' || f.type === 'input_request')).map((f) => f.tool ?? '?');

async function sendPrompt(page: Page, prompt: string): Promise<Turn> {
  test.info().setTimeout(test.info().timeout + TURN_MS);
  const turn: Turn = {prompt, startedAt: Date.now(), from: tapOf(page).length, before: await panel(page), answer: '', confirmations: []};
  turnsOf(page).push(turn);
  const box = page.locator(AI_PROMPT_BOX).filter({visible: true}).first();
  await box.fill(prompt);
  await box.press('Enter');
  await expect.poll(async () => {
    turn.native = (await panel(page)).handledNatively > turn.before.handledNatively;
    return turn.native || framesOf(page, turn).length > 0;
  }, {message: 'the panel handed the prompt to the platform\'s interpreter or to the runtime', timeout: 30000}).toBe(true);
  return turn;
}

function requireRuntime(turn: Turn): Turn {
  if (turn.native)
    throw new Error('the platform\'s own prompt interpreter handled the prompt, so the runtime was never called; reword it');
  return turn;
}

async function turnStatus(page: Page, turn: Turn): Promise<'pending' | 'question' | 'confirm' | 'done' | 'lost'> {
  const s = await panel(page);
  if (s.question)
    return 'question';
  if (s.confirmCards - turn.before.confirmCards > turn.confirmations.length)
    return 'confirm';
  const end = framesOf(page, turn).find((f) => f.dir === 'in' && TERMINAL.includes(f.type));
  if (end && !s.busy) {
    turn.end = end.type;
    turn.error = end.error;
    turn.answer = s.responses > turn.before.responses || s.lastText !== turn.before.lastText ? s.lastText : '';
    return 'done';
  }
  return !end && tapOf(page).slice(turn.from).some((f) => f.type === 'close' && f.socket === turn.socket) ? 'lost' : 'pending';
}

interface Answers {
  question?: () => Promise<void>;
  confirm?: () => Promise<void>;
}

async function waitTurn(page: Page, turn: Turn, answers: Answers = {}): Promise<void> {
  while (!turn.native) {
    let status = 'pending';
    await expect.poll(async () => (status = await turnStatus(page, turn)),
      {message: `the AI turn "${turn.prompt.slice(0, 60)}" ended`, timeout: TURN_MS}).not.toBe('pending');
    if (status === 'lost')
      throw new Error('the runtime WebSocket closed before the turn ended: the runtime drops the turns of a closed connection');
    if (status === 'done' && turn.end === 'error')
      throw new Error(`the AI turn ended in an error: ${turn.error}`);
    if (status === 'done')
      return;
    const answer = status === 'question' ? answers.question : answers.confirm;
    if (!answer)
      throw new Error(`the AI asked for ${status === 'question' ? 'an answer (AskUserQuestion form)' : 'a confirmation (confirm card)'} that this step does not give`);
    await answer();
    if (status === 'question')
      answers.question = undefined;
  }
}

/** Every question of the form is answered: the option matching the pattern where one does, else its first. */
async function answerQuestion(page: Page, turn: Turn, option: RegExp): Promise<void> {
  const form = page.locator(`${AI_PANEL} .d4-ai-panel-output-area`).filter({visible: true}).first();
  const offered: string[][] = [];
  for (const select of await form.locator('select').filter({visible: true}).all()) {
    const options = (await select.locator('option').allTextContents()).map((o) => o.trim());
    offered.push(options);
    const picked = options.find((o) => option.test(o));
    turn.question ??= picked ? {options, picked} : undefined;
    await select.selectOption({label: picked ?? options[0]});
  }
  if (!turn.question)
    throw new Error(`the AI's questions offer no option matching ${option}; they offer: ${offered.map((o) => o.join(', ')).join(' | ')}`);
  await form.locator('button', {hasText: /^\s*submit\s*$/i}).filter({visible: true}).last().click();
  await expect.poll(async () => (await panel(page)).question, {message: 'the question form was submitted'}).toBe(false);
}

async function answerConfirmation(page: Page, turn: Turn, button: string): Promise<void> {
  const card = page.locator(`${AI_PANEL} .grokky-confirm-card`).filter({visible: true}).last();
  turn.confirmations.push((await card.innerText()).trim());
  await card.locator('button', {hasText: exactText(button)}).click();
}

function answerOf(turn: Turn): string {
  if (turn.native)
    throw new Error('the platform\'s own prompt interpreter took the prompt: the AI never answered it');
  expect(turn.end, 'how the turn ended').toBe('final');
  return turn.answer;
}

/** Some number of the answer within the tolerance: a model shows its working even when asked for the number alone. */
function expectNumberIn(answer: string, want: number, tolerance: number, what: string): void {
  const numbers = (answer.match(/-?\d+(?:\.\d+)?/g) ?? []).map(Number);
  expect(numbers.some((n) => Math.abs(n - want) < tolerance),
    `${what}, ${want.toFixed(3)}, among the numbers of the answer: ${numbers.join(', ') || 'none'}`).toBe(true);
}

function median(xs: number[]): number {
  xs.sort((a, b) => a - b);
  const n = xs.length;
  return n % 2 ? xs[(n - 1) / 2] : (xs[n / 2 - 1] + xs[n / 2]) / 2;
}

/** Missing values left out. */
const columnMedian = async (page: Page, column: string): Promise<number> => median(await page.evaluate((name) => {
  const c = grok.shell.t?.col(name);
  if (!c)
    throw new Error(`the current table has no "${name}" column`);
  const values: number[] = [];
  for (let i = 0; i < c.length; i++) {
    if (!c.isNone(i))
      values.push(c.get(i));
  }
  return values;
}, column));

/** The stand's runtime container in this machine's Docker: the one the platform runs as the service the panel uses. */
async function runtimeContainer(page: Page): Promise<string> {
  const known = containers.get(page);
  if (known)
    return known;
  const id: string = await page.evaluate(async () =>
    (await grok.dapi.docker.dockerContainers.filter('name = "grokky-claude-runtime"').list())[0]?.id ?? '');
  if (!id)
    throw new Error('no grokky-claude-runtime container on the stand');
  const {serviceName} = await (await serverRequests(page)).get<{serviceName: string}>(`/docker/containers/${id}`);
  const name = execFileSync('docker', ['ps', '--filter', `label=service=${serviceName}`, '--format', '{{.Names}}'],
    {encoding: 'utf8'}).trim();
  if (!name)
    throw new Error(`no container labelled service=${serviceName} in this machine's Docker: the stand's runtime runs elsewhere`);
  containers.set(page, name);
  return name;
}

const inRuntime = (container: string, script: string): string =>
  execFileSync('docker', ['exec', container, 'sh', '-c', script], {encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe']});

function runtimeProcesses(container: string): Process[] {
  const listing = inRuntime(container,
    'for p in /proc/[0-9]*; do printf "%s %s\\n" "${p#/proc/}" "$(tr "\\0" " " 2>/dev/null < "$p/cmdline")"; done');
  return listing.split('\n').flatMap((line) => {
    const m = /^(\d+) (.+)$/.exec(line.trim());
    return m && !m[2].includes('/proc/[0-9]') ? [{pid: Number(m[1]), cmd: m[2].trim()}] : [];
  });
}

const sandboxFiles = (container: string, dir: string): string[] =>
  inRuntime(container, `find ${dir} -path ${dir}/agents -prune -o -type f -print 2>/dev/null || true`).split('\n').filter((f) => f !== '');

function removeFiles(container: string, files: string[]): void {
  if (files.length)
    execFileSync('docker', ['exec', container, 'rm', '-f', '--', ...files]);
}

/** Key-shaped text and the session's own token, on the page or in its inputs; the token never leaves the page. */
const secretsOnPage = (page: Page): Promise<string[]> => page.evaluate(([re, sel]) => {
  let text = document.body.innerText + '\n' + (document.querySelector(sel)?.textContent ?? '');
  for (const el of Array.from(document.querySelectorAll('textarea, input')) as HTMLInputElement[])
    text += '\n' + (el.value || '');
  const token = String(grok.dapi.token ?? '');
  return [new RegExp(re).test(text) ? 'key-shaped text' : '', token.length > 8 && text.includes(token) ? 'the session token' : '']
    .filter((s) => s !== '');
}, [SECRET_RE.source, AI_PANEL] as [string, string]);

function canaryOf(page: Page): Canary {
  const canary = canaries.get(page);
  if (!canary)
    throw new Error('no note was planted: "Given another user\'s folder in the AI\'s sandbox holds a private note" comes first');
  return canary;
}

// --- the chat ---------------------------------------------------------------------------------------

export const trafficRecorded = Given('the AI runtime traffic is recorded', async (page: Page) => {
  if (taps.has(page))
    return;
  const frames: Frame[] = [];
  taps.set(page, frames);
  const turnLog = process.env.GROKKY_BDD_TURNS;
  const scenario = test.info().title;
  let sockets = 0;
  page.on('websocket', (ws: WebSocket) => {
    if (!RUNTIME_WS.test(ws.url()))
      return;
    const socket = ++sockets;
    const record = (dir: 'in' | 'out') => (data: {payload: string | Buffer}) => {
      let d: any;
      try {
        d = JSON.parse(String(data.payload));
      }
      catch {
        return;
      }
      if (!d?.type)
        return;
      frames.push({dir, type: d.type, socket, sessionId: d.sessionId, tool: d.name ?? d.toolName,
        summary: d.type === 'tool_activity' ? String(d.summary ?? '').replace(SECRET_TEXT, '$1<redacted>') : undefined,
        error: d.type === 'error' ? String(d.message ?? '').replace(SECRET_TEXT, '$1<redacted>') : undefined,
        mode: d.type === 'user_message' ? String(d.systemPromptMode ?? 'full') : undefined,
        restored: d.type === 'user_message' ? String(d.message ?? '').startsWith(RESTORED_TRANSCRIPT) : undefined});
      if (turnLog && ['models', 'user_message', 'tool_activity', 'final', 'error'].includes(d.type)) {
        const {tool, summary} = frames[frames.length - 1];
        fs.appendFileSync(turnLog, JSON.stringify({scenario, dir, type: d.type, sessionId: d.sessionId, tool, summary,
          model: d.type === 'models' ? {default: d.defaultModel, fast: d.fastModel} : d.model, metrics: d.metrics,
          error: frames[frames.length - 1].error, revision: d.revision,
          answer: d.type === 'final' ? String(d.content ?? '').replace(SECRET_TEXT, '$1<redacted>').slice(0, 2000) :
            undefined}) + '\n');
      }
    };
    ws.on('framereceived', record('in'));
    ws.on('framesent', record('out'));
    ws.on('close', () => frames.push({dir: 'in', type: 'close', socket}));
  });
}, {tier: 'api', description: 'records the frames of every claude-runtime WebSocket and its close, across ' +
  'reloads; must precede the first page load that may open one; with GROKKY_BDD_TURNS=<file> it also appends ' +
  'the models, prompts, tool calls, answers and errors there as JSON lines (model, cost, tokens)'});

export const liveModelGate = Given('live AI model calls are allowed', async (_page: Page) => {
  test.skip(process.env.GROKKY_BDD_LIVE !== '1',
    'set GROKKY_BDD_LIVE=1: the rest of this test sends prompts through the claude-runtime container to a paid model');
}, {description: 'a gate: every scenario after it calls the model provider, so a plain suite run skips them'});

export const runtimeHealthy = Given('the AI runtime container is healthy', async (page: Page) => {
  await expect.poll(() => page.evaluate(async () => {
    try {
      const c = (await grok.dapi.docker.dockerContainers.filter('name = "grokky-claude-runtime"').list())[0];
      if (!c)
        return 'no grokky-claude-runtime container on the stand';
      const r = await grok.dapi.docker.dockerContainers.fetchProxy(c.id, '/health');
      return r.ok ? String((await r.json()).status) : `HTTP ${r.status}`;
    }
    catch (e) {
      return String(e);
    }
  }), {message: 'the runtime\'s /health through the stand\'s docker proxy', timeout: pollMs(120000)}).toBe('ok');
}, {tier: 'api', description: '/health of the grokky-claude-runtime container through grok.dapi.docker, polled while an on-demand container starts'});

export const aiPanelNewChat = Given('the AI panel is open on a new chat', async (page: Page) => {
  await page.evaluate(() => { grok.shell.windows.showAI = true; });
  await expect(page.locator(AI_PROMPT_BOX).filter({visible: true}).first(), 'the AI panel\'s prompt box').toBeVisible({timeout: pollMs(30000)});
  const root = page.locator(AI_PANEL).filter({visible: true}).first();
  const messages = root.locator('.d4-ai-assistant-response-container, .d4-ai-user-prompt-container');
  if (await messages.count() > 0)
    await root.locator('i.fa-plus').filter({visible: true}).first().click();
  await expect(messages, 'the messages of the panel\'s chat').toHaveCount(0);
  turns.set(page, []);
}, {tier: 'ui', description: 'shows the AI panel and starts a new chat with its + icon when the current one has messages'});

export const askAI = When('user asks the AI panel {string}', async (page: Page, prompt: string) => {
  await waitTurn(page, requireRuntime(await sendPrompt(page, prompt)));
}, {tier: 'ui', description: 'types the prompt, presses Enter and waits for the runtime\'s final, error or aborted event for the panel\'s session and the panel to go idle'});

export const askAIAndAnswer = When('user asks the AI panel {string} and picks the option matching {string} in its question',
  async (page: Page, prompt: string, option: string) => {
    const turn = requireRuntime(await sendPrompt(page, prompt));
    await waitTurn(page, turn, {question: () => answerQuestion(page, turn, new RegExp(option, 'i'))});
    if (!turn.question)
      throw new Error(`the AI ended the turn without asking a question; it answered: "${turn.answer.slice(0, 200)}"`);
  }, {tier: 'ui', description: 'answers every question of the AskUserQuestion form the turn shows (the option matching the pattern, case-insensitive, where one does), then Submit'});

export const askAIAndConfirm = When('user asks the AI panel {string} and clicks {string} on its confirmation card',
  async (page: Page, prompt: string, button: string) => {
    const turn = requireRuntime(await sendPrompt(page, prompt));
    await waitTurn(page, turn, {confirm: () => answerConfirmation(page, turn, button)});
    if (turn.confirmations.length === 0)
      throw new Error(`the AI ended the turn without asking for a confirmation; it answered: "${turn.answer.slice(0, 200)}"`);
  }, {tier: 'ui', description: 'clicks the button (OK or CANCEL) on every confirm card the turn shows, and waits for the turn to end'});

export const sendToAI = When('user sends {string} to the AI panel', async (page: Page, prompt: string) => {
  requireRuntime(await sendPrompt(page, prompt));
}, {tier: 'ui', description: 'types the prompt and presses Enter; done once the runtime has the prompt, while the turn runs'});

export const submitToAI = When('user submits {string} in the AI panel', async (page: Page, prompt: string) => {
  await waitTurn(page, await sendPrompt(page, prompt));
}, {tier: 'ui', description: 'whoever takes it, the platform\'s own prompt interpreter or the runtime; done when that one is — for a claim on the outcome alone'});

export const reopenChat = When('user reopens the AI chat from the chat history', async (page: Page) => {
  const first = turnsOf(page)[0];
  if (!first)
    throw new Error('no chat was held in this scenario to reopen');
  await page.evaluate(() => { grok.shell.windows.showAI = true; });
  await page.locator(`${AI_PANEL} .d4-ai-panel-run-controls .fa-history`).filter({visible: true}).first().click();
  await page.locator('.d4-menu-popup .d4-menu-item').filter({visible: true, hasText: first.prompt.slice(0, 40)}).first().click();
  await expect(page.locator(`${AI_PANEL} .d4-ai-user-prompt-container`).filter({visible: true}).first(),
    'the reopened chat\'s first prompt').toContainText(first.prompt.slice(0, 40));
}, {tier: 'ui', description: 'the panel\'s Chat History icon, then the entry the chat\'s first prompt names'});

// --- what the turns did -----------------------------------------------------------------------------

export const turnStopped = Then('the AI turn should have been stopped', async (page: Page) => {
  const turn = lastTurn(page);
  await expect.poll(async () => {
    const ends = framesOf(page, turn).filter((f) => f.dir === 'in' && TERMINAL.includes(f.type)).map((f) => f.type);
    return `${ends.join(', ') || 'no terminal event'}; panel ${(await panel(page)).busy ? 'busy' : 'idle'}`;
  }, {message: 'how the runtime ended the stopped turn, and the panel', timeout: pollMs(15000)}).toBe('aborted; panel idle');
  turn.end = 'aborted';
}, {description: 'the runtime\'s aborted event for the turn\'s session, and the panel idle again'});

export const nothingAfterStop = Then('the stopped AI turn should have sent nothing after it was aborted', async (page: Page) => {
  const list = turnsOf(page);
  const stopped = list.findIndex((t) => t.end === 'aborted');
  if (stopped < 0 || stopped === list.length - 1)
    throw new Error('claim it after a stopped turn and the turn that followed it');
  const all = tapOf(page);
  const {from, sessionId} = list[stopped];
  const abortedAt = all.findIndex((f, i) => i >= from && f.sessionId === sessionId && f.type === 'aborted');
  const late = all.slice(abortedAt + 1, list[stopped + 1].from)
    .filter((f) => f.dir === 'in' && f.sessionId === sessionId && (f.type === 'chunk' || f.type === 'final'));
  expect(late.length, 'chunk or final events for the session after its aborted event, until the next prompt').toBe(0);
}, {description: 'read after the next turn: no chunk or final arrived for the stopped turn after its aborted event'});

export const toolRan = Then('the AI should have used the {string} tool', async (page: Page, tool: string) => {
  expect(toolsOf(page, lastTurn(page)), 'the tools of the last turn, as tool_activity and input_request events name them').toContain(tool);
}, {description: 'a tool_activity or input_request of the last turn names the tool (bare name)'});

export const oneOfToolsRan = Then('the AI should have used one of the tools {string}', async (page: Page, list: string) => {
  const used = toolsOf(page, lastTurn(page));
  expect(list.split(/\s*,\s*/).filter((t) => used.includes(t)), `the tools of the last turn: ${used.join(', ')}`).not.toEqual([]);
}, {description: 'comma-separated bare tool names; the last turn used at least one'});

export const terminalRan = Then('a terminal command of the AI should have run {word}', async (page: Page, program: string) => {
  const turn = lastTurn(page);
  const container = await runtimeContainer(page);
  const since = Number(inRuntime(container, 'date +%s.%N')) - (Date.now() - turn.startedAt) / 1000 - 1;
  const commands = execFileSync('docker', ['exec', container, 'python3', '-c', TURN_COMMANDS, String(since), turn.prompt.slice(0, 60)],
    {encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe']}).split('\n').filter((l) => l !== '').map((l) => String(JSON.parse(l)));
  expect(commands.some((c) => new RegExp(`\\b${escapeRegExp(program)}\\b`).test(c)),
    `the terminal commands of the last turn, from its agent's session store: ${commands.join(' | ') || 'none'}`).toBe(true);
}, {description: 'one of the terminal commands the last turn\'s agent ran names the program, read from its Hermes session store in the container (the activity summary does not carry the command)'});

export const confirmationAsked = Then('the AI should have asked to confirm with a card naming {string}', async (page: Page, name: string) => {
  const turn = lastTurn(page);
  expect(toolsOf(page, turn), 'the tools of the last turn').toContain('datagrok_confirm');
  expect(turn.confirmations.some((c) => c.includes(name)), `the confirm cards the turn showed: ${turn.confirmations.join(' | ')}`).toBe(true);
}, {description: 'a datagrok_confirm round trip of the last turn, whose card names the entity'});

export const promptMode = Then('the last prompt should have reached the runtime in {word} mode', async (page: Page, mode: string) => {
  const modes = framesOf(page, lastTurn(page)).filter((f) => f.dir === 'out' && f.type === 'user_message').map((f) => f.mode);
  expect([...new Set(modes)], 'the systemPromptMode of the last turn\'s user_message frames (full when none)').toEqual([mode]);
}, {description: 'bash, none or full (no systemPromptMode): the mode the panel sent the last prompt in'});

export const questionAsked = Then('the AI should have asked a question offering {int} or more options', async (page: Page, n: number) => {
  const q = lastTurn(page).question;
  expect(q, 'the question the AI asked in the panel').toBeTruthy();
  expect(q!.options.length, `the options it offered: ${q!.options.join(', ')}`).toBeGreaterThanOrEqual(n);
}, {description: 'the AskUserQuestion question whose option was picked, as answered'});

export const systemMessage = Then('the AI panel should show the message {string}', async (page: Page, text: string) => {
  await expect.poll(async () => (await panel(page)).systemMessages.join(' | '),
    {message: `the panel's system messages include "${text}"`}).toContain(text);
}, {description: 'one of the panel\'s system messages (grokky-system-message) contains the text'});

export const earlierTurnsShown = Then('the AI panel should show the earlier turns of the chat', async (page: Page) => {
  const shown = (await page.locator(AI_PANEL).filter({visible: true}).first()
    .locator('.d4-ai-user-prompt-container, .d4-ai-assistant-response-container').allInnerTexts()).map((t) => t.trim());
  for (const turn of turnsOf(page)) {
    const answer = answerOf(turn);
    expect(shown.some((t) => t.includes(turn.prompt)), `the prompt "${turn.prompt.slice(0, 60)}" among: ${shown.join(' | ').slice(0, 300)}`).toBe(true);
    expect(answer !== '' && shown.some((t) => t.includes(answer)), `the answer "${answer.slice(0, 60)}"`).toBe(true);
  }
}, {description: 'every prompt the chat held so far and its answer, as the panel shows them'});

export const sameSession = Then('every turn of the chat should have used the same AI session', async (page: Page) => {
  const ids = turnsOf(page).map((t) => t.sessionId);
  expect(ids.length, 'turns in the chat').toBeGreaterThanOrEqual(2);
  expect(new Set(ids).size, `the session ids of the turns: ${ids.join(', ')}`).toBe(1);
}, {description: 'the user_message frames of all the scenario\'s turns carry one sessionId, so the runtime continued it'});

export const historyKept = Then('the AI runtime should have kept the chat\'s history itself', async (page: Page) => {
  const frames = framesOf(page, lastTurn(page));
  const resets = frames.filter((f) => f.type === 'session_reset').length;
  const replays = frames.filter((f) => f.dir === 'out' && f.type === 'user_message' && f.restored).length;
  expect(`${resets} session_reset, ${replays} prompt(s) replaying the panel's transcript`,
    'how the last turn got the chat\'s history').toBe('0 session_reset, 0 prompt(s) replaying the panel\'s transcript');
}, {description: 'the last turn ran on the history the runtime holds: no session_reset, no prompt carrying the panel\'s saved transcript'});

// --- answers against the data -----------------------------------------------------------------------

export const highestMeanNamed = Then('the AI answer should name the {string} with the highest mean {string} where {string} is {string}',
  async (page: Page, group: string, value: string, filterColumn: string, filterValue: string) => {
    const answer = answerOf(lastTurn(page));
    const means: Record<string, number> = await page.evaluate(([g, v, fc, fv]) => {
      const df = grok.shell.t;
      const sums: Record<string, [number, number]> = {};
      for (let i = 0; i < df.rowCount; i++) {
        if (String(df.get(fc, i)) !== fv || df.col(v).isNone(i))
          continue;
        const key = String(df.get(g, i));
        sums[key] ??= [0, 0];
        sums[key][0] += df.get(v, i);
        sums[key][1]++;
      }
      return Object.fromEntries(Object.entries(sums).map(([k, [sum, n]]) => [k, sum / n]));
    }, [group, value, filterColumn, filterValue]);
    const [top, mean] = Object.entries(means).sort((a, b) => b[1] - a[1])[0];
    const named = Object.keys(means).filter((k) => new RegExp(`\\b${escapeRegExp(k)}\\b`).test(answer));
    expect(named, `the ${group} values the answer names: "${answer.slice(0, 200)}"`).toEqual([top]);
    expectNumberIn(answer, mean, 0.051, `the mean ${value} of ${top}`);
  }, {description: 'computed from the current table: the answer names the top group alone and gives its mean to one decimal'});

export const medianGiven = Then('the AI answer should be the median of {string}', async (page: Page, column: string) => {
  expectNumberIn(answerOf(lastTurn(page)), await columnMedian(page, column), 0.051, `the median ${column}`);
}, {description: 'a number of the answer, to one decimal, against the median computed from the current table'});

export const medianDivided = Then('the AI answer should be the median of {string} divided by {float}', async (page: Page, column: string, divisor: number) => {
  expectNumberIn(answerOf(lastTurn(page)), await columnMedian(page, column) / divisor, 0.06, `the median ${column} over ${divisor}`);
}, {description: 'a number of the answer, to one decimal, against the median computed from the current table over the divisor'});

export const groupMedians = Then('the AI answer should give the median {string} of each {string} in {string}',
  async (page: Page, value: string, group: string, file: string) => {
    const answer = answerOf(lastTurn(page));
    const groups: Record<string, number[]> = await page.evaluate(async ([v, g, f]) => {
      const df = await grok.dapi.files.readCsv(f);
      const byGroup: Record<string, number[]> = {};
      for (let i = 0; i < df.rowCount; i++) {
        if (!df.col(v).isNone(i))
          (byGroup[String(df.get(g, i))] ??= []).push(df.get(v, i));
      }
      return byGroup;
    }, [value, group, file]);
    for (const [key, xs] of Object.entries(groups)) {
      const want = median(xs);
      const after = new RegExp(`\\b${escapeRegExp(key)}\\b\\D*?(-?\\d+(?:\\.\\d+)?)`, 'i').exec(answer);
      if (!after)
        throw new Error(`the answer gives no number after "${key}": "${answer.slice(0, 300)}"`);
      expect(Math.abs(Number(after[1]) - want), `the number after "${key}" against its median ${want}`).toBeLessThan(0.051);
    }
  }, {description: 'computed from the file: for every group, the first number after its name is its median to one decimal'});

export const numericColumnsNamed = Then('the AI answer should name every numeric column of the table', async (page: Page) => {
  const answer = answerOf(lastTurn(page));
  const numeric: string[] = await page.evaluate(() => (grok.shell.t as PageTable).columns.toList()
    .filter((c) => c.isNumerical && c.type !== 'datetime').map((c) => c.name));
  expect(numeric.filter((c) => !new RegExp(`\\b${escapeRegExp(c)}\\b`, 'i').test(answer)),
    `numeric columns of the table the answer leaves out: "${answer.slice(0, 300)}"`).toEqual([]);
}, {description: 'every numerical column of the current table, dates aside (core counts them numerical), by name in any case, in the last answer'});

export const hostNameAnswered = Then('the AI answer should be exactly the host name of the AI runtime container', async (page: Page) => {
  const host = inRuntime(await runtimeContainer(page), 'hostname').trim();
  expect(answerOf(lastTurn(page)), 'the answer, against `hostname` run in the container').toBe(host);
}, {description: 'the raw stdout of the command and nothing else, compared with `hostname` run through docker exec'});

export const helpPageCited = Then('the AI answer should link to the help page {string} it read', async (page: Page, helpPage: string) => {
  const turn = lastTurn(page);
  const read = framesOf(page, turn).filter((f) => f.tool === 'read_file').map((f) => f.summary ?? '');
  expect(read, 'the files the turn read').toContain(`Read workspace/help/${helpPage}.md`);
  const hrefs = await page.locator(AI_PANEL).filter({visible: true}).first().locator('.d4-ai-assistant-response-container').last()
    .locator('a').evaluateAll((links) => links.map((a) => a.getAttribute('href') ?? ''));
  expect(hrefs.map((h) => h.replace(/(\.md)?\/?([?#].*)?$/, '')), 'the links of the answer').toContain(`https://datagrok.ai/help/${helpPage}`);
}, {description: 'the turn read workspace/help/<page>.md (its read_file activity) and the last answer links to https://datagrok.ai/help/<page>'});

// --- the table and the grid -------------------------------------------------------------------------

export const bmiValues = Then('every value of {string} column should be the BMI of {string} in kg and {string} in cm',
  async (page: Page, bmi: string, weight: string, height: string) => {
    await expect.poll(() => page.evaluate(([b, w, h]) => {
      const df = grok.shell.t;
      const [cb, cw, ch] = [df.col(b), df.col(w), df.col(h)];
      if (!cb)
        return `no "${b}" column; the table has: ${df.columns.names().join(', ')}`;
      const wrong: string[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        const want = cw.isNone(i) || ch.isNone(i) ? null : cw.get(i) / (ch.get(i) / 100) ** 2;
        const got = cb.isNone(i) ? null : cb.get(i);
        if (want === null ? got !== null : got === null || Math.abs(got - want) > 0.051)
          wrong.push(`row ${i + 1}: ${got} for ${want?.toFixed(2)}`);
      }
      return wrong.length ? `${wrong.length} of ${df.rowCount} rows off, e.g. ${wrong.slice(0, 3).join('; ')}` : 'every row';
    }, [bmi, weight, height]), {message: `${bmi} against ${weight} / (${height} / 100)², to one decimal`}).toBe('every row');
  }, {tier: 'api', description: 'every row of the current table, to one decimal; empty where the weight or the height is'});

export const groupSummary = Then('the current table should hold the count and the mean {string} of every {string} group of table {string}',
  async (page: Page, values: string, group: string, source: string) => {
    await expect.poll(() => page.evaluate(([vs, g, src]) => {
      const from = (grok.shell.tables as PageTable[]).find((t) => t.name === src);
      const df: PageTable = grok.shell.t;
      if (!from)
        return [`no table "${src}" is open`];
      if (df.dart === from.dart || grok.shell.tv?.dataFrame?.dart !== df.dart)
        return [`the current view shows ${grok.shell.tv?.dataFrame?.name ?? 'no table'}, not a new summary table`];
      const means = vs.split(/\s*,\s*/);
      const truth: Record<string, number[]> = {};
      for (let i = 0; i < from.rowCount; i++) {
        const sums = truth[String(from.get(g, i))] ??= new Array(means.length + 1).fill(0);
        sums[0]++;
        for (const [j, m] of means.entries())
          sums[j + 1] += from.get(m, i);
      }
      const groups = Object.keys(truth).sort().join(', ');
      const columns = df.columns.toList();
      const numbers = columns.filter((c) => c.isNumerical && c.type !== 'datetime');
      const groupCol = columns.find((c) => !c.isNumerical &&
        Array.from({length: df.rowCount}, (_, i) => String(c.get(i))).sort().join(', ') === groups);
      if (!groupCol)
        return [`no column of ${df.name} lists the ${g} groups ${groups} once each; its columns: ${df.columns.names().join(', ')}`];
      const problems: string[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        const sums = truth[String(groupCol.get(i))];
        const row = numbers.map((c) => Number(c.get(i)));
        for (const [j, want] of [sums[0], ...means.map((_, k) => sums[k + 1] / sums[0])].entries()) {
          if (!row.some((x) => Math.abs(x - want) <= (j === 0 ? 0 : 0.051)))
            problems.push(`${groupCol.get(i)}: no ${j === 0 ? 'count' : `mean ${means[j - 1]}`} ${want.toFixed(j === 0 ? 0 : 2)} among ${row.join(', ')}`);
        }
      }
      return problems;
    }, [values, group, source] as [string, string, string]), {message: `the current table against the ${group} groups of ${source}`}).toEqual([]);
  }, {tier: 'api', description: 'a new table, current in its own view: one row per group of the source table, each holding the group\'s row count and the means (to one decimal) of the comma-separated columns, computed from the source'});

export const onlyGrid = Then('the current table view should show nothing but its grid', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => Array.from(grok.shell.tv.viewers as Iterable<{type: string}>).map((v) => v.type)),
    {message: 'the viewers of the current table view'}).toEqual(['Grid']);
}, {tier: 'api', description: 'grok.shell.tv.viewers holds the grid alone'});

export const hiddenInGrid = Then('the {string} column should be hidden in the grid', async (page: Page, column: string) => {
  await expect.poll(() => page.evaluate((c) => {
    const gc = grok.shell.tv.grid.col(c);
    return gc ? (gc.visible ? 'shown' : 'hidden') : `the grid has no "${c}" column`;
  }, column), {message: `the "${column}" column of the current table view's grid`}).toBe('hidden');
}, {tier: 'api', description: 'the grid column of the current table view is not visible; the table keeps the column'});

// --- the server -------------------------------------------------------------------------------------

export const fileCards = Then('the AI panel should show the files of {string} as cards', async (page: Page, folder: string) => {
  await expect.poll(() => page.evaluate(async ([sel, f]) => {
    const listed: {name: string; isDirectory: boolean}[] = await grok.dapi.files.list(f);
    const root = Array.from(document.querySelectorAll(sel)).find((p) => (p as HTMLElement).offsetParent !== null) as HTMLElement;
    const lists = root.querySelectorAll('.grokky-entity-cards');
    const items = Array.from(lists[lists.length - 1]?.children ?? []) as HTMLElement[];
    const cards = items.map((e) => e.querySelector('.grok-file-info label')?.textContent?.trim() ?? '').filter((n) => n !== '');
    const shown = items.map((e) => e.querySelector('.grok-file-info label')?.textContent?.trim() ?? e.innerText.trim());
    const files = listed.filter((x) => !x.isDirectory).map((x) => x.name);
    const entries = listed.map((x) => x.name);
    return [
      ...files.filter((n) => !cards.includes(n)).map((n) => `no card for ${n}`),
      ...shown.filter((n) => !entries.includes(n)).map((n) => `"${n}" is not in ${f}`),
      ...shown.filter((n, i) => shown.indexOf(n) !== i).map((n) => `${n} shown twice`),
    ];
  }, [AI_PANEL, folder] as [string, string]), {message: `the entity cards of the last answer against the server's listing of ${folder}`})
    .toEqual([]);
}, {tier: 'api', description: 'every file the server lists in the folder has a rendered file card in the last cards block, and nothing else shown is missing from the folder (its subfolders may show)'});

// --- the runtime container --------------------------------------------------------------------------

export const runsHermes = Then('the AI runtime container should run the Hermes agent and not the Claude Agent SDK', async (page: Page) => {
  const container = await runtimeContainer(page);
  const sdk = inRuntime(container, 'ls -d /app/node_modules/@anthropic-ai/claude-agent-sdk 2>/dev/null || echo none').trim();
  expect(runtimeProcesses(container).filter((p) => AGENT_PROCESS.test(p.cmd)).length,
    `Hermes agent processes (python -m grokky_agent) in ${container}`).toBeGreaterThanOrEqual(1);
  expect(sdk, `@anthropic-ai/claude-agent-sdk in ${container}'s /app/node_modules`).toBe('none');
}, {description: 'docker exec into the stand\'s runtime container (same machine): a Hermes agent process runs and the SDK is not installed'});

export const runningInRuntime = Then('the AI runtime container should be running {string}', async (page: Page, command: string) => {
  const container = await runtimeContainer(page);
  await expect.poll(() => runtimeProcesses(container).some((p) => p.cmd.includes(command)),
    {message: `a process "${command}" in ${container}`, timeout: pollMs(TURN_MS)}).toBe(true);
}, {description: 'polled over docker exec: a process whose command line contains the text'});

export const notRunningInRuntime = Then('the AI runtime container should not be running {string}', async (page: Page, command: string) => {
  const container = await runtimeContainer(page);
  await expect.poll(() => runtimeProcesses(container).filter((p) => p.cmd.includes(command)).map((p) => `${p.pid} ${p.cmd}`),
    {message: `processes "${command}" left in ${container}`, timeout: pollMs(15000)}).toEqual([]);
}, {description: 'polled over docker exec until no process\'s command line contains the text'});

export const sandboxSwept = Given('no file named {string} is in the AI\'s sandbox', async (page: Page, names: string) => {
  const container = await runtimeContainer(page);
  const dir = `/users/${await page.evaluate(() => String(grok.shell.user.id))}`;
  const wanted = names.split(/\s*,\s*/);
  const named = (): string[] => sandboxFiles(container, dir).filter((f) => wanted.includes(path.basename(f)));
  removeFiles(container, named());
  expect(named(), `files named ${names} in ${dir}`).toEqual([]);
  const before = new Set(sandboxFiles(container, dir));
  atFeatureEnd(page, async () => {
    removeFiles(container, sandboxFiles(container, dir).filter((f) => !before.has(f)));
    expect(sandboxFiles(container, dir).filter((f) => !before.has(f)), `files the AI left in ${dir}`).toEqual([]);
  });
}, {tier: 'api', description: 'the named files go from the user\'s runtime directory now, and every file added there after this step goes when the feature ends, read back gone; the synced agents/ folder is left alone'});

export const otherUsersNote = Given('another user\'s folder in the AI\'s sandbox holds a private note', async (page: Page) => {
  const container = await runtimeContainer(page);
  const sweep = (): string => inRuntime(container, `rm -rf ${CANARY_DIRS}*; ls -d ${CANARY_DIRS}* 2>/dev/null || true`).trim();
  expect(sweep(), 'canary folders an earlier run left').toBe('');
  const run = randomUUID();
  const dir = `${CANARY_DIRS}${run.slice(-12)}`;
  const note = `private-note-${run}`;
  // no chown -R: its own listing of the folder could keep relatime from recording the agent's
  inRuntime(container, `mkdir ${dir} && printf %s ${note} > ${dir}/notes.txt && chown --reference=/users ${dir} ${dir}/notes.txt`);
  canaries.set(page, {note, dir, accessed: inRuntime(container, `stat -c %x ${dir} ${dir}/notes.txt`).trim()});
  atFeatureEnd(page, async () => {
    canaries.delete(page);
    expect(sweep(), 'canary folders left in the runtime').toBe('');
  });
}, {tier: 'api', description: 'a folder under the runtime\'s /users, owned as the users\' folders are, holding a note with a random text; its access times are kept; removed at feature end and read back gone, a killed run\'s swept first'});

export const agentsKilled = When('the AI runtime\'s agent processes are killed', async (page: Page) => {
  const container = await runtimeContainer(page);
  const pids = runtimeProcesses(container).filter((p) => AGENT_PROCESS.test(p.cmd)).map((p) => p.pid);
  expect(pids.length, `agent processes (python -m grokky_agent) in ${container}`).toBeGreaterThanOrEqual(1);
  inRuntime(container, `kill -9 ${pids.join(' ')}`);
  await expect.poll(() => runtimeProcesses(container).filter((p) => pids.includes(p.pid)).length,
    {message: `the killed agent processes ${pids.join(', ')} still in ${container}`}).toBe(0);
}, {tier: 'api', description: 'docker exec … kill -9 on every python -m grokky_agent of the stand\'s runtime container (the chat\'s own and the warm spares), then waits until they are gone'});

export const canaryUnopened = Then('the other user\'s folder should not have been opened', async (page: Page) => {
  const {dir, accessed} = canaryOf(page);
  expect(inRuntime(await runtimeContainer(page), `stat -c %x ${dir} ${dir}/notes.txt`).trim(),
    `the access times of ${dir} and its note, against when it was planted (relatime moves them on the first listing or read)`).toBe(accessed);
}, {description: 'the folder was never listed and the note never read, by any process: their access times are as planted'});

// --- what reaches the page --------------------------------------------------------------------------

export const noSecrets = Then('no key-shaped text or session token should be on the page', async (page: Page) => {
  expect(await secretsOnPage(page), 'secrets on the page, in the AI panel or in an input').toEqual([]);
}, {description: 'API-key- and bearer-shaped text, and the session\'s own token, over the page text, the whole AI panel and every input'});

export const noteNotShown = Then('the other user\'s note should not be on the page', async (page: Page) => {
  const {note} = canaryOf(page);
  expect(await page.evaluate(([n, sel]) => (document.body.innerText + (document.querySelector(sel)?.textContent ?? '')).includes(n),
    [note, AI_PANEL] as [string, string]), 'the note\'s text on the page or in the AI panel').toBe(false);
}, {description: 'the random text of the note the canary folder holds, over the page text and the whole AI panel'});

export const screenshot = Then('a screenshot {string} of the AI panel is saved', async (page: Page, name: string) => {
  const secrets = await secretsOnPage(page);
  if (secrets.length)
    throw new Error(`${secrets.join(' and ')} on the page: screenshot ${name} not written`);
  fs.mkdirSync(SHOTS, {recursive: true});
  const full = path.join(SHOTS, `${name}.png`);
  const crop = path.join(SHOTS, `${name}-panel.png`);
  await page.screenshot({path: full});
  await page.locator(AI_PANEL).filter({visible: true}).first().screenshot({path: crop});
  await test.info().attach(`${name}.png`, {path: full, contentType: 'image/png'});
  await test.info().attach(`${name}-panel.png`, {path: crop, contentType: 'image/png'});
}, {description: 'the page and a crop of the AI panel to GROKKY_BDD_SHOTS (default bdd/test-results/screenshots), refused while a secret is shown'});
