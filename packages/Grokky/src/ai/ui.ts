/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

import {CombinedAISearchAssistant} from './search/combined-search';
import {fireAIPanelToggleEvent, getAIAbortSubscription, fireBeforeUserPromptEvent,
  fireAfterUserPromptEvent, UserPromptEventArgs, createStyledMarkdown, tryMarkdown, plainTextEl, isEnterKey,
  SHORTCUT_HINT, copyToClipboard} from '../utils';
import {AIPanel, StreamingPanel} from './panel';
import {AIWindowManager} from './ai-window';
import {ClaudeRuntimeClient, ErrorEvent, FinalEvent, ToolActivityEvent} from '../claude/runtime-client';
import {executeSingleBlock, runVerification, renderEntityRefList} from '../claude/exec-blocks';
import {sendChatTurn} from '../claude/queue-task';
import {UsageLimiter} from './usage-limiter';
import {viewFunctionTools, NO_VIEW_TOOLS} from './view-tools';
import {resolveScopes, showSuggestionsMenu} from './prompt-suggestions';
import {AiFeature} from '../claude/ai-feature';

interface ExecutionPlan {
  plan: string[];
  code: string;
}

const executionPlanSchema = {
  type: 'object',
  properties: {
    plan: {
      type: 'array',
      items: {type: 'string'},
      description: 'Steps to achieve the goal. Each item is a single step, plain text, without any leading number, bullet, or punctuation prefix.',
    },
    code: {
      type: 'string',
      description: 'Datagrok JavaScript code to execute. Uses globals: grok, ui, DG, view (current view), t (DataFrame if in TableView). Use await freely. Return an HTMLElement to display a result.',
    },
  },
  required: ['plan', 'code'],
  additionalProperties: false,
} as const;

const MAX_ACTIVITY_LINES = 5;

interface StreamingOpts {
  sessionPrefix: AiFeature;
  sessionId?: string;
  outputSchema?: object;
  onFinal: (evt: FinalEvent, host: HTMLElement) => void | Promise<void>;
}


export const RENDERED_EVENT = 'grokky-rendered';

async function streamingWidget(prompt: string, opts: StreamingOpts): Promise<DG.Widget> {
  const client = ClaudeRuntimeClient.getInstance();
  await client.ensureConnected();
  const sessionId = opts.sessionId ?? `${opts.sessionPrefix}-${crypto.randomUUID()}`;

  const activityHost = ui.divV([], 'grokky-stream-activity');
  const contentHost = ui.div([], 'grokky-stream-content');
  const container = ui.divV([activityHost, contentHost], 'grokky-stream');

  let onFirstEvent: () => void = () => {};
  const firstEvent = new Promise<void>((res) => { onFirstEvent = res; });

  const widget = DG.Widget.fromRoot(ui.wait(async () => {
    await firstEvent;
    return container;
  }));

  const signalRendered = () => grok.events.fireCustomEvent(RENDERED_EVENT, sessionId);

  const sessionSubs: rxjs.Subscription[] = [];
  const forSession = <T extends {sessionId: string}>(
    src: rxjs.Observable<T>, handler: (evt: T) => void,
  ) => {
    const sub = src.subscribe((evt) => {
      if (evt.sessionId === sessionId)
        handler(evt);
    });
    widget.subs.push(sub);
    sessionSubs.push(sub);
  };
  const stopListening = () => sessionSubs.forEach((s) => s.unsubscribe());

  forSession<ToolActivityEvent>(client.onToolActivity, (evt) => {
    activityHost.appendChild(ui.divText(evt.summary, 'grokky-stream-activity-line'));
    while (activityHost.children.length > MAX_ACTIVITY_LINES)
      activityHost.removeChild(activityHost.firstChild!);
    onFirstEvent();
  });

  forSession<FinalEvent>(client.onFinal, async (evt) => {
    activityHost.remove();
    onFirstEvent();
    try {
      await opts.onFinal(evt, contentHost);
    } finally {
      signalRendered();
      stopListening();
    }
  });

  forSession<ErrorEvent>(client.onError, (evt) => {
    activityHost.remove();
    contentHost.appendChild(ui.info(`Error: ${evt.message}`));
    onFirstEvent();
    signalRendered();
    stopListening();
  });

  client.send(sessionId, prompt, {logPrompt: prompt, outputSchema: opts.outputSchema});
  return widget;
}

export async function smartExecution(prompt: string, sessionId?: string): Promise<DG.Widget> {
  return streamingWidget(prompt, {
    sessionId,
    sessionPrefix: AiFeature.Execute,
    outputSchema: executionPlanSchema,
    onFinal: async (evt, host) => {
      const result = evt.structured_output as ExecutionPlan | undefined;
      if (!result) {
        host.appendChild(ui.divText(evt.content || 'No plan returned.'));
        return;
      }

      host.appendChild(ui.markdown(result.plan.map((s, i) => `${i + 1}. ${s}`).join('\n'), {sanitize: true}));
      if (!result.code)
        return;

      const codeAcc = ui.accordion();
      const fencedCode = '```datagrok-exec\n' + result.code + '\n```';
      codeAcc.addPane('Code', () => createStyledMarkdown(fencedCode), false);
      host.appendChild(codeAcc.root);

      let resolveExec: () => void = () => {};
      const execDone = new Promise<void>((res) => { resolveExec = res; });
      host.appendChild(ui.wait(async () => {
        try {
          const {element, error} = await executeSingleBlock(result.code, grok.shell.v, 0);
          if (error)
            return ui.info(`Execution error: ${error.error}`);
          return element ?? ui.divText('');
        } finally {
          resolveExec();
        }
      }));
      await execDone;
    },
  });
}

export async function askWiki(question: string, sessionId?: string): Promise<DG.Widget> {
  return streamingWidget(question, {
    sessionId,
    sessionPrefix: AiFeature.Help,
    onFinal: (evt, host) => {
      const text = evt.content || '';
      const md = tryMarkdown(text) ?? plainTextEl(text);
      md.style.userSelect = 'text';
      host.appendChild(md);
    },
  });
}

// sets up the ui button for the input
export function setupSearchUI() {
  if (!ClaudeRuntimeClient.getInstance().available) {
    console.warn('Claude runtime container not found. Search UI will not have AI assistance.');
    return;
  }

  const maxRetries = 10;
  function searchForSearchBox() {
    return document.querySelector<HTMLInputElement>('.power-search-search-everywhere-input');
  }

  function initInput(searchInput: HTMLInputElement) {
    const parent: HTMLElement = searchInput.parentElement!;
    const wandIcon = ui.iconFA('magic', async (e) => {
      const scopes = await resolveScopes('powerSearch');
      showSuggestionsMenu(scopes, (s) => {
        const prompt = s.prompt ?? '';
        searchInput.value = prompt;
        searchInput.dispatchEvent(new Event('input', {bubbles: true}));
        // PowerPack's input handler debounces 500ms and then calls ui.empty(host) on power-pack-search-host;
        // wait past that so the AI loader/result isn't wiped immediately after we prepend it.
        setTimeout(() => aiCombinedSearch(prompt), 600);
      }, e);
    }, 'Prompt suggestions');
    wandIcon.classList.add('grokky-search-wand');
    parent.insertBefore(wandIcon, searchInput);

    const searchHelpDiv = document.getElementsByClassName('power-search-help-text-container')[0] as HTMLDivElement;
    if (searchHelpDiv) {
      searchHelpDiv.innerText = `Press Enter to grok. ${searchHelpDiv.innerText}`;
      searchHelpDiv.style.display = 'none';
    }

    searchInput.addEventListener('keyup', (event: KeyboardEvent) => {
      if (isEnterKey(event) && searchInput.value?.trim()) {
        event.preventDefault();
        setTimeout(() => aiCombinedSearch(searchInput.value), 400); // timeout needed to allow other enter handlers to run first
      }
      if (searchHelpDiv)
        searchHelpDiv.style.display = searchInput.value?.trim().split(' ').length > 1 ? '' : 'none';
    });
  }

  let retries = 0;
  const intervalId = setInterval(() => {
    const searchInput = searchForSearchBox();
    if (searchInput) {
      clearInterval(intervalId);
      initInput(searchInput);
    } else if (++retries >= maxRetries) {
      clearInterval(intervalId);
      console.warn('Search input box not found after maximum retries.');
    }
  }, 1000);
}

export async function aiCombinedSearch(prompt: string) {
  // hide the menu
  document.querySelector('.d4-menu-popup')?.remove();
  await CombinedAISearchAssistant.instance.searchUI(prompt);
}

// Called by the core query editor to decide whether to show its AI toggle icon.
// The query view needs no dedicated panel anymore: its functions (getQueryInfo,
// setQueryAndRun + the meta.viewType-registered SQL schema functions) are reached
// by the singleton panel through the view-function meta-tools.
export async function setupAIQueryEditorUI(_v: DG.ViewBase, _connectionID: string, _queryEditorRoot: HTMLElement, _setAndRunFunc: (query: string) => void): Promise<boolean> {
  if (!ClaudeRuntimeClient.getInstance().available)
    return false;
  initAIWindow();
  return true;
}

async function withinLimit(panel: StreamingPanel, quotaCategory: string): Promise<boolean> {
  try {
    await UsageLimiter.getInstance().checkAndIncrement(quotaCategory);
    return true;
  } catch (e: any) {
    panel.showError(e.message);
    return false;
  }
}

export async function runPromptWithLifecycle(
  panel: StreamingPanel,
  prompt: string,
  view: DG.ViewBase,
  quotaCategory: string,
  displayPrompt?: string,
  isAnswer = false,
): Promise<void> {
  if (prompt.startsWith('!!'))
    prompt = prompt.slice(1); // !!cmd → !cmd, falls through to normal pipeline
  else if (prompt.startsWith('!')) {
    const command = prompt.slice(1).trimStart();
    if (!await withinLimit(panel, quotaCategory))
      return;
    await runClaudeStreaming(panel, command, view, 'bash');
    return;
  }
  let session: ReturnType<StreamingPanel['startChatSession']> | undefined;
  if (!panel.rawRender) {
    const args: UserPromptEventArgs = {prompt, context: view, handled: false};
    fireBeforeUserPromptEvent(args);
    if (args.handled)
      return;

    // Echo the user message before routing so it never disappears, even when
    // grok.ai.processPrompt's built-in handler claims the prompt.
    session = panel.startChatSession();
    if (displayPrompt) {
      session.session.addEngineMessage({role: 'user', content: [{type: 'text', text: prompt}]});
      session.session.addUiMessage(displayPrompt, false, {system: true});
    } else
      session.session.addUserMessage({role: 'user', content: [{type: 'text', text: prompt}]}, prompt);

    if (!isAnswer && await grok.ai.processPrompt(prompt)) {
      // Built-in handler claimed the prompt — no AI response, just a green check next to it.
      session.session.markHandledNatively();
      session.endSession();
      panel.pushNativeContext(prompt);
      return;
    }
  }
  if (!await withinLimit(panel, quotaCategory)) {
    session?.endSession();
    return;
  }
  await runClaudeStreaming(panel, prompt, view, undefined, session);
  if (!panel.rawRender)
    fireAfterUserPromptEvent({prompt, context: view, handled: false});
}

function buildConfirmCard(input: any, resolve: (confirmed: boolean) => void): HTMLElement {
  const answer = (confirmed: boolean) => {
    const [style, icon, text] = confirmed ? ['success', 'check', 'Confirmed'] : ['muted', 'times', 'Cancelled'];
    buttons.remove();
    strip.classList.replace('grokky-strip-danger', `grokky-strip-${style}`);
    strip.append(ui.divH([ui.iconFA(icon), ui.divText(text)], 'grokky-strip-answer'));
    resolve(confirmed);
  };
  const strip = ui.divH([ui.divText(input.action ?? 'Confirm action')], 'grokky-strip grokky-strip-danger');
  const buttons = ui.divH([ui.bigButton('OK', () => answer(true)), ui.button('CANCEL', () => answer(false))]);
  return ui.divV([
    strip,
    ui.divV([
      ui.tableFromMap({
        ...(input.name ? {Name: input.name} : {}),
        ...(input.owner ? {Owner: input.owner} : {}),
        ...(input.created ? {Created: input.created} : {}),
        ...(input.details ? {Details: input.details} : {}),
      }),
      buttons,
    ], 'grokky-confirm-body'),
  ], 'grokky-confirm-card');
}

async function runClaudeStreaming(
  panel: StreamingPanel, userPrompt: string, view: DG.ViewBase,
  systemPromptMode?: string,
  existingSession?: ReturnType<StreamingPanel['startChatSession']>,
): Promise<void> {
  const session = existingSession ?? panel.startChatSession();
  if (!existingSession)
    session.session.addUserMessage({role: 'user', content: [{type: 'text', text: userPrompt}]}, userPrompt);
  try {
    await streamOnce(panel, userPrompt, view, systemPromptMode, session);
  } finally {
    session.endSession();
  }
}

enum Status {
  None = '',
  Tool = 'tool',
  Writing = 'writing',
  Revising = 'revising',
}

async function streamOnce(
  panel: StreamingPanel, userPrompt: string, view: DG.ViewBase,
  systemPromptMode: string | undefined,
  chatSession: ReturnType<StreamingPanel['startChatSession']>,
): Promise<void> {
  return new Promise<void>(async (resolve) => {
    const sessionId = panel.sessionId;
    let accumulated = '';
    let segmentStart = 0;
    let status = Status.None;
    let nextBlockIndex = 0;
    const execCodes: string[] = [];
    let preRevisionExecCodes: string[] | null = null;
    const subs: {unsubscribe: () => void}[] = [];
    const cleanup = () => subs.forEach((s) => s.unsubscribe());
    let resentAfterReset = false;

    const forSession = <T extends {sessionId: string}>(
      source: {subscribe: (cb: (evt: T) => void) => {unsubscribe: () => void}},
      handler: (evt: T) => void,
    ) => subs.push(source.subscribe((evt) => {
        if (evt.sessionId === sessionId)
          handler(evt);
      }));

    const endWithError = (msg: string) => {
      panel.clearStreaming();
      panel.showError(msg);
      cleanup();
      resolve();
    };

    try {
      const client = ClaudeRuntimeClient.getInstance();
      const nativeCtx = panel.flushNativeContext();
      const enrichedUserPrompt = nativeCtx ? nativeCtx + userPrompt : userPrompt;
      const prompt = panel.rawRender ? enrichedUserPrompt : panel.prependViewContext(panel.prependEntityContext(enrichedUserPrompt), view);

      // Three static meta-tools let Claude search and invoke the current view's functions
      // (view.getFunctions()) without declaring them all. Defs are identical every turn
      // (prompt-cache friendly); the runners resolve the live grok.shell.v at call time.
      const fullMode = !systemPromptMode && !panel.noPrompt;
      const viewTools = fullMode ? viewFunctionTools() : NO_VIEW_TOOLS;

      await client.ensureConnected();

      // Text before and after a tool call arrives as separate assistant segments; without a
      // separator they concatenate mid-sentence ("…check those.No, there are…").
      let pendingSegmentBreak = false;
      forSession(client.onChunk, (evt) => {
        if (pendingSegmentBreak) {
          pendingSegmentBreak = false;
          if (accumulated && !/\n\s*$/.test(accumulated))
            accumulated += '\n\n';
        }
        accumulated += evt.content;
        if (status !== Status.Writing) {
          status = Status.Writing;
          panel.setStatus('Writing…');
        }
        panel.updateStreaming(accumulated.slice(segmentStart));
      });

      forSession(client.onToolActivity, (evt) => {
        pendingSegmentBreak = true;
        status = Status.Tool;
        panel.setStatus(evt.summary);
        chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: `[tool-activity] ${evt.summary}`}]});
      });

      let revisingTimer: number | null = null;
      forSession(client.onRevisionStart, () => {
        // A gate demanded a revision. The visible answer stays put; the revision streams hidden
        // (the runtime suppresses its chunks) and `final.revision` decides whether it replaces the
        // original. The status appears only if revising takes noticeable time — an immediate
        // NO_REVISION resolves before it ever shows.
        revisingTimer = window.setTimeout(() => {
          if (status === Status.Tool) return;
          status = Status.Revising;
          panel.setStatus('Revising…');
        }, 1200);
        preRevisionExecCodes = [...execCodes];
      });


      forSession(client.onFinal, async (evt) => {
        // Always end the turn: otherwise the panel stays busy and the chat is not saved.
        try {
          if (revisingTimer != null) {
            clearTimeout(revisingTimer);
            revisingTimer = null;
          }
          if (evt.revision === 'replaced') {
            // The revision supersedes the visible answer — render it in place of the original.
            accumulated = evt.content;
            segmentStart = 0;
          }
          const fullContent = accumulated || evt.content;
          const segmentContent = accumulated ? accumulated.slice(segmentStart) : fullContent;
          chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: fullContent}]});
          const finalExecCodes = evt.revision === 'kept' && preRevisionExecCodes ? preRevisionExecCodes : execCodes;
          await panel.finalizeStreaming(segmentContent, finalExecCodes, view);
          if (evt.unverified) {
            const warn = 'Not verified — the assistant could not confirm this action took effect.';
            panel.appendStreamedElement(ui.divText(warn, 'grokky-unverified-warning'));
          }
        } finally {
          cleanup();
          resolve();
        }
      });

      forSession(client.onError, (evt) => {
        endWithError(evt.message);
      });

      forSession(client.onAborted, async () => {
        panel.clearStreaming();
        chatSession.session.addUiMessage('Processing aborted by user', false, {system: true});
        cleanup();
        resolve();
      });

      forSession(client.onInputRequest, async (evt) => {
        // datagrok_exec: run the JS here, return the outcome so Claude responds AFTER knowing it.
        if (evt.toolName === 'datagrok_exec') {
          // Record the code that actually ran — this is what makes a history-restored
          // conversation reproducible ("do that again") instead of a text-only memory.
          chatSession.session.addEngineMessage({role: 'assistant',
            content: [{type: 'text', text: `[executed datagrok_exec]\n${(evt.input.code ?? '').slice(0, 1500)}`}]});
          const {element, value, error} = await executeSingleBlock(evt.input.code ?? '', view, nextBlockIndex++);
          if (evt.input.code && !error)
            execCodes.push(evt.input.code);
          if (element) {
            panel.appendStreamedElement(element);
            segmentStart = accumulated.length;
          }
          let result: any = error ?
            {success: false, error: error.error} :
            {success: true, ...(value != null ? {returnValue: value} : {})};
          // In-round-trip verification: run the provided assertion right after the action code, so a
          // passing action needs no separate datagrok_verify round-trip (see verify.ts in the runtime).
          if (!error && evt.input.verify?.assertion) {
            const v = await runVerification(evt.input.verify.assertion, view);
            result = {...result, verified: {passed: v.passed,
              ...(v.observed !== undefined ? {observed: v.observed} : {}), ...(v.error ? {error: v.error} : {})}};
          }
          client.respondToInput(sessionId, evt.requestId, result);
          return;
        }
        if (evt.toolName === 'datagrok_verify') {
          const {passed, observed, error} = await runVerification(evt.input.assertion ?? '', view);
          client.respondToInput(sessionId, evt.requestId, {
            passed,
            ...(observed !== undefined ? {observed} : {}),
            ...(error ? {error} : {}),
          });
          return;
        }
        // datagrok_confirm: render the approval card and hold the round-trip until a button is clicked.
        if (evt.toolName === 'datagrok_confirm') {
          const confirmed = await new Promise<boolean>((resolve) => {
            panel.appendStreamedElement(buildConfirmCard(evt.input ?? {}, resolve));
          });
          segmentStart = accumulated.length;
          client.respondToInput(sessionId, evt.requestId, {confirmed});
          return;
        }
        // datagrok_show_entities: render entity cards immediately, no user interaction needed.
        if (evt.toolName === 'datagrok_show_entities') {
          panel.appendStreamedElement(renderEntityRefList(evt.input.entities ?? []));
          segmentStart = accumulated.length;
          client.respondToInput(sessionId, evt.requestId, {success: true});
          return;
        }
        // View-declared tool call — run the tool's implementation from the current view.
        const runner = viewTools.runners.get(evt.toolName);
        if (runner) {
          try {
            const result = await runner(evt.input ?? {});
            client.respondToInput(sessionId, evt.requestId, result === undefined ? {success: true} : result);
          } catch (e: any) {
            client.respondToInput(sessionId, evt.requestId, {success: false, error: e.message});
          }
          return;
        }
        // Default: AskUserQuestion
        accumulated = '';
        segmentStart = 0;
        panel.clearStreaming();
        panel.showInputRequest(evt.input);
      });

      subs.push(client.onClose.subscribe(() =>
        endWithError('Connection lost, the request was stopped. Send it again.')));

      subs.push(getAIAbortSubscription().subscribe(() => {
        client.abort(sessionId);
      }));

      const resolvedMode = systemPromptMode ?? (panel.noPrompt ? 'none' : undefined);
      const sendPrompt = () => {
        const transcript = resolvedMode === 'bash' || client.isResumable(sessionId) ? '' : panel.restoredTranscript();
        const message = transcript ? transcript + '\n---\n\n' + prompt : prompt;
        const options = {
          ...(resolvedMode ? {systemPromptMode: resolvedMode} : {}),
          ...(viewTools.defs.length ? {clientTools: viewTools.defs} : {}),
          logPrompt: userPrompt,
        };
        // Queued-task admission (queue-task.ts): the queued call only holds the turn's
        // admission slot; the turn itself streams over this socket as usual.
        sendChatTurn(client, sessionId, message, options)
          .catch((e: any) => endWithError(e?.message ?? String(e)));
      };

      forSession(client.onSessionReset, () => {
        if (resentAfterReset)
          return endWithError('The session was lost and could not be restored. Start a new chat.');
        resentAfterReset = true;
        sendPrompt();
      });

      sendPrompt();
    } catch (e: any) {
      console.error('Claude runtime error:', e);
      endWithError(e.message);
    }
  });
}

let _shellAIPanel: AIPanel | null = null;

export function initAIWindow(): AIPanel | null {
  if (!ClaudeRuntimeClient.getInstance().available)
    return null;
  if (!_shellAIPanel) {
    _shellAIPanel = new AIPanel('shell-ai-panel', null as any);
    _shellAIPanel.onRunRequest.subscribe(async (args) => {
      if (args.currentPrompt.prompt.trim() === '/noprompt') {
        _shellAIPanel!.enableNoPrompt();
        return;
      }
      await runPromptWithLifecycle(_shellAIPanel!, args.currentPrompt.prompt, grok.shell.v, 'shell-ai', undefined, args.currentPrompt.isAnswer);
    });
    AIWindowManager.instance.init(_shellAIPanel);
  }
  return _shellAIPanel;
}

export function setupShellAIPanelUI(): void {
  if (!initAIWindow())
    return;
  AIWindowManager.instance.show();
}

const AI_ICON_SELECTOR = 'i[data-name="ai"]';

export async function setupTableViewAIPanelUI() {
  if (!ClaudeRuntimeClient.getInstance().available)
    return;
  const handleView = (tableView: DG.TableView) => {
    if (tableView.root?.parentElement?.querySelector(AI_ICON_SELECTOR) != null)
      return;
    // Ribbon icon toggles the shared singleton panel; table context is picked up at prompt time.
    const iconFse = ui.iconFA('user-robot', () => fireAIPanelToggleEvent(tableView), `Ask AI \n ${SHORTCUT_HINT}`);
    iconFse.style.width = iconFse.style.height = '18px';
    tableView.setRibbonPanels([...tableView.getRibbonPanels(), [iconFse]]);
  };
  // also handle already opened views
  Array.from(grok.shell.tableViews).filter((v) => v.dataFrame != null).forEach((view) => {
    handleView(view as DG.TableView);
  });

  grok.events.onViewAdded.subscribe((view) => {
    if (view.type === DG.VIEW_TYPE.TABLE_VIEW)
      handleView(view as DG.TableView);
  });
}

export async function setupScriptsAIPanelUI() {
  const handleView = (scriptView: DG.ScriptView) => {
    if (scriptView.root?.parentElement?.querySelector(AI_ICON_SELECTOR) != null)
      return;
    // The singleton panel picks up the script view's AI tools (get/set code) at prompt time.
    const iconFse = ui.iconSvg('ai.svg', () => fireAIPanelToggleEvent(scriptView), `Ask AI \n ${SHORTCUT_HINT}`);
    iconFse.style.width = iconFse.style.height = '18px';
    scriptView.setRibbonPanels([...scriptView.getRibbonPanels(), [iconFse]]);
  };

  grok.events.onViewAdded.subscribe((view) => {
    if (view.type === 'ScriptView') {
      const scriptView = view as DG.ScriptView;
      setTimeout(() => {
        handleView(scriptView);
        scriptView.subs.push(scriptView.tabs.onTabChanged.subscribe(() => handleView(scriptView)));
      }, 500);
    }
  });
}

// Adds a "Run" button to the ribbon of file views opened from MyFiles/agents/scripts/.
export function setupAgentScriptsUI(): void {
  grok.events.onViewAdded.subscribe((view) => {
    try {
      const basePath = view.basePath ?? view.path ?? '';
      if (!basePath.includes('/agents/scripts/'))
        return;
      const {name} = view;
      const runIcon = ui.iconFA('play', () => runAgentScript(name), `Run ${name}`);
      runIcon.classList.add('fas');
      view.setRibbonPanels([...view.getRibbonPanels(), [runIcon]]);
    } catch (e: any) {
      console.warn('Grokky: failed to add run button:', e.message);
    }
  });
}

async function runAgentScript(name: string): Promise<void> {
  try {
    const shell = initAIWindow();
    if (!shell)
      return;
    AIWindowManager.instance.showPanel(shell);
    shell.resetSession();
    const conn = await grok.dapi.connections.filter('name = "My files"').first();
    const workflow = await grok.dapi.files.readAsText(`${conn.nqName}/agents/scripts/${name}.md`);
    const prompt =
      `Execute the following workflow. After each step, post a one-line status update to chat.\n\n` +
      `---\n${workflow}\n---`;
    const displayPrompt = `▶ Running workflow: ${name}`;
    await runPromptWithLifecycle(shell, prompt, grok.shell.v, 'shell-ai', displayPrompt);
  } catch (e: any) {
    grok.shell.error(`Failed to run ${name}: ${e.message}`);
  }
}
