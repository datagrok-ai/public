/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

import {CombinedAISearchAssistant} from './search/combined-search';
import {fireAIPanelToggleEvent, getAIAbortSubscription, fireBeforeUserPromptEvent,
  fireAfterUserPromptEvent, UserPromptEventArgs, createStyledMarkdown, isEnterKey, SHORTCUT_HINT,
  copyToClipboard} from '../utils';
import {AIPanel, StreamingPanel} from './panel';
import {AIWindowManager} from './ai-window';
import {ClaudeRuntimeClient, ClaudeModel, ErrorEvent, FinalEvent, ToolActivityEvent, AuthUrlEvent, AuthErrorEvent} from '../claude/runtime-client';
import {executeSingleBlock, runVerification, renderEntityRefList} from '../claude/exec-blocks';
import {UsageLimiter} from './usage-limiter';
import {collectViewAITools, NO_VIEW_TOOLS} from './view-tools';
import {resolveScopes, showSuggestionsMenu} from './prompt-suggestions';
import {_package} from '../package';

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
  sessionPrefix?: string;
  sessionId?: string;
  outputSchema?: object;
  onFinal: (evt: FinalEvent, host: HTMLElement) => void | Promise<void>;
}


export const RENDERED_EVENT = 'grokky-rendered';

async function streamingWidget(prompt: string, opts: StreamingOpts): Promise<DG.Widget> {
  const client = ClaudeRuntimeClient.getInstance();
  await client.ensureConnected();
  const sessionId = opts.sessionId ?? `${opts.sessionPrefix ?? 'stream'}-${crypto.randomUUID()}`;

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

  client.send(sessionId, prompt, {model: ClaudeModel.Sonnet, ...(opts.outputSchema ? {outputSchema: opts.outputSchema} : {})});
  return widget;
}

export async function smartExecution(prompt: string, sessionId?: string): Promise<DG.Widget> {
  return streamingWidget(prompt, {
    sessionId,
    sessionPrefix: 'execute',
    outputSchema: executionPlanSchema,
    onFinal: async (evt, host) => {
      const result = evt.structured_output as ExecutionPlan | undefined;
      if (!result) {
        host.appendChild(ui.divText(evt.content || 'No plan returned.'));
        return;
      }

      host.appendChild(ui.markdown(result.plan.map((s, i) => `${i + 1}. ${s}`).join('\n')));
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
    sessionPrefix: 'help',
    onFinal: (evt, host) => {
      const md = ui.markdown(evt.content || '');
      md.style.userSelect = 'text';
      host.appendChild(md);
    },
  });
}

// sets up the ui button for the input
export function setupSearchUI() {
  if (!grok.ai.config.configured) {
    console.warn('LLM API key is not set up. Search UI will not have AI assistance.');
    return;
  }

  const maxRetries = 10;
  function searchForSearchBox() {
    const searchBoxContainer = document.getElementsByClassName('power-pack-welcome-view')[0];
    if (!searchBoxContainer)
      return null;
    return Array.from(searchBoxContainer.getElementsByClassName('power-search-search-everywhere-input'))[0] as HTMLInputElement;
  }

  function initInput(searchInput: HTMLInputElement) {
    const parent: HTMLElement = searchInput.parentElement!;
    const wandIcon = ui.iconFA('magic', async (e) => {
      const scopes = await resolveScopes('powerSearch');
      showSuggestionsMenu(scopes, (prompt) => {
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

  const retries = 0;
  const intervalId = setInterval(() => {
    const searchInput = searchForSearchBox();
    if (searchInput) {
      clearInterval(intervalId);
      initInput(searchInput);
    } else if (retries >= maxRetries) {
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
// The query view needs no dedicated panel anymore: its AI tools (get_query_info,
// set_query_and_run + the SQL schema tools from the viewAIToolsProvider) are collected
// by the singleton panel at prompt time.
export async function setupAIQueryEditorUI(_v: DG.ViewBase, _connectionID: string, _queryEditorRoot: HTMLElement, _setAndRunFunc: (query: string) => void): Promise<boolean> {
  if (!grok.ai.config.configured)
    return false;
  initAIWindow();
  return true;
}

export async function runPromptWithLifecycle(
  panel: StreamingPanel,
  prompt: string,
  view: DG.ViewBase,
  quotaCategory: string,
  displayPrompt?: string,
): Promise<void> {
  if (prompt.startsWith('!!'))
    prompt = prompt.slice(1); // !!cmd → !cmd, falls through to normal pipeline
  else if (prompt.startsWith('!')) {
    const command = prompt.slice(1).trimStart();
    if (!await UsageLimiter.getInstance().tryCheckAndIncrement(quotaCategory, command))
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

    if (await grok.ai.processPrompt(prompt)) {
      // Built-in handler claimed the prompt — no AI response, just a green check next to it.
      session.session.markHandledNatively();
      session.endSession();
      panel.pushNativeContext(prompt);
      return;
    }
  }
  if (!await UsageLimiter.getInstance().tryCheckAndIncrement(quotaCategory, prompt)) {
    session?.endSession();
    return;
  }
  await runClaudeStreaming(panel, prompt, view, undefined, session);
  if (!panel.rawRender)
    fireAfterUserPromptEvent({prompt, context: view, handled: false});
}

function buildAuthRenewalWidget(client: ClaudeRuntimeClient): HTMLElement {
  const errorDiv = ui.divText('', 'grokky-auth-error');

  const submitBtn = ui.button('Submit', () => submitCode()) as HTMLButtonElement;
  submitBtn.disabled = true;

  const codeInput = ui.input.string('Authorization code', {
    placeholder: 'Paste code here',
    onValueChanged: () => {
      submitBtn.disabled = !codeInput.value?.trim();
      errorDiv.classList.add('grokky-auth-error');
    },
  });

  const openLink = ui.link('Open authorization page', () => {
    const sub = client.onAuthUrl.subscribe((evt) => {
      window.open(evt.url, '_blank');
      sub.unsubscribe();
    });
    client.startAuth();
  });

  const pendingStrip = ui.divV([
    ui.divText('Session expired'),
    ui.divText('Open the auth page and paste the code below'),
  ], 'grokky-auth-strip-pending');

  const successStrip = ui.divV([
    ui.divText('Session renewed'),
    ui.divText('Re-send your message to continue'),
  ], 'grokky-auth-strip-success');
  successStrip.style.display = 'none';

  const pendingBody = ui.divV([openLink, codeInput.root, errorDiv, submitBtn], 'grokky-auth-body');

  const widget = ui.div([pendingStrip, pendingBody, successStrip], 'grokky-auth-widget');

  const subs: {unsubscribe: () => void}[] = [];
  const cleanupSubs = () => subs.forEach((s) => s.unsubscribe());

  subs.push(
    client.onAuthDone.subscribe(() => {
      if (!widget.isConnected) {
        cleanupSubs();
        return;
      }
      pendingStrip.remove();
      pendingBody.remove();
      successStrip.style.display = '';
      cleanupSubs();
    }),
    client.onAuthError.subscribe((evt) => {
      if (!widget.isConnected) {
        cleanupSubs();
        return;
      }
      errorDiv.textContent = evt.message;
      errorDiv.classList.remove('grokky-auth-error');
      submitBtn.disabled = false;
      submitBtn.textContent = 'Submit';
      codeInput.input.focus();
    }),
  );

  function submitCode(): void {
    const code = codeInput.value?.trim();
    if (!code)
      return;
    submitBtn.disabled = true;
    submitBtn.textContent = 'Verifying…';
    client.sendAuthCode(code);
  }

  return widget;
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

async function streamOnce(
  panel: StreamingPanel, userPrompt: string, view: DG.ViewBase,
  systemPromptMode: string | undefined,
  chatSession: ReturnType<StreamingPanel['startChatSession']>,
): Promise<void> {
  return new Promise<void>(async (resolve) => {
    const sessionId = panel.sessionId;
    let accumulated = '';
    let segmentStart = 0;
    let toolStatus = '';
    let nextBlockIndex = 0;
    const subs: {unsubscribe: () => void}[] = [];
    const cleanup = () => subs.forEach((s) => s.unsubscribe());

    const forSession = <T extends {sessionId: string}>(
      source: {subscribe: (cb: (evt: T) => void) => {unsubscribe: () => void}},
      handler: (evt: T) => void,
    ) => subs.push(source.subscribe((evt) => {
        if (evt.sessionId === sessionId)
          handler(evt);
      }));

    const endWithError = (msg: string) => {
      panel.clearStreaming();
      grok.shell.error(msg);
      cleanup();
      resolve();
    };

    try {
      const client = ClaudeRuntimeClient.getInstance();
      const nativeCtx = panel.flushNativeContext();
      const enrichedUserPrompt = nativeCtx ? nativeCtx + userPrompt : userPrompt;
      const prompt = panel.rawRender ? enrichedUserPrompt : panel.prependViewContext(panel.prependEntityContext(enrichedUserPrompt), view);

      // Fresh per prompt: the current view may ship its own AI tools (natively or via
      // viewAIToolsProvider functions). Their defs go to the runtime; calls come back as input_request.
      const fullMode = !systemPromptMode && !panel.noPrompt;
      const viewTools = fullMode ? await collectViewAITools(grok.shell.v ?? view) : NO_VIEW_TOOLS;

      await client.ensureConnected();

      forSession(client.onChunk, (evt) => {
        accumulated += evt.content;
        toolStatus = '';
        panel.updateStreaming(accumulated.slice(segmentStart), chatSession.loader);
      });

      forSession(client.onToolActivity, (evt) => {
        toolStatus = `\n\n---\n**${evt.summary}**`;
        panel.updateStreaming(accumulated.slice(segmentStart) + toolStatus, chatSession.loader);
        chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: `[tool-activity] ${evt.summary}`}]});
      });


      forSession(client.onFinal, async (evt) => {
        panel.cancelInputRequest();
        const fullContent = accumulated || evt.content;
        if (/Failed to authenticate.*API Error: 401|authentication_error|\/login/i.test(fullContent)) {
          panel.clearStreaming();
          panel.appendStreamedElement(buildAuthRenewalWidget(client));
          cleanup();
          resolve();
          return;
        }
        const segmentContent = accumulated ? accumulated.slice(segmentStart) : fullContent;
        chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: fullContent}]});
        await panel.finalizeStreaming(segmentContent, fullContent, view);
        if (evt.unverified) {
          const warn = 'Not verified — the assistant could not confirm this action took effect.';
          panel.appendStreamedElement(ui.divText(warn, 'grokky-unverified-warning'));
        }
        cleanup();
        resolve();
      });

      forSession(client.onError, (evt) => {
        panel.cancelInputRequest();
        if (/401|authentication|credentials|\/login/i.test(evt.message)) {
          panel.clearStreaming();
          panel.appendStreamedElement(buildAuthRenewalWidget(client));
          cleanup();
          resolve();
          return;
        }
        endWithError(`Claude: ${evt.message}`);
      });

      forSession(client.onAborted, async () => {
        panel.cancelInputRequest();
        panel.clearStreaming();
        chatSession.session.addUiMessage('**Processing aborted by user**', false, {system: true});
        cleanup();
        resolve();
      });

      forSession(client.onInputRequest, async (evt) => {
        // datagrok_exec: run the JS here, return the outcome so Claude responds AFTER knowing it.
        if (evt.toolName === 'datagrok_exec') {
          const {element, value, error} = await executeSingleBlock(evt.input.code ?? '', view, nextBlockIndex++);
          if (element) {
            panel.appendStreamedElement(element);
            segmentStart = accumulated.length;
            toolStatus = '';
          }
          const result = error ?
            {success: false, error: error.error} :
            {success: true, ...(value != null ? {returnValue: value} : {})};
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
        // datagrok_show_entities: render entity cards immediately, no user interaction needed.
        if (evt.toolName === 'datagrok_show_entities') {
          panel.appendStreamedElement(renderEntityRefList(evt.input.entities ?? []));
          segmentStart = accumulated.length;
          toolStatus = '';
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
        toolStatus = '';
        panel.clearStreaming();
        const response = await panel.showInputRequest(evt.input);
        if (response) {
          client.respondToInput(sessionId, evt.requestId, response);
          const answerText = Object.values(response.answers).join(', ');
          chatSession.session.addEngineMessage({role: 'user', content: [{type: 'text', text: answerText}]});
        }
      });

      subs.push(client.onClose.subscribe(() => endWithError('Claude: connection lost')));

      subs.push(getAIAbortSubscription().subscribe(() => {
        client.abort(sessionId);
      }));

      const resolvedMode = systemPromptMode ?? (panel.noPrompt ? 'none' : undefined);
      client.send(sessionId, prompt, {
        ...(resolvedMode ? {systemPromptMode: resolvedMode} : {}),
        ...(viewTools.defs.length ? {clientTools: viewTools.defs} : {}),
      });
    } catch (e: any) {
      panel.clearStreaming();
      grok.shell.error(`Claude runtime: ${e.message}`);
      console.error('Claude runtime error:', e);
      cleanup();
      resolve();
    }
  });
}

let _shellAIPanel: AIPanel | null = null;

export function initAIWindow(): AIPanel | null {
  if (!grok.ai.config.configured)
    return null;
  if (!_shellAIPanel) {
    _shellAIPanel = new AIPanel('shell-ai-panel', null as any);
    _shellAIPanel.onRunRequest.subscribe(async (args) => {
      if (args.currentPrompt.prompt.trim() === '/noprompt') {
        _shellAIPanel!.enableNoPrompt();
        return;
      }
      await runPromptWithLifecycle(_shellAIPanel!, args.currentPrompt.prompt, grok.shell.v, 'shell-ai');
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
  if (!grok.ai.config.configured)
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
    if (view.type === 'ScriptView')
      setTimeout(() => handleView(view as DG.ScriptView), 500);
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
    const workflow = await _package.files.readAsText(`scripts/${name}.md`);
    const prompt =
      `Execute the following workflow. After each step, post a one-line status update to chat.\n\n` +
      `---\n${workflow}\n---`;
    const displayPrompt = `▶ Running workflow: ${name}`;
    await runPromptWithLifecycle(shell, prompt, grok.shell.v, 'shell-ai', displayPrompt);
  } catch (e: any) {
    grok.shell.error(`Failed to run ${name}: ${e.message}`);
  }
}
