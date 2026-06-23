/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

import {CombinedAISearchAssistant} from './search/combined-search';
import {fireAIPanelToggleEvent, getAIAbortSubscription, fireBeforeUserPromptEvent,
  fireAfterUserPromptEvent, UserPromptEventArgs, createStyledMarkdown, isEnterKey} from '../utils';
import {BuiltinDBInfoMeta} from '../db/query-meta-utils';
import {DBAIPanel, ScriptingAIPanel, ShellAIPanel, StreamingPanel, TVAIPanel} from './panel';
import {ClaudeRuntimeClient, ClaudeModel, ErrorEvent, FinalEvent, ToolActivityEvent} from '../claude/runtime-client';
import {executeSingleBlock, renderEntityBlocks} from '../claude/exec-blocks';
import {UsageLimiter} from './usage-limiter';
import {SQLGenerationContext} from '../db/sql-tools';
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
      renderEntityBlocks(contentHost);
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

export async function setupAIQueryEditorUI(v: DG.ViewBase, connectionID: string, queryEditorRoot: HTMLElement, setAndRunFunc: (query: string) => void): Promise<boolean> {
  if (!grok.ai.config.configured)
    return false;
  const connection = await grok.dapi.connections.find(connectionID);
  if (!connection) {
    grok.shell.error(`Connection with ID ${connectionID} not found.`);
    return false;
  }

  const allDbInfos = await BuiltinDBInfoMeta.allFromConnection(connection);
  const catalogs = allDbInfos.map((d) => d.name);
  const defaultCatalog = connection.parameters?.['catalog'] ?? connection.parameters?.['db'] ?? catalogs[0] ?? '';

  const panel = new DBAIPanel(catalogs, defaultCatalog, connectionID, v, setAndRunFunc);
  panel.show();

  let sqlContext: SQLGenerationContext | null = null;

  panel.onRunRequest.subscribe(async (args) => {
    if (!sqlContext)
      sqlContext = await SQLGenerationContext.create(connectionID, args.currentPrompt.catalogName);
    await runPromptWithLifecycle(panel, args.currentPrompt.prompt, v, 'db-query', (toolName, input) => sqlContext!.handleToolCall(toolName, input));
  });
  return true;
}

export async function runPromptWithLifecycle(
  panel: StreamingPanel,
  prompt: string,
  view: DG.ViewBase,
  quotaCategory: string,
  clientToolHandler?: (toolName: string, input: any) => Promise<string>,
  displayPrompt?: string,
): Promise<void> {
  if (prompt.startsWith('!!'))
    prompt = prompt.slice(1); // !!cmd → !cmd, falls through to normal pipeline
  else if (prompt.startsWith('!')) {
    const command = prompt.slice(1).trimStart();
    if (!await UsageLimiter.getInstance().tryCheckAndIncrement(quotaCategory, command))
      return;
    await runClaudeStreaming(panel, command, view, clientToolHandler, 'bash');
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
  await runClaudeStreaming(panel, prompt, view, clientToolHandler, undefined, session);
  if (!panel.rawRender)
    fireAfterUserPromptEvent({prompt, context: view, handled: false});
}

async function runClaudeStreaming(
  panel: StreamingPanel, userPrompt: string, view: DG.ViewBase,
  clientToolHandler?: (toolName: string, input: any) => Promise<string>,
  systemPromptMode?: string,
  existingSession?: ReturnType<StreamingPanel['startChatSession']>,
): Promise<void> {
  const session = existingSession ?? panel.startChatSession();
  if (!existingSession)
    session.session.addUserMessage({role: 'user', content: [{type: 'text', text: userPrompt}]}, userPrompt);
  try {
    await streamOnce(panel, userPrompt, view, clientToolHandler, systemPromptMode, session);
  } finally {
    session.endSession();
  }
}

async function streamOnce(
  panel: StreamingPanel, userPrompt: string, view: DG.ViewBase,
  clientToolHandler: ((toolName: string, input: any) => Promise<string>) | undefined,
  systemPromptMode: string | undefined,
  chatSession: ReturnType<StreamingPanel['startChatSession']>,
): Promise<void> {
  return new Promise<void>(async (resolve) => {
    const sessionId = panel.sessionId;
    let accumulated = ''; // full markdown — kept for session history
    let displayBuffer = ''; // shown during streaming (entity chunks excluded from display, kept in finalBuffer)
    let finalBuffer = ''; // entity chunks stay so renderEntityBlocks finds them at finalize
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

      await client.ensureConnected();

      forSession(client.onChunk, (evt) => {
        accumulated += evt.content;
        if (evt.kind !== 'entity') {
          displayBuffer += evt.content;
          toolStatus = '';
          panel.updateStreaming(displayBuffer, chatSession.loader);
        }
        finalBuffer += evt.content;
      });

      forSession(client.onToolActivity, (evt) => {
        toolStatus = `\n\n---\n**${evt.summary}**`;
        panel.updateStreaming(displayBuffer + toolStatus, chatSession.loader);
        chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: `[tool-activity] ${evt.summary}`}]});
      });

      forSession(client.onToolResult, (evt) => {
        // datagrok_exec results are internal — Claude's prose is the user-facing feedback.
        if (evt.toolName?.endsWith('datagrok_exec'))
          return;
        toolStatus = `\n\n---\n\`\`\`\n${evt.content}\n\`\`\``;
        panel.updateStreaming(displayBuffer + toolStatus, chatSession.loader);
        chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: `[tool-result] ${evt.content}`}]});
      });

      forSession(client.onFinal, async (evt) => {
        panel.cancelInputRequest();
        const exec = accumulated || evt.content;
        const display = finalBuffer || evt.content;
        chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: exec}]});
        await panel.finalizeStreaming(display, exec, view);
        cleanup();
        resolve();
      });

      forSession(client.onError, (evt) => {
        panel.cancelInputRequest();
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
          if (element)
            panel.appendStreamedElement(element);
          const result = error ?
            {success: false, error: error.error} :
            {success: true, ...(value != null ? {returnValue: value} : {})};
          client.respondToInput(sessionId, evt.requestId, result);
          return;
        }
        // Dispatch client-side DB tool calls (arrive as mcp__datagrok__<name>)
        const mcpToolName = evt.toolName.replace(/^mcp__datagrok__/, '');
        if (clientToolHandler && mcpToolName !== evt.toolName) {
          try {
            const result = await clientToolHandler(mcpToolName, evt.input);
            client.respondToInput(sessionId, evt.requestId, result);
          } catch (e: any) {
            client.respondToInput(sessionId, evt.requestId, `Error: ${e.message}`);
          }
          return;
        }
        // Default: AskUserQuestion
        accumulated = '';
        displayBuffer = '';
        finalBuffer = '';
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

let _shellAIPanel: ShellAIPanel | null = null;

export function setupShellAIPanelUI(): void {
  if (!grok.ai.config.configured) return;
  if (!_shellAIPanel) {
    _shellAIPanel = new ShellAIPanel();
    _shellAIPanel.onRunRequest.subscribe(async (args) => {
      if (args.currentPrompt.prompt.trim() === '/noprompt') {
        _shellAIPanel!.enableNoPrompt();
        return;
      }
      await runPromptWithLifecycle(_shellAIPanel!, args.currentPrompt.prompt, grok.shell.v, 'shell-ai');
    });
  }

  if (grok.shell.windows.ai.childElementCount === 0 || _shellAIPanel.isFrontPanel)
    _shellAIPanel.show(true);
}

const AI_ICON_SELECTOR = 'i[data-name="ai"]';

export async function setupTableViewAIPanelUI() {
  if (!grok.ai.config.configured)
    return;
  const handleView = (tableView: DG.TableView) => {
    if (tableView.root?.parentElement?.querySelector(AI_ICON_SELECTOR) != null)
      return;
    // setup ribbon panel icon
    const iconFse = ui.iconFA('user-robot', () => fireAIPanelToggleEvent(tableView), 'Ask AI \n Ctrl+I');
    iconFse.style.width = iconFse.style.height = '18px';
    tableView.setRibbonPanels([...tableView.getRibbonPanels(), [iconFse]]);
    // setup the panel itself
    const panel = new TVAIPanel(tableView);
    panel.hide();

    // Setup request handler
    panel.onRunRequest.subscribe(async (args) => {
      const prompt = args.currentPrompt.prompt;
      await runPromptWithLifecycle(panel, prompt, tableView, 'tableview');
    });
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

// TODO: rewrite to use Claude engine instead of deprecated script-tools
export async function setupScriptsAIPanelUI() {
  const handleView = (scriptView: DG.ScriptView) => {
    if (scriptView.root?.parentElement?.querySelector(AI_ICON_SELECTOR) != null)
      return;
    const iconFse = ui.iconSvg('ai.svg', () => fireAIPanelToggleEvent(scriptView), 'Ask AI \n Ctrl+I');
    iconFse.style.width = iconFse.style.height = '18px';
    scriptView.setRibbonPanels([...scriptView.getRibbonPanels(), [iconFse]]);
    const panel = new ScriptingAIPanel(scriptView);
    panel.hide();
    panel.onRunRequest.subscribe(async (args) => {
      await runPromptWithLifecycle(panel, args.currentPrompt.prompt, scriptView, 'scripting');
    });
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
    setupShellAIPanelUI();
    _shellAIPanel!.resetSession();
    const workflow = await _package.files.readAsText(`scripts/${name}.md`);
    const prompt =
      `Execute the following workflow. After each step, post a one-line status update to chat.\n\n` +
      `---\n${workflow}\n---`;
    const displayPrompt = `▶ Running workflow: ${name}`;
    await runPromptWithLifecycle(_shellAIPanel!, prompt, grok.shell.v, 'shell-ai', undefined, displayPrompt);
  } catch (e: any) {
    grok.shell.error(`Failed to run ${name}: ${e.message}`);
  }
}
