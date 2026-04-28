/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

import {CombinedAISearchAssistant} from './search/combined-search';
import {fireAIPanelToggleEvent, getAIAbortSubscription, fireBeforeUserPromptEvent,
  fireAfterUserPromptEvent, UserPromptEventArgs, createStyledMarkdown} from '../utils';
import {BuiltinDBInfoMeta} from '../db/query-meta-utils';
import {DBAIPanel, ScriptingAIPanel, ShellAIPanel, StreamingPanel, TVAIPanel} from './panel';
import {ClaudeRuntimeClient, ErrorEvent, FinalEvent, ToolActivityEvent} from '../claude/runtime-client';
import {executeDatagrokBlocks} from '../claude/exec-blocks';
import {UsageLimiter} from './usage-limiter';
import {SQLGenerationContext} from '../db/sql-tools';

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

export async function smartExecution(prompt: string): Promise<DG.Widget> {
  const client = ClaudeRuntimeClient.getInstance();
  await client.ensureConnected();
  const sessionId = `execute-${crypto.randomUUID()}`;

  const activityHost = ui.divV([], 'grokky-execute-activity');
  const container = ui.divV([activityHost], 'grokky-execute');

  let onFirstEvent: () => void = () => {};
  const firstEvent = new Promise<void>((res) => { onFirstEvent = res; });

  const widget = DG.Widget.fromRoot(ui.wait(async () => {
    await firstEvent;
    return container;
  }));

  const forSession = <T extends {sessionId: string}>(
    src: rxjs.Observable<T>, handler: (evt: T) => void,
  ) => widget.subs.push(src.subscribe((evt) => {
      if (evt.sessionId === sessionId)
        handler(evt);
    }));

  forSession<ToolActivityEvent>(client.onToolActivity, (evt) => {
    activityHost.appendChild(ui.divText(evt.summary, 'grokky-execute-activity-line'));
    while (activityHost.children.length > MAX_ACTIVITY_LINES)
      activityHost.removeChild(activityHost.firstChild!);
    onFirstEvent();
  });

  forSession<FinalEvent>(client.onFinal, async (evt) => {
    activityHost.remove();
    onFirstEvent();
    const result = evt.structured_output as ExecutionPlan | undefined;
    if (!result) {
      container.appendChild(ui.divText(evt.content || 'No plan returned.'));
      return;
    }

    container.appendChild(ui.markdown(result.plan.map((s, i) => `${i + 1}. ${s}`).join('\n')));
    if (!result.code)
      return;

    const codeAcc = ui.accordion();
    const code = '```datagrok-exec\n' + result.code + '\n```';
    codeAcc.addPane('Code', () => createStyledMarkdown(code), false);
    container.appendChild(codeAcc.root);

    const spinner = ui.loader();
    container.appendChild(spinner);
    try {
      const execResults = await executeDatagrokBlocks(code, grok.shell.v);
      container.append(...execResults);
    } catch (e: any) {
      container.appendChild(ui.info(`Execution error: ${e?.message ?? e}`));
    } finally {
      spinner.remove();
    }
  });

  forSession<ErrorEvent>(client.onError, (evt) => {
    activityHost.remove();
    container.appendChild(ui.info(`Error: ${evt.message}`));
    onFirstEvent();
  });

  client.send(sessionId, prompt, {outputSchema: executionPlanSchema});
  return widget;
}

export async function askWiki(question: string) {
  try {
    const res = await ClaudeRuntimeClient.getInstance().query(question);
    const markdown = ui.markdown(res);
    markdown.style.userSelect = 'text';
    return DG.Widget.fromRoot(markdown);
  } catch (error: any) {
    console.error('Error during AI help:', error);
    return DG.Widget.fromRoot(ui.divText(`Error during AI help: ${error.message}`));
  }
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
    const aiIcon = ui.iconFA('magic', () => { aiCombinedSearch(searchInput.value); }, 'Ask AI');
    aiIcon.style.cursor = 'pointer';
    aiIcon.style.color = 'var(--blue-1)';
    aiIcon.style.marginLeft = '4px';
    aiIcon.style.marginRight = '4px';
    parent.appendChild(aiIcon);

    const onSearchChanged = () => {
      if (!searchInput.value?.trim())
        aiIcon.style.display = 'none';
      else
        aiIcon.style.display = 'flex';
    };
    searchInput.addEventListener('input', onSearchChanged);
    // set up the enter key listener and modify suggestion
    const searchHelpDiv = document.getElementsByClassName('power-search-help-text-container')[0] as HTMLDivElement;
    if (searchHelpDiv) {
      searchHelpDiv.innerText = `Press Enter to grok. ${searchHelpDiv.innerText}`;
      searchHelpDiv.style.visibility = 'hidden';
    }

    searchInput.addEventListener('keyup', (event: KeyboardEvent) => {
      if (event.key === 'Enter' && searchInput.value?.trim()) {
        event.preventDefault();
        setTimeout(() => aiCombinedSearch(searchInput.value), 400); // timeout needed to allow other enter handlers to run first
      }
      if (searchHelpDiv)
        searchHelpDiv.style.visibility = searchInput.value?.trim().split(' ').length > 1 ? 'visible' : 'hidden';
    });

    onSearchChanged();
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

async function aiCombinedSearch(prompt: string) {
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

async function runPromptWithLifecycle(
  panel: StreamingPanel,
  prompt: string,
  view: DG.ViewBase,
  quotaCategory: string,
  clientToolHandler?: (toolName: string, input: any) => Promise<string>,
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
  if (!panel.rawRender) {
    const args: UserPromptEventArgs = {prompt, context: view, handled: false};
    fireBeforeUserPromptEvent(args);
    if (args.handled)
      return;
    if (await grok.ai.processPrompt(prompt))
      return;
  }
  if (!await UsageLimiter.getInstance().tryCheckAndIncrement(quotaCategory, prompt))
    return;
  await runClaudeStreaming(panel, prompt, view, clientToolHandler);
  if (!panel.rawRender)
    fireAfterUserPromptEvent({prompt, context: view, handled: false});
}

async function runClaudeStreaming(panel: StreamingPanel, userPrompt: string, view: DG.ViewBase, clientToolHandler?: (toolName: string, input: any) => Promise<string>, systemPromptMode?: string) {
  const chatSession = panel.startChatSession();
  const sessionId = panel.sessionId;
  let accumulated = '';
  let toolStatus = '';
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
    chatSession.endSession();
    cleanup();
  };

  try {
    const client = ClaudeRuntimeClient.getInstance();
    const prompt = panel.rawRender ? userPrompt : panel.prependViewContext(userPrompt, view);

    await client.ensureConnected();

    chatSession.session.addUserMessage({role: 'user', content: [{type: 'text', text: userPrompt}]}, userPrompt);

    forSession(client.onChunk, (evt) => {
      accumulated += evt.content;
      toolStatus = '';
      panel.updateStreaming(accumulated, chatSession.loader);
    });

    forSession(client.onToolActivity, (evt) => {
      toolStatus = `\n\n---\n**${evt.summary}**`;
      panel.updateStreaming(accumulated + toolStatus, chatSession.loader);
      chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: `[tool-activity] ${evt.summary}`}]});
    });

    forSession(client.onToolResult, (evt) => {
      toolStatus = `\n\n---\n\`\`\`\n${evt.content}\n\`\`\``;
      panel.updateStreaming(accumulated + toolStatus, chatSession.loader);
      chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: `[tool-result] ${evt.content}`}]});
    });

    forSession(client.onFinal, (evt) => {
      panel.cancelInputRequest();
      chatSession.session.addEngineMessage({role: 'assistant', content: [{type: 'text', text: evt.content}]});
      panel.finalizeStreaming(evt.content, view);
      chatSession.endSession();
      cleanup();
    });

    forSession(client.onError, (evt) => {
      panel.cancelInputRequest();
      endWithError(`Claude: ${evt.message}`);
    });

    forSession(client.onAborted, () => {
      panel.cancelInputRequest();
      panel.clearStreaming();
      chatSession.session.addUiMessage('**Processing aborted by user**', false);
      chatSession.endSession();
      cleanup();
    });

    forSession(client.onInputRequest, async (evt) => {
      // Dispatch client-side DB tool calls (arrive as mcp__datagrok__<name>)
      const mcpToolName = evt.toolName.replace(/^mcp__datagrok__/, '');
      if (clientToolHandler && mcpToolName !== evt.toolName) {
        try {
          const result = await clientToolHandler(mcpToolName, evt.input);
          client.respondToInput(sessionId, result);
        } catch (e: any) {
          client.respondToInput(sessionId, `Error: ${e.message}`);
        }
        return;
      }
      // Default: AskUserQuestion
      accumulated = '';
      toolStatus = '';
      panel.clearStreaming();
      const response = await panel.showInputRequest(evt.input);
      if (response) {
        client.respondToInput(sessionId, response);
        const answerText = Object.values(response.answers).join(', ');
        chatSession.session.addEngineMessage({role: 'user', content: [{type: 'text', text: answerText}]});
      }
    });

    subs.push(client.onClose.subscribe(() => endWithError('Claude: connection lost')));

    subs.push(getAIAbortSubscription().subscribe(() => {
      client.abort(sessionId);
    }));

    const resolvedMode = systemPromptMode ?? (panel.noPrompt ? 'none' : undefined);
    client.send(sessionId, prompt, resolvedMode ? {systemPromptMode: resolvedMode} : undefined);
  } catch (e: any) {
    panel.clearStreaming();
    grok.shell.error(`Claude runtime: ${e.message}`);
    console.error('Claude runtime error:', e);
    chatSession.endSession();
    cleanup();
  }
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
  _shellAIPanel.show();
}

export async function setupTableViewAIPanelUI() {
  if (!grok.ai.config.configured)
    return;
  const handleView = (tableView: DG.TableView) => {
    // setup ribbon panel icon
    const iconFse = ui.iconSvg('ai.svg', () => fireAIPanelToggleEvent(tableView), 'Ask AI \n Ctrl+I');
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
