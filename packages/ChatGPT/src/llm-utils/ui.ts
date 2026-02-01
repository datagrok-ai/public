/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {askDeepWiki} from './deepwikiclient';
import {askOpenAIHelp, ModelType} from './openAI-client';
import {LLMCredsManager} from './creds';
import {CombinedAISearchAssistant} from './combined-search';
import {ChatGPTPromptEngine} from '../prompt-engine/prompt-engine';
import {ChatGptAssistant} from '../prompt-engine/chatgpt-assistant';
import * as api from '../package-api';
import {Plan} from '../prompt-engine/interfaces';
import {AssistantRenderer} from '../prompt-engine/rendering-tools';
import {fireAIAbortEvent, fireAIPanelToggleEvent} from '../utils';
import {generateAISqlQueryWithTools} from './sql-tools';
import {processTableViewAIRequest} from './tableview-tools';
import {DBAIPanel, ScriptingAIPanel, TVAIPanel} from './panel';
import {generateDatagrokScript} from './script-tools';
import {ChatModel} from 'openai/resources/shared';


export async function askWiki(question: string, useOpenAI: boolean = true) {
  try {
    const res = !useOpenAI ? (await askDeepWiki(question)) : (await askOpenAIHelp(question));
    const markdown = ui.markdown(res);
    markdown.style.userSelect = 'text';
    return DG.Widget.fromRoot(markdown);
  } catch (error: any) {
    console.error('Error during DeepGROK ask:', error);
    return DG.Widget.fromRoot(ui.divText(`Error during DeepGROK ask: ${error.message}`));
  }
}

export async function smartExecution(prompt: string, modelName: ChatModel) {
  const gptEngine = ChatGPTPromptEngine.getInstance(modelName);
  const gptAssistant = new ChatGptAssistant(gptEngine);

  const mainWaitDiv = ui.wait(async () => {
    const resDiv = ui.divV([], 'chatgpt-ask-ai-result');
    const plan = JSON.parse(await api.funcs.getExecutionPlan(prompt)) as Plan;
    const planDiv = AssistantRenderer.renderPlan(plan);
    resDiv.appendChild(planDiv);
    resDiv.appendChild(ui.wait(async () => {
      const result = await gptAssistant.execute(plan);
      const resRen = AssistantRenderer.renderResult(result);
      return resRen;
    }));
    return resDiv;
  });
  return DG.Widget.fromRoot(mainWaitDiv);
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
    const searchInput = Array.from(searchBoxContainer.getElementsByClassName('power-search-search-everywhere-input'))[0] as HTMLInputElement;
    return searchInput;
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
    if (searchHelpDiv)
      searchHelpDiv.innerText = `Press Enter to grok. ${searchHelpDiv.innerText}`;

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
  if (!connection) { // should not happen but just in case
    grok.shell.error(`Connection with ID ${connectionID} not found.`);
    return false;
  }

  const schemas = await grok.dapi.connections.getSchemas(connection);
  const defaultSchema = schemas.includes('public') ? 'public' : schemas[0];

  const panel = new DBAIPanel(schemas, defaultSchema, connectionID, v);
  panel.show();

  panel.onRunRequest.subscribe(async (args) => {
    ui.setUpdateIndicator(queryEditorRoot, true, 'Grokking Query...', () => { fireAIAbortEvent(); });
    const session = panel.startChatSession();
    try {
      const sqlQuery = await generateAISqlQueryWithTools(args.currentPrompt.prompt, connectionID, args.currentPrompt.schemaName!, {
        oldMessages: args.prevMessages,
        aiPanel: session.session, modelType: panel.getCurrentInputs().model,
      });
      if (sqlQuery && typeof sqlQuery === 'string' && sqlQuery.trim().length > 0)
        setAndRunFunc(sqlQuery);
    } catch (error: any) {
      ui.setUpdateIndicator(queryEditorRoot, false);
      grok.shell.error(`Error during AI query generation`);
      console.error('Error during AI query generation:', error);
    }
    ui.setUpdateIndicator(queryEditorRoot, false);
    session.endSession();
  });
  return true;
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
      const session = panel.startChatSession();
      try {
        await processTableViewAIRequest(
          args.currentPrompt.prompt,
          tableView,
          {
            oldMessages: args.prevMessages,
            aiPanel: session.session,
            modelType: panel.getCurrentInputs().model,
          }
        );
      } catch (error: any) {
        grok.shell.error('Error during AI table view processing');
        console.error('Error during AI table view processing:', error);
      }
      session.endSession();
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

export async function setupScriptsAIPanelUI() {
  const handleView = (scriptView: DG.ScriptView) => {
    // setup ribbon panel icon
    const iconFse = ui.iconSvg('ai.svg', () => fireAIPanelToggleEvent(scriptView), 'Ask AI \n Ctrl+I');
    iconFse.style.width = iconFse.style.height = '18px';
    scriptView.setRibbonPanels([...scriptView.getRibbonPanels(), [iconFse]]);
    // Setup request handler
    // setup the panel itself
    const panel = new ScriptingAIPanel(scriptView);
    panel.hide();
    panel.onRunRequest.subscribe(async (args) => {
      ui.setUpdateIndicator(scriptView.root, true, 'Vibe-Grokking Script...', () => { fireAIAbortEvent(); });
      const indicator = scriptView.root.querySelector('.d4-update-shadow') as HTMLElement;
      if (indicator)
        indicator.style.zIndex = '1000';
      const session = panel.startChatSession();
      try {
        const res = await generateDatagrokScript(
          args.currentPrompt.prompt,
          panel.getCurrentInputs().language,
          {
            oldMessages: args.prevMessages,
            aiPanel: session.session
          }
        );
        if (res && res.trim().length > 0)
          scriptView.code = res;
        console.log('Generated script:', res);
      } catch (error: any) {
        grok.shell.error('Error during AI table view processing');
        console.error('Error during AI table view processing:', error);
      }
      ui.setUpdateIndicator(scriptView.root, false);
      session.endSession();
    });
  };

  grok.events.onViewAdded.subscribe((view) => {
    if (view.type === 'ScriptView')
      setTimeout(() => handleView(view as DG.ScriptView), 500);
  });
}
