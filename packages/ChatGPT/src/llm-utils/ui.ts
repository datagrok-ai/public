/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {askDeepWiki} from './deepwikiclient';
import {askOpenAIHelp} from './openAI-client';
import {LLMCredsManager} from './creds';
import {CombinedAISearchAssistant} from './combined-search';
import {ChatGPTPromptEngine} from '../prompt-engine/prompt-engine';
import {ChatGptAssistant} from '../prompt-engine/chatgpt-assistant';
import * as api from '../package-api';
import {Plan} from '../prompt-engine/interfaces';
import {AssistantRenderer} from '../prompt-engine/rendering-tools';
import {AI_SQL_QUERY_ABORT_EVENT, dartLike} from '../utils';
import {generateAISqlQueryWithTools} from './sql-tools';
import {DBAIPanel, ModelType} from './panel';


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

export async function smartExecution(prompt: string, modelName: string) {
  const gptEngine = ChatGPTPromptEngine.getInstance(LLMCredsManager.getApiKey(), modelName);
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
  if (!LLMCredsManager.getApiKey()) {
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

export async function setupAIQueryEditorUI(connectionID: string, aiElement: HTMLElement, queryEditorRoot: HTMLElement, setAndRunFunc: (query: string) => void): Promise<void> {
  const connection = await grok.dapi.connections.find(connectionID);
  if (!connection) {
    grok.shell.error(`Connection with ID ${connectionID} not found.`);
    return;
  }

  const schemas = await grok.dapi.connections.getSchemas(connection);
  const defaultSchema = schemas.includes('public') ? 'public' : schemas[0];

  const panel = new DBAIPanel(schemas, defaultSchema);
  panel.show();

  const v = grok.shell.v?.root.contains(aiElement) ? grok.shell.v : null;

  if (v) {
    // do some subscriptions
    const sub = grok.events.onCurrentViewChanged.subscribe(() => {
      if (grok.shell.v != v)
        panel.hide();
      else
        panel.show();
    });
    const closeSub = grok.events.onViewRemoved.subscribe((view) => {
      if (view == v) {
        sub.unsubscribe();
        closeSub.unsubscribe();
        panel.hide();
        panel.dispose();
      }
    });
  }


  // temp hack
  aiElement.remove();
  const mutationObserver = new MutationObserver((_d) => {
    panel.show();
  });
  mutationObserver.observe(aiElement, {attributes: true, attributeFilter: ['style']});

  panel.onRunRequest.subscribe(async (args) => {
    ui.setUpdateIndicator(queryEditorRoot, true, 'Grokking Query...', () => { grok.events.fireCustomEvent(AI_SQL_QUERY_ABORT_EVENT, null); });
    const session = panel.startChatSession();
    try {
      const sqlQuery = await generateAISqlQueryWithTools(args.currentPrompt.prompt, connectionID, args.currentPrompt.schemaName!, {
        oldMessages: args.prevMessages,
        aiPanel: {
          addAIMessage: (aiMessage, title, content) => {
            session.addAIMessage(aiMessage, title, content);
          },
          addEngineMessage: (aiMessage) => session.addAIMessage(aiMessage, '', '', true),
          addUiMessage: (msg: string, fromUser: boolean) => session.addUIMessage(msg, fromUser),
          addUserMessage: (aiMsg, content) => session.addUserMessage(aiMsg, content),
        }, modelName: ModelType[panel.getCurrentInputs().model],
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
}
