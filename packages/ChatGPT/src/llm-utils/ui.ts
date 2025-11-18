/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {askDeepGrok} from './deepwikiclient';
import {askOpenAIHelp} from './openAI-client';
import {LLMCredsManager} from './creds';

let searching = false;

export async function askDeepWiki(question: string, apiKey: string, vectorStoreId: string, useOpenAI: boolean = true) {
  if (searching)
    return;
  searching = true;
  try {
    const searchResultHost = document.getElementsByClassName('power-pack-search-host')[0];
    if (!searchResultHost)
      return;
    const spinner = ui.icons.spinner();
    spinner.style.color = 'var(--blue-1)';
    const searchingText = ui.divText('Grokking your question...', {style: {marginLeft: '8px'}});
    const loader = ui.divH([spinner, searchingText], {style: {alignItems: 'center', marginTop: '8px'}});
    const widgetDiv = ui.divV([loader], {style: {marginTop: '8px', marginBottom: '8px', alignItems: 'center', justifyContent: 'center'}});
    searchResultHost.prepend(widgetDiv);

    const res = !useOpenAI ? (await askDeepGrok(question)) : (await askOpenAIHelp(question, apiKey, vectorStoreId));
    loader.remove();
    const answerDiv = ui.markdown(res);
    answerDiv.style.width;
    answerDiv.style.userSelect = 'text';
    answerDiv.style.webkitUserSelect = 'text';
    answerDiv.style.width = '100%';
    widgetDiv.appendChild(answerDiv);
  } finally {
    searching = false;
  }
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
      searchHelpDiv.innerText = `Click Enter to search using AI. ${searchHelpDiv.innerText}`;

    searchInput.addEventListener('keydown', (event: KeyboardEvent) => {
      if (event.key === 'Enter' && searchInput.value?.trim()) {
        event.preventDefault();
        aiCombinedSearch(searchInput.value);
      }
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

function aiCombinedSearch(prompt: string) {
  console.log(`AI combined search for: ${prompt}`);
}
