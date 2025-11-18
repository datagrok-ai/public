/* eslint-disable max-len */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {askDeepGrok} from './deepwikiclient';
import {askOpenAIHelp} from './openAI-client';

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

export function setupSearchUI(getApiKey: () => string, getVectorStoreId: () => string) {
  const maxRetries = 10;
  function searchForSearchBox() {
    const searchBoxContainer = document.getElementsByClassName('power-pack-welcome-view')[0];
    if (!searchBoxContainer)
      return null;
    const searchInput = Array.from(searchBoxContainer.getElementsByTagName('input')).filter((el) => el.placeholder?.toLowerCase()?.includes('search everywhere'))[0];
    return searchInput;
  }

  function initInput(searchInput: HTMLInputElement) {
    const parent: HTMLElement = searchInput.parentElement!;
    const aiIcon = ui.iconFA('magic', () => { askDeepWiki(searchInput.value, getApiKey(), getVectorStoreId()); }, 'Ask AI');
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
