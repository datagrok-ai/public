/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChatGptAssistant} from './prompt-engine/chatgpt-assistant';
import {ChatGPTPromptEngine} from './prompt-engine/prompt-engine';
import {AssistantRenderer} from './prompt-engine/rendering-tools';
import {getAiPanelVisibility, initAiPanel, setAiPanelVisibility} from './ai-panel';
import {findBestFunction} from './prompts/find-best-function';
import {askDeepGrok} from './llm-utils/deepwikiclient';
import {askDeepWiki} from './llm-utils/ui';
import {ChatGptFuncParams, Plan} from './prompt-engine/interfaces';
import {OpenAIHelpClient} from './llm-utils/openAI-client';

export * from './package.g';
export const _package = new DG.Package();

export let apiKey: string = '';
export let vectorStoreId: string = '';
export const modelName: string = 'gpt-4.1-mini-2025-04-14';

export class PackageFunctions {
  @grok.decorators.init()
  static async init() {
    apiKey = _package.settings['apiKey'];
    vectorStoreId = _package.settings['vectorStoreId'];
  }

  @grok.decorators.func()
  static async deepDemo(@grok.decorators.param({name: 'question', type: 'string'}) question: string) {
    // if (!apiKey || apiKey.length === 0) {
    //   grok.shell.error('Please set API key in ChatGPT package settings');
    //   return;
    // }
    return await askDeepGrok(question);
  }


  @grok.decorators.autostart()
  static autostart() {
    // setupSearchUI(() => apiKey, () => vectorStoreId);
    // grok.shell.info('started')
    //
    grok.events.onViewAdded.subscribe((view) => {
      if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
        const tableView = view as DG.TableView;
        const iconFse = ui.iconSvg('ai.svg', () => setAiPanelVisibility(true), 'Ask AI');
        tableView.setRibbonPanels([...tableView.getRibbonPanels(), [iconFse]]);
      }
    });

    initAiPanel();
    // Add keyboard shortcut for toggling AI panel
    document.addEventListener('keydown', (event) => {
      // Check for Ctrl+I (Ctrl key + I key)
      if (event.ctrlKey && event.key === 'i') {
        event.preventDefault(); // Prevent default browser behavior
        const isVisible = getAiPanelVisibility();
        setAiPanelVisibility(!isVisible);
      }
    });
  }

  @grok.decorators.func({tags: ['searchProvider']})
  static askHelpLLMProvider(): DG.SearchProvider {
    const provider: DG.SearchProvider = {
      'home': {
        name: 'Ask DeepGROK AI',
        search: async (_query: string) => null,
        getSuggestions: (_query: string) => [
          {
            priority: 9,
            suggestionText: 'Ask Documentation'
          }
        ],
        isApplicable: (query: string) => query.length >= 5,
        onValueEnter: async (query: string) => {
          await askDeepWiki(query, apiKey, vectorStoreId);
        },
        description: 'Get answers from DeepGROK AI assistant based on Datagrok documentation.'

      }
    };
    return provider;
  }

 @grok.decorators.func({tags: ['searchProvider']})
  static smartChainExecutionProvider(): DG.SearchProvider {
    return {
      'home': {
        name: 'LLM Smart chain execution',
        description: 'Plans and executes function chains based on user goals.',
        isApplicable: (s: string) => s.length > 5,
        returnType: 'widget',
        options: {
          widgetHeight: 500,
        },
        search: async (_s, _v) => {
          return null;
        },
        onValueEnter: async (s: string, v) => {
          const searchResultHost = document.getElementsByClassName('power-pack-search-host')[0];
          if (!searchResultHost)
            return;
          const gptEngine = new ChatGPTPromptEngine(apiKey, modelName);
          const gptAssistant = new ChatGptAssistant(gptEngine);

          const resultWait = ui.wait(async () => {
            const resDiv = ui.divV([], 'chatgpt-ask-ai-result');
            const plan = JSON.parse(await grok.functions.call('ChatGPT:getExecutionPlan', {userGoal: s})) as Plan;
            const planDiv = AssistantRenderer.renderPlan(plan);
            resDiv.appendChild(planDiv);
            const result = await gptAssistant.execute(plan);

            resDiv.appendChild(AssistantRenderer.renderResult(result));
            return resDiv;
          });
          searchResultHost.prepend(resultWait);
        },
        getSuggestions: (s: string) => s.length > 5 ? [
          {
            suggestionText: 'Smart Grok',
            priority: 11,
          }
        ] : null,
      }
    };
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
 static async getExecutionPlan(userGoal: string): Promise<string> {
   const gptEngine = new ChatGPTPromptEngine(apiKey, modelName);
   const gptAssistant = new ChatGptAssistant(gptEngine);
   const plan: Plan = await gptAssistant.plan(userGoal);
   // Cache only works with scalar values, so we serialize the plan to a string
   return JSON.stringify(plan);
 }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async fuzzyMatch(prompt: string, searchPatterns: string[]): Promise<string> {
    const queryMatchResult = await findBestFunction(prompt, searchPatterns);
    return JSON.stringify(queryMatchResult);
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async askDocumentationCached(prompt: string): Promise<string> {
    const client = OpenAIHelpClient.getInstance(apiKey, vectorStoreId);
    return await client.getHelpAnswer(prompt);
  }
}
