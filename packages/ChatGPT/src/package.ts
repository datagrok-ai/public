/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChatGptAssistant} from './prompt-engine/chatgpt-assistant';
import {ChatGPTPromptEngine} from './prompt-engine/prompt-engine';
import {getAiPanelVisibility, initAiPanel, setAiPanelVisibility} from './ai-panel';
import {findBestFunction, tableQueriesFunctionsSearchLlm} from './prompts/find-best-function';
import {askWiki, setupSearchUI, smartExecution} from './llm-utils/ui';
import {Plan} from './prompt-engine/interfaces';
import {OpenAIHelpClient} from './llm-utils/openAI-client';
import {LLMCredsManager} from './llm-utils/creds';
import {CombinedAISearchAssistant} from './llm-utils/combined-search';
import {JsonSchema} from './prompt-engine/interfaces';
import {generateAISqlQuery} from './llm-utils/sql-utils';

export * from './package.g';
export const _package = new DG.Package();

export const modelName: string = 'gpt-4o-mini';

export class PackageFunctions {
  @grok.decorators.init()
  static async init() {
    LLMCredsManager.init(_package);
    setupSearchUI();
  }


  @grok.decorators.autostart()
  static autostart() {
    try {
      grok.events.onViewAdded.subscribe((view) => {
        if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
          const tableView = view as DG.TableView;
          const iconFse = ui.iconSvg('ai.svg', () => setAiPanelVisibility(true), 'Ask AI \n Ctrl+I');
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
    } catch (e) {
      console.log('AI autostart failed.');
      console.log(e);
    }
  }

  @grok.decorators.func({tags: ['searchProvider']})
  static combinedLLMSearchProvider(): DG.SearchProvider {
    return {
      'home': {
        name: 'Ask AI Assistant',
        search: async (_query: string) => {
          CombinedAISearchAssistant.instance.resetSearchUI();
          return null;
        },
        getSuggestions: (_query: string) => [],
        isApplicable: (query: string) => LLMCredsManager.getApiKey() ? query?.trim().length >= 2 : false,
        description: 'Get answers form AI assistant',
        onValueEnter: async (query) => {
          LLMCredsManager.getApiKey() && query?.trim().length >= 2 &&
            await CombinedAISearchAssistant.instance.searchUI(query);
        }
      }
    } satisfies DG.SearchProvider;
  }

  @grok.decorators.func({meta: {
    role: 'aiSearchProvider',
    useWhen: 'If the user is asking questions about how to do something, how to write the code on platform, how to execute tasks, or any other questions related to Datagrok platform functionalities and capabilities. The tone of the prompt should generally sound like "how do I do this" / "what is this". for example, "what sequence notations are supported?'
  }, name: 'Help',
  description: 'Get answers from DeepGROK AI assistant based on Datagrok documentation and public code.', result: {type: 'widget', name: 'result'}})
  static async askHelpLLMProvider(@grok.decorators.param({type: 'string'})prompt: string): Promise<DG.Widget | null> {
    return await askWiki(prompt);
  }

 @grok.decorators.func({meta: {
   role: 'aiSearchProvider',
   useWhen: 'If the prompt looks like a user has a goal to achieve something with concrete input(s), and wants the system to plan and execute a series of steps/functions to achieve that goal. for example, adme properties of CHEMBL1234, enumerate some peptide, etc.. . Also, if the tone of the prompt sounds like "Do something", use this function'
 }, name: 'Execute',
 description: 'Plans and executes function steps to achieve needed results', result: {type: 'widget', name: 'result'}})
  static async smartChainExecutionProvider(@grok.decorators.param({type: 'string'})prompt: string): Promise<DG.Widget | null> {
    return await smartExecution(prompt, modelName);
  }

  @grok.decorators.func({meta: {
    role: 'aiSearchProvider',
    useWhen: 'if the prompt suggest that the user is looking for a data table result and the prompt resembles a query pattern. for example, "bioactivity data for shigella" or "compounds similar to aspirin" or first 100 chembl compounds. there should be some parts of user prompt that could match parameters in some query, like shigella, aspirin, first 100 etc.'
  }, name: 'Query',
  description: 'Tries to find a query which has the similar pattern as the prompt user entered and executes it', result: {type: 'widget', name: 'result'}})
 static async llmSearchQueryProvider(@grok.decorators.param({type: 'string'})prompt: string): Promise<DG.Widget | null> {
   return await tableQueriesFunctionsSearchLlm(prompt);
 }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async askAIGeneralCached(
    model: string,
    systemPrompt: string,
    prompt: string,
    schema?: JsonSchema
  ): Promise<string> {
    const client = OpenAIHelpClient.getInstance();
    // this is used only here to provide caching
    return await client.generalPrompt(model, systemPrompt, prompt, schema);
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async ask(
    question: string,
  ): Promise<string> {
    const client = OpenAIHelpClient.getInstance();
    // this is used only here to provide caching
    return await client.generalPrompt(modelName, 'You are a helpful assistant.', question);
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async getExecutionPlan(userGoal: string): Promise<string> {
    const gptEngine = ChatGPTPromptEngine.getInstance(LLMCredsManager.getApiKey(), modelName);
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
  static async fuzzyMatch(prompt: string, searchPatterns: string[], descriptions: string[]): Promise<string> {
    const queryMatchResult = await findBestFunction(prompt, searchPatterns, descriptions);
    return JSON.stringify(queryMatchResult);
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async askDocumentationCached(prompt: string): Promise<string> {
    const client = OpenAIHelpClient.getInstance();
    return await client.getHelpAnswer(prompt);
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async generateSqlQuery(prompt: string, connectionID: string, schemaName: string): Promise<string> {
    console.log(prompt);
    console.log(connectionID);
    console.log(schemaName);

    return await generateAISqlQuery(prompt, connectionID, schemaName);
  }
}
