/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ChatGptAssistant} from './prompt-engine/chatgpt-assistant';
import {ChatGPTPromptEngine} from './prompt-engine/prompt-engine';
import {findBestMatchingQuery, tableQueriesFunctionsSearchLlm} from './llm-utils/query-matching';
import {askWiki, setupAIQueryEditorUI, setupScriptsAIPanelUI, setupSearchUI, setupTableViewAIPanelUI, smartExecution} from './llm-utils/ui';
import {Plan} from './prompt-engine/interfaces';
import {initModelTypeNames, ModelType, OpenAIClient} from './llm-utils/openAI-client';
import {LLMCredsManager} from './llm-utils/creds';
import {CombinedAISearchAssistant} from './llm-utils/combined-search';
import {JsonSchema} from './prompt-engine/interfaces';
import {genDBConnectionMeta, moveDBMetaToStickyMetaOhCoolItEvenRhymes} from './llm-utils/db-index-tools';
import {biologicsIndex} from './llm-utils/indexes/biologics-index';
import {chemblIndex} from './llm-utils/indexes/chembl-index';
import {uploadFilesToVectorStroreOneByOne} from './llm-utils/indexes/dg-repository-index';

export * from './package.g';
export const _package = new DG.Package();


export class PackageFunctions {
  @grok.decorators.init()
  static async init() {
    initModelTypeNames();
    LLMCredsManager.init(_package);
    setupSearchUI();
    setupTableViewAIPanelUI();
    setupScriptsAIPanelUI();

    // @ts-ignore
    window.openAI = OpenAIClient.getInstance().openai;

    // AIProvider.setProvider(_package.settings.APIName === 'openai chat completions' ?
    //   new OpenAIChatCompletionsProvider(OpenAIClient.getInstance().openai) : new OpenAIResponsesProvider(OpenAIClient.getInstance().openai)
    // );
  }


  @grok.decorators.autostart()
  static autostart() {
    // try {
    //   grok.events.onViewAdded.subscribe((view) => {
    //     if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
    //       const tableView = view as DG.TableView;
    //       const iconFse = ui.iconSvg('ai.svg', () => setAiPanelVisibility(true), 'Ask AI \n Ctrl+I');
    //       tableView.setRibbonPanels([...tableView.getRibbonPanels(), [iconFse]]);
    //     }
    //   });

    //   initAiPanel();
    //   // Add keyboard shortcut for toggling AI panel
    //   document.addEventListener('keydown', (event) => {
    //     // Check for Ctrl+I (Ctrl key + I key)
    //     if (event.ctrlKey && event.key === 'i') {
    //       event.preventDefault(); // Prevent default browser behavior
    //       const isVisible = getAiPanelVisibility();
    //       setAiPanelVisibility(!isVisible);
    //     }
    //   });
    // } catch (e) {
    //   console.log('AI autostart failed.');
    //   console.log(e);
    // }
  }

  @grok.decorators.func({
    meta: {role: 'searchProvider'},
  })
  static combinedLLMSearchProvider(): DG.SearchProvider {
    const isAiConfigured = grok.ai.config.configured;
    return {
      'home': {
        name: 'Ask AI Assistant',
        search: async (_query: string) => {
          CombinedAISearchAssistant.instance.resetSearchUI();
          return null;
        },
        getSuggestions: (_query: string) => [],
        isApplicable: (query: string) => isAiConfigured ? query?.trim().length >= 2 : false,
        description: 'Get answers form AI assistant',
        onValueEnter: async (query) => {
          isAiConfigured && query?.trim().length >= 2 &&
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
   useWhen: 'If the prompt looks like a user has a goal to achieve something with concrete input(s), and wants the system to plan and execute a series of steps/functions to achieve that goal. This relates to functions that analyse or mutate data, not get it. for example, adme properties of CHEMBL1234, enumerate some peptide, etc... Also, if the tone of the prompt sounds like "Do something to something", use this function'
 }, name: 'Execute',
 description: 'Plans and executes function steps to achieve needed results', result: {type: 'widget', name: 'result'}})
  static async smartChainExecutionProvider(@grok.decorators.param({type: 'string'})prompt: string): Promise<DG.Widget | null> {
    return await smartExecution(prompt, ModelType.Fast);
  }

  @grok.decorators.func({meta: {
    role: 'aiSearchProvider',
    useWhen: 'if the prompt suggest that the user is looking for a data table result and the prompt resembles a query pattern. for example, "bioactivity data for shigella" or "compounds similar to aspirin" or first 100 chembl compounds. there should be some parts of user prompt that could match parameters in some query, like shigella, aspirin, first 100 etc. Always use this function when user wants to get the data without any further processing or calculating'
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
    const client = OpenAIClient.getInstance();
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
    const client = OpenAIClient.getInstance();
    // this is used only here to provide caching
    return await client.generalPrompt(ModelType.Fast, 'You are a helpful assistant.', question);
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async getExecutionPlan(userGoal: string): Promise<string> {
    const gptEngine = ChatGPTPromptEngine.getInstance(ModelType['Deep Research']);
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
  static async findMatchingPatternQuery(prompt: string): Promise<string> {
    const queryMatchResult = await findBestMatchingQuery(prompt);
    return JSON.stringify(queryMatchResult);
  }

  @grok.decorators.func({
    'meta': {
      'cache': 'all',
      'cache.invalidateOn': '0 0 1 * *'
    }
  })
  static async askDocumentationCached(prompt: string): Promise<string> {
    const client = OpenAIClient.getInstance();
    return await client.getHelpAnswer(prompt);
  }

  // @grok.decorators.func({
  //   'meta': {
  //     'cache': 'all',
  //     'cache.invalidateOn': '0 0 1 * *'
  //   }
  // })
  // static async generateSqlQuery(prompt: string, connectionID: string, schemaName: string): Promise<string> {
  //   return await generateAISqlQuery(prompt, connectionID, schemaName);
  // }
  @grok.decorators.func({})
  static async setupAIQueryEditor(view: DG.ViewBase, connectionID: string, queryEditorRoot: HTMLElement, @grok.decorators.param({type: 'dynamic'}) setAndRunFunc: Function): Promise<boolean> {
    return setupAIQueryEditorUI(view, connectionID, queryEditorRoot, setAndRunFunc as (query: string) => void);
  }

  @grok.decorators.func({})
  static async moveMetaToDB(@grok.decorators.param({type: 'string', options: {choices: ['biologics', 'chembl']}})dbName: string): Promise<void> {
    const meta = dbName === 'biologics' ? biologicsIndex : chemblIndex;
    await moveDBMetaToStickyMetaOhCoolItEvenRhymes(meta);
  }

  @grok.decorators.func({})
  static async setupVectorStore() {
    const keyWordInput = ui.input.string('PassKey', {tooltipText: 'Enter the passkey for this action', placeholder: 'Type passkey here'});
    ui.dialog('Index DG Repository Files to Vector Store')
      .add(ui.divText('This will upload all files from DG Repository to the vector store for AI search and assistance. This may take up to an hour depending on the number of files in the repository.'))
      .add(ui.divText('This action will delete the previous files and update all files. make sure to LEAVE THE PLATFORM OPEN until the process is finished.'))
      .add(ui.divText('Make sure that you have set up the OpenAI API key and vector store ID in the package settings before proceeding.'))
      .add(keyWordInput)
      .onOK(async () => {
        if (keyWordInput.value !== 'DG-index-admin') {
          grok.shell.error('Incorrect passkey. Operation aborted.');
          return;
        }
        await
        uploadFilesToVectorStroreOneByOne();
      })
      .show();
  }

  @grok.decorators.func({})
  static async searchForSomething() {
    const query = 'anything';
    const openAIClient = OpenAIClient.getInstance().openai;
    const f = await openAIClient.vectorStores.search('vs_696fbf9985f8819191dc6f71e04b6593', {query: query, max_num_results: 5, filters: {
      key: 'firstParentFolder',
      type: 'eq',
      value: 'datagrok-celery-task'
    }});
    console.log(f.data);
  }


  @grok.decorators.func({})
  static async indexDatabaseSchema() {
    grok.shell.info('Do not touch this function if you are not a developer. :D :D :D');
    return;
    const connections = await grok.dapi.connections.list();
    const connectionsInput = ui.input.choice('Connection', {items: connections.map((c) => c.nqName), value: connections[0].nqName, nullable: false});
    const getConnectionByName = (): DG.DataConnection => {
      return connections.find((c) => c.nqName === connectionsInput.value)!;
    };
    const schemaInputRoot = ui.div([]);
    let schemaInput: DG.ChoiceInput<string>;
    const createSchemaInput = async () => {
      ui.empty(schemaInputRoot);
      const connection = getConnectionByName();
      const schemas = await grok.dapi.connections.getSchemas(connection);
      schemaInput = ui.input.choice('Schema', {items: schemas, value: schemas[0], nullable: false}) as DG.ChoiceInput<string>;
      schemaInputRoot.append(schemaInput.root);
    };
    connectionsInput.onChanged.subscribe(() => createSchemaInput());
    await createSchemaInput();
    ui.dialog('Index Database Schema')
      .add(connectionsInput)
      .add(schemaInputRoot)
      .onOK(async () => {
        const res = await genDBConnectionMeta(getConnectionByName(), [schemaInput.value]);
        grok.shell.info('Database schema indexed successfully.');
        console.log(res);
        DG.Utils.download(`${connectionsInput.value}_${schemaInput.value}_db_index.json`, JSON.stringify(res, (_, v) => typeof v === 'bigint' ? Number(v) : v, 2));
      }).show();
  }
}
