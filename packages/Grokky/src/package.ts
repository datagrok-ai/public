/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {interval} from 'rxjs';
import {findBestMatchingQuery, tableQueriesFunctionsSearchLlm} from './ai/search/query-matching';
import {askWiki, smartExecution, setupAIQueryEditorUI, setupScriptsAIPanelUI, setupSearchUI, setupShellAIPanelUI, setupTableViewAIPanelUI} from './ai/ui';
import {CombinedAISearchAssistant} from './ai/search/combined-search';
import {UsageLimiter} from './ai/usage-limiter';
import {ClaudeRuntimeClient} from './claude/runtime-client';
import {genDBConnectionMeta, moveDBMetaToStickyMetaOhCoolItEvenRhymes} from './db/db-index-tools';
import {biologicsIndex} from './db/indexes/biologics-index';
import {chemblIndex} from './db/indexes/chembl-index';
export * from './package.g';

export class ChatGPTPackage extends DG.Package {
}

export const _package = new ChatGPTPackage();


export class PackageFunctions {
  @grok.decorators.init({tags: ['init']})
  static async init() {
    await UsageLimiter.getInstance().init();
    setupSearchUI();
    setupTableViewAIPanelUI();
    setupScriptsAIPanelUI();
    PackageFunctions.ensureAgentsFolder();
    PackageFunctions.subscribeToSyncEvents();
  }

  // Creates agents/ folder in My Files if it doesn't exist yet.
  static async ensureAgentsFolder(): Promise<void> {
    try {
      const conn = await grok.dapi.connections.filter('name = "My files"').first();
      if (!conn)
        return;
      const agentsPath = `${conn.nqName}/agents`;
      const exists = await grok.dapi.files.exists(agentsPath);
      if (!exists) {
        await grok.dapi.files.writeAsText(`${agentsPath}/README.md`,
          'Place your personal knowledge files here. Claude will use them as context.');
        console.log('Grokky: created agents/ folder');
      }
    } catch (e: any) {
      console.warn('Grokky: failed to ensure agents folder:', e.message);
    }
  }

  static isAgentsFile(fi: DG.FileInfo): boolean {
    return (fi.fullPath ?? fi.path ?? fi.name ?? '').includes('agents');
  }

  // Subscribes to platform events that should trigger file sync.
  static subscribeToSyncEvents(): void {
    const sync = (...args: Parameters<ClaudeRuntimeClient['syncUserFiles']>) =>
      ClaudeRuntimeClient.getInstance().syncUserFiles(...args);

    // MyFiles agents: file operations (create, upload, delete, rename, move)
    grok.events.onEvent('d4-file-event').subscribe((eventData: any) => {
      const dartFiles = eventData?.dart?.files;
      if (!dartFiles)
        return;
      const files: DG.FileInfo[] = Array.from({length: dartFiles.length}, (_: any, i: number) => DG.toJs(dartFiles[i]));
      if (files.some(PackageFunctions.isAgentsFile))
        sync('user-files');
    });

    // MyFiles agents: in-place file edits (save)
    grok.events.onFileEdited.subscribe((fi: DG.FileInfo) => {
      if (PackageFunctions.isAgentsFile(fi))
        sync('user-files');
    });

    // Packages: when a JS bundle is loaded
    grok.events.onPackageLoaded.subscribe((pkg: DG.Package) => {
      sync('packages', pkg.name);
    });

    // Poll for shared connections and package updates every 10 minutes.
    // No reliable push events exist for sharing or other users' publishes.
    // TODO: think about more efficient strategies here.
    interval(15 * 60 * 1000).subscribe(() => {
      sync('shared');
      sync('packages');
    });
  }

  @grok.decorators.autostart({tags: ['autostart']})
  static autostart() {
    if (grok.shell.windows.showAI)
      setupShellAIPanelUI();
    grok.shell.windows.onPanelVisibilityChanged.subscribe(() => {
      if (grok.shell.windows.showAI)
        setupShellAIPanelUI();
    });
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
  description: 'Get answers from AI assistant based on Datagrok documentation and public code.', result: {type: 'widget', name: 'result'}})
  static async askHelpLLMProvider(@grok.decorators.param({type: 'string'})prompt: string): Promise<DG.Widget | null> {
    return await askWiki(prompt);
  }

  @grok.decorators.func({meta: {
    role: 'aiSearchProvider',
    useWhen: 'If the prompt looks like a user has a goal to achieve something with concrete input(s), and wants the system to plan and execute a series of steps/functions to achieve that goal. This relates to functions that analyse or mutate data, not get it. for example, adme properties of CHEMBL1234, enumerate some peptide, etc... Also, if the tone of the prompt sounds like "Do something to something", use this function'
  }, name: 'Execute',
  description: 'Plans and executes function steps to achieve needed results', result: {type: 'widget', name: 'result'}})
  static async smartChainExecutionProvider(@grok.decorators.param({type: 'string'})prompt: string): Promise<DG.Widget | null> {
    return await smartExecution(prompt);
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
  static async findMatchingPatternQuery(prompt: string): Promise<string> {
    const queryMatchResult = await findBestMatchingQuery(prompt);
    return JSON.stringify(queryMatchResult);
  }

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
