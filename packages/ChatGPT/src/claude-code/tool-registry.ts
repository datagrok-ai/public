import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {TableViewContext} from '../llm-utils/tableview-tools';
import {SQLGenerationContext} from '../llm-utils/sql-tools';
import {ScriptGenerationContext} from '../llm-utils/script-tools';
import {BuiltinDBInfoMeta} from '../llm-utils/query-meta-utils';
import {LLMCredsManager} from '../llm-utils/creds';

export class BrowserToolExecutor {
  private tableViewContext: TableViewContext | null;
  private sqlContext: SQLGenerationContext | null = null;
  private scriptContext: ScriptGenerationContext | null = null;

  constructor(
    private tableView: DG.TableView | null,
    private connectionId: string | null,
  ) {
    this.tableViewContext = tableView ? new TableViewContext(tableView) : null;
  }

  private async getSqlContext(): Promise<SQLGenerationContext> {
    if (!this.sqlContext) {
      if (!this.connectionId)
        throw new Error('No database connection available for SQL tools');
      const connection = await grok.dapi.connections.find(this.connectionId);
      if (!connection)
        throw new Error(`Connection with ID ${this.connectionId} not found`);
      const allDbInfos = await BuiltinDBInfoMeta.allFromConnection(connection);
      const defaultCatalog = connection.parameters?.['catalog'] ?? connection.parameters?.['db'] ?? allDbInfos[0]?.name ?? '';
      this.sqlContext = new SQLGenerationContext(this.connectionId, defaultCatalog, allDbInfos, connection);
    }
    return this.sqlContext;
  }

  private getScriptContext(): ScriptGenerationContext {
    if (!this.scriptContext) {
      const vectorStoreId = LLMCredsManager.getVectorStoreId();
      this.scriptContext = new ScriptGenerationContext(vectorStoreId);
    }
    return this.scriptContext;
  }

  async execute(msg: {callId: string, tool: string, args: any}): Promise<string> {
    const {tool, args} = msg;
    switch (tool) {
    // Table View tools -> TableViewContext
    case 'add_viewer':
      if (!this.tableViewContext) return 'Error: no active table view';
      return this.tableViewContext.addViewer(args.viewerType, args.viewerProperties);
    case 'describe_viewer':
      if (!this.tableViewContext) return 'Error: no active table view';
      return this.tableViewContext.describeViewer(args.viewerType);
    case 'list_current_viewers':
      if (!this.tableViewContext) return 'Error: no active table view';
      return this.tableViewContext.listCurrentViewers();
    case 'adjust_viewer':
      if (!this.tableViewContext) return 'Error: no active table view';
      return this.tableViewContext.adjustViewer(args.viewerId, args.viewerProperties);
    case 'find_viewers_by_type':
      if (!this.tableViewContext) return 'Error: no active table view';
      return this.tableViewContext.findViewersByType(args.viewerType);
    case 'get_available_viewers':
      if (!this.tableViewContext) return 'Error: no active table view';
      return this.tableViewContext.getAvailableViewers();
    case 'highlight_element':
      if (!this.tableViewContext) return 'Error: no active table view';
      return (await this.tableViewContext.highlightElement(args.elementId, args.description)) ?? 'Element highlighted';
    case 'list_ui_elements':
      if (!this.tableViewContext) return 'Error: no active table view';
      return this.tableViewContext.listUIElements();

    // SQL tools -> SQLGenerationContext
    case 'list_tables_in_schema': {
      const ctx = await this.getSqlContext();
      return (await ctx.listTables(false, args.schemaName, args.catalogName)).description;
    }
    case 'describe_tables': {
      const ctx = await this.getSqlContext();
      return ctx.describeTables(args.tables);
    }
    case 'list_joins': {
      const ctx = await this.getSqlContext();
      return ctx.listJoins(args.tables);
    }
    case 'try_sql': {
      const ctx = await this.getSqlContext();
      return ctx.trySql(args.sql, args.description);
    }
    case 'list_all_schemas': {
      const ctx = await this.getSqlContext();
      return (await ctx.listAllSchemas(args.catalogName)).join(', ');
    }
    case 'list_catalogs': {
      const ctx = await this.getSqlContext();
      return (await ctx.listCatalogs()).join(', ');
    }
    case 'find_similar_queries': {
      if (!this.connectionId)
        return 'No database connection available for query search';
      const prompt = args.prompt ?? '';
      const entities = await grok.ai.searchEntities(prompt, 0.3, 10, ['DataQuery']);
      const queries = entities
        .filter((e): e is DG.DataQuery => !!e && e instanceof DG.DataQuery && e.connection?.id === this.connectionId)
        .filter((_, i) => i < 4);
      if (queries.length === 0)
        return 'No similar queries found for this connection';
      let result = `Found ${queries.length} similar queries:\n\n`;
      for (const q of queries) {
        result += `**${q.name}**\n`;
        if (q.description)
          result += `Description: ${q.description}\n`;
        result += `\`\`\`sql\n${q.query}\n\`\`\`\n\n`;
      }
      return result;
    }

    // Script tools -> ScriptGenerationContext
    case 'find_similar_script_samples':
      return this.getScriptContext().findSimilarScriptSamples(args.description, args.language);
    case 'search_documentation':
      return this.getScriptContext().searchDocumentation(args.query, args.maxResults);
    case 'list_js_members':
      return this.getScriptContext().listJsMembers(args.expression);
    case 'search_js_api_sources':
      return this.getScriptContext().searchJsApiSources(args.query, args.maxResults);

    default:
      return `Error: unknown tool "${tool}"`;
    }
  }
}
