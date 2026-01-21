/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBExplorerConfig, ExplicitReference, QueryJoinOptions} from './types';
import * as rxjs from 'rxjs';
import {renderExampleCard} from './renderer';

type EntryPoint = DBExplorerConfig['entryPoints'][string];

type JoinOption = QueryJoinOptions

type SchemaCacheEntry = {
  schemaInfo: any[];
  tableNames: string[];
  columnNamesByTable: {[tableName: string]: string[]};
};

type CustomRenderer = NonNullable<DBExplorerConfig['customRenderers']>[number];

interface DBExplorerState {
  connectionName: string | null;
  connectionNqName: string | null;
  connectionProvider: string | null;
  schemaName: string | null;
  entryPoints: {[semanticType: string]: EntryPoint};
  joinOptions: JoinOption[];
  headerNames: {[tableName: string]: string};
  uniqueColumns: {[tableName: string]: string};
  customSelectedColumns: {[tableName: string]: string[]};
  customRenderers?: CustomRenderer[];
  explicitReferences?: ExplicitReference[];
}

export class DBExplorerEditor {
  private state: DBExplorerState = {
    connectionName: null,
    connectionNqName: null,
    connectionProvider: null,
    schemaName: null,
    entryPoints: {},
    joinOptions: [],
    headerNames: {},
    uniqueColumns: {},
    customSelectedColumns: {},
    customRenderers: []
  };

  private connections: DG.DataConnection[] = [];
  private schemaInfo: any[] = [];
  private tableNames: string[] = [];
  private columnNamesByTable: {[tableName: string]: string[]} = {};
  private schemaInfoCache: {[schemaName: string]: SchemaCacheEntry} = {};
  private validationErrors: string[] = [];
  private onSave?: (config: any) => void;

  private rootDiv!: HTMLDivElement;
  private connectionInput!: DG.ChoiceInput<string | null>;
  private schemaInput!: DG.ChoiceInput<string | null>;
  private accordionDiv!: HTMLDivElement;

  constructor(onSave?: (config: any) => void, private showConnectionSection: boolean = true) {
    this.onSave = onSave;
  }

  private createActionButtons(onEdit: () => void, onDelete: () => void): HTMLDivElement {
    const btnDiv = ui.divH([]);
    btnDiv.style.gap = '8px';
    const editBtn = ui.icons.edit(onEdit, 'Edit');
    const deleteBtn = this.getDeleteIcon(onDelete, 'Delete');
    btnDiv.appendChild(editBtn);
    btnDiv.appendChild(deleteBtn);
    return btnDiv;
  }

  private getDeleteIcon(handler: Function, tooltip?: string) {
    const icon = ui.icons.delete(() => handler(), tooltip);
    icon.style.color = 'var(--red-3)';
    return icon;
  }

  private get dataConnection(): DG.DataConnection | null {
    if (!this.state.connectionName || !this.state.connectionNqName)
      return null;
    return this.connections.find((c) => c.nqName?.toLowerCase() === this.state.connectionNqName?.toLowerCase()) || null;
  }

  private getAvailableColumnsForTable(tableName: string): string[] {
    // Start with the base table columns
    const baseColumns = this.columnNamesByTable[tableName] || [];

    // Find all join options that start from this table
    const joinColumns: string[] = [];
    for (const joinOption of this.state.joinOptions) {
      const joinSchema = joinOption.fromSchema ?? this.state.schemaName;
      if (joinOption && joinOption.fromTable === tableName && joinSchema === this.state.schemaName) {
        // Add the aliased column names from join selects
        if (joinOption.select) {
          for (const selectStr of joinOption.select) {
            if (selectStr) {
              // Parse "column as alias" or just "column"
              const parts = selectStr.trim().split(/\s+as\s+/i);
              const alias = parts.length > 1 ? parts[1].trim() : parts[0].trim();
              if (alias)
                joinColumns.push(alias);
            }
          }
        }
      }
    }

    return [...baseColumns, ...joinColumns];
  }

  public async getUI(): Promise<HTMLDivElement> {
    this.rootDiv = ui.divV([], 'db-explorer-editor');

    // Create basic connection info section
    await this.createConnectionSection();

    // Create accordion for advanced options
    this.accordionDiv = ui.div([]);
    this.rootDiv.appendChild(this.accordionDiv);

    // Create action buttons
    this.createButtonBar();

    return this.rootDiv;
  }

  private async createConnectionSection(): Promise<void> {
    const connectionForm = ui.divV([], 'db-explorer-connection-form ui-form');
    connectionForm.style.maxHeight = 'fit-content';

    // Get all connections
    this.connections = await grok.dapi.connections.list();
    const connectionBqNames = this.connections.map((c) => c.nqName).sort();

    this.connectionInput = ui.input.choice('Connection', {
      value: this.state.connectionNqName,
      nullable: true,
      items: connectionBqNames,
      tooltipText: 'Select a connection where your identifiers live'
    });
    ui.tooltip.bind(this.connectionInput.input, 'Database connection to use for exploration');

    this.schemaInput = ui.input.choice('Schema', {
      value: this.state.schemaName,
      nullable: true,
      items: [] as string[],
      tooltipText: 'Select a schema'
    });
    ui.tooltip.bind(this.schemaInput.input, 'Database schema name where the identifiers are located');

    // Handle connection change
    this.connectionInput.onChanged.subscribe(async () => {
      const selectedConnection = this.connectionInput.value ?
        this.connections.find((c) => c.nqName === this.connectionInput.value) : null;

      if (!selectedConnection) {
        this.schemaInput.items = [];
        this.schemaInput.value = null;
        this.state.connectionName = null;
        this.state.connectionNqName = null;
        this.state.connectionProvider = null;
        this.schemaInfoCache = {};
        this.clearAdvancedSections();
        return;
      }

      this.state.connectionNqName = this.connectionInput.value;
      this.state.connectionName = selectedConnection.name || null;
      this.state.connectionProvider = selectedConnection.dataSource || null;
      this.schemaInfoCache = {};
      const schemas = await grok.dapi.connections.getSchemas(selectedConnection);
      const oldSchemaValue = this.schemaInput.value;
      this.schemaInput.items = schemas.sort();
      if (oldSchemaValue && schemas.includes(oldSchemaValue))
        this.schemaInput.value = oldSchemaValue;
      else
        this.schemaInput.value = null;
    });

    // Handle schema change
    this.schemaInput.onChanged.subscribe(async () => {
      this.state.schemaName = this.schemaInput.value;

      if (!this.state.connectionName || !this.state.schemaName) {
        this.clearAdvancedSections();
        return;
      }

      await this.loadSchemaInfo();
      await this.buildAdvancedSections();
    });

    // Handle data source change

    connectionForm.appendChild(this.connectionInput.root);
    connectionForm.appendChild(this.schemaInput.root);

    this.rootDiv.appendChild(connectionForm);
    if (!this.showConnectionSection)
      connectionForm.style.display = 'none';
  }

  private async loadSchemaInfo(): Promise<void> {
    const entry = await this.ensureSchemaLoaded(this.state.schemaName);
    if (!entry)
      return;
    this.schemaInfo = entry.schemaInfo;
    this.tableNames = entry.tableNames;
    this.columnNamesByTable = entry.columnNamesByTable;
  }

  private async ensureSchemaLoaded(schemaName: string | null): Promise<SchemaCacheEntry | null> {
    if (!schemaName)
      return null;
    if (this.schemaInfoCache[schemaName])
      return this.schemaInfoCache[schemaName];
    const connection = this.dataConnection;
    if (!connection)
      return null;
    const schemaInfo = await grok.dapi.connections.getSchema(connection, schemaName);
    const tableNames = schemaInfo.map((t) => t.friendlyName).sort();
    const columnNamesByTable: {[tableName: string]: string[]} = {};
    for (const table of schemaInfo)
      columnNamesByTable[table.friendlyName] = table.columns.map((c: any) => c.name).sort();
    const entry: SchemaCacheEntry = {schemaInfo, tableNames, columnNamesByTable};
    this.schemaInfoCache[schemaName] = entry;
    return entry;
  }

  private getSchemaTableNames(schemaName: string | null): string[] {
    if (!schemaName)
      return [];
    return this.schemaInfoCache[schemaName]?.tableNames ?? [];
  }

  private getSchemaColumnNames(schemaName: string | null, tableName: string | null): string[] {
    if (!schemaName || !tableName)
      return [];
    return this.schemaInfoCache[schemaName]?.columnNamesByTable?.[tableName] ?? [];
  }

  private clearAdvancedSections(): void {
    this.accordionDiv.innerHTML = '';
    this.schemaInfo = [];
    this.tableNames = [];
    this.columnNamesByTable = {};
  }

  private async buildAdvancedSections(): Promise<void> {
    this.accordionDiv.innerHTML = '';

    const accordion = ui.accordion('db-explorer-advanced-accordion');
    const tryGetPaneHeaderRoot = (paneRoot: HTMLElement) => paneRoot.getElementsByClassName('d4-accordion-pane-header')?.item(0) as HTMLElement ?? paneRoot;
    // Entry Points Section
    const entryPane = accordion.addPane('Identifiers', () => this.createEntryPointsSection(), true, null, false);
    ui.tooltip.bind(tryGetPaneHeaderRoot(entryPane.root), 'Define identifiers that trigger exploration');

    // Join Options Section
    const joinPane = accordion.addPane('Joins', () => this.createJoinOptionsSection(), true, null, false);
    ui.tooltip.bind(tryGetPaneHeaderRoot(joinPane.root), 'Define how to include data from related tables');

    // Explicit References Section
    const explicitRefsPane = accordion.addPane('Explicit References', () => this.createExplicitReferencesSection(), false, null, false);
    ui.tooltip.bind(tryGetPaneHeaderRoot(explicitRefsPane.root), 'Define references between schema-table-columns that are not in schema metadata');

    // Header Names Section
    const headerPane = accordion.addPane('Header Names', () => this.createHeaderNamesSection(), false, null, false);
    ui.tooltip.bind(tryGetPaneHeaderRoot(headerPane.root), 'Customize column headers in the exploration view. Define a column for tables which will be used as a header for multiple records.');

    // Unique Columns Section
    const uniquePane = accordion.addPane('Unique Columns', () => this.createUniqueColumnsSection(), false, null, false);
    ui.tooltip.bind(tryGetPaneHeaderRoot(uniquePane.root), 'Specify which column uniquely identifies rows in tables');

    // Custom Selected Columns Section
    const customPane = accordion.addPane('Custom Selected Columns', () => this.createCustomSelectedColumnsSection(), false, null, false);
    ui.tooltip.bind(tryGetPaneHeaderRoot(customPane.root), 'Override which columns to display for specific tables');

    // Custom Renderers Section
    const renderersPane = accordion.addPane('Renderers', () => this.createCustomRenderersSection(), false, null, false);
    ui.tooltip.bind(tryGetPaneHeaderRoot(renderersPane.root), 'Define custom renderers for specific table columns');


    this.accordionDiv.appendChild(accordion.root);
  }

  private createEntryPointsSection(): HTMLDivElement {
    const section = ui.divV([], 'entry-points-section');

    const infoText = ui.divText('Define identifiers that trigger exploration', 'db-explorer-info-text');
    ui.tooltip.bind(infoText, 'Identifiers map to semantic types which map to table columns where identifiers are located');

    const tableContainer = ui.divV([], 'entry-points-table-container');
    this.renderEntryPointsTable(tableContainer);

    const addButtonContainer = ui.divH([]);
    addButtonContainer.style.marginTop = '10px';
    const addButton = this.getAddIcon(() => this.showEntryPointDialog(null, tableContainer), 'Add a new identifier');
    addButtonContainer.appendChild(addButton);

    section.appendChild(infoText);
    section.appendChild(tableContainer);
    section.appendChild(addButtonContainer);

    return section;
  }

  private getAddIcon(handler: Function, tooltip?: string) {
    const icon = ui.icons.add(() => handler(), tooltip);
    icon.style.color = 'var(--blue-1)';
    return icon;
  }

  private renderEntryPointsTable(container: HTMLDivElement): void {
    container.innerHTML = '';

    if (Object.keys(this.state.entryPoints).length === 0) {
      container.appendChild(ui.divText('No Identifiers defined yet', 'db-explorer-empty-message'));
      return;
    }

    const table = ui.table(
      Object.keys(this.state.entryPoints),
      (semanticType) => {
        const ep = this.state.entryPoints[semanticType];
        const regexpInfo = ep.regexpExample ?
          `${ep.regexpExample.example} (${ep.regexpExample.nonVariablePart})` :
          '-';
        return [
          semanticType,
          ep.table,
          ep.column,
          ep.matchRegexp || '-',
          regexpInfo,
          this.createActionButtons(
            () => this.showEntryPointDialog(semanticType, container),
            () => {
              delete this.state.entryPoints[semanticType];
              this.renderEntryPointsTable(container);
            }
          )
        ];
      },
      ['Semantic Type', 'Table', 'Column', 'Match Regexp', 'Regexp Example', 'Actions']
    );
    table.style.width = '100%';
    container.appendChild(table);
  }

  private showEntryPointDialog(semanticType: string | null, tableContainer: HTMLDivElement): void {
    const isEdit = semanticType !== null;
    const entryPoint = isEdit ? this.state.entryPoints[semanticType!] : {
      table: this.tableNames[0] || '',
      column: '',
      matchRegexp: undefined,
      regexpExample: undefined
    };

    const dialog = ui.dialog(isEdit ? `Edit Identifier: ${semanticType}` : 'Add Identifier');

    const semanticTypeInput = ui.input.string('Semantic Type', {value: semanticType || ''});
    ui.tooltip.bind(semanticTypeInput.input, 'The semantic type name for this identifier (e.g., CHEMBL_ID)');

    const tableInput = ui.input.choice('Table', {value: entryPoint.table, items: this.tableNames});
    ui.tooltip.bind(tableInput.input, 'Table containing this identifier');

    const columnInput = ui.input.choice('Column', {
      value: entryPoint.column,
      items: this.columnNamesByTable[entryPoint.table] || []
    });
    ui.tooltip.bind(columnInput.input, 'Column containing this identifier');

    tableInput.onChanged.subscribe(() => {
      if (tableInput.value)
        columnInput.items = this.columnNamesByTable[tableInput.value] || [];
    });

    const matchRegexpInput = ui.input.string('Match Regexp', {
      value: entryPoint.matchRegexp || '',
      nullable: true
    });
    ui.tooltip.bind(matchRegexpInput.input, 'JavaScript RegExp pattern to match identifiers (optional). Leave empty if the given semantic type already exists in Datagrok. E.g., CHEMBL\\d+');

    let runTestFuncRef: (() => void) | null = null;
    matchRegexpInput.onChanged.subscribe(() => {
      if (runTestFuncRef)
        runTestFuncRef();
    });
    // Add test option for regexp
    const testLink = ui.link('test', () => {
      if (!matchRegexpInput.value?.trim()) {
        grok.shell.error('Please enter a regular expression first');
        return;
      }

      // Create test dialog
      const testDialog = ui.dialog('Test Regular Expression');

      const testInput = ui.input.string('Test Identifier', {value: '', placeholder: 'Enter identifier to test'});
      ui.tooltip.bind(testInput.input, 'Enter an identifier to test against the regular expression');

      const resultDiv = ui.divV([], 'regexp-test-result');
      resultDiv.style.marginTop = '10px';
      resultDiv.style.padding = '10px';
      resultDiv.style.borderRadius = '4px';

      const runTest = () => {
        const regexpStr = matchRegexpInput.value?.trim();

        if (!regexpStr) {
          grok.shell.error('Please enter a regular expression first');
          return;
        }
        resultDiv.innerHTML = '';
        const testValue = testInput.value?.trim();

        if (!testValue) {
          resultDiv.style.backgroundColor = 'var(--grey-2)';
          resultDiv.appendChild(ui.divText('Enter an identifier to test', 'regexp-test-message'));
          return;
        }

        try {
          const regex = new RegExp(regexpStr);
          const matchResult = testValue.match(regex);

          if (matchResult) {
            resultDiv.style.backgroundColor = 'var(--green-1)';
            resultDiv.appendChild(ui.divText('✓ Match successful', 'regexp-test-success'));
            const matchedText = ui.divText(`Matched: "${matchResult[0]}"`, 'regexp-test-matched');
            matchedText.style.marginTop = '5px';
            matchedText.style.fontWeight = 'bold';
            resultDiv.appendChild(matchedText);
          } else {
            resultDiv.style.backgroundColor = 'var(--red-1)';
            resultDiv.appendChild(ui.divText('✗ No match', 'regexp-test-failure'));
          }
        } catch (e: any) {
          resultDiv.style.backgroundColor = 'var(--orange-1)';
          resultDiv.appendChild(ui.divText(`Error: ${e.message}`, 'regexp-test-error'));
        }
      };
      runTestFuncRef = runTest;

      testInput.onChanged.subscribe(runTest);

      testDialog.add(ui.divV([
        testInput.root,
        resultDiv
      ]));

      testDialog.show();
      dialog.onClose.subscribe(() => {
        runTestFuncRef = null;
      });
    }, 'Test the regular expression with sample identifiers');
    testLink.style.marginLeft = '5px';
    testLink.style.fontSize = '12px';

    matchRegexpInput.addOptions(testLink);

    const hasRegexpExample = !!entryPoint.regexpExample;
    const regexpExampleCheckbox = ui.input.bool('Include Regexp Example', {value: hasRegexpExample});
    ui.tooltip.bind(regexpExampleCheckbox.input, 'Add example for documentation and autosuggestions');

    const regexpExampleDiv = ui.divV([], 'ui-form');

    const exampleInput = ui.input.string('Example', {value: entryPoint.regexpExample?.example || ''});
    ui.tooltip.bind(exampleInput.input, 'Example identifier (e.g., CHEMBL1234)');

    const nonVarInput = ui.input.string('Non-Variable Part', {value: entryPoint.regexpExample?.nonVariablePart || ''});
    ui.tooltip.bind(nonVarInput.input, 'Fixed prefix (e.g., CHEMBL). Used for autosuggestions');

    const markupInput = ui.input.string('Regexp Markup', {value: entryPoint.regexpExample?.regexpMarkup || ''});
    ui.tooltip.bind(markupInput.input, 'Pattern for display suggestion (e.g., CHEMBL[0-9]+)');

    const updateRegexpExampleSection = () => {
      regexpExampleDiv.innerHTML = '';
      if (regexpExampleCheckbox.value) {
        regexpExampleDiv.appendChild(exampleInput.root);
        regexpExampleDiv.appendChild(nonVarInput.root);
        regexpExampleDiv.appendChild(markupInput.root);
      }
    };

    regexpExampleCheckbox.onChanged.subscribe(updateRegexpExampleSection);
    updateRegexpExampleSection();

    // add a preview card of how this thing will look like in the UI
    const previewCardDiv = ui.divV([], {classes: 'db-explorer-entrypoint-preview-card', style: {minWidth: '350px', flexGrow: '1'}});

    const updatePreviewCard = () => {
      ui.empty(previewCardDiv);
      if (!this.dataConnection || !this.state.schemaName || !tableInput.value) {
        previewCardDiv.appendChild(ui.divText('Please select the table to see the preview', 'db-explorer-preview-message'));
        return;
      }
      const card = renderExampleCard(
        this.dataConnection, this.state.schemaName, tableInput.value, {joinOptions: this.state.joinOptions, customRenderers: this.state.customRenderers,
          uniqueColumns: this.state.uniqueColumns, customSelectedColumns: this.state.customSelectedColumns
        });
      previewCardDiv.appendChild(card);
    };
    updatePreviewCard();

    const onAnythingChanged = DG.debounce(rxjs.merge(
      tableInput.onChanged
    //   columnInput.onChanged, // these changes are actually doing nothing :D
    //   matchRegexpInput.onChanged
    ), 500).subscribe(() => {
      updatePreviewCard();
    });

    dialog.add(ui.divH([ui.divV([
      semanticTypeInput.root,
      tableInput.root,
      columnInput.root,
      matchRegexpInput.root,
      regexpExampleCheckbox.root,
      regexpExampleDiv
    ], {classes: 'ui-form', style: {minWidth: '300px'}}), previewCardDiv]));

    dialog.addButton(isEdit ? 'Save' : 'Add', () => {
      if (!semanticTypeInput.value) {
        grok.shell.error('Semantic type is required');
        return;
      }
      if (!tableInput.value || !columnInput.value) {
        grok.shell.error('Both table and column are required');
        return;
      }

      if (regexpExampleCheckbox.value && (
        !exampleInput.value || !nonVarInput.value || !markupInput.value)) {
        grok.shell.error('All regexp example fields are required');
        return;
      }

      // If editing and name changed, delete old entry
      if (isEdit && semanticTypeInput.value !== semanticType)
        delete this.state.entryPoints[semanticType!];

      this.state.entryPoints[semanticTypeInput.value] = {
        table: tableInput.value || '',
        column: columnInput.value || '',
        matchRegexp: matchRegexpInput.value || undefined,
        regexpExample: regexpExampleCheckbox.value ? {
          example: exampleInput.value,
          nonVariablePart: nonVarInput.value,
          regexpMarkup: markupInput.value
        } : undefined
      };

      this.renderEntryPointsTable(tableContainer);
      dialog.close();
    });

    const closeSub = dialog.onClose.subscribe(() => {
      onAnythingChanged.unsubscribe();
      closeSub.unsubscribe();
    });

    dialog.show({resizable: true});
  }

  private createJoinOptionsSection(): HTMLDivElement {
    const section = ui.divV([], 'join-options-section');

    const infoText = ui.divText('Define how to include data from related tables', 'db-explorer-info-text');
    ui.tooltip.bind(infoText, 'Join options are useful for constructing cards that show up in tooltips and panels');

    const tableContainer = ui.divV([], 'join-options-table-container');
    this.renderJoinOptionsTable(tableContainer);

    const addButtonContainer = ui.divH([]);
    addButtonContainer.style.marginTop = '10px';
    const addButton = this.getAddIcon(() => this.showJoinOptionDialog(null, tableContainer), 'Add a new join option');
    addButtonContainer.appendChild(addButton);

    section.appendChild(infoText);
    section.appendChild(tableContainer);
    section.appendChild(addButtonContainer);

    return section;
  }

  private renderJoinOptionsTable(container: HTMLDivElement): void {
    container.innerHTML = '';

    if (this.state.joinOptions.length === 0) {
      container.appendChild(ui.divText('No join options defined yet', 'db-explorer-empty-message'));
      return;
    }

    const table = ui.table(
      this.state.joinOptions,
      (joinOption, index) => { // index is starting from 1, not 0, like... WTF
        const selectDisplay = joinOption.select.join(', ');
        return [
          joinOption.fromSchema ?? this.state.schemaName ?? '-',
          joinOption.fromTable,
          joinOption.columnName,
          joinOption.onSchema ?? this.state.schemaName ?? '-',
          joinOption.tableName,
          joinOption.onColumn,
          selectDisplay,
          this.createActionButtons(
            () => this.showJoinOptionDialog(index - 1, container),
            () => {
              this.state.joinOptions.splice(index - 1, 1);
              this.renderJoinOptionsTable(container);
            }
          )
        ];
      },
      ['From Schema', 'From Table', 'Column Name', 'To Schema', 'Join Table', 'On Column', 'Select Columns', 'Actions']
    );
    table.style.width = '100%';
    container.appendChild(table);
  }

  private showJoinOptionDialog(index: number | null, tableContainer: HTMLDivElement): void {
    const isEdit = index !== null;
    let onAnythingChangedRef: () => void = () => {};
    const joinOption = isEdit ? this.state.joinOptions[index!] : {
      fromSchema: this.state.schemaName || undefined,
      fromTable: this.tableNames[0] || '',
      columnName: '',
      onSchema: this.state.schemaName || undefined,
      tableName: this.tableNames[0] || '',
      onColumn: '',
      select: []
    };

    const dialog = ui.dialog(isEdit ? `Edit Join Option ${index! + 1}` : 'Add Join Option');

    const schemaItems = (this.schemaInput?.items ?? []) as string[];
    const fromSchemaInput = ui.input.choice('From Schema', {
      value: joinOption.fromSchema ?? this.state.schemaName ?? null,
      items: schemaItems,
      onValueChanged: onAnythingChangedRef
    });
    ui.tooltip.bind(fromSchemaInput.input, 'Schema of the starting table (the schema of the table that will be joined by the related table)');

    const fromTableInput = ui.input.choice('From Table', {value: joinOption.fromTable, items: this.tableNames, onValueChanged: onAnythingChangedRef});
    ui.tooltip.bind(fromTableInput.input, 'Starting table for the join (the table that will be joined by the related table)');

    const columnNameInput = ui.input.choice('Column Name', {
      value: joinOption.columnName,
      items: this.getSchemaColumnNames(joinOption.fromSchema ?? this.state.schemaName, joinOption.fromTable),
      onValueChanged: onAnythingChangedRef,
    });
    ui.tooltip.bind(columnNameInput.input, 'Foreign key column in starting table');

    const onSchemaInput = ui.input.choice('To Schema', {
      value: joinOption.onSchema ?? this.state.schemaName ?? null,
      items: schemaItems,
      onValueChanged: onAnythingChangedRef
    });
    ui.tooltip.bind(onSchemaInput.input, 'Schema of the related table (one that will join the starting table)');

    const tableNameInput = ui.input.choice('Table Name (to join)', {value: joinOption.tableName, items: this.tableNames, onValueChanged: onAnythingChangedRef});
    ui.tooltip.bind(tableNameInput.input, 'Related table to join (the table that will join the starting table)');

    const onColumnInput = ui.input.choice('On Column', {
      value: joinOption.onColumn,
      items: this.getSchemaColumnNames(joinOption.onSchema ?? this.state.schemaName, joinOption.tableName),
      onValueChanged: onAnythingChangedRef,
    });
    ui.tooltip.bind(onColumnInput.input, 'Key column in related table');

    this.ensureSchemaLoaded(fromSchemaInput.value ?? this.state.schemaName).then((entry) => {
      if (!entry)
        return;
      fromTableInput.items = entry.tableNames;
      if (fromTableInput.value && !fromTableInput.items.includes(fromTableInput.value))
        fromTableInput.value = fromTableInput.items[0] ?? null;
      columnNameInput.items = this.getSchemaColumnNames(fromSchemaInput.value ?? this.state.schemaName, fromTableInput.value);
    });

    this.ensureSchemaLoaded(onSchemaInput.value ?? this.state.schemaName).then((entry) => {
      if (!entry)
        return;
      tableNameInput.items = entry.tableNames;
      if (tableNameInput.value && !tableNameInput.items.includes(tableNameInput.value))
        tableNameInput.value = tableNameInput.items[0] ?? null;
      onColumnInput.items = this.getSchemaColumnNames(onSchemaInput.value ?? this.state.schemaName, tableNameInput.value);
    });

    fromSchemaInput.onChanged.subscribe(async () => {
      if (!fromSchemaInput.value)
        return;
      const entry = await this.ensureSchemaLoaded(fromSchemaInput.value);
      fromTableInput.items = entry?.tableNames ?? [];
      if (fromTableInput.value && !fromTableInput.items.includes(fromTableInput.value))
        fromTableInput.value = fromTableInput.items[0] ?? null;
      columnNameInput.items = this.getSchemaColumnNames(fromSchemaInput.value, fromTableInput.value);
    });

    fromTableInput.onChanged.subscribe(() => {
      columnNameInput.items = this.getSchemaColumnNames(fromSchemaInput.value ?? this.state.schemaName, fromTableInput.value);
    });

    onSchemaInput.onChanged.subscribe(async () => {
      if (!onSchemaInput.value)
        return;
      const entry = await this.ensureSchemaLoaded(onSchemaInput.value);
      tableNameInput.items = entry?.tableNames ?? [];
      if (tableNameInput.value && !tableNameInput.items.includes(tableNameInput.value))
        tableNameInput.value = tableNameInput.items[0] ?? null;
      onColumnInput.items = this.getSchemaColumnNames(onSchemaInput.value, tableNameInput.value);
      availableColumns = this.getSchemaColumnNames(onSchemaInput.value, tableNameInput.value);
      renderSelectColumns();
    });

    let availableColumns: string[] = this.getSchemaColumnNames(onSchemaInput.value ?? this.state.schemaName, joinOption.tableName);
    tableNameInput.onChanged.subscribe(() => {
      onColumnInput.items = this.getSchemaColumnNames(onSchemaInput.value ?? this.state.schemaName, tableNameInput.value);
      availableColumns = this.getSchemaColumnNames(onSchemaInput.value ?? this.state.schemaName, tableNameInput.value);
      renderSelectColumns();
    });

    // Select columns section
    const selectColumnsLabel = ui.label('Select Columns');
    ui.tooltip.bind(selectColumnsLabel, 'Columns to include from related table (format: "column" or "column as alias")');
    const selectedColumnsContainer = ui.divV([], 'selected-columns-container');

    const renderSelectColumns = () => {
      selectedColumnsContainer.innerHTML = '';
      joinOption.select.forEach((selectStr, colIndex) => {
        const colDiv = ui.divH([]);
        colDiv.style.marginBottom = '5px';
        colDiv.style.gap = '8px';
        colDiv.style.alignItems = 'center';

        // Parse existing string to extract column and alias
        const parts = selectStr.trim().split(/\s+as\s+/i);
        const column = parts[0].trim();
        const alias = parts.length > 1 ? parts[1].trim() : '';

        const columnChoice = ui.input.choice('', {
          value: column,
          items: availableColumns,
          onValueChanged: (newVal) => {
            if (newVal) {
              // If alias was same as old column or empty, update it
              const newAlias = (alias === column || !alias) ? newVal : alias;
              joinOption.select[colIndex] = newAlias ? `${newVal} as ${newAlias}` : newVal;
              renderSelectColumns();
            }
          }
        });
        columnChoice.root.style.flex = '1';

        const asLabel = ui.divText('as', 'db-explorer-as-label');
        asLabel.style.padding = '0 5px';
        asLabel.style.color = 'var(--grey-6)';

        const asInput = ui.input.string('', {
          value: alias,
          onValueChanged: (newVal) => {
            joinOption.select[colIndex] = newVal ? `${column} as ${newVal}` : column;
            onAnythingChangedRef();
          }
        });
        asInput.root.style.flex = '1';
        (asInput.input as HTMLInputElement).placeholder = 'alias (optional)';

        const deleteBtn = this.getDeleteIcon(() => {
          joinOption.select.splice(colIndex, 1);
          renderSelectColumns();
        }, 'Remove column');

        colDiv.appendChild(columnChoice.root);
        colDiv.appendChild(asLabel);
        colDiv.appendChild(asInput.root);
        colDiv.appendChild(deleteBtn);
        selectedColumnsContainer.appendChild(colDiv);
      });

      // Add plus icon at the bottom
      if (availableColumns.length > 0) {
        const addColumnDiv = ui.divH([]);
        addColumnDiv.style.marginTop = '5px';
        const addIcon = this.getAddIcon(() => {
          const firstCol = availableColumns[0] || '';
          joinOption.select.push(firstCol);
          renderSelectColumns();
        }, 'Add column');
        addColumnDiv.appendChild(addIcon);
        selectedColumnsContainer.appendChild(addColumnDiv);
      }
      onAnythingChangedRef();
    };

    const previewCardDiv = ui.divV([], {classes: 'db-explorer-entrypoint-preview-card', style: {minWidth: '350px', flexGrow: '1'}});
    let prevTimer: any = null;
    const updatePreviewCard = () => {
      if (prevTimer)
        clearTimeout(prevTimer);
      prevTimer = setTimeout(() => {
        ui.empty(previewCardDiv);
        if (!this.dataConnection || !this.state.schemaName || !fromTableInput.value || !tableNameInput.value || !columnNameInput.value || !onColumnInput.value) {
          previewCardDiv.appendChild(ui.divText('Please fill in all join option fields to see preview', 'db-explorer-preview-message'));
          return;
        }
        // build the joinOptions array such that it includes the
        const filteredJoins = this.state.joinOptions.filter((_, i) => i !== index);
        const tempJoinOptions = [...filteredJoins];
        const tempJoinOption: JoinOption = {
          fromSchema: fromSchemaInput.value ?? this.state.schemaName ?? undefined,
          fromTable: fromTableInput.value!,
          columnName: columnNameInput.value!,
          onSchema: onSchemaInput.value ?? this.state.schemaName ?? undefined,
          tableName: tableNameInput.value!,
          onColumn: onColumnInput.value!,
          select: joinOption.select
        };
        tempJoinOptions.push(tempJoinOption);
        const previewSchema = fromSchemaInput.value ?? this.state.schemaName!;
        const card = renderExampleCard(
          this.dataConnection, previewSchema, fromTableInput.value!, {joinOptions: tempJoinOptions, customRenderers: this.state.customRenderers,
            uniqueColumns: this.state.uniqueColumns, // note: specifically leaving out custom selected columns to show all columns in preview
          });
        previewCardDiv.appendChild(card);
      }, 500);
    };
    onAnythingChangedRef = () => {
      updatePreviewCard();
    };

    renderSelectColumns();


    dialog.add(ui.divH([ui.divV([
      fromSchemaInput.root,
      fromTableInput.root,
      columnNameInput.root,
      onSchemaInput.root,
      tableNameInput.root,
      onColumnInput.root,
      selectColumnsLabel,
      selectedColumnsContainer
    ], {classes: 'ui-form', style: {minWidth: '300px'}}), previewCardDiv]));

    dialog.addButton(isEdit ? 'Save' : 'Add', () => {
      if (!fromTableInput.value || !columnNameInput.value || !tableNameInput.value || !onColumnInput.value) {
        grok.shell.error('All fields are required');
        return;
      }

      if (joinOption.select.length === 0) {
        grok.shell.error('At least one select column is required');
        return;
      }

      const newJoinOption: JoinOption = {
        fromSchema: fromSchemaInput.value ?? this.state.schemaName ?? undefined,
        fromTable: fromTableInput.value,
        columnName: columnNameInput.value,
        onSchema: onSchemaInput.value ?? this.state.schemaName ?? undefined,
        tableName: tableNameInput.value,
        onColumn: onColumnInput.value,
        select: joinOption.select
      };

      if (isEdit)
        this.state.joinOptions[index!] = newJoinOption;
      else
        this.state.joinOptions.push(newJoinOption);

      this.renderJoinOptionsTable(tableContainer);
      dialog.close();
    });

    dialog.show({resizable: true});
  }

  private createHeaderNamesSection(): HTMLDivElement {
    const section = ui.divV([], 'header-names-section');

    const infoText = ui.divText('Customize column headers in the exploration view', 'db-explorer-info-text');
    ui.tooltip.bind(infoText, 'Each key is a table name, and the value is the column to use as the header if there are multiple entries');

    const tableContainer = ui.divV([], 'header-names-table-container');
    this.renderHeaderNamesTable(tableContainer);

    const addButtonContainer = ui.divH([]);
    addButtonContainer.style.marginTop = '10px';
    const addButton = this.getAddIcon(() => this.showHeaderNameDialog(null, tableContainer), 'Add a custom header name mapping');
    addButtonContainer.appendChild(addButton);

    section.appendChild(infoText);
    section.appendChild(tableContainer);
    section.appendChild(addButtonContainer);

    return section;
  }

  private renderHeaderNamesTable(container: HTMLDivElement): void {
    container.innerHTML = '';

    if (Object.keys(this.state.headerNames).length === 0) {
      container.appendChild(ui.divText('No header names defined yet', 'db-explorer-empty-message'));
      return;
    }

    const table = ui.table(
      Object.keys(this.state.headerNames),
      (tableName) => {
        const columnName = this.state.headerNames[tableName];
        return [
          tableName,
          columnName,
          this.createActionButtons(
            () => this.showHeaderNameDialog(tableName, container),
            () => {
              delete this.state.headerNames[tableName];
              this.renderHeaderNamesTable(container);
            }
          )
        ];
      },
      ['Table', 'Column', 'Actions']
    );
    table.style.width = '100%';
    container.appendChild(table);
  }

  private showHeaderNameDialog(tableName: string | null, tableContainer: HTMLDivElement): void {
    const isEdit = tableName !== null;
    const columnName = isEdit ? this.state.headerNames[tableName!] : '';

    const dialog = ui.dialog(isEdit ? `Edit Header Name: ${tableName}` : 'Add Header Name');

    const tableInput = ui.input.choice('Table', {value: tableName || this.tableNames[0] || '', items: this.tableNames});
    ui.tooltip.bind(tableInput.input, 'Table name');

    const columnInput = ui.input.choice('Column', {
      value: columnName,
      items: this.columnNamesByTable[tableName || this.tableNames[0] || ''] || []
    });
    ui.tooltip.bind(columnInput.input, 'Column to use as header');

    tableInput.onChanged.subscribe(() => {
      if (tableInput.value)
        columnInput.items = this.columnNamesByTable[tableInput.value] || [];
    });

    dialog.add(ui.divV([
      tableInput.root,
      columnInput.root
    ], 'ui-form'));

    dialog.addButton(isEdit ? 'Save' : 'Add', () => {
      if (!tableInput.value || !columnInput.value) {
        grok.shell.error('Both table and column are required');
        return;
      }

      // If editing and name changed, delete old entry
      if (isEdit && tableInput.value !== tableName)
        delete this.state.headerNames[tableName!];

      this.state.headerNames[tableInput.value] = columnInput.value;
      this.renderHeaderNamesTable(tableContainer);
      dialog.close();
    });

    dialog.show({resizable: true});
  }

  private createUniqueColumnsSection(): HTMLDivElement {
    const section = ui.divV([], 'unique-columns-section');

    const infoText = ui.divText('Specify which column uniquely identifies rows in tables', 'db-explorer-info-text');
    ui.tooltip.bind(infoText, 'Fallback if database does not contain primary key information');

    const tableContainer = ui.divV([], 'unique-columns-table-container');
    this.renderUniqueColumnsTable(tableContainer);

    const addButtonContainer = ui.divH([]);
    addButtonContainer.style.marginTop = '10px';
    const addButton = this.getAddIcon(() => this.showUniqueColumnDialog(null, tableContainer), 'Add a unique column mapping');
    addButtonContainer.appendChild(addButton);

    section.appendChild(infoText);
    section.appendChild(tableContainer);
    section.appendChild(addButtonContainer);

    return section;
  }

  private renderUniqueColumnsTable(container: HTMLDivElement): void {
    container.innerHTML = '';

    if (Object.keys(this.state.uniqueColumns).length === 0) {
      container.appendChild(ui.divText('No unique columns defined yet', 'db-explorer-empty-message'));
      return;
    }

    const table = ui.table(
      Object.keys(this.state.uniqueColumns),
      (tableName) => {
        const columnName = this.state.uniqueColumns[tableName];
        return [
          tableName,
          columnName,
          this.createActionButtons(
            () => this.showUniqueColumnDialog(tableName, container),
            () => {
              delete this.state.uniqueColumns[tableName];
              this.renderUniqueColumnsTable(container);
            }
          )
        ];
      },
      ['Table', 'Column', 'Actions']
    );
    table.style.width = '100%';
    container.appendChild(table);
  }

  private showUniqueColumnDialog(tableName: string | null, tableContainer: HTMLDivElement): void {
    const isEdit = tableName !== null;
    const columnName = isEdit ? this.state.uniqueColumns[tableName!] : '';

    const dialog = ui.dialog(isEdit ? `Edit Unique Column: ${tableName}` : 'Add Unique Column');

    const tableInput = ui.input.choice('Table', {value: tableName || this.tableNames[0] || '', items: this.tableNames});
    ui.tooltip.bind(tableInput.input, 'Table name');

    const columnInput = ui.input.choice('Column', {
      value: columnName,
      items: this.columnNamesByTable[tableName || this.tableNames[0] || ''] || []
    });
    ui.tooltip.bind(columnInput.input, 'Column that uniquely identifies rows');

    tableInput.onChanged.subscribe(() => {
      if (tableInput.value)
        columnInput.items = this.columnNamesByTable[tableInput.value] || [];
    });

    dialog.add(ui.divV([
      tableInput.root,
      columnInput.root
    ], 'ui-form'));

    dialog.addButton(isEdit ? 'Save' : 'Add', () => {
      if (!tableInput.value || !columnInput.value) {
        grok.shell.error('Both table and column are required');
        return;
      }

      // If editing and name changed, delete old entry
      if (isEdit && tableInput.value !== tableName)
        delete this.state.uniqueColumns[tableName!];

      this.state.uniqueColumns[tableInput.value] = columnInput.value;
      this.renderUniqueColumnsTable(tableContainer);
      dialog.close();
    });

    dialog.show({resizable: true});
  }

  private createCustomSelectedColumnsSection(): HTMLDivElement {
    const section = ui.divV([], 'custom-selected-columns-section');

    const infoText = ui.divText('Override which columns to display for specific tables', 'db-explorer-info-text');
    ui.tooltip.bind(infoText, 'This allows showing more relevant information in tooltips and cards');

    const tableContainer = ui.divV([], 'custom-selected-columns-table-container');
    this.renderCustomSelectedColumnsTable(tableContainer);

    const addButtonContainer = ui.divH([]);
    addButtonContainer.style.marginTop = '10px';
    const addButton = this.getAddIcon(() => this.showCustomSelectedColumnsDialog(null, tableContainer), 'Add custom column selection for a table');
    addButtonContainer.appendChild(addButton);

    section.appendChild(infoText);
    section.appendChild(tableContainer);
    section.appendChild(addButtonContainer);

    return section;
  }

  private renderCustomSelectedColumnsTable(container: HTMLDivElement): void {
    container.innerHTML = '';

    if (Object.keys(this.state.customSelectedColumns).length === 0) {
      container.appendChild(ui.divText('No custom selected columns defined yet', 'db-explorer-empty-message'));
      return;
    }

    const table = ui.table(
      Object.keys(this.state.customSelectedColumns),
      (tableName) => {
        const columns = this.state.customSelectedColumns[tableName];
        return [
          tableName,
          columns.join(', '),
          this.createActionButtons(
            () => this.showCustomSelectedColumnsDialog(tableName, container),
            () => {
              delete this.state.customSelectedColumns[tableName];
              this.renderCustomSelectedColumnsTable(container);
            }
          )
        ];
      },
      ['Table', 'Columns', 'Actions']
    );
    table.style.width = '100%';
    container.appendChild(table);
  }

  private showCustomSelectedColumnsDialog(tableName: string | null, tableContainer: HTMLDivElement): void {
    const isEdit = tableName !== null;
    const columns = isEdit ? [...this.state.customSelectedColumns[tableName!]] : [];

    let onAnythingChangedRef: () => void = () => {};

    const dialog = ui.dialog(isEdit ? `Edit Custom Columns: ${tableName}` : 'Add Custom Columns');

    const tableInput = ui.input.choice('Table', {value: tableName || this.tableNames[0] || '', items: this.tableNames,
      onValueChanged: () => onAnythingChangedRef()
    });
    ui.tooltip.bind(tableInput.input, 'Table name');

    let availableColumns = this.getAvailableColumnsForTable(tableInput.value || '');
    if (columns.length === 0 && availableColumns.length > 0) {
      // by default, select all columns
      availableColumns.forEach((c) => columns.push(c));
    }

    // Columns selector section
    const columnsLabel = ui.label('Columns');
    ui.tooltip.bind(columnsLabel, 'Select columns to display for this table (includes base columns and join columns)');
    const columnsContainer = ui.divV([], 'columns-container');

    const renderColumns = () => {
      columnsContainer.innerHTML = '';
      columns.forEach((col, colIndex) => {
        const colDiv = ui.divH([]);
        colDiv.style.marginBottom = '5px';
        colDiv.style.gap = '8px';
        colDiv.style.alignItems = 'center';

        const columnChoice = ui.input.choice('', {
          value: col,
          items: availableColumns,
          onValueChanged: (newVal) => {
            if (newVal)
              columns[colIndex] = newVal;
            onAnythingChangedRef();
          }
        });
        columnChoice.root.style.flex = '1';

        const deleteBtn = this.getDeleteIcon(() => {
          columns.splice(colIndex, 1);
          renderColumns();
        }, 'Remove column');

        colDiv.appendChild(columnChoice.root);
        colDiv.appendChild(deleteBtn);
        columnsContainer.appendChild(colDiv);
      });

      // Add plus icon at the bottom
      if (availableColumns.length > 0) {
        const addColumnDiv = ui.divH([]);
        addColumnDiv.style.marginTop = '5px';
        const addIcon = this.getAddIcon(() => {
          const firstCol = availableColumns[0] || '';
          columns.push(firstCol);
          renderColumns();
        }, 'Add column');
        addColumnDiv.appendChild(addIcon);
        columnsContainer.appendChild(addColumnDiv);
      }
      onAnythingChangedRef();
    };

    tableInput.onChanged.subscribe(() => {
      if (tableInput.value) {
        availableColumns = this.getAvailableColumnsForTable(tableInput.value);
        // when table changes, reset selected columns
        if (columns.some((c) => !availableColumns.includes(c)) || columns.length === 0) {
          columns.length = 0;
          availableColumns.forEach((c) => columns.push(c));
        }

        renderColumns();
      }
    });


    const previewCardDiv = ui.divV([], {classes: 'db-explorer-entrypoint-preview-card', style: {minWidth: '350px', flexGrow: '1'}});
    let prevTimer: any = null;
    const updatePreviewCard = () => {
      if (prevTimer)
        clearTimeout(prevTimer);
      prevTimer = setTimeout(() => {
        ui.empty(previewCardDiv);
        if (!this.dataConnection || !this.state.schemaName || !tableInput.value) {
          previewCardDiv.appendChild(ui.divText('Please select a table to see the preview', 'db-explorer-preview-message'));
          return;
        }
        const modifiedSelection = {...this.state.customSelectedColumns, [tableInput.value]: columns};
        if (columns.length === 0)
          delete modifiedSelection[tableInput.value];
        const card = renderExampleCard(
          this.dataConnection, this.state.schemaName, tableInput.value, {joinOptions: this.state.joinOptions, customRenderers: this.state.customRenderers,
            uniqueColumns: this.state.uniqueColumns, customSelectedColumns: modifiedSelection,
          });
        previewCardDiv.appendChild(card);
      }, 500);
    };
    onAnythingChangedRef = () => {
      updatePreviewCard();
    };

    renderColumns();

    dialog.add(ui.divH([ui.divV([
      tableInput.root,
      columnsLabel,
      columnsContainer
    ], {classes: 'ui-form', style: {minWidth: '20px'}}), previewCardDiv]));

    dialog.addButton(isEdit ? 'Save' : 'Add', () => {
      if (!tableInput.value) {
        grok.shell.error('Table is required');
        return;
      }

      // If editing and name changed, delete old entry
      if (isEdit && tableInput.value !== tableName)
        delete this.state.customSelectedColumns[tableName!];

      this.state.customSelectedColumns[tableInput.value] = columns;
      this.renderCustomSelectedColumnsTable(tableContainer);
      dialog.close();
    });

    dialog.show({resizable: true});
  }

  private createExplicitReferencesSection(): HTMLDivElement {
    const section = ui.divV([], 'explicit-references-section');

    const infoText = ui.divText('Define explicit references between tables and schemas', 'db-explorer-info-text');
    ui.tooltip.bind(infoText, 'Explicit references allow linking schema-table-columns that are not in DB schema metadata');

    const tableContainer = ui.divV([], 'explicit-references-table-container');
    this.renderExplicitReferencesTable(tableContainer);

    const addButtonContainer = ui.divH([]);
    addButtonContainer.style.marginTop = '10px';
    const addButton = this.getAddIcon(() => this.showExplicitReferenceDialog(null, tableContainer), 'Add an explicit reference');
    addButtonContainer.appendChild(addButton);

    section.appendChild(infoText);
    section.appendChild(tableContainer);
    section.appendChild(addButtonContainer);

    return section;
  }

  private renderExplicitReferencesTable(container: HTMLDivElement): void {
    container.innerHTML = '';

    if (!this.state.explicitReferences || this.state.explicitReferences.length === 0) {
      container.appendChild(ui.divText('No explicit references defined yet', 'db-explorer-empty-message'));
      return;
    }

    const table = ui.table(
      this.state.explicitReferences,
      (ref, index) => {
        return [
          ref.schema ?? this.state.schemaName ?? '-',
          ref.table,
          ref.column,
          ref.refSchema ?? this.state.schemaName ?? '-',
          ref.refTable,
          ref.refColumn,
          this.createActionButtons(
            () => this.showExplicitReferenceDialog(index - 1, container),
            () => {
              this.state.explicitReferences!.splice(index - 1, 1);
              this.renderExplicitReferencesTable(container);
            }
          )
        ];
      },
      ['Schema', 'Table', 'Column', 'Ref Schema', 'Ref Table', 'Ref Column', 'Actions']
    );
    table.style.width = '100%';
    container.appendChild(table);
  }

  private showExplicitReferenceDialog(index: number | null, tableContainer: HTMLDivElement): void {
    const isEdit = index !== null;
    const ref = isEdit ? this.state.explicitReferences![index!] : {
      schema: this.state.schemaName ?? undefined,
      table: this.tableNames[0] || '',
      column: '',
      refSchema: this.state.schemaName ?? undefined,
      refTable: this.tableNames[0] || '',
      refColumn: ''
    };

    const dialog = ui.dialog(isEdit ? `Edit Explicit Reference ${index! + 1}` : 'Add Explicit Reference');
    const schemaItems = (this.schemaInput?.items ?? []) as string[];

    const schemaInput = ui.input.choice('Schema', {value: ref.schema ?? this.state.schemaName ?? null, items: schemaItems});
    ui.tooltip.bind(schemaInput.input, 'Schema of the source table');

    const tableInput = ui.input.choice('Table', {value: ref.table, items: this.getSchemaTableNames(schemaInput.value ?? this.state.schemaName)});
    ui.tooltip.bind(tableInput.input, 'Source table');

    const columnInput = ui.input.choice('Column', {
      value: ref.column,
      items: this.getSchemaColumnNames(schemaInput.value ?? this.state.schemaName, tableInput.value)
    });
    ui.tooltip.bind(columnInput.input, 'Source column');

    const refSchemaInput = ui.input.choice('Ref Schema', {value: ref.refSchema ?? this.state.schemaName ?? null, items: schemaItems});
    ui.tooltip.bind(refSchemaInput.input, 'Schema of the referenced table');

    const refTableInput = ui.input.choice('Ref Table', {value: ref.refTable, items: this.getSchemaTableNames(refSchemaInput.value ?? this.state.schemaName)});
    ui.tooltip.bind(refTableInput.input, 'Referenced table');

    const refColumnInput = ui.input.choice('Ref Column', {
      value: ref.refColumn,
      items: this.getSchemaColumnNames(refSchemaInput.value ?? this.state.schemaName, refTableInput.value)
    });
    ui.tooltip.bind(refColumnInput.input, 'Referenced column');

    this.ensureSchemaLoaded(schemaInput.value ?? this.state.schemaName).then((entry) => {
      if (!entry)
        return;
      tableInput.items = entry.tableNames;
      if (tableInput.value && !tableInput.items.includes(tableInput.value))
        tableInput.value = tableInput.items[0] ?? null;
      columnInput.items = this.getSchemaColumnNames(schemaInput.value ?? this.state.schemaName, tableInput.value);
    });

    this.ensureSchemaLoaded(refSchemaInput.value ?? this.state.schemaName).then((entry) => {
      if (!entry)
        return;
      refTableInput.items = entry.tableNames;
      if (refTableInput.value && !refTableInput.items.includes(refTableInput.value))
        refTableInput.value = refTableInput.items[0] ?? null;
      refColumnInput.items = this.getSchemaColumnNames(refSchemaInput.value ?? this.state.schemaName, refTableInput.value);
    });

    schemaInput.onChanged.subscribe(async () => {
      if (!schemaInput.value)
        return;
      const entry = await this.ensureSchemaLoaded(schemaInput.value);
      tableInput.items = entry?.tableNames ?? [];
      if (tableInput.value && !tableInput.items.includes(tableInput.value))
        tableInput.value = tableInput.items[0] ?? null;
      columnInput.items = this.getSchemaColumnNames(schemaInput.value, tableInput.value);
    });

    tableInput.onChanged.subscribe(() => {
      columnInput.items = this.getSchemaColumnNames(schemaInput.value ?? this.state.schemaName, tableInput.value);
    });

    refSchemaInput.onChanged.subscribe(async () => {
      if (!refSchemaInput.value)
        return;
      const entry = await this.ensureSchemaLoaded(refSchemaInput.value);
      refTableInput.items = entry?.tableNames ?? [];
      if (refTableInput.value && !refTableInput.items.includes(refTableInput.value))
        refTableInput.value = refTableInput.items[0] ?? null;
      refColumnInput.items = this.getSchemaColumnNames(refSchemaInput.value, refTableInput.value);
    });

    refTableInput.onChanged.subscribe(() => {
      refColumnInput.items = this.getSchemaColumnNames(refSchemaInput.value ?? this.state.schemaName, refTableInput.value);
    });

    dialog.add(ui.divV([
      schemaInput.root,
      tableInput.root,
      columnInput.root,
      refSchemaInput.root,
      refTableInput.root,
      refColumnInput.root
    ], 'ui-form'));

    dialog.addButton(isEdit ? 'Save' : 'Add', () => {
      if (!schemaInput.value || !tableInput.value || !columnInput.value || !refSchemaInput.value || !refTableInput.value || !refColumnInput.value) {
        grok.shell.error('All fields are required');
        return;
      }

      const newRef: ExplicitReference = {
        schema: schemaInput.value ?? undefined,
        table: tableInput.value,
        column: columnInput.value,
        refSchema: refSchemaInput.value ?? undefined,
        refTable: refTableInput.value,
        refColumn: refColumnInput.value
      };

      if (!this.state.explicitReferences)
        this.state.explicitReferences = [];

      if (isEdit)
        this.state.explicitReferences[index!] = newRef;
      else
        this.state.explicitReferences.push(newRef);

      this.renderExplicitReferencesTable(tableContainer);
      dialog.close();
    });

    dialog.show({resizable: true});
  }

  private createCustomRenderersSection(): HTMLDivElement {
    const section = ui.divV([], 'custom-renderers-section');

    const infoText = ui.divText('Define custom renderers for specific table columns', 'db-explorer-info-text');
    ui.tooltip.bind(infoText, 'Custom renderers allow you to display column data in special formats (images, molecules, etc.)');

    const tableContainer = ui.divV([], 'custom-renderers-table-container');
    this.renderCustomRenderersTable(tableContainer);

    const addButtonContainer = ui.divH([]);
    addButtonContainer.style.marginTop = '10px';
    const addButton = this.getAddIcon(() => this.showCustomRendererDialog(null, tableContainer), 'Add a custom renderer');
    addButtonContainer.appendChild(addButton);

    section.appendChild(infoText);
    section.appendChild(tableContainer);
    section.appendChild(addButtonContainer);

    return section;
  }

  private renderCustomRenderersTable(container: HTMLDivElement): void {
    container.innerHTML = '';

    if (!this.state.customRenderers || this.state.customRenderers.length === 0) {
      container.appendChild(ui.divText('No custom renderers defined yet', 'db-explorer-empty-message'));
      return;
    }

    const table = ui.table(
      this.state.customRenderers,
      (renderer, index) => {
        return [
          renderer.table,
          renderer.column,
          renderer.renderer,
          this.createActionButtons(
            () => this.showCustomRendererDialog(index - 1, container),
            () => {
              this.state.customRenderers!.splice(index - 1, 1);
              this.renderCustomRenderersTable(container);
            }
          )
        ];
      },
      ['Table', 'Column', 'Renderer', 'Actions']
    );
    table.style.width = '100%';
    container.appendChild(table);
  }

  private showCustomRendererDialog(index: number | null, tableContainer: HTMLDivElement): void {
    const isEdit = index !== null;
    const renderer = isEdit ? this.state.customRenderers![index!] : {
      table: this.tableNames[0] || '',
      column: '',
      renderer: 'molecule' as const
    };

    const dialog = ui.dialog(isEdit ? `Edit Custom Renderer ${index! + 1}` : 'Add Custom Renderer');

    const tableInput = ui.input.choice('Table', {value: renderer.table, items: this.tableNames});
    ui.tooltip.bind(tableInput.input, 'Table containing the column to render');

    const columnInput = ui.input.choice('Column', {
      value: renderer.column,
      items: this.columnNamesByTable[renderer.table] || []
    });
    ui.tooltip.bind(columnInput.input, 'Column to apply custom renderer to');

    const rendererInput = ui.input.choice('Renderer', {
      value: renderer.renderer,
      items: ['rawImage', 'imageURL', 'molecule', 'helm']
    });
    ui.tooltip.bind(rendererInput.input, 'Type of renderer to use');

    tableInput.onChanged.subscribe(() => {
      if (tableInput.value)
        columnInput.items = this.columnNamesByTable[tableInput.value] || [];
    });

    dialog.add(ui.divV([
      tableInput.root,
      columnInput.root,
      rendererInput.root
    ], 'ui-form'));

    dialog.addButton(isEdit ? 'Save' : 'Add', () => {
      if (!tableInput.value || !columnInput.value || !rendererInput.value) {
        grok.shell.error('All fields are required');
        return;
      }

      const newRenderer: CustomRenderer = {
        table: tableInput.value,
        column: columnInput.value,
        renderer: rendererInput.value as 'rawImage' | 'imageURL' | 'molecule' | 'helm'
      };

      if (!this.state.customRenderers)
        this.state.customRenderers = [];

      if (isEdit)
        this.state.customRenderers[index!] = newRenderer;
      else
        this.state.customRenderers.push(newRenderer);

      this.renderCustomRenderersTable(tableContainer);
      dialog.close();
    });

    dialog.show({resizable: true});
  }

  private createButtonBar(): void {
    const buttonBar = ui.divH([], 'db-explorer-action-bar');
    buttonBar.style.marginTop = '20px';
    buttonBar.style.gap = '10px';

    //const validateButton = ui.button('Validate', () => this.validateConfig());
    const exportButton = ui.button('Export JSON', () => this.exportConfig());
    const importButton = ui.button('Import JSON', () => this.importConfig());

    //buttonBar.appendChild(validateButton);
    buttonBar.appendChild(importButton);
    buttonBar.appendChild(exportButton);

    if (this.onSave) {
      const saveButton = ui.bigButton('Save', () => {
        const validation = this.validateConfig();
        if (validation.valid && this.onSave) {
          const config = this.getConfig();
          this.onSave(config);
        }
      });
      buttonBar.appendChild(saveButton);
    }

    this.rootDiv.appendChild(buttonBar);
  }

  public validateConfig(): {valid: boolean; errors: string[]} {
    this.validationErrors = [];

    // Validate basic fields
    if (!this.state.connectionNqName)
      this.validationErrors.push('Connection name is required');
    if (!this.state.schemaName)
      this.validationErrors.push('Schema name is required');

    // Validate entry points
    if (Object.keys(this.state.entryPoints).length === 0)
      this.validationErrors.push('At least one identifier is required');

    Object.keys(this.state.entryPoints).forEach((semanticType) => {
      const ep = this.state.entryPoints[semanticType];
      if (!ep.table)
        this.validationErrors.push(`Identifier '${semanticType}': table is required`);
      if (!ep.column)
        this.validationErrors.push(`Identifier '${semanticType}': column is required`);

      if (ep.regexpExample) {
        if (!ep.regexpExample.example)
          this.validationErrors.push(`Identifier '${semanticType}': regexp example is incomplete`);
        if (!ep.regexpExample.nonVariablePart)
          this.validationErrors.push(`Identifier '${semanticType}': regexp non-variable part is required`);
        if (!ep.regexpExample.regexpMarkup)
          this.validationErrors.push(`Identifier '${semanticType}': regexp markup is required`);
      }
    });

    // Validate join options
    this.state.joinOptions.forEach((jo, index) => {
      if (!jo.fromTable)
        this.validationErrors.push(`Join option ${index + 1}: fromTable is required`);
      if (!jo.columnName)
        this.validationErrors.push(`Join option ${index + 1}: columnName is required`);
      if (!jo.tableName)
        this.validationErrors.push(`Join option ${index + 1}: tableName is required`);
      if (!jo.onColumn)
        this.validationErrors.push(`Join option ${index + 1}: onColumn is required`);
      if (jo.select.length === 0)
        this.validationErrors.push(`Join option ${index + 1}: at least one select column is required`);

      jo.select.forEach((sel, selIndex) => {
        if (!sel || sel.trim() === '')
          this.validationErrors.push(`Join option ${index + 1}, select ${selIndex + 1}: column is required`);
      });
    });

    // Validate explicit references
    (this.state.explicitReferences ?? []).forEach((ref, index) => {
      if (!ref.schema)
        this.validationErrors.push(`Explicit reference ${index + 1}: schema is required`);
      if (!ref.table)
        this.validationErrors.push(`Explicit reference ${index + 1}: table is required`);
      if (!ref.column)
        this.validationErrors.push(`Explicit reference ${index + 1}: column is required`);
      if (!ref.refSchema)
        this.validationErrors.push(`Explicit reference ${index + 1}: refSchema is required`);
      if (!ref.refTable)
        this.validationErrors.push(`Explicit reference ${index + 1}: refTable is required`);
      if (!ref.refColumn)
        this.validationErrors.push(`Explicit reference ${index + 1}: refColumn is required`);
    });

    // Display validation results
    if (this.validationErrors.length > 0) {
      grok.shell.error(`Validation failed with ${this.validationErrors.length} error(s)`);
      const errorList = ui.divV(this.validationErrors.map((err) => ui.divText(`• ${err}`)));
      grok.shell.info(errorList);
      return {valid: false, errors: this.validationErrors};
    }
    return {valid: true, errors: []};
  }

  public getConfig(): DBExplorerConfig {
    const config: DBExplorerConfig = {
      connectionName: this.state.connectionName!,
      nqName: this.state.connectionNqName ?? undefined,
      schemaName: this.state.schemaName!,
      dataSourceName: this.dataConnection?.dataSource || undefined,
      entryPoints: this.state.entryPoints,
      joinOptions: this.state.joinOptions,
      headerNames: this.state.headerNames,
      uniqueColumns: this.state.uniqueColumns,
      customSelectedColumns: this.state.customSelectedColumns,
    };

    // Only include customRenderers if it has values
    if (this.state.customRenderers && this.state.customRenderers.length > 0)
      config.customRenderers = this.state.customRenderers;

    if (this.state.explicitReferences && this.state.explicitReferences.length > 0)
      config.explicitReferences = this.state.explicitReferences;

    return config;
  }

  private exportConfig(): void {
    const validation = this.validateConfig();
    if (!validation.valid)
      return;

    const config = this.getConfig();
    const jsonStr = JSON.stringify(config, null, 2);

    // Copy to clipboard
    navigator.clipboard.writeText(jsonStr).then(() => {
      grok.shell.info('Configuration copied to clipboard');
    });

    // Also show in a dialog
    const dialog = ui.dialog('DB Explorer Configuration');
    const textArea = ui.element('textarea') as HTMLTextAreaElement;
    textArea.value = jsonStr;
    textArea.style.width = '600px';
    textArea.style.height = '400px';
    textArea.style.fontFamily = 'monospace';
    dialog.add(textArea);
    dialog.show({resizable: true});
  }


  /**
   * BE WARE: getUI method needs to be called before setConfig, otherwise UI elements and connections won't be initialized
   */
  public async setConfig(config: Partial<DBExplorerConfig>, warnIfNotFound = true): Promise<void> {
    // Update state
    this.state = {
      connectionName: config.connectionName || null,
      connectionNqName: config.nqName || null,
      connectionProvider: config.dataSourceName || null,
      schemaName: config.schemaName || null,
      entryPoints: config.entryPoints || {},
      joinOptions: config.joinOptions || [],
      headerNames: config.headerNames || {},
      uniqueColumns: config.uniqueColumns || {},
      customSelectedColumns: config.customSelectedColumns || {},
      customRenderers: config.customRenderers || [],
      explicitReferences: config.explicitReferences || []
    };

    // Update UI inputs
    if (this.connectionInput && this.state.connectionName) {
      if (this.dataConnection) { this.connectionInput.value = this.dataConnection.nqName; } else {
        // if the nqName is not set (which can be, trust me), try to find it from connection name
        const connections = this.connections;
        const con = connections.find((c) => c.name.toLowerCase() === this.state.connectionName!.toLowerCase() &&
                (!this.state.connectionProvider || c.dataSource === this.state.connectionProvider));
        if (con) {
          this.state.connectionNqName = con.nqName;
          this.state.connectionName = con.name;
          this.state.connectionProvider = con.dataSource;
          this.connectionInput.value = con.nqName;
        }
      }
    }

    // Trigger schema load if connection and schema are set
    if (this.state.connectionName && this.state.schemaName && this.dataConnection) {
      const connection = this.dataConnection;
      const schemas = await grok.dapi.connections.getSchemas(connection);
      this.schemaInput.items = schemas.sort();
      this.schemaInput.value = this.state.schemaName ?? config.schemaName ?? null;
      await this.loadSchemaInfo();
      await this.buildAdvancedSections();
    } else {
      if (warnIfNotFound)
        grok.shell.warning('Error loading configuration: Connection or schema not found.');
    }
  }

  private importConfig(): void {
    const dialog = ui.dialog('Import Configuration');
    const textArea = ui.element('textarea') as HTMLTextAreaElement;
    textArea.placeholder = 'Paste your JSON configuration here...';
    textArea.style.width = '600px';
    textArea.style.height = '400px';
    textArea.style.fontFamily = 'monospace';

    dialog.add(textArea);
    dialog.addButton('Import', async () => {
      try {
        const config = JSON.parse(textArea.value);
        await this.setConfig(config);

        grok.shell.info('Configuration imported successfully');
        dialog.close();
      } catch (e: any) {
        grok.shell.error(`Failed to import configuration: ${e.message}`);
      }
    });

    dialog.show({resizable: true});
  }
}

