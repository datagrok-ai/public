import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/**
 * Class AddNewColumnDialog is a useful method to add a new column to the table.
 *
 * It uses user-friendly and functional Dialog Window that allows to use
 * formulas, functions and other columns to create a new column and immediately
 * see the preview result.
 */
export class AddNewColumnDialog {
  title: string = 'Add New Column Ext';
  helpUrl: string = '/help/transform/add-new-column.md';
  visibleTags: string[] = ['math', 'text', 'date', 'timespan', 'binning', 'logic', 'stats'];
  supportedTypes: string[] = [DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.INT,
      DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.DATE_TIME, DG.COLUMN_TYPE.BOOL];
  defaultType: string = DG.COLUMN_TYPE.STRING;
  autoType: string = 'auto';
  plainTextType: string = 'plain text';
  placeholderName: string = 'Name';
  placeholderType: string = 'Type';  // Used only for uniformity when saving inputs history.
  placeholderExpression: string = 'Expression';
  maxAutoNameLength: number = 25;
  maxPreviwRowCount: number = 100;
  newColumnBgColor: number = 0xFFFDFFE7;  // The same bg-color as the bg-color of tooltips.
  colNamePattern: RegExp = /\${(.+?)}/g;
  tooltips = {
    name: 'New unique column name.',
    type: 'Column type. When set to "auto", type is determined based on the expression.',
    expression: `Formula for calculating new column values.<br>Columns and functions can be drag-n-dropped into this field.`,
    preview: 'Preview result columns.'
  };

  sourceDf: DG.DataFrame | null = null;
  previwDf?: DG.DataFrame;

  inputName?: DG.InputBase;
  inputType?: DG.InputBase;
  inputExpression?: DG.InputBase;
  uiPreview?: HTMLDivElement;
  widgetColumns?: HTMLDivElement;
  widgetFunctions?: HTMLDivElement;
  previewGrid?: DG.Grid;
  uiDialog?: DG.Dialog;
  resultColumnType?: string;

  constructor() {
    this.sourceDf = grok.shell.t;  // An attempt to access the active table.
    if (this.sourceDf)
      this.init();
    else
      grok.shell.error('No open tables found');
  }

  /** Initializes all parameters and opens a Dialog Window. */
  async init(): Promise<void> {
    this.uiDialog = ui.dialog({ title: this.title, helpUrl: this.helpUrl })
        .add(await this.initUiLayout())
        .onOK(async () => await this.addNewColumnAction()).show();

    this.uiDialog.history(
      () => this.saveInputHistory(),
      (x) => this.loadInputHistory(x)
    );

    await this.updatePreview();
  }

  /** Returns values of the input fields to be saved in history. */
  saveInputHistory(): any {
    let typeForHistory =
        this.getSelectedType()[0] == this.autoType
        ? this.resultColumnType
        : this.inputType!.value;

    return Object.fromEntries([
      [this.placeholderName, this.inputName!.value],
      [this.placeholderType, typeForHistory],
      [this.placeholderExpression, this.inputExpression!.value]
    ]);
  }

  /** Loads values of the input fields from the history and applies them. */
  loadInputHistory(history: any): void {
    this.inputName!.value = history[this.placeholderName];
    this.inputType!.value = history[this.placeholderType];
    this.inputExpression!.value = history[this.placeholderExpression];
  }

  /** Creates and initializes the "Column Name" input field. */
  initInputName(): DG.InputBase {
    let control = ui.stringInput('', '');
        control.onInput(async () => await this.updatePreview());
        control.setTooltip(this.tooltips['name']);

    let input = control.input as HTMLInputElement;
        input.classList.add('ui-input-addnewcolumn-name');
        input.placeholder = this.placeholderName;

    return control;
  }

  /** Creates and initializes the "Column Type" input field. */
  initInputType(): DG.InputBase {
    let defaultChoise = `${this.autoType} (${this.defaultType})`;
    this.supportedTypes.unshift(defaultChoise);    // The first item of the ChoiceBox will be "Auto".
    this.supportedTypes.push(this.plainTextType);  // The last item of the ChoiceBox will be "Treat As String".

    let control = ui.choiceInput('', defaultChoise, this.supportedTypes);
        control.onInput(async () => await this.updatePreview());
        control.setTooltip(this.tooltips['type']);

    let input = control.input as HTMLInputElement;
        input.classList.add('ui-input-addnewcolumn-type');
        input.insertBefore(ui.element('hr'),       // Separator before the last item in the ChoiceBox.
            (input as any)[(input as any).length - 1]);

    return control;
  }

  /** Creates and initializes the "Expression" input field. */
  initInputExpression(): DG.InputBase {
    let control = ui.textInput('', '', async () => {
      // The first characters of the Expression become the default name of the new column:
      (this.inputName!.input as HTMLInputElement).placeholder =
          (!control.value || (control.value.length > this.maxAutoNameLength))
          ? this.placeholderName
          : control.value;
      await this.updatePreview();
    });
    control.setTooltip(this.tooltips['expression']);

    let input = control.input as HTMLInputElement;
        input.classList.add('ui-input-addnewcolumn-expression');
        input.placeholder = this.placeholderExpression;

    // Columns and functions can be drag-n-dropped into the Expression field:
    ui.makeDroppable(input, {
      acceptDrop: (dragObject) => this.typeOf(dragObject, DG.Column, DG.Func),
      doDrop: (dragObject, _) => this.insertIntoExpression(dragObject)
    });

    return control;
  }

  /** Creates and initializes the Preview Grid. */
  initUiPreviewGrid(): HTMLDivElement {
    // Limiting the number of rows in the Preview Grid:
    let previewRowCount = Math.min(this.sourceDf!.rowCount, this.maxPreviwRowCount);
    this.previwDf = this.sourceDf!.clone(DG.BitSet.create(previewRowCount, (idx) => idx < previewRowCount));
    this.previewGrid = DG.Viewer.grid(this.previwDf!);
    this.previewGrid.root.classList.add('ui-grid-with-thin-scrollbars');

    // Making the Preview Grid less interactive:
    let props = this.previewGrid.props;
        props.showCellTooltip = false;
        props.showCurrentCellOutline = false;
        props.showDefaultPopupMenu = false;
        props.allowDynamicMenus = false;
        props.autoScrollRowIntoView = false;
        props.autoScrollColumnIntoView = false;
        props.allowRowDragging = false;
        props.allowColReordering = false;
        props.allowEdit = false;
        props.allowColSelection = false;
        props.allowRowSelection = false;
        props.allowBlockSelection = false;
        props.colHeaderFont = props.defaultCellFont;

    let control = ui.div();
        control.append(this.previewGrid.root);
        control.classList.add('ui-addnewcolumn-preview');
    ui.tooltip.bind(control, this.tooltips['preview']);

    return control;
  }

  /** Creates and initializes the "Column List Widget". */
  async initWidgetColumns(): Promise<HTMLDivElement> {
    let widget = await DG.Func.byName('ColumnGridWidget').apply({ df: this.sourceDf });

    let control = ui.box();
        control.append(widget.root);
        control.classList.add('ui-widget-addnewcolumn-columns');

    return control;
  }

  /** Creates and initializes the "Function List Widget". */
  async initWidgetFunctions(): Promise<HTMLDivElement> {
    let widget = await DG.Func.byName('FunctionsWidget').apply();
        widget.props.visibleTags = this.visibleTags.join(',');
        widget.props.showSignature = true;

    let control = ui.box();
        control.append(widget.root);
        control.classList.add('ui-widget-addnewcolumn-functions');

    return control;
  }

  /** Creates the Global UI Layout for the Dialog Window. */
  async initUiLayout(): Promise<HTMLDivElement> {
    this.inputName = this.initInputName();
    this.inputType = this.initInputType();
    this.inputExpression = this.initInputExpression();
    this.uiPreview = this.initUiPreviewGrid();
    this.widgetColumns = await this.initWidgetColumns();
    this.widgetFunctions = await this.initWidgetFunctions();

    let layout =
        ui.div([
            ui.block50([
                ui.block([this.inputName.root], { style: { width: '65%' } }),
                ui.block([this.inputType.root], { style: { width: '35%' } }),
                ui.block([this.inputExpression.root]),
                ui.block([this.uiPreview])
            ], { style: { paddingRight: '20px' } }),

            ui.block25([
                ui.block([this.widgetColumns])
            ], { style: { paddingRight: '20px' } }),

            ui.block25([
                ui.block([this.widgetFunctions])
            ])
        ]);
    layout.classList.add('ui-addnewcolumn-layout');

    return layout;
  }

  /** Updates the Preview Grid. Executed every time the controls are changed. */
  async updatePreview(): Promise<void> {
    // Making the Column List for the Preview Grid:
    let columnIds = this.findUniqueColumnNamesInExpression();
    let colName = this.getResultColumnName();
    columnIds.push(colName);

    // Making the Preview Grid:
    await this.previwDf!.columns.addNewCalculated(
        colName,
        this.inputExpression!.value,
        ...this.getSelectedType()
    );
    this.previewGrid!.dataFrame = this.previwDf!.clone(null, columnIds);
    this.previewGrid!.col(colName)!.backColor = this.newColumnBgColor;
    this.resultColumnType = this.previwDf!.col(colName)!.type;
    this.previwDf!.columns.remove(colName);

    this.setAutoType();  // Adding (or removing) the column auto-type caption to "Auto" item in the ChoiceBox.
  }

  /** Finds all unique column names used in the Expression input field. */
  findUniqueColumnNamesInExpression(): string[] {
    let names = this.inputExpression!.value.match(this.colNamePattern) as string[];

    if (!names)
      return [];

    return names
        .filter((v, i, a) => a.indexOf(v) === i)                     // Selecting only unique names.
        .map((name) => name.slice(2, -1))                            // Removing the naming syntax: ${}.
        .filter((v) => this.sourceDf!.columns.names().includes(v));  // Selecting only real column names.
  }

  /** Inserts drag-n-dropping object into the Expression input field. */
  insertIntoExpression(x: any): void {
    let snippet: string = '';

    if (this.typeOf(x, DG.Column))
      snippet = `\${${x.name}}`;
    else if (this.typeOf(x, DG.Func))
      snippet = `${x.name}()`;
    else
      return;

    let input = this.inputExpression!.input as HTMLInputElement;
    let text = input.value;
    let cursorPos = input.selectionStart!;
    input.value = text.slice(0, cursorPos) + snippet + text.slice(cursorPos);

    this.inputExpression!.fireChanged();
  }

  /** Creates new unique Column Name. */
  getResultColumnName(): string {
    let input = this.inputName!.input as HTMLInputElement;

    return DG.utils.getUniqueName(
        input.value
        ? input.value
        : input.placeholder,
        this.sourceDf!.columns.names()
    );
  }

  /** Detects selected item in the ChoiceBox. */
  getSelectedType(): any[] {
    let selectedType: any[] = [];

    if (this.inputType!.value == this.plainTextType)
      selectedType.push(DG.COLUMN_TYPE.STRING, true);  // Prepare for AddNewColumn with the treatAsString param.
    else
      selectedType.push(
          this.inputType!.value
          ? this.inputType!.value.split(' ')[0]
          : this.autoType
      );

    return selectedType;
  }

  /** Adds (or removes) the column auto-type caption to "Auto" item in the ChoiceBox. */
  setAutoType(): void {
    let input = (this.inputType!.input as any)[0];

    if (this.getSelectedType()[0] == this.autoType)
      input.innerHTML = `${this.autoType} (${this.resultColumnType})`;  // Adding auto-type caption if "Auto" item is selected.
    else
      input.innerHTML = this.autoType;  // Removing auto-type caption from "Auto" item if any other items are selected.
  }

  /** Checks if an object is an instance of one of the listed classes. */
  typeOf(x: any, ...types: any[]): boolean {
    return types.some((t) => x instanceof t);
  }

  /** Adds a New Column to the source table. */
  async addNewColumnAction(): Promise<void> {
    await this.sourceDf!.columns.addNewCalculated(
        this.getResultColumnName(),
        this.inputExpression!.value,
        ...this.getSelectedType()
    );
  }

  /** Closes Add New Column Dialog Window. */
  close(): void {
    this.uiDialog!.close();
  }

  /** Opens Add New Column Dialog Window. */
  open(): void {
    this.uiDialog!.show();
  }
}
