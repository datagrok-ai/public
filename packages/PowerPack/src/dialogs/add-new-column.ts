/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { EditorState, Extension } from '@codemirror/state'
import { EditorView, ViewUpdate, hoverTooltip } from '@codemirror/view'
import {Completion, CompletionContext, CompletionResult, autocompletion, completeFromList} from "@codemirror/autocomplete"

/**
 * Class AddNewColumnDialog is a useful method to add a new column to the table
 * or to edit formula for an existing column.
 *
 * It uses user-friendly and functional Dialog Window that allows to use
 * formulas, functions and other columns to create a new column and immediately
 * see the preview result.
 */

type PropInfo = {
  propName: string,
  propType: string
}

const VALIDATION_TYPES_MAPPING: { [key: string]: string[] } = {
  'num': ['number', 'int', 'double', 'float'],
  'number': ['num', 'int', 'double', 'float'],
  'double': ['int', 'float', 'number', 'num'],
  'float': ['int']
};

const COLUMN_FUNCTION_NAME = 'GetCurrentRowField';
const GET_VAR_FUNCTION_NAME = 'GetVar';
const RESERVED_FUNC_NAMES_AND_TYPES: {[key: string]: string} = {
  'GetRowIndex': DG.TYPE.INT,
}

export class AddNewColumnDialog {
  addColumnTitle: string = 'Add New Column';
  editColumnTitle: string = 'Edit Column Formula';
  helpUrl: string = '/help/transform/add-new-column.md';
  visibleTags: string[] = ['math', 'text', 'date', 'timespan', 'binning', 'logic', 'stats'];
  supportedTypes: string[] = [DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.INT,
    DG.COLUMN_TYPE.STRING, DG.COLUMN_TYPE.DATE_TIME, DG.COLUMN_TYPE.BOOL, DG.COLUMN_TYPE.QNUM];
  defaultType: string = DG.COLUMN_TYPE.STRING;
  autoType: string = 'auto';
  plainTextType: string = 'plain text';
  placeholderName: string = 'Name';
  placeholderType: string = 'Type'; // Used only for uniformity when saving inputs history.
  placeholderExpression: string = 'Expression';
  maxAutoNameLength: number = 50;
  maxPreviwRowCount: number = 100;
  newColumnBgColor: number = 0xFFFDFFE7; // The same bg-color as the bg-color of tooltips.
  colNamePattern: RegExp = /\${(.+?)}/g;
  tooltips = {
    name: 'Сolumn name.',
    type: 'Column type. When set to "auto", type is determined based on the expression.',
    expression: `Formula for calculating column values.<br>
      Columns and functions can be drag-n-dropped into this field.`,
    preview: 'Preview result columns.',
  };

  sourceDf: DG.DataFrame | null = null; // Represents Source Table.
  previwDf?: DG.DataFrame; // Represents Preview Table.
  columnsDf?: DG.DataFrame; // Represents Columns Widget.
  gridPreview?: DG.Grid;
  widgetColumns?: DG.Widget;
  widgetFunctions?: DG.Widget;
  resultColumnType?: string;
  dialogTitle: string = '';
  call: DG.FuncCall | null = null;
  edit: boolean = false;

  inputName?: DG.InputBase;
  inputType?: DG.InputBase;
  uiPreview?: HTMLDivElement;
  uiColumns?: HTMLDivElement;
  uiFunctions?: HTMLDivElement;
  uiDialog?: DG.Dialog;
  codeMirror?: EditorView;
  codeMirrorDiv?: HTMLDivElement;
  errorDiv?: HTMLDivElement;
  columnNames: string[] = [];
  coreFunctionsNames: string[] = [];
  coreFunctionsParams: {[key: string]: PropInfo[]} = {};
  packageFunctionsNames: {[key: string]: string[]} = {};
  packageFunctionsParams: {[key: string]: PropInfo[]} = {};
  packageNames: string[] = [];

  constructor(call: DG.FuncCall | null = null) {
    const table = call?.getParamValue('table');

    if (table) {
      this.call = call;
      this.sourceDf = table;
      this.edit = table.columns.names().includes(call?.getParamValue('name'));
    } else
      this.sourceDf = grok.shell.t;


    if (this.sourceDf) {
      this.columnNames = this.sourceDf.columns.names();
      this.init();
    }
    else
      grok.shell.error('Table not found');
  }

  /** Initializes all parameters and opens a Dialog Window. */
  async init(): Promise<void> {
    this.dialogTitle = this.edit ? this.editColumnTitle : this.addColumnTitle;

    this.uiDialog = ui.dialog({title: this.dialogTitle, helpUrl: this.helpUrl});

    this.inputName = this.initInputName();
    this.inputType = this.initInputType();
    const codeMirrorRes = this.initCodeMirror();
    this.codeMirror = codeMirrorRes.codeMirror;
    this.codeMirrorDiv = codeMirrorRes.coreMirrorDiv;
    this.errorDiv = codeMirrorRes.errorDiv;

    // Not necessary, but if the Dialog knows about inputs, then it can implement extra-logic:
    this.uiDialog
      .add(this.inputName)
      .add(this.inputType)

    this.uiDialog
      .add(await this.initUiLayout())
      .onOK(async () => await this.addNewColumnAction())
      .show({resizable: true, width: 750, height: 500});

    this.uiDialog.history(
      () => this.saveInputHistory(),
      (x) => this.loadInputHistory(x),
    );

    this.prepareForSeleniumTests();
    await this.updatePreview(this.codeMirror!.state.doc.toString());
    this.prepareFunctionsListForAutocomplete();
    
  }

  prepareFunctionsListForAutocomplete() {
    //filter functions with one input (multiple inputs or functions returning void are not included)
    const allFunctionsList = DG.Func.find().filter((it) => it.outputs.length === 1);
    for (const func of allFunctionsList) {
      const params: PropInfo[] = func.inputs.map((it) => {
        return {propName: it.name, propType: it.propertyType};
      });
      try {
        const packageName = func.package.name;
        if (!this.packageFunctionsNames[packageName]) {
          this.packageNames.push(packageName);
          this.packageFunctionsNames[packageName] = [];
        }
        this.packageFunctionsNames[packageName].push(func.name);
        this.packageFunctionsParams[`${packageName}:${func.name}`] = params;
      } catch { //in case of core functions calling func.package throws an exception
        this.coreFunctionsNames.push(func.name);
        this.coreFunctionsParams[func.name] = params;
      }
    }
  }

  /** Prepares form inputs for Selenium tests by given them unique names. */
  prepareForSeleniumTests() {
    (this.inputName!.input as HTMLInputElement).name = 'input-Add-New-Column---Name';
    (this.inputType!.input as HTMLInputElement).name = 'input-Add-New-Column---Type';
    this.uiDialog!.getButton('OK').name = 'button-Add-New-Column---OK';
    this.uiDialog!.getButton('CANCEL').name = 'button-Add-New-Column---CANCEL';
  }

  /** Returns values of the input fields to be saved in history. */
  saveInputHistory(): any {
    const typeForHistory =
        this.getSelectedType()[0] == this.autoType ?
          this.resultColumnType :
        this.inputType!.value;

    return Object.fromEntries([
      [this.placeholderName, this.inputName!.value],
      [this.placeholderType, typeForHistory],
      [this.placeholderExpression, this.codeMirror!.state.doc.toString()],
    ]);
  }

  /** Loads values of the input fields from the history and applies them. */
  loadInputHistory(history: any): void {
    this.inputName!.value = history[this.placeholderName];
    this.inputType!.value = history[this.placeholderType];
    this.codeMirror!.dispatch({changes: {
      from: 0,
      to: this.codeMirror!.state.doc.length,
      insert: history[this.placeholderExpression]
    }});
  }

  /** Creates and initializes the "Column Name" input field. */
  initInputName(): DG.InputBase {
    const control = ui.stringInput('', '');
    control.onInput(async () => await this.updatePreview(this.codeMirror!.state.doc.toString()));
    control.setTooltip(this.tooltips['name']);

    const input = control.input as HTMLInputElement;
    input.classList.add('ui-input-addnewcolumn-name');
    input.placeholder = this.placeholderName;
    if (this.call)
      input.value = this.call.getParamValue('name');

    return control;
  }

  /** Creates and initializes the "Column Type" input field. */
  initInputType(): DG.InputBase {
    const defaultChoise = `${this.autoType} (${this.defaultType})`;
    this.supportedTypes.unshift(defaultChoise); // The first item of the ChoiceBox will be "Auto".
    this.supportedTypes.push(this.plainTextType); // The last item of the ChoiceBox will be "Treat As String".

    const control = ui.choiceInput('', this.call ?
      this.call.getParamValue('type') : defaultChoise, this.supportedTypes);
    control.onInput(async () => await this.updatePreview(this.codeMirror!.state.doc.toString()));
    control.setTooltip(this.tooltips['type']);

    const input = control.input as HTMLInputElement;
    input.classList.add('ui-input-addnewcolumn-type');
    input.insertBefore(ui.element('hr'), // Separator before the last item in the ChoiceBox.
      (input as any)[(input as any).length - 1]);

    return control;
  }


  initCodeMirror(): {codeMirror: EditorView, coreMirrorDiv: HTMLDivElement, errorDiv: HTMLDivElement} {
    const editorDiv = ui.div('', { style: { height: '140px' } });
    const formulaErrorDiv = ui.div();
    //let errorMsg = '';
    editorDiv.onclick = () => {
      cm.focus();
    }
    // Columns and functions can be drag-n-dropped into the Expression field:
    ui.makeDroppable(editorDiv, {
      acceptDrop: (dragObject) => this.typeOf(dragObject, DG.Column, DG.Func),
      doDrop: (dragObject, _) => {
        cm.focus();
        this.insertIntoCodeMirror(dragObject, cm);
      },
    });

    const autocomplete = autocompletion({
      activateOnTyping: true,
      override: [this.functionsCompletions(this.columnNames, this.packageNames, this.coreFunctionsNames,
        this.packageFunctionsNames, this.packageFunctionsParams, this.coreFunctionsParams)],
    });

    const wordHover = this.hoverTooltipCustom(this.packageFunctionsParams, this.coreFunctionsParams);

    const cm = new EditorView({
      parent: editorDiv,
      state: EditorState.create({
        doc: '',
        extensions: [
          EditorView.lineWrapping,
          autocomplete,
          wordHover,
          EditorView.updateListener.of(async (e: ViewUpdate) => {
            if (!e.docChanged)
              return;
            const cmValue = cm.state.doc.toString();
            (this.inputName!.input as HTMLInputElement).placeholder =
              ((!cmValue || (cmValue.length > this.maxAutoNameLength)) ? this.placeholderName : cmValue).trim();
            const error = cmValue ? this.validateFormula(cmValue) : '';
            ui.empty(formulaErrorDiv);
            if (!error) {
              await this.updatePreview(cmValue);
            }
            else {
              formulaErrorDiv.append(ui.divText(error, 'cm-error-div'));
            }
          }),
        ],
      }),
    });

    if (this.call)
      this.codeMirror!.dispatch({changes: {
        from: 0,
        to: this.codeMirror!.state.doc.length,
        insert: this.call.getParamValue('expression')
      }});

    return {codeMirror: cm, coreMirrorDiv: editorDiv, errorDiv: formulaErrorDiv};
  }

  validateFormula(formula: string): string {
    //check unmatched columns
    const matchesAll = [...formula.matchAll(this.colNamePattern)];
    const unmatchedCols: string[] = [];
    if (matchesAll.length) {
      for (const match of matchesAll)
        if (!this.columnNames.includes(match[1]))
          unmatchedCols.push(match[1]);
    }
    if (unmatchedCols.length)
      return unmatchedCols.length > 1 ? `Columns ${unmatchedCols.join(',')} are missing` :
        `Column ${unmatchedCols[0]} is missing`;
    //check syntax errors
    try {
      const funcCall = grok.functions.parse(formula, false);
      console.log(funcCall);
      this.validateFuncCallTypes(funcCall);
    } catch (e: any) {
      return e.message;
    }
    return '';
  }

  validateFuncCallTypes(funcCall: DG.FuncCall) {

    const innerFuncCalls: string[] = [];
    const actualInputParamTypes: { [key: string]: string } = {};

    //collect actual input parameter types
    for (const key of Object.keys(funcCall.inputs)) {
      if (funcCall.inputs[key] instanceof DG.FuncCall) {
        //do not allow functions with multiple inputs 
        if (funcCall.inputs[key].func.outputs.length > 1)
          throw new Error(`Function ${funcCall.inputs[key].func} returns multiple values`);
        //treat $COLUMN_NAME as scalar
        if (funcCall.inputs[key].func.name !== COLUMN_FUNCTION_NAME)
          innerFuncCalls.push(key);
        actualInputParamTypes[key] = funcCall.inputs[key].func.outputs[0].propertyType;
      } else
        actualInputParamTypes[key] = typeof funcCall.inputs[key];
    }

    //validate types for current function
    for (const property of funcCall.func.inputs) {
      let actualInputType = actualInputParamTypes[property.name];
      const input = funcCall.inputs[property.name];
      //check for variables missing in the context
      //TODO: preview with variables doesn't work
      if (input.func && input.func.name === GET_VAR_FUNCTION_NAME) {
        try {
          const res = input.callSync();
          actualInputType = typeof res.getOutputParamValue();
        } catch (e: any) {
          //throw new Error(`Variable ${funcCall.inputs[property.name].inputs['variableName']} is not declared`);
          throw e;
        }
      }
      //handling dynamic type in actual input
      if (actualInputType === DG.TYPE.DYNAMIC) {
        //extract type from column in case $COLUMN_NAME was passed to formula
        if (input.func && input.func.name === COLUMN_FUNCTION_NAME)
          actualInputType = this.sourceDf!.col(input.inputs['field'])!.type;
        //handling 'row' function
        //TODO: Handle other similar functions
        for (const reservedFunc of Object.keys(RESERVED_FUNC_NAMES_AND_TYPES)) {
          if (input.func.name === reservedFunc) {
            actualInputType = RESERVED_FUNC_NAMES_AND_TYPES[reservedFunc];
            break;
          }
        }
      }
      //dynamic allows any type
      if (property.propertyType === DG.TYPE.DYNAMIC)
        continue;
      //check for exact match
      if (property.propertyType === actualInputType)
        continue;
      //check for type match in mapping
      if (VALIDATION_TYPES_MAPPING[property.propertyType] && VALIDATION_TYPES_MAPPING[property.propertyType].includes(actualInputType))
        continue;
      throw new Error(`Function ${funcCall.func.name} '${property.name}' param should be ${property.propertyType} type instead of ${actualInputType}`);
    }
    //validate inner func calls recursively
    if (innerFuncCalls.length) {
      for (const param of innerFuncCalls)
        this.validateFuncCallTypes(funcCall.inputs[param]);
    }
  }

  hoverTooltipCustom(packageFunctionsParams: { [key: string]: PropInfo[] },
    coreFunctionsParams: { [key: string]: PropInfo[] }): Extension {
    return hoverTooltip((view: EditorView, pos: number, side: number) => {
      let {from, to, text} = view.state.doc.lineAt(pos)
      let start = pos, end = pos
      while (start > from && /\w|:/.test(text[start - from - 1]))
        start--;
      while (end < to && /\w|:/.test(text[end - from]))
        end++;
      if (start == pos && side < 0 || end == pos && side > 0)
        return null;
      return {
        pos: start,
        end,
        above: true,
        create(view) {
          let dom = document.createElement("div");
          const funcName = text.slice(start - from, end - from);
          const funcParams = funcName.includes(':') ? packageFunctionsParams[funcName] : coreFunctionsParams[funcName]; 
          dom.textContent = `${funcName}${funcParams ? `(${funcParams.map((it) => `${it.propName}:${it.propType}`).join(',')})` : ''}`;
          return {dom}
        }
      }
    });
  }



  /** Creates and initializes the Preview Grid. */
  initUiPreview(): HTMLDivElement {
    // Limiting the number of rows in the Preview Grid:
    const previewRowCount = Math.min(this.sourceDf!.rowCount, this.maxPreviwRowCount);
    this.previwDf = this.sourceDf!.clone(DG.BitSet.create(previewRowCount, (idx) => idx < previewRowCount));
    this.gridPreview = DG.Viewer.grid(this.previwDf!);
    this.gridPreview.root.classList.add('ui-grid-with-thin-scrollbars');

    // Making the Preview Grid less interactive:
    const props = this.gridPreview.props;
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

    const previewRoot = this.gridPreview.root;
    previewRoot.setAttribute('style', 'height: -webkit-fill-available !important;; min-height:225px;');

    const control = ui.div(previewRoot);
    control.append(this.gridPreview.root);
    control.classList.add('ui-addnewcolumn-preview');
    control.style.flexGrow = '1';
    ui.tooltip.bind(control, this.tooltips['preview']);

    return control;
  }

  /** Creates and initializes the "Column List Widget". */
  async initUiColumns(): Promise<HTMLDivElement> {
    this.widgetColumns = await DG.Func.byName('ColumnGridWidget').apply({df: this.sourceDf});

    this.columnsDf = DG.toJs(this.widgetColumns!.props.dfColumns);

    const control = ui.box();
    control.append(this.widgetColumns!.root);
    control.classList.add('ui-widget-addnewcolumn-columns');

    return control;
  }

  /** Creates and initializes the "Function List Widget". */
  async initUiFunctions(): Promise<HTMLDivElement> {
    this.widgetFunctions = await DG.Func.byName('FunctionsWidget').apply();
    this.widgetFunctions!.props.visibleTags = this.visibleTags.join(',');
    this.widgetFunctions!.props.showSignature = true;
    const control = ui.box();
    control.append(this.widgetFunctions!.root);
    control.classList.add('ui-widget-addnewcolumn-functions');
    control.style.height = 'inherit';
    return control;
  }

  /** Creates the Global UI Layout for the Dialog Window. */
  async initUiLayout(): Promise<HTMLDivElement> {
    this.uiPreview = this.initUiPreview();
    this.uiColumns = await this.initUiColumns();
    this.uiFunctions = await this.initUiFunctions();

    const flexStyle = {display: 'flex', flexGrow: '1'};
    const layout =
        ui.div([
          ui.block50([
            ui.block([
              ui.block([this.inputName!.root], {style: {width: '65%'}}),
              ui.block([this.inputType!.root], {style: {width: '35%'}}),
            ]),
            ui.block([this.codeMirrorDiv!, this.errorDiv!]),
            ui.block([this.uiPreview], {style: flexStyle}),
          ], {style: Object.assign({}, {paddingRight: '20px', flexDirection: 'column'}, flexStyle)}),

          ui.block25([
            ui.block([this.uiColumns], {style: Object.assign({}, {flexDirection: 'column'}, flexStyle)}),
          ], {style: Object.assign({}, {paddingRight: '20px'}, flexStyle)}),

          ui.block25([
            ui.block([this.uiFunctions], {style: flexStyle}),
          ], {style: flexStyle}),
        ]);
    layout.classList.add('ui-addnewcolumn-layout');
    return layout;
  }

  /** Updates the Preview Grid. Executed every time the controls are changed. */
  async updatePreview(expression: string): Promise<void> {
    // Making the Column List for the Preview Grid:
    const columnIds = this.findUniqueColumnNamesInExpression(expression);
    const colName = this.getResultColumnName();
    columnIds.push(colName);

    const type = this.getSelectedType()[0];
    // Making the Preview Grid:
    const call = (DG.Func.find({name: 'AddNewColumn'})[0]).prepare({table: this.previwDf!,
      name: colName, expression: expression, type: type});
    await call.call(false, undefined, {processed: true, report: false});
    /*    await this.previwDf!.columns.addNewCalculated(
        colName,
        this.inputExpression!.value,
        ...this.getSelectedType()
    );*/
    this.gridPreview!.dataFrame = this.previwDf!.clone(null, columnIds);
    this.gridPreview!.col(colName)!.backColor = this.newColumnBgColor;
    this.resultColumnType = this.previwDf!.col(colName)!.type;
    this.previwDf!.columns.remove(colName);

    this.setAutoType(); // Adding (or removing) the column auto-type caption to "Auto" item in the ChoiceBox.
  }

  /** Finds all unique column names used in the Expression input field. */
  findUniqueColumnNamesInExpression(expression: string): string[] {
    const expr = expression;
    const names = expr.match(this.colNamePattern) as string[];

    if (!names)
      return [];

    const lcSourceColumnNames = (this.sourceDf!.columns.names() as string[]).map((name) => name.toLowerCase());

    return names
      .filter((v, i, a) => a.indexOf(v) === i) // Selecting only unique names.
      .map((name) => name.slice(2, -1)) // Removing the naming syntax: ${}.
      .filter((v) => lcSourceColumnNames.includes(v.toLowerCase())); // Selecting only real column names.
  }

  /** Inserts drag-n-dropping object into the Expression input field. */
  insertIntoCodeMirror(x: any, cm: EditorView): void {
    let snippet: string = '';

    if (this.typeOf(x, DG.Column))
      snippet = `\${${x.name}}`;
    else if (this.typeOf(x, DG.Func)) {
      const paramsStr = (x as DG.Func).inputs.map((it) =>it.propertyType).join(', ');
      snippet = `${x.name}(${paramsStr})`;
    }
    else
      return;
    const value = cm.state.doc.toString();
    const cursorPos = cm.state.selection.main.head;
    const newValue = value.slice(0, cursorPos) + snippet + value.slice(cursorPos);
    cm.dispatch({
      changes: {
        from: 0,
        to: cm.state.doc.length,
        insert: newValue
      }
    })
  }

  /** Creates new unique Column Name. */
  getResultColumnName(): string {
    const input = this.inputName!.input as HTMLInputElement;
    let value: string = input.value;
    if ((value ?? '') == '')
      value = input.placeholder;
    return this.edit ?
      value :
      this.sourceDf!.columns.getUnusedName(value);
  }

  /** Detects selected item in the ChoiceBox. */
  getSelectedType(): any[] {
    const selectedType: any[] = [];

    if (this.inputType!.value == this.plainTextType)
      selectedType.push(DG.COLUMN_TYPE.STRING, true); // Prepare for AddNewColumn with the treatAsString param.
    else {
      selectedType.push(
          this.inputType!.value ?
          this.inputType!.value.split(' ')[0] :
            this.autoType,
          false);
    }

    return selectedType;
  }

  /** Adds (or removes) the column auto-type caption to "Auto" item in the ChoiceBox. */
  setAutoType(): void {
    const input = (this.inputType!.input as any)[0];

    if (this.getSelectedType()[0] == this.autoType) {
      // Adding auto-type caption if "Auto" item is selected.
      input.innerHTML = `${this.autoType} (${this.resultColumnType})`;
    } else
      input.innerHTML = this.autoType; // Removing auto-type caption from "Auto" item if any other items are selected.
  }

  /** Checks if an object is an instance of one of the listed classes. */
  typeOf(x: any, ...types: any[]): boolean {
    return types.some((t) => x instanceof t);
  }

  // Returns types and semtypes of selected items in the Columns Widget.
  getTypesOfSelectedColumns(): Array<string[]> {
    const selectedIndexes: Int32Array = this.columnsDf!.selection.getSelectedIndexes();
    const selectedColumns: DG.Column[] = (this.sourceDf!.columns as DG.ColumnList)
      .toList().filter((v, i) => selectedIndexes.includes(i));
    let types: string[] = [];
    let semTypes: string[] = [];
    selectedColumns.forEach((v) => {types.push(v.type), semTypes.push(v.semType);});
    types = types.filter((v, i, a) => a.indexOf(v) === i);
    semTypes = semTypes.filter((v, i, a) => a.indexOf(v) === i && v != null);
    return [types, semTypes];
  }

  /** Adds a New Column to the source table or edit formula for an existing column. */
  async addNewColumnAction(): Promise<void> {
    if (this.edit) {
      this.call!.setParamValue('name', this.inputName!.value);
      this.call!.setParamValue('expression', this.codeMirror!.state.doc.toString());
      this.call!.setParamValue('type', this.getSelectedType()[0]);
      this.call!.setParamValue('treatAsString', this.getSelectedType()[1]);
      await this.call!.call();
    } else {
      await this.sourceDf!.columns.addNewCalculated(
        this.getResultColumnName(),
          this.codeMirror!.state.doc.toString().trim(),
          ...this.getSelectedType(),
      );
    }
  }

  /** Closes Add New Column Dialog Window. */
  close(): void {
    this.uiDialog!.close();
  }

  /** Opens Add New Column Dialog Window. */
  open(): void {
    this.uiDialog!.show({resizable: true});
  }

  functionsCompletions(colNames: string[], packageNames: string[], 
    coreFunctionsNames: string[], packageFunctionsNames: {[key: string]: string[]},
    packageFunctionsParams: {[key: string]: PropInfo[]}, coreFunctionsParams:  {[key: string]: PropInfo[]}) {
    return (context: CompletionContext) => {
      let word = context.matchBefore(/[\w|:|$]*/);
      if (!word || word?.from === word?.to && !context.explicit)
        return null;
      let options: Completion[] = [];
      let index = word!.from;
      let filter = true;
      if (word.text.includes(':')) {
        const colonIdx = word.text.indexOf(":");
        const packName = word.text.substring(0, word.text.indexOf(":"));
        if (packageFunctionsNames[packName])
          packageFunctionsNames[packName].forEach((name: string, idx: number) => {
            options.push({ label: name, type: "variable",
              apply: `${name}(${packageFunctionsParams[`${packName}:${name}`].map((it)=> it.propType).join(',')})` })
          });
        index = word!.from + colonIdx + 1;
        filter = !word.text.endsWith(':');
      } else if (word.text.includes('$')) {
        colNames.forEach((name: string) => options.push({ label: name, type: "variable", apply: `{${name}}` }));
        index = word!.from + word.text.indexOf("$") + 1;
        filter = !word.text.endsWith('$');
      } else
        coreFunctionsNames.concat(packageNames)
          .forEach((name: string, idx: number) => options.push({
            label: name, type: "variable",
            apply: idx < coreFunctionsNames.length ? `${name}(${coreFunctionsParams[name].map((it)=> it.propType).join(',')})` : name,
            detail: idx < coreFunctionsNames.length ? '' : 'package'
          }));;
      return {
        from: index,
        options: options,
        filter: filter
      } as CompletionResult
    };
  }
}


