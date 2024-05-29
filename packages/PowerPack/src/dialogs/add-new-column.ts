/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { EditorState } from '@codemirror/state'
import { EditorView, lineNumbers } from '@codemirror/view'
import { basicSetup } from 'codemirror'
import {gutter, GutterMarker} from '@codemirror/view';
import {Completion, CompletionContext, CompletionResult, autocompletion} from "@codemirror/autocomplete"

/**
 * Class AddNewColumnDialog is a useful method to add a new column to the table
 * or to edit formula for an existing column.
 *
 * It uses user-friendly and functional Dialog Window that allows to use
 * formulas, functions and other columns to create a new column and immediately
 * see the preview result.
 */
let coreFunctionsNames: string[] = [];
let coreFunctionsParams: string[] = [];
let packageFunctionsNames: {[key: string]: string[]} = {};
let packageFunctionsParams: {[key: string]: string[]} = {};
let packageNames: string[] = [];
let columnNames: string[] = [];
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
  maxAutoNameLength: number = 25;
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
  //inputExpression?: DG.InputBase;
  uiPreview?: HTMLDivElement;
  uiColumns?: HTMLDivElement;
  uiFunctions?: HTMLDivElement;
  uiDialog?: DG.Dialog;
  codeMirror?: EditorView;
  codeMirrorDiv?: HTMLDivElement;

  constructor(call: DG.FuncCall | null = null) {
    const table = call?.getParamValue('table');

    if (table) {
      this.call = call;
      this.sourceDf = table;
      this.edit = table.columns.names().includes(call?.getParamValue('name'));
    } else
      this.sourceDf = grok.shell.t;


    if (this.sourceDf)
      this.init();
    else
      grok.shell.error('Table not found');
  }

  /** Initializes all parameters and opens a Dialog Window. */
  async init(): Promise<void> {
    this.dialogTitle = this.edit ? this.editColumnTitle : this.addColumnTitle;

    this.uiDialog = ui.dialog({title: this.dialogTitle, helpUrl: this.helpUrl});

    this.inputName = this.initInputName();
    this.inputType = this.initInputType();
    //this.inputExpression = this.initInputExpression();
    const codeMirrorRes = this.initCodeMirror();
    this.codeMirror = codeMirrorRes.codeMirror;
    this.codeMirrorDiv = codeMirrorRes.coreMirrorDiv;

    // Not necessary, but if the Dialog knows about inputs, then it can implement extra-logic:
    this.uiDialog
      .add(this.inputName)
      .add(this.inputType)
     // .add(this.inputExpression);

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
    const allFunctionsList = DG.Func.find();
    coreFunctionsNames = [];
    coreFunctionsParams = [];
    packageFunctionsNames = {};
    packageFunctionsParams = {};
    packageNames = [];
    columnNames = [];
    for (const func of allFunctionsList) {
      const paramsStr = func.inputs.map((it) => `${it.name}:${it.propertyType}`).join(', ');
      try {
        const packageName = func.package.name;
        if (!packageFunctionsNames[packageName]) {
          packageFunctionsNames[packageName] = [];
          packageFunctionsParams[packageName] = [];
        }
        packageFunctionsNames[packageName].push(func.name);
        packageFunctionsParams[packageName].push(paramsStr);
      } catch { //in case of core functions calling func.package throws an exception
        coreFunctionsNames.push(func.name);
        coreFunctionsParams.push(paramsStr);
      }
    }
    packageNames = Object.keys(packageFunctionsNames);
    columnNames = this.sourceDf!.columns.names();
  }

  /** Prepares form inputs for Selenium tests by given them unique names. */
  prepareForSeleniumTests() {
    (this.inputName!.input as HTMLInputElement).name = 'input-Add-New-Column---Name';
    (this.inputType!.input as HTMLInputElement).name = 'input-Add-New-Column---Type';
   // (this.inputExpression!.input as HTMLInputElement).name = 'input-Add-New-Column---Expression';
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
    //this.inputExpression!.value = history[this.placeholderExpression];
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

  // /** Creates and initializes the "Expression" input field. */
  // initInputExpression(): DG.InputBase {
  //   const control = ui.textInput('', '', async () => {
  //     // The first characters of the Expression become the default name of the new column:
  //     (this.inputName!.input as HTMLInputElement).placeholder =
  //         ((!control.value || (control.value.length > this.maxAutoNameLength)) ?
  //           this.placeholderName :
  //           control.value).trim();
  //     await this.updatePreview(this.codeMirror!.state.doc.toString());
  //   });
  //   control.setTooltip(this.tooltips['expression']);

  //   const input = control.input as HTMLInputElement;
  //   input.classList.add('ui-input-addnewcolumn-expression');
  //   input.placeholder = this.placeholderExpression;
  //   if (this.call)
  //     input.value = this.call.getParamValue('expression');

  //   // Columns and functions can be drag-n-dropped into the Expression field:
  //   ui.makeDroppable(input, {
  //     acceptDrop: (dragObject) => this.typeOf(dragObject, DG.Column, DG.Func),
  //     doDrop: (dragObject, _) => this.insertIntoExpression(dragObject),
  //   });

  //   ui.tools.initFormulaAccelerators(control, this.sourceDf!);

  //   return control;
  // }

  initCodeMirror(): {codeMirror: EditorView, coreMirrorDiv: HTMLDivElement} {
    const editorDiv = ui.div('', { style: { height: '140px' } });
    let errorMsg = '';
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

    const errorMarker = new class extends GutterMarker {
      toDOM() {
        const el = document.createElement('span');
        el.innerText = "⬤";
        el.style.color = 'red';
        el.style.paddingRight = '5px';
        ui.tooltip.bind(el, errorMsg);
        return el;
      }
    };

    const errorGutter = gutter({
      lineMarker(view, line) {
        errorMsg = '';
        const regex = /\${(\w*)}/g;
        const matchesAll = [...view.state.doc.toString().matchAll(regex)];
        const unmatchedCols: string[] = [];
        if (matchesAll.length) {
          for(const match of matchesAll)
            if(!columnNames.includes(match[1]))
              unmatchedCols.push(match[1]);
        }
        if (unmatchedCols.length)
          errorMsg = unmatchedCols.length > 1 ? `Columns ${unmatchedCols.join(',')} are missing` :
            `Column ${unmatchedCols[0]} is missing`;
        
        return errorMsg ? errorMarker : null;
      },
      initialSpacer: () => errorMarker,
    });

    const autocomplete = autocompletion({
      activateOnTyping: true,
      override: [this.myCompletions],
    });

    const cm = new EditorView({
      parent: editorDiv,
      state: EditorState.create({
        doc: '',
        extensions: [
          basicSetup,
          lineNumbers(),
          errorGutter,
          EditorView.lineWrapping,
          autocomplete,
          EditorView.updateListener.of(async () => {
            const cmValue = cm.state.doc.toString();
            (this.inputName!.input as HTMLInputElement).placeholder =
              ((!cmValue || (cmValue.length > this.maxAutoNameLength)) ? this.placeholderName : cmValue).trim();
            await this.updatePreview(cmValue);
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

    //ui.tools.initFormulaAccelerators(control, this.sourceDf!); ????????????????????

    return {codeMirror: cm, coreMirrorDiv: editorDiv};
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
            ui.block([this.codeMirrorDiv!]),
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

  // /** Inserts drag-n-dropping object into the Expression input field. */
  // insertIntoExpression(x: any): void {
  //   let snippet: string = '';

  //   if (this.typeOf(x, DG.Column))
  //     snippet = `\${${x.name}}`;
  //   else if (this.typeOf(x, DG.Func))
  //     snippet = `${x.name}()`;
  //   else
  //     return;

  //   const input = this.inputExpression!.input as HTMLInputElement;
  //   const text = input.value;
  //   const cursorPos = input.selectionStart!;
  //   input.value = text.slice(0, cursorPos) + snippet + text.slice(cursorPos);

  //   this.inputExpression!.fireChanged();
  // }

    /** Inserts drag-n-dropping object into the Expression input field. */
    insertIntoCodeMirror(x: any, cm: EditorView): void {
      let snippet: string = '';
  
      if (this.typeOf(x, DG.Column))
        snippet = `\${${x.name}}`;
      else if (this.typeOf(x, DG.Func)) {
        const paramsStr = (x as DG.Func).inputs.map((it) => `${it.name}:${it.propertyType}`).join(', ');
        snippet = `${x.name}(${paramsStr})`;
      }
      else
        return;
      const value = cm.state.doc.toString();
      const cursorPos = cm.state.selection.main.head;
      const newValue = value.slice(0, cursorPos) + snippet + value.slice(cursorPos);
      cm.dispatch({changes: {
        from: 0,
        to: cm.state.doc.length,
        insert: newValue
      }})
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


  myCompletions(context: CompletionContext): CompletionResult | null {
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
          options.push({ label: name, type: "variable", apply: `${name}(${packageFunctionsParams[packName][idx]})` })
      });
      index = word!.from + colonIdx + 1;
      filter = !word.text.endsWith(':');
    } else if (word.text.includes('$')) {
      columnNames.forEach((name: string) => options.push({ label: name, type: "variable", apply: `{${name}}` }));
      index = word!.from + word.text.indexOf("$") + 1;
      filter = !word.text.endsWith('$');
    } else
      coreFunctionsNames.concat(packageNames)
        .forEach((name: string, idx: number) => options.push({ label: name, type: "variable",
          apply: idx < coreFunctionsNames.length ? `${name}(${coreFunctionsParams[idx]})`: name,
          detail: idx < coreFunctionsNames.length ? '' : 'package'}));;
    return {
      from: index,
      options: options,
      filter: filter
    }
  }
}


