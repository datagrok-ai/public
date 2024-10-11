/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { EditorSelection, EditorState, Extension, StateEffectType, Prec } from '@codemirror/state'
import { DecorationSet, EditorView, ViewUpdate, hoverTooltip, keymap } from '@codemirror/view'
import {RegExpCursor} from "@codemirror/search"
import {Completion, CompletionContext, CompletionResult, autocompletion, startCompletion, completeFromList} from "@codemirror/autocomplete"
import {StateEffect, StateField} from "@codemirror/state"
import {Decoration} from "@codemirror/view"
import { minimalSetup } from 'codemirror';
import {bracketMatching} from "@codemirror/language"
import { Subject } from 'rxjs';

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

type UpdatePreviewParams = {
  expression: string,
  changeName: boolean,
}

const VALIDATION_TYPES_MAPPING: { [key: string]: string[] } = {
  'num': ['number', 'int', 'double', 'float', 'qnum'],
  'number': ['num', 'int', 'double', 'float', 'qnum'],
  'double': ['int', 'float', 'number', 'num', 'qnum'],
  'float': ['int', 'qnum'],
  'int': ['num', 'number', 'qnum'],
  'bool': ['boolean'],
  'boolean': ['bool']
};

const FLOATING_POINT_TYPES = ['float', 'double'];
const ALLOWED_OUTPUT_TYPES = ['dynamic', DG.TYPE.DATE_TIME, DG.TYPE.QNUM];
const PACKAGES_TO_EXCLUDE = ['ApiTests', 'CvmTests'];
const TAGS_TO_EXCLUDE = ['internal'];

const COLUMN_FUNCTION_NAME = 'GetCurrentRowField';
const GET_VAR_FUNCTION_NAME = 'GetVar';
const RESERVED_FUNC_NAMES_AND_TYPES: {[key: string]: string} = {
  'GetRowIndex': DG.TYPE.INT,
}

const DEFAULT_HINT = `Type '$' to select a column or press 'Ctrl + Space' to select a function`;
const FUNC_OUTPUT_TYPE = 'output';

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
  maxPreviewRowCount: number = 20;
  newColumnBgColor: number = 0xFFFDFFE7; // The same bg-color as the bg-color of tooltips.
  colNamePattern: RegExp = /\${(.+?)(?<!\\)}|\$\[(.+?)(?<!\\)\]/g;
  colNamePatternWithoutDollar: RegExp = /{(.+?)}|\[(.+?)\]/g;
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
  currentCalculatedColName = '';

  inputName?: DG.InputBase;
  inputType?: DG.InputBase;
  uiPreview?: HTMLDivElement;
  uiColumns?: HTMLDivElement;
  uiFunctions?: HTMLDivElement;
  uiDialog?: DG.Dialog;
  codeMirror?: EditorView;
  codeMirrorDiv = ui.div('', { style: { height: '140px' } });
  errorDiv = ui.div('', 'cm-errort-div cm-hint-div');
  hintDiv = ui.div('', 'cm-hint-div');
  columnNames: string[] = [];
  columnNamesLowerCase: string[] = [];
  coreFunctionsNames: string[] = [];
  coreFunctionsParams: {[key: string]: PropInfo[]} = {};
  packageFunctionsNames: {[key: string]: string[]} = {};
  packageFunctionsParams: {[key: string]: PropInfo[]} = {};
  packageNames: string[] = [];
  packageAutocomplete = false;
  functionAutocomplete = false;
  autocompleteEnter = false;
  selectedColumn: DG.Column | null = null;
  error = '';
  mutationObserver: MutationObserver | null = null;
  mouseDownOnCm = false;
  updatePreviewEvent = new Subject<UpdatePreviewParams>();

  constructor(call: DG.FuncCall | null = null) {
    const table = call?.getParamValue('table');

    DG.debounce(this.updatePreviewEvent, 1000).subscribe(async (params: UpdatePreviewParams) => {
      await this.updatePreview(params.expression, params.changeName);
    });

    if (table) {
      this.call = call;
      this.sourceDf = table;
      this.edit = table.columns.names().includes(call?.getParamValue('name'));
    } else
      this.sourceDf = grok.shell.t;


    if (this.sourceDf) {
      this.columnNames = this.sourceDf.columns.names();
      this.columnNamesLowerCase = this.sourceDf.columns.names().map((it) => it.toLowerCase());
      this.hintDiv.append(ui.divText(DEFAULT_HINT));
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

    // Not necessary, but if the Dialog knows about inputs, then it can implement extra-logic:
    this.uiDialog
      .add(this.inputName)
      .add(this.inputType)

    this.uiDialog
      .add(await this.initUiLayout())
      .onOK(async () => {
        this.codeMirror?.destroy();
        await this.addNewColumnAction();
      })
      .onCancel(async () => {
        this.codeMirror?.destroy();
      })
      .show({resizable: true, width: 750, height: 500});

    this.uiDialog.onClose.subscribe((_) => {
      this.mutationObserver?.disconnect(); 
    });

    this.uiDialog.history(
      () => this.saveInputHistory(),
      (x) => this.loadInputHistory(x),
    );

    this.codeMirror = this.initCodeMirror();
    this.codeMirrorDiv.onkeydown = (e: KeyboardEvent) => {
      if (e.key === 'Enter' && this.autocompleteEnter) { //do not close the dialog when autocompleting using Enter button
        e.stopImmediatePropagation();
        this.autocompleteEnter = false;
      }
      if (e.key === 'Escape') //do not close the dialog if press Ecs over the codeMirror
        e.stopImmediatePropagation();
    }

    this.prepareForSeleniumTests();
    if (!this.call)
      await this.updatePreview(this.codeMirror!.state.doc.toString(), false);
    this.prepareFunctionsListForAutocomplete();
    //set initial focus on code mirror
    ui.tools.waitForElementInDom(this.codeMirrorDiv).then(() => setTimeout(() => this.codeMirror?.focus(), 50));
  }

  prepareFunctionsListForAutocomplete() {
    //filter functions with one input (multiple inputs or functions returning void are not included)
    //also filter functions returning scalar param unless it is a vector function
    const returnTypeCond = (it: DG.Func) => {
      return (DG.TYPES_SCALAR.has(it.outputs[0].propertyType) || ALLOWED_OUTPUT_TYPES.includes(it.outputs[0].propertyType))
        || it.options['vectorFunc'];
    }
    const allFunctionsList = DG.Func.find()
      .filter((it) => TAGS_TO_EXCLUDE.every((tag) => !it.hasTag(tag)) && it.outputs.length === 1 
      && returnTypeCond(it));
    for (const func of allFunctionsList) {
      const params: PropInfo[] = func.inputs.map((it) => {
        return {propName: it.name, propType: it.semType ?? it.propertyType};
      });
      //the last param in the list is return value
      params.push({propName: FUNC_OUTPUT_TYPE, propType: func.outputs[0].semType ?? func.outputs[0].propertyType});
      try {
        const packageName = func.package.name;
        if (PACKAGES_TO_EXCLUDE.includes(packageName))
          continue;
        if (!this.packageFunctionsNames[packageName]) {
          this.packageNames.push(packageName);
          this.packageFunctionsNames[packageName] = [];
        }
        this.packageFunctionsNames[packageName].push(func.name);
        this.packageFunctionsParams[`${packageName}:${func.name}`] = params;
      } catch { //in case of core functions calling func.package throws an exception
        const funcName = func.nqName.startsWith('core:') ? func.name : func.nqName;
        this.coreFunctionsNames.push(funcName);
        this.coreFunctionsParams[funcName] = params;
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
    const control = ui.input.string('', {value: ''});
    control.onInput.subscribe(async () => await this.updatePreview(this.codeMirror!.state.doc.toString(), true));
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
    const defaultChoice = `${this.autoType} (${this.defaultType})`;
    this.supportedTypes.unshift(defaultChoice); // The first item of the ChoiceBox will be "Auto".
    this.supportedTypes.push(this.plainTextType); // The last item of the ChoiceBox will be "Treat As String".

    const control = ui.input.choice('', {value: this.call ?
      this.call.getParamValue('type') : defaultChoice, items: this.supportedTypes});
    control.onInput.subscribe(async () => await this.updatePreview(this.codeMirror!.state.doc.toString(), false));
    control.setTooltip(this.tooltips['type']);

    const input = control.input as HTMLInputElement;
    input.classList.add('ui-input-addnewcolumn-type');
    input.insertBefore(ui.element('hr'), // Separator before the last item in the ChoiceBox.
      (input as any)[(input as any).length - 1]);

    return control;
  }


  initCodeMirror(): EditorView {
    this.codeMirrorDiv!.onclick = () => {
      cm.focus();
    };

    this.codeMirrorDiv!.onmousedown = () => {
      this.mouseDownOnCm = true;
    };

    this.uiDialog!.root.onmouseleave = () => {
      if(this.mouseDownOnCm) {
        this.mouseDownOnCm = false;
        cm.focus();
      }
    };

    this.uiDialog!.root.onclick = () => {
      setTimeout(() => {
        this.setCodeMirrorFocus(cm);
      }, 100);
    };

    // Columns and functions can be drag-n-dropped into the Expression field:
    ui.makeDroppable(this.codeMirrorDiv!, {
      acceptDrop: (dragObject) => this.typeOf(dragObject, DG.Column, DG.Func),
      doDrop: (dragObject, _) => {
        cm.focus();
        this.insertIntoCodeMirror(dragObject, cm);
      },
    });

    //autocompletion extension
    const autocomplete = autocompletion({
      override: [this.functionsCompletions(this.columnNames, this.packageNames, this.coreFunctionsNames,
        this.packageFunctionsNames, this.packageFunctionsParams, this.coreFunctionsParams)],
      activateOnCompletion: ({ apply }) => {
        this.autocompleteEnter = true;
        if (typeof apply === 'string' && apply.match(this.colNamePatternWithoutDollar))
          return true;
        this.packageAutocomplete = typeof apply === 'string' && apply.slice(-1) === ':';
        this.functionAutocomplete = !this.packageAutocomplete;
        return this.packageAutocomplete;
      }
    });

    //functions tooltip extension
    const wordHover = this.hoverTooltipCustom(this.packageFunctionsParams, this.coreFunctionsParams);

    //highlight column names
    const addColHighlight = StateEffect.define<{from: number, to: number}>({
      map: ({from, to}, change) => ({from: change.mapPos(from), to: change.mapPos(to)})
    });

    //highlight unmatched parentheses
    const addUnmatchedParentheses = StateEffect.define<{ from: number, to: number }>({
      map: ({ from, to }, change) => ({ from: change.mapPos(from), to: change.mapPos(to) })
    });

    //highlight text within quotes
    const addTextWithinQuotes = StateEffect.define<{ from: number, to: number }>({
      map: ({ from, to }, change) => ({ from: change.mapPos(from), to: change.mapPos(to) })
    });

    //remove all highlights
    const removeHighlight = StateEffect.define<{ from: number, to: number }>({
      map: ({ from, to }, change) => ({ from: change.mapPos(from), to: change.mapPos(to) })
    });

    const highlightMark = Decoration.mark({class: "cm-column-name"});
    const unmatchedParenthesesMark = Decoration.mark({class: "cm-unmatched-bracket"});
    const withinQuotesMark = Decoration.mark({class: "cm-within-quotes"});
    const highlightTheme = EditorView.baseTheme({
      ".cm-column-name": { 
        'color': 'var(--blue-2)',
        'font-weight': 'bold' 
       },
      ".cm-unmatched-bracket": {
        'color': 'red',
        'font-weight': 'bold' 
      },
      ".cm-within-quotes": {
        'color': '#c27706',
      },
    });

    const highlightField = StateField.define<DecorationSet>({
      create() {
        return Decoration.none;
      },
      update(underlines, tr) {
        underlines = underlines.map(tr.changes)
        for (let e of tr.effects)
          if (e.is(addColHighlight))
            underlines = underlines.update({
              add: [highlightMark.range(e.value.from, e.value.to)]
            });
          else if (e.is(addUnmatchedParentheses))
            underlines = underlines.update({
              add: [unmatchedParenthesesMark.range(e.value.from, e.value.to)]
            })
          else if (e.is(addTextWithinQuotes))
            underlines = underlines.update({
              add: [withinQuotesMark.range(e.value.from, e.value.to)]
            })
          else if (e.is(removeHighlight))
            underlines = Decoration.none;
        return underlines;
      },
      provide: f => EditorView.decorations.from(f)
    });

    const setSelection = (view: EditorView, selections: any[], stateEffect: StateEffectType<unknown>) => {
      let effects: StateEffect<unknown>[] = selections.map(({from, to}) => stateEffect.of({from, to}));
      if (!effects.length)
        return false;
    
      if (!view.state.field(highlightField, false))
        effects.push(StateEffect.appendConfig.of([highlightField, highlightTheme]));
      view.dispatch({effects});
      return true;
    }

    const cutSelection = (v: EditorView) => {
      const selectionContents = v.state.sliceDoc(v.state.selection.main.from, v.state.selection.main.to);
      navigator.clipboard.writeText(selectionContents);
      v.dispatch({ changes: { from: v.state.selection.main.from, to: v.state.selection.main.to, insert: '' } })
      return true;
    }

    const addRegexpSelection = (regexp: string, stateEffect: StateEffectType<{ from: number; to: number; }>) => {
      const cursor = new RegExpCursor(cm.state.doc, regexp);

      const selections = [];
      while (!cursor.done) {
        cursor.next();
        if (cursor.value.from !== -1 && cursor.value.to !== -1)
          selections.push({from: cursor.value.from, to: cursor.value.to});
      }
      if (selections.length)
        setSelection(cm, selections, stateEffect);
    }

    //create code mirror
    const cm = new EditorView({
      parent: this.codeMirrorDiv!,
      state: EditorState.create({
        doc: '',
        extensions: [
          EditorView.lineWrapping,
          autocomplete,
          wordHover,
          minimalSetup,
          highlightTheme,
          highlightField,
          bracketMatching({brackets : "()[]{}"}),
          keymap.of([
            {
              key: 'Shift-Delete',
              run: cutSelection
            }
          ]),
          EditorView.updateListener.of(async (e: ViewUpdate) => {

            //update hint
            ui.empty(this.hintDiv);
            const resFunc = this.getFunctionNameAtPosition(cm, cm.state.selection.main.head, -1,
              this.packageFunctionsParams, this.coreFunctionsParams);
            const fullFuncName = resFunc?.funcName;
            this.hintDiv.append(ui.divText(resFunc?.signature ?? DEFAULT_HINT));

            //return in case formula hasn't been changed
            if (!e.docChanged)
              return;

            this.setCodeMirrorFocus(cm);

            const cmValue = cm.state.doc.toString();

            //remove highlight
            setSelection(cm, [{from: 0, to: cmValue.length}], removeHighlight);

            //add column highlight
            addRegexpSelection('\\$\\{(.+?)(?<!\\\\)\\}|\\$\\[(.+?)(?<!\\\\)\\]', addColHighlight);
            
            //add text in quotes highlight
            addRegexpSelection(`".*?"|'.*?'`, addTextWithinQuotes);

            //add unmatched parentheses highlight
            const openBrackets: number[] = [];
            const closeBrackets: number[] = [];
          
            for (let i = 0; i < cmValue.length; i++) {
              if (cmValue[i] === '(') {
                openBrackets.push(i);
              }
              if (cmValue[i] === ')') {
                if (!openBrackets.length)
                  closeBrackets.push(i)
                else
                  openBrackets.pop();
              }
            }
            const unmatchedBracketsSelections: any[] = [];
            openBrackets.concat(closeBrackets).forEach((it) => {
              unmatchedBracketsSelections.push({from: it, to: it + 1});
            });

            if (unmatchedBracketsSelections.length)
              setSelection(cm, unmatchedBracketsSelections, addUnmatchedParentheses);

            (this.inputName!.input as HTMLInputElement).placeholder =
              ((!cmValue || (cmValue.length > this.maxAutoNameLength)) ? this.placeholderName : cmValue).trim();
            this.error = '';
            if (cmValue) {
              if (this.packageAutocomplete)
                setTimeout(() => {
                  startCompletion(cm);
                }, 100);
              else if (fullFuncName?.includes(':')) {
                const packAndFuncNames = fullFuncName.split(':');
                if (!this.packageNames.includes(packAndFuncNames[0]))
                  this.error = `Package ${packAndFuncNames[0]} not found`;
                else if (!packAndFuncNames[1])
                  this.error = `Start typing to see ${packAndFuncNames[0]} package functions`;
                else if (!this.packageFunctionsNames[packAndFuncNames[0]].includes(packAndFuncNames[1]))
                  this.error = `Function ${packAndFuncNames[1]} not found in ${packAndFuncNames[0]} package`;
                else
                  this.error = this.validateFormula(cmValue);
              } else {
                if (this.functionAutocomplete)
                  this.setSelection(cm.state.selection.main.head, true);
                this.error = this.validateFormula(cmValue);
              }
            }
            this.packageAutocomplete = false;
            this.functionAutocomplete = false;
            ui.empty(this.errorDiv);
            if (this.error)
              this.errorDiv.append(ui.divText(this.error, 'cm-error-div'));
            this.updatePreviewEvent.next({expression: cmValue, changeName: false});
          }),
        ],
      }),
    });

    //remove error in case autocomplete is open
    this.mutationObserver = new MutationObserver((mutationsList, observer) => {
      mutationsList.forEach((m) => {
        if (Array.from(m.removedNodes).filter((it) => (it as HTMLElement).classList.contains('cm-tooltip-autocomplete')).length) {
          ui.empty(this.errorDiv);
          this.errorDiv.append(ui.divText(this.error, 'cm-error-div'));
          return;
        }
        if (Array.from(m.addedNodes).filter((it) => (it as HTMLElement).classList.contains('cm-tooltip-autocomplete')).length)
          ui.empty(this.errorDiv);
      });
    });
    this.mutationObserver.observe(cm.dom, {attributes: true, childList: true});

    if (this.call)
      cm!.dispatch({changes: {
        from: 0,
        to: cm.state.doc.length,
        insert: this.call.getParamValue('expression')
      }});

    return cm;
  }

  setCodeMirrorFocus(cm: EditorView) {
    if (!this.inputName!.root.contains(document.activeElement)
      && !this.uiColumns!.contains(document.activeElement)
      && !this.uiFunctions!.contains(document.activeElement)
      && !this.inputType!.root.contains(document.activeElement)) {
        cm!.focus();
    }
  }

  validateFormula(formula: string): string {
    const matchesAll = formula.match(this.colNamePattern) as string[];
    const unmatchedCols: string[] = [];
    if (matchesAll?.length) {
      for (const match of matchesAll) {
        const matchCol = match.substring(2, match.length - 1);
        const unescapedMatch = grok.functions.handleOuterBracketsInColName(matchCol, false);
        if (!this.columnNamesLowerCase.includes(unescapedMatch.toLowerCase()))
          unmatchedCols.push(matchCol);
      }
    }
    if (unmatchedCols.length)
      return unmatchedCols.length > 1 ? `Columns ${unmatchedCols.join(',')} are missing` :
        `Column ${unmatchedCols[0]} is missing`;
    //check cases when only one column is entered
    const singleColumnPattern = /^\${(.+?)}$|^\$\[(.+?)\]$/;
    const found = formula.trim().match(singleColumnPattern);
    if (found)
      return '';
    //check syntax errors
    try {
      const funcCall = grok.functions.parse(formula, false);
      this.validateFuncCallTypes(funcCall);
    } catch (e: any) {
      return e.message.endsWith(': end of input expected]') ? 'Possible syntax error' : e.message;
    }
    return '';
  }


  validateFuncCallTypes(funcCall: DG.FuncCall) {
    const innerFuncCalls: string[] = [];
    const actualInputParamTypes: { [key: string]: string } = {};
    const actualInputSemTypes: { [key: string]: string } = {};

    //collect actual input parameter types
    for (const key of Object.keys(funcCall.inputs)) {
      if (funcCall.inputs[key] instanceof DG.FuncCall) {
        //do not allow functions with multiple inputs 
        if (funcCall.inputs[key].func.outputs.length > 1)
          throw new Error(`Function ${funcCall.inputs[key].func} returns multiple values`);
        //treat $COLUMN_NAME as scalar
        if (funcCall.inputs[key].func.name !== COLUMN_FUNCTION_NAME)
          innerFuncCalls.push(key);
        if (funcCall.inputs[key].func.outputs.length) {
          actualInputParamTypes[key] = funcCall.inputs[key].func.outputs[0].propertyType;
          actualInputSemTypes[key] = funcCall.inputs[key].func.outputs[0].semType;
        } else
          actualInputParamTypes[key] = 'dynamic';
      } else
        actualInputParamTypes[key] = typeof funcCall.inputs[key];
    }

    //validate types for current function
    for (const property of funcCall.func.inputs) {
      let actualInputType = actualInputParamTypes[property.name];
      let actualSemType = actualInputSemTypes[property.name];
      const input = funcCall.inputs[property.name];
      //check for variables missing in the context
      //TODO: preview with variables doesn't work
      if (input?.func && input.func.name === GET_VAR_FUNCTION_NAME) {
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
        if (input?.func && input.func.name === COLUMN_FUNCTION_NAME) {
          actualInputType = this.sourceDf!.col(input.inputs['field'])!.type;
          actualSemType = this.sourceDf!.col(input.inputs['field'])!.semType;
        }
        //handling 'row' function
        //TODO: Handle other similar functions
        for (const reservedFunc of Object.keys(RESERVED_FUNC_NAMES_AND_TYPES)) {
          if (input?.func.name === reservedFunc) {
            actualInputType = RESERVED_FUNC_NAMES_AND_TYPES[reservedFunc];
            break;
          }
        }
      }
      //dynamic allows any type
      if (property.propertyType === DG.TYPE.DYNAMIC || actualInputType === DG.TYPE.DYNAMIC)
        continue;
      //check for semType match
      if (property.semType && actualSemType && property.semType !== actualSemType)
        throw new Error(`Function ${funcCall.func.name} '${property.name}' param should be ${property.semType} type instead of ${actualSemType}`);
      //check for column and list types
      if (property.propertyType === DG.TYPE.COLUMN || actualInputType === DG.TYPE.LIST) {
        if (property.propertySubType && property.propertySubType !== actualInputType)
          throw new Error(`Function ${funcCall.func.name} '${property.name}' param should be ${property.propertySubType} type instead of ${actualInputType}`);
      } else {
        //check for type match
        if (property.propertyType !== actualInputType) {
          //check for type match in mapping
          const mappingMatch = VALIDATION_TYPES_MAPPING[property.propertyType] && VALIDATION_TYPES_MAPPING[property.propertyType].includes(actualInputType);
          if (!mappingMatch)
            throw new Error(`Function ${funcCall.func.name} '${property.name}' param should be ${property.propertyType} type instead of ${actualInputType}`);
        }
      }
    }
    //validate inner func calls recursively
    if (innerFuncCalls.length) {
      for (const param of innerFuncCalls)
        this.validateFuncCallTypes(funcCall.inputs[param]);
    }
  }

  getFunctionNameAtPosition(view: EditorView, pos: number, side: number,
    packageFunctionsParams: { [key: string]: PropInfo[] }, coreFunctionsParams: { [key: string]: PropInfo[] },
    withoutSignature?: boolean): { funcName: string, signature?: string, start: number, end: number } | null {
    let { from, to, text } = view.state.doc.lineAt(pos);
    let start = pos, end = pos;
    while (start > from && /\w|:/.test(text[start - from - 1]))
      start--;
    while (end < to && /\w|:/.test(text[end - from]))
      end++;
    if (start == pos && side < 0 || end == pos && side > 0)
      return null;
    const funcName = text.slice(start - from, end - from);
    if (!packageFunctionsParams[funcName] && !coreFunctionsParams[funcName])
      return null;
    if (withoutSignature)
      return {funcName: funcName, start: start, end: end};
    const funcParams = funcName.includes(':') ? packageFunctionsParams[funcName] : coreFunctionsParams[funcName];
    if (!funcParams)
      return null;
    const funcInputs = funcParams.filter((it) => it.propName !== FUNC_OUTPUT_TYPE);
    const funcOutputs = funcParams.filter((it) => it.propName === FUNC_OUTPUT_TYPE);
    let funcOutputType = '';
    if (funcOutputs.length)
      funcOutputType = funcOutputs[0].propType;
    return {
      signature: `${funcName}${funcInputs.length ? `(${funcInputs.map((it) => `${it.propName}:${it.propType}`).join(', ')})` : ''}: ${funcOutputType}`,
      funcName: funcName,
      start: start,
      end: end
    }
  }

  hoverTooltipCustom(packageFunctionsParams: { [key: string]: PropInfo[] },
    coreFunctionsParams: { [key: string]: PropInfo[] }): Extension {
    return hoverTooltip((view: EditorView, pos: number, side: number) => {
      const res = this.getFunctionNameAtPosition(view, pos, side, packageFunctionsParams, coreFunctionsParams);
      if (!res)
        return null;
      return {
        pos: res.start,
        end: res.end,
        above: true,
        create(view) {
          let dom = document.createElement("div");          
          dom.textContent = res.signature!;
          return {dom}
        }
      }
    });
  }

  setSelection(cursorPos: number, fromEnd?: boolean) {
    const openParenthesis = fromEnd ? this.codeMirror!.state.doc.toString().lastIndexOf('(', cursorPos) :
      this.codeMirror!.state.doc.toString().indexOf('(', cursorPos);
    if (openParenthesis === -1)
      return;
    const closeParenthesisIdx = this.codeMirror!.state.doc.toString().indexOf(')', openParenthesis);
    if (closeParenthesisIdx === -1)
      return;
    const commaIdx = this.codeMirror!.state.doc.toString().indexOf(',', openParenthesis);
    let firstParamEnd = commaIdx;
    if (commaIdx === -1 || commaIdx > closeParenthesisIdx)
      firstParamEnd = closeParenthesisIdx;
    setTimeout(() => this.codeMirror!.focus(), 100);
    this.codeMirror!.dispatch({
      selection: EditorSelection.create([
        EditorSelection.range(openParenthesis + 1, firstParamEnd),
      ])
    });
  }



  /** Creates and initializes the Preview Grid. */
  initUiPreview(): HTMLDivElement {
    // Limiting the number of rows in the Preview Grid:
    const previewRowCount = Math.min(this.sourceDf!.rowCount, this.maxPreviewRowCount);
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
    previewRoot.setAttribute('style', 'height: -webkit-fill-available !important;');

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

    if (this.widgetColumns!.getProperties().filter((it => it.name === 'grid')).length) { //added check for grid property for backward compatibility
      const columnsGrid = DG.toJs(this.widgetColumns!.props.grid) as DG.Grid;
      columnsGrid.autoSize(350, 345, undefined, undefined, true);
      columnsGrid.root.classList.add('add-new-column-columns-grid');
      ui.onSizeChanged(this.uiDialog!.root).subscribe(() => {
        const newHeight = this.uiDialog!.root.getElementsByClassName('add-new-column-columns-grid')[0].clientHeight - 5;
        columnsGrid.autoSize(350, newHeight, undefined, undefined, true);
      });
    }
    
    this.columnsDf = DG.toJs(this.widgetColumns!.props.dfColumns);
    this.columnsDf?.onCurrentRowChanged.subscribe(() => {
      if (this.columnsDf && this.columnsDf!.currentRowIdx !== -1) {
        const colName = this.columnsDf!.get('name', this.columnsDf!.currentRowIdx);
        if (this.sourceDf) {
          this.selectedColumn = this.sourceDf.col(colName);
          this.widgetFunctions!.props.sortByColType = this.selectedColumn;
        } else {
          this.selectedColumn = null;
          this.widgetFunctions!.props.sortByColType = null;
        }
      }
    })

    const control = ui.box();
    control.append(this.widgetColumns!.root);
    control.classList.add('ui-widget-addnewcolumn-columns');

    return control;
  }

  /** Creates and initializes the "Function List Widget". */
  async initUiFunctions(): Promise<HTMLDivElement> {
    this.widgetFunctions = await DG.Func.byName('FunctionsWidget').apply({
      scalarOnly: true,
      plusIconOnHover: true,
      includeVectorFuncs: true,
      ignoreTags: TAGS_TO_EXCLUDE,
      ignorePackages: PACKAGES_TO_EXCLUDE
    });
    this.widgetFunctions!.props.visibleCategories = this.visibleTags.join(',');
    this.widgetFunctions!.props.showSignature = true;
    (this.widgetFunctions as DG.FunctionsWidget)!.onActionPlusIconClicked.subscribe((e: DG.Func) => {
      if (this.codeMirror)
        this.insertIntoCodeMirror(e, this.codeMirror);
    });
    (this.widgetFunctions as DG.FunctionsWidget)!.onActionClicked.subscribe((e: DG.Func) => {
      if (this.codeMirror)
        this.insertIntoCodeMirror(e, this.codeMirror);
    });
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
            ui.block([this.codeMirrorDiv!, this.hintDiv, this.errorDiv]),
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
  async updatePreview(expression: string, changeName: boolean): Promise<void> {
    //get result column name
    const colName = this.getResultColumnName();

    //in case name was changed in nameInput, do not recalculate preview
    if (changeName) {
      this.gridPreview!.dataFrame.col(this.currentCalculatedColName)!.name = colName;
      this.gridPreview!.invalidate();
      this.currentCalculatedColName = colName;
      return;
    }

    this.currentCalculatedColName = colName;
    // Making the Column List for the Preview Grid:
    const columnIds = this.findUniqueColumnNamesInExpression(expression);
    columnIds.push(colName);

    const type = this.getSelectedType()[0];
    // Making the Preview Grid:
    const call = (DG.Func.find({name: 'AddNewColumn'})[0]).prepare({table: this.previwDf!,
      name: colName, expression: expression, type: type});
    ui.setUpdateIndicator(this.gridPreview!.root, true);
    await call.call(false, undefined, {processed: true, report: false});
    /*    await this.previwDf!.columns.addNewCalculated(
        colName,
        this.inputExpression!.value,
        ...this.getSelectedType()
    );*/
    ui.setUpdateIndicator(this.gridPreview!.root, false);
    //temporary fix to activate macromolecule cell renderer
    await grok.functions.call('Bio:detectMacromolecule', {col: this.previwDf!.col(colName)});

    this.gridPreview!.dataFrame = this.previwDf!.clone(null, columnIds);
    this.gridPreview!.col(colName)!.backColor = this.newColumnBgColor;
    this.resultColumnType = this.previwDf!.col(colName)!.type;
    this.previwDf!.columns.remove(colName);

    if (FLOATING_POINT_TYPES.includes(this.resultColumnType))
      this.gridPreview!.dataFrame.col(colName)!.tags[DG.TAGS.FORMAT] = '#.00000';

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

    if (this.typeOf(x, DG.Column)) {
      const cursorPos = cm.state.selection.main.head;
      let parenthesesPos = cursorPos;
      while (parenthesesPos > 0) {
        if (cm.state.doc.toString()[parenthesesPos] === '(')
          break;
        parenthesesPos--;
      }
      const funcName = this.getFunctionNameAtPosition(cm, parenthesesPos, -1, this.packageFunctionsParams, this.coreFunctionsParams, true)?.funcName;
      const isAggr = funcName ? Object.values(DG.AGG).includes(funcName!.toLocaleLowerCase() as DG.AGG) : false;
      const escapedColName = grok.functions.handleOuterBracketsInColName(x.name, true);
      snippet = isAggr ? `\$[${escapedColName}]` : `\${${escapedColName}}`;
    }
    else if (this.typeOf(x, DG.Func)) {
      const params = (x as DG.Func).inputs.map((it) => it.semType ?? it.propertyType);
      const colPos = this.findColumnTypeMatchingParam(x);
      if (colPos !== -1) {
        const isAggr = Object.values(DG.AGG).includes((x as DG.Func).name.toLocaleLowerCase() as DG.AGG);
        const escapedColName = grok.functions.handleOuterBracketsInColName(this.selectedColumn!.name, true);
        params[colPos] = isAggr ? `\$[${escapedColName}]` : `\${${escapedColName}}`;
      }
      const paramsStr = params.join(', ');
      const funcName = (x as DG.Func).nqName.startsWith('core:') ? (x as DG.Func).name : (x as DG.Func).nqName;
      snippet = `${funcName}(${paramsStr})`;
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
    });
    this.setSelection(cursorPos);
  }

  findColumnTypeMatchingParam(f: DG.Func): number {
    if (!this.selectedColumn)
      return -1;
    let finalPos = -1;
    let bestTypePosFound = false;
    let bestDynamicTypePosFound = false;
    for (let i = 0; i < f.inputs.length; i++) {
      const mappedTypes = VALIDATION_TYPES_MAPPING[f.inputs[i].propertyType] ?? [];
      if (this.selectedColumn.semType && f.inputs[i].semType === this.selectedColumn.semType) {
        finalPos = i;
        break;
      } else if ((this.selectedColumn.type === f.inputs[i].propertyType ||
        mappedTypes.includes(this.selectedColumn.type)) && f.inputs[i].semType == null && !bestTypePosFound) {
        bestTypePosFound = true;
        finalPos = i;
        if (!this.selectedColumn.semType)
          break;
      } else if (f.inputs[i].propertyType === 'dynamic' && !bestDynamicTypePosFound) {
        bestDynamicTypePosFound = true;
        finalPos = i;
      }
    }
    return finalPos;
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
      const name = this.getResultColumnName();
      await this.sourceDf!.columns.addNewCalculated(
        name,
          this.codeMirror!.state.doc.toString().trim(),
          ...this.getSelectedType(),
      );
      //temporary fix to activate macromolecule cell renderer
      const semType = await grok.functions.call('Bio:detectMacromolecule', {col: this.sourceDf?.col(name)});
      if (semType)
        this.sourceDf!.col(name)!.semType = semType;
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
      let word = context.matchBefore(/[\w|:|$|${|$\[|'|"]*/);
      if (!word || word?.from === word?.to && !context.explicit)
        return null;

      //check if word is inside quotes and if yes - do not show autocomplete
      const beforeWord = context.state.doc.toString().substring(0, word.to)
      const quoteSym = beforeWord.includes('\'') ? '\'' : beforeWord.includes('\"') ? '\"' : '';
      if (quoteSym) {
        const closingQuote = context.state.doc.toString().substring(word.to).includes(quoteSym);
        if (closingQuote)
          return null;
      }

      let options: Completion[] = [];
      let index = word!.from;
      let filter = true;
      if (word.text.includes(':')) {
        const colonIdx = word.text.indexOf(":");
        const packName = word.text.substring(0, word.text.indexOf(":"));
        if (packageFunctionsNames[packName])
          packageFunctionsNames[packName].forEach((name: string, idx: number) => {
            options.push({ label: name, type: "variable",
              apply: `${name}(${packageFunctionsParams[`${packName}:${name}`]
                .filter((it) => it.propName !== FUNC_OUTPUT_TYPE)
                .map((it)=> it.propType).join(',')})` })
          });
        index = word!.from + colonIdx + 1;
        filter = !word.text.endsWith(':');
      } else if (word.text.includes('$') || word.text.includes('${') || word.text.includes('$[')) {
        const openingSym = word.text.includes('$[') ? '[' : '{';
        const closingSym = openingSym === '{' ? '}' : ']';
        const openingBracketIdx = word.text.indexOf(openingSym);
        const closingBracket = context.state.doc.length > word.text.length ? context.state.doc.toString().at(word.to) === openingSym : false;
        colNames.forEach((name: string) => options.push({ label: name, type: "variable",
          apply: openingBracketIdx !== -1 ? closingBracket ? `${grok.functions.handleOuterBracketsInColName(name, true)}` : 
            `${grok.functions.handleOuterBracketsInColName(name, true)}${closingSym}` :
              closingBracket ? `${openingSym}${grok.functions.handleOuterBracketsInColName(name, true)}` : 
                `${openingSym}${grok.functions.handleOuterBracketsInColName(name, true)}${closingSym}`}));
        index = word!.from + (openingBracketIdx === -1 ? word.text.indexOf("$") + 1 : openingBracketIdx + 1);
        filter = !word.text.endsWith('$') && !word.text.endsWith(openingSym);
      } else
        coreFunctionsNames.concat(packageNames)
          .forEach((name: string, idx: number) => options.push({
            label: name, type: "variable",
            apply: idx < coreFunctionsNames.length ? `${name}(${coreFunctionsParams[name]
              .filter((it) => it.propName !== FUNC_OUTPUT_TYPE)
              .map((it)=> it.propType).join(',')})` : `${name}:`,
            detail: idx < coreFunctionsNames.length ? '' : 'package',
            section: idx < coreFunctionsNames.length ? '' : 'packages'
          }));
      return {
        from: index,
        options: options,
        filter: filter
      } as CompletionResult
    };
  }
}