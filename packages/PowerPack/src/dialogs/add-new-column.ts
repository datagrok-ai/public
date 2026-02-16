/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {EditorSelection, EditorState, Extension, StateEffectType} from '@codemirror/state';
import {DecorationSet, EditorView, ViewUpdate, hoverTooltip, keymap} from '@codemirror/view';
import {Completion, CompletionContext, CompletionResult,
  autocompletion, startCompletion} from '@codemirror/autocomplete';
import {StateEffect, StateField} from '@codemirror/state';
import {Decoration} from '@codemirror/view';
import {minimalSetup} from 'codemirror';
import {bracketMatching} from '@codemirror/language';
import {Subject} from 'rxjs';

/**
 * Class AddNewColumnDialog is a useful method to add a new column to the table
 * or to edit formula for an existing column.
 *
 * It uses user-friendly and functional Dialog Window that allows to use
 * formulas, functions and other columns to create a new column and immediately
 * see the preview result.
 */

const SYNTAX_ERROR = 'Possible syntax error';

type PropInfo = {
  propName: string,
  propType: string
}

type FuncInfo = {
  params: PropInfo[],
  isVectorFunc: boolean,
}

type UpdatePreviewParams = {
  expression: string,
  changeName: boolean
}

type ColumnNamesAndSelections = {
  quotesSelection: {from: number, to: number}[],
  columnSelections: {from: number, to: number}[],
  unmatchedBracketsSelections: {from: number, to: number}[],
  columnNames: string[],
  isSingleCol: boolean,
}

const VALIDATION_TYPES_MAPPING: { [key: string]: string[] } = {
  'num': ['number', 'int', 'double', 'float', 'qnum'],
  'number': ['num', 'int', 'double', 'float', 'qnum'],
  'double': ['int', 'float', 'number', 'num', 'qnum'],
  'float': ['int', 'qnum'],
  'int': ['num', 'number', 'qnum'],
  'bool': ['boolean'],
  'boolean': ['bool'],
  [DG.TYPE.LIST]: [DG.TYPE.OBJECT],
  [DG.TYPE.COLUMN_LIST]: [DG.TYPE.OBJECT],
  [DG.TYPE.DATA_FRAME_LIST]: [DG.TYPE.OBJECT],
  ['num_list']: [DG.TYPE.OBJECT],
};

const TYPED_LISTS: {[key: string]: string[]} = {
  [DG.TYPE.COLUMN_LIST]: [DG.TYPE.COLUMN],
  [DG.TYPE.DATA_FRAME_LIST]: [DG.TYPE.DATA_FRAME],
  [DG.TYPE.STRING_LIST]: [DG.TYPE.STRING],
  ['num_list']: [DG.TYPE.QNUM, DG.TYPE.FLOAT, DG.TYPE.INT, DG.TYPE.NUM, DG.TYPE.BIG_INT, 'number', 'float'],
};


const FLOATING_POINT_TYPES = ['float', 'double'];
const ALLOWED_OUTPUT_TYPES = ['dynamic', DG.TYPE.DATE_TIME, DG.TYPE.QNUM];
const PACKAGES_TO_EXCLUDE = ['ApiTests', 'CvmTests'];
const TAGS_TO_EXCLUDE = ['internal'];

const COLUMN_FUNCTION_NAME = 'GetCurrentRowField';
const GET_VAR_FUNCTION_NAME = 'GetVar';
const RESERVED_FUNC_NAMES_AND_TYPES: {[key: string]: string} = {
  'GetRowIndex': DG.TYPE.INT,
};

const DEFAULT_HINT = `Type '$' to select a column or press 'Ctrl + Space' to select a function`;
const FUNC_OUTPUT_TYPE = 'output';

const isNumerical = (type: string) => type == DG.TYPE.INT || type == DG.TYPE.FLOAT ||
  type == DG.TYPE.NUM || type == DG.TYPE.QNUM || type == DG.TYPE.BIG_INT;

const numericCast: {[key: string]: {[key: string]: string}} = {
  [DG.TYPE.INT]: {[DG.TYPE.FLOAT]: DG.TYPE.FLOAT, [DG.TYPE.NUM]: DG.TYPE.FLOAT,
    [DG.TYPE.QNUM]: DG.TYPE.QNUM, [DG.TYPE.BIG_INT]: ''},
  [DG.TYPE.FLOAT]: {[DG.TYPE.INT]: DG.TYPE.FLOAT, [DG.TYPE.NUM]: DG.TYPE.FLOAT,
    [DG.TYPE.QNUM]: '', [DG.TYPE.BIG_INT]: ''},
  [DG.TYPE.NUM]: {[DG.TYPE.INT]: DG.TYPE.FLOAT, [DG.TYPE.FLOAT]: DG.TYPE.FLOAT,
    [DG.TYPE.QNUM]: '', [DG.TYPE.BIG_INT]: ''},
  'number': {[DG.TYPE.INT]: DG.TYPE.FLOAT, [DG.TYPE.FLOAT]: DG.TYPE.FLOAT, [DG.TYPE.QNUM]: '', [DG.TYPE.BIG_INT]: ''},
  [DG.TYPE.QNUM]: {[DG.TYPE.INT]: DG.TYPE.QNUM, [DG.TYPE.FLOAT]: '', [DG.TYPE.NUM]: '', [DG.TYPE.BIG_INT]: ''},
  [DG.TYPE.BIG_INT]: {[DG.TYPE.INT]: '', [DG.TYPE.FLOAT]: '', [DG.TYPE.NUM]: '', [DG.TYPE.QNUM]: ''},
};

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
  widgetColumns?: DG.ColumnGrid;
  widgetFunctions?: DG.Widget;
  resultColumnType?: string;
  dialogTitle: string = '';
  call: DG.FuncCall;
  edit: boolean = false;
  currentCalculatedColName = '';

  inputName?: DG.InputBase;
  inputType?: DG.InputBase;
  uiPreview?: HTMLDivElement;
  uiColumns?: HTMLDivElement;
  uiFunctions?: HTMLDivElement;
  uiDialog?: DG.Dialog;
  codeMirror?: EditorView;
  codeMirrorDiv = ui.div('', {style: {border: 'dotted 1px var(--grey-3)'}});
  errorDiv = ui.div('', 'cm-errort-div cm-hint-div');
  hintDiv = ui.div('', 'cm-hint-div');
  columnNames: string[] = [];
  columnNamesLowerCase: string[] = [];
  coreFunctionsNames: string[] = [];
  coreFunctionsParams: {[key: string]: FuncInfo} = {};
  packageFunctionsNames: {[key: string]: string[]} = {};
  fullPackageFunctionNames: string[] = [];
  packageFunctionsParams: {[key: string]: FuncInfo} = {};
  packageNames: string[] = [];
  packageAutocomplete = false;
  functionAutocomplete = false;
  autocompleteEnter = false;
  selectedColumn: DG.Column | null = null;
  error = '';
  changedType = false;
  multipleColsOutput = false;
  mutationObserver: MutationObserver | null = null;
  mouseDownOnCm = false;
  updatePreviewEvent = new Subject<UpdatePreviewParams>();
  widget?: DG.Widget;
  applyFormulaButton?: HTMLButtonElement;
  colNameWidget = '';
  colTypeWidget = '';
  autocompleteOpened = false;

  private get isFilterFormulaEditor(): boolean {
    return this.call.aux['filterFormulaEditor'] == true;
  }

  constructor(call: DG.FuncCall, widget?: DG.Widget) {
    this.codeMirrorDiv.classList.add(this.widget ? 'add-new-column-widget-cm-div' : 'add-new-column-dialog-cm-div');
    this.call = call;
    const table = call.getParamValue('table');

    DG.debounce(this.updatePreviewEvent, 1000).subscribe(async (params: UpdatePreviewParams) => {
      await this.updatePreview(params.expression, params.changeName);
    });

    if (table) {
      this.sourceDf = table;
      this.edit = table.columns.names().includes(call?.getParamValue('name'));
    } else
      this.sourceDf = grok.shell.t;

    if (widget)
      this.widget = widget;
    else {
      this.dialogTitle = this.isFilterFormulaEditor ?
        'Edit Formula' : (this.edit ? this.editColumnTitle : this.addColumnTitle);
      this.uiDialog = ui.dialog({title: this.dialogTitle, helpUrl: this.helpUrl});
      this.uiDialog.root.classList.add('add-new-column-dialog-root');
    }

    if (this.sourceDf) {
      if (!this.sourceDf.rowCount) {
        grok.shell.error('Column can not be added to empty dataframe');
        return;
      }
      this.columnNames = this.sourceDf.columns.names();
      this.columnNamesLowerCase = this.sourceDf.columns.names().map((it) => it.toLowerCase());
      this.hintDiv.append(ui.divText(DEFAULT_HINT));
      this.init();
    } else
      grok.shell.error('Table not found');
  }

  /** Initializes all parameters and opens a Dialog Window. */
  async init(): Promise<void> {
    // preview df will be used for both widget and dialog
    const previewRowCount = Math.min(this.sourceDf!.rowCount, this.maxPreviewRowCount);
    this.previwDf = this.sourceDf!.clone(DG.BitSet.create(previewRowCount, (idx) => idx < previewRowCount));

    if (!this.widget) {
      this.inputName = this.initInputName();
      this.inputType = this.initInputType();

      // Not necessary, but if the Dialog knows about inputs, then it can implement extra-logic:
      this.uiDialog!
        .add(this.inputName)
        .add(this.inputType);

      if (this.isFilterFormulaEditor) {
        // if the editor is created for editing viewer filter formula, we do not need name and type inputs, and they must be fixed
        this.inputName!.root.style.display = 'none';
        this.inputType!.root.style.display = 'none';
        this.inputName!.value = 'Filter Formula'; // does not matter, not used anywhere
        this.inputType!.value = DG.COLUMN_TYPE.BOOL; // needs to be bool only to assist with validation
      }

      this.uiDialog!
        .add(await this.initUiLayout())
        .onOK(async () => {
          await this.addNewColumnAction();
        })
        .show({resizable: true, width: 750, height: 500});

      this.uiDialog!.onClose.subscribe((_) => {
        this.codeMirror?.destroy();
        this.mutationObserver?.disconnect();
      });

      this.uiDialog!.history(
        () => this.saveInputHistory(),
        (x) => this.loadInputHistory(x),
      );
      this.prepareForSeleniumTests();
    } else {
      const editButton = ui.button('Edit in dialog', () => {
        this.call.setParamValue('expression', this.codeMirror!.state.doc.toString());
        this.applyFormulaButton!.disabled = true;
        editButton.disabled = true;
        const editDlg = new AddNewColumnDialog(this.call);
        editDlg.uiDialog?.onCancel(() => {
          this.applyFormulaButton!.disabled = false;
          editButton.disabled = false;
        });
      }, 'Edit with preview');
      this.applyFormulaButton = ui.button('Apply', async () => {
        this.applyFormulaButton!.disabled = true;
        await this.addNewColumnAction();
        this.applyFormulaButton!.disabled = false;
      }, 'Apply to the column');
      this.widget.root.append(ui.divV([
        this.codeMirrorDiv!,
        this.errorDiv,
        ui.divH([editButton, this.applyFormulaButton], {style: {justifyContent: 'right'}}),
      ], {style: {maxWidth: '400px'}}));
      this.widget.detach = () => {
        this.mutationObserver?.disconnect();
        this.codeMirror?.destroy();
      };
    }

    this.codeMirror = this.initCodeMirror();
    this.codeMirrorDiv.addEventListener('keydown', (e: KeyboardEvent) => {
      //do not close the dialog when autocompleting using Enter button
      if (e.code === 'Enter' && this.autocompleteEnter) {
        e.stopPropagation();
        this.autocompleteEnter = false;
      } else if (e.key === 'Escape') //do not close the dialog if press Ecs over the codeMirror
        e.stopPropagation();
      else if (e.code === 'KeyA' && e.ctrlKey) {
        e.stopPropagation();
        this.codeMirror?.dispatch({
          selection: {
            anchor: 0,
            head: this.codeMirror?.state.doc.length,
          },
        });
      }
    });
    //set initial focus on code mirror
    ui.tools.waitForElementInDom(this.codeMirrorDiv).then(() => setTimeout(() => this.codeMirror?.focus(), 50));

    if (!this.call.getParamValue('expression')) {
      const columnsAndSelections = this.getColumnNamesAndSelections(this.codeMirror!.state.doc.toString());
      this.error = this.validateFormula(this.codeMirror!.state.doc.toString(), columnsAndSelections.columnNames,
        columnsAndSelections.isSingleCol);
      await this.updatePreview(this.codeMirror!.state.doc.toString(), false);
    }
    this.prepareFunctionsListForAutocomplete();
  }

  prepareFunctionsListForAutocomplete() {
    //filter functions with one output (multiple outputs or functions returning void are not included)
    //also filter functions returning scalar param unless it is a vector function
    const returnTypeCond = (it: DG.Func) => {
      return (DG.TYPES_SCALAR.has(it.outputs[0].propertyType) ||
       ALLOWED_OUTPUT_TYPES.includes(it.outputs[0].propertyType)) || it.options['vectorFunc'];
    };
    const allFunctionsList = DG.Func.find()
      .filter((it) => TAGS_TO_EXCLUDE.every((tag) => !it.hasTag(tag)) && it.outputs.length === 1 &&
      returnTypeCond(it));
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
        this.fullPackageFunctionNames.push(`${packageName}:${func.name}`);
        this.packageFunctionsParams[`${packageName}:${func.name}`] =
          {params: params, isVectorFunc: func.options['vectorFunc']};
      } catch { //in case of core functions calling func.package throws an exception
        const funcName = func.nqName.startsWith('core:') ? func.name : func.nqName;
        this.coreFunctionsNames.push(funcName);
        this.coreFunctionsParams[funcName] = {params: params, isVectorFunc: func.options['vectorFunc']};
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
      insert: history[this.placeholderExpression],
    }});
  }

  /** Creates and initializes the "Column Name" input field. */
  initInputName(): DG.InputBase {
    const control = ui.input.string('', {value: ''});
    control.onInput.subscribe(async () => {
      await this.updatePreview(this.codeMirror!.state.doc.toString(), true);
    });
    control.setTooltip(this.tooltips['name']);

    const input = control.input as HTMLInputElement;
    input.classList.add('ui-input-addnewcolumn-name');
    input.placeholder = this.placeholderName;
    if (this.call.getParamValue('name'))
      input.value = this.call.getParamValue('name');

    return control;
  }

  /** Creates and initializes the "Column Type" input field. */
  initInputType(): DG.InputBase {
    const defaultChoice = `${this.autoType} (${this.defaultType})`;
    this.supportedTypes.unshift(defaultChoice); // The first item of the ChoiceBox will be "Auto".
    this.supportedTypes.push(this.plainTextType); // The last item of the ChoiceBox will be "Treat As String".

    const control = ui.input.choice('', {value: this.call.getParamValue('table') && this.call.getParamValue('type') ?
      this.call.getParamValue('type') : defaultChoice, items: this.supportedTypes});
    control.onInput.subscribe(async () => {
      this.changedType = true;
      await this.updatePreview(this.codeMirror!.state.doc.toString(), false);
    });
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

    if (!this.widget) {
      this.uiDialog!.root.onmouseleave = () => {
        if (this.mouseDownOnCm) {
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
    }

    //autocompletion extension
    const autocomplete = autocompletion({
      override: [this.functionsCompletions(this.columnNames, this.packageNames, this.fullPackageFunctionNames,
        this.coreFunctionsNames, this.packageFunctionsNames, this.packageFunctionsParams, this.coreFunctionsParams)],
      activateOnCompletion: ({apply}) => {
        this.autocompleteEnter = true;
        //check for column autocompletion
        if (typeof apply === 'string' && (apply.startsWith('{') ||
          apply.endsWith('}') || apply.startsWith('[') || apply.endsWith(']')))
          return true;
        //in case this is not a column name - we have either package or function
        this.packageAutocomplete = typeof apply === 'string' && apply.slice(-1) === ':';
        this.functionAutocomplete = !this.packageAutocomplete;
        return this.packageAutocomplete;
      },
    });

    //functions tooltip extension
    const wordHover = this.hoverTooltipCustom(this.packageFunctionsParams, this.coreFunctionsParams);

    //highlight column names
    const addColHighlight = StateEffect.define<{from: number, to: number}>({
      map: ({from, to}, change) => ({from: change.mapPos(from), to: change.mapPos(to)}),
    });

    //highlight unmatched parentheses
    const addUnmatchedParentheses = StateEffect.define<{ from: number, to: number }>({
      map: ({from, to}, change) => ({from: change.mapPos(from), to: change.mapPos(to)}),
    });

    //highlight text within quotes
    const addTextWithinQuotes = StateEffect.define<{ from: number, to: number }>({
      map: ({from, to}, change) => ({from: change.mapPos(from), to: change.mapPos(to)}),
    });

    //remove all highlights
    const removeHighlight = StateEffect.define<{ from: number, to: number }>({
      map: ({from, to}, change) => ({from: change.mapPos(from), to: change.mapPos(to)}),
    });

    const highlightMark = Decoration.mark({class: 'cm-column-name'});
    const unmatchedParenthesesMark = Decoration.mark({class: 'cm-unmatched-bracket'});
    const withinQuotesMark = Decoration.mark({class: 'cm-within-quotes'});
    const highlightTheme = EditorView.baseTheme({
      '.cm-column-name': {
        'color': 'var(--blue-2)',
        'font-weight': 'bold',
      },
      '.cm-unmatched-bracket': {
        'color': 'red',
        'font-weight': 'bold',
      },
      '.cm-within-quotes': {
        'color': '#c27706',
      },
    });

    const highlightField = StateField.define<DecorationSet>({
      create() {
        return Decoration.none;
      },
      update(underlines, tr) {
        underlines = underlines.map(tr.changes);
        for (const e of tr.effects) {
          if (e.is(addColHighlight)) {
            underlines = underlines.update({
              add: [highlightMark.range(e.value.from, e.value.to)],
            });
          } else if (e.is(addUnmatchedParentheses)) {
            underlines = underlines.update({
              add: [unmatchedParenthesesMark.range(e.value.from, e.value.to)],
            });
          } else if (e.is(addTextWithinQuotes)) {
            underlines = underlines.update({
              add: [withinQuotesMark.range(e.value.from, e.value.to)],
            });
          } else if (e.is(removeHighlight))
            underlines = Decoration.none;
        }
        return underlines;
      },
      provide: (f) => EditorView.decorations.from(f),
    });

    const setSelection = (view: EditorView, selections: any[], stateEffect: StateEffectType<unknown>) => {
      const effects: StateEffect<unknown>[] = selections.map(({from, to}) => stateEffect.of({from, to}));
      if (!effects.length)
        return false;

      if (!view.state.field(highlightField, false))
        effects.push(StateEffect.appendConfig.of([highlightField, highlightTheme]));
      view.dispatch({effects});
      return true;
    };

    const cutSelection = (v: EditorView) => {
      const selectionContents = v.state.sliceDoc(v.state.selection.main.from, v.state.selection.main.to);
      navigator.clipboard.writeText(selectionContents);
      v.dispatch({changes: {from: v.state.selection.main.from, to: v.state.selection.main.to, insert: ''}});
      return true;
    };

    // const addRegexpSelection = (regexp: string, stateEffect: StateEffectType<{ from: number; to: number; }>) => {
    //   const cursor = new RegExpCursor(cm.state.doc, regexp);

    //   const selections = [];
    //   while (!cursor.done) {
    //     cursor.next();
    //     if (cursor.value.from !== -1 && cursor.value.to !== -1)
    //       selections.push({from: cursor.value.from, to: cursor.value.to});
    //   }
    //   if (selections.length)
    //     setSelection(cm, selections, stateEffect);
    // };

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
          bracketMatching({brackets: '()[]{}'}),
          keymap.of([
            {
              key: 'Shift-Delete',
              run: cutSelection,
            },
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

            const columnsAndSelections = this.getColumnNamesAndSelections(cmValue);

            //add column highlight
            setSelection(cm, columnsAndSelections.columnSelections, addColHighlight);

            //add text in quotes addTextWithinQuotes
            setSelection(cm, columnsAndSelections.quotesSelection, addTextWithinQuotes);

            setSelection(cm, columnsAndSelections.unmatchedBracketsSelections, addUnmatchedParentheses);

            if (!this.widget) {
              (this.inputName!.input as HTMLInputElement).placeholder =
                ((!cmValue || (cmValue.length > this.maxAutoNameLength)) ? this.placeholderName : cmValue).trim();
            }
            this.error = '';
            if (cmValue) {
              if (this.packageAutocomplete) {
                setTimeout(() => {
                  startCompletion(cm);
                }, 100);
              } else if (fullFuncName?.includes(':')) {
                const packAndFuncNames = fullFuncName.split(':');
                const packName: string | undefined =
                  getKeyCaseInsensitive(this.packageFunctionsNames, packAndFuncNames[0]);
                if (!packName || !this.packageNames.includes(packName))
                  this.error = `Package ${packAndFuncNames[0]} not found`;
                else if (packAndFuncNames[1] && packName &&
                    !this.packageFunctionsNames[packName]
                      .some((it) => it.toLowerCase() === packAndFuncNames[1].toLowerCase()))
                  this.error = `Function ${packAndFuncNames[1]} not found in ${packAndFuncNames[0]} package`;
                else {
                  this.error = this.validateFormula(cmValue, columnsAndSelections.columnNames,
                    columnsAndSelections.isSingleCol);
                }
              } else {
                if (this.functionAutocomplete)
                  this.setSelection(cm.state.selection.main.head, true);
                this.error = this.validateFormula(cmValue, columnsAndSelections.columnNames,
                  columnsAndSelections.isSingleCol);
              }
            }
            this.packageAutocomplete = false;
            this.functionAutocomplete = false;
            if (!this.autocompleteOpened)
              this.updateError();
            //in case of syntax error we try to run expression to save string interpolation functionality
            this.updatePreviewEvent
              .next({expression: cmValue, changeName: false});
          }),
        ],
      }),
    });

    //remove error in case autocomplete is open
    this.mutationObserver = new MutationObserver((mutationsList, observer) => {
      mutationsList.forEach((m) => {
        if (Array.from(m.removedNodes)
          .filter((it) => (it as HTMLElement).classList.contains('cm-tooltip-autocomplete')).length) {
          this.autocompleteOpened = false;
          ui.empty(this.errorDiv);
          this.errorDiv.append(ui.divText(this.error, 'cm-error-div'));
          return;
        }
        if (Array.from(m.addedNodes)
          .filter((it) => (it as HTMLElement).classList.contains('cm-tooltip-autocomplete')).length) {
          this.autocompleteOpened = true;
          ui.empty(this.errorDiv);
        }
      });
    });
    this.mutationObserver.observe(cm.dom, {attributes: true, childList: true});

    if (this.call.getParamValue('expression')) {
      cm!.dispatch({
        changes: {
          from: 0,
          to: cm.state.doc.length,
          insert: this.call.getParamValue('expression'),
        },
      });
    }

    return cm;
  }

  updateError() {
    ui.empty(this.errorDiv);
    if (this.error)
      this.errorDiv.append(ui.divText(this.error, 'cm-error-div'));
    const buttonToDisable = this.widget ? this.applyFormulaButton : this.uiDialog!.getButton('OK');
    buttonToDisable!.disabled = !!this.error && this.error !== SYNTAX_ERROR;
  }

  setCodeMirrorFocus(cm: EditorView) {
    if (!this.widget) {
      if (!this.inputName!.root.contains(document.activeElement) &&
        !this.uiColumns!.contains(document.activeElement) &&
        !this.uiFunctions!.contains(document.activeElement) &&
        !this.inputType!.root.contains(document.activeElement))
        cm!.focus();
    }
  }

  getIntervalsWithinAndOutsideQuotes(formula: string):
    {quotesSelection: {from: number, to: number}[], intervalsWithoutQuotes: [number, number][]} {
    const re = /".*?"|'.*?'/gm;
    let match = null;
    const quotesSelection: {from: number, to: number}[] = [];
    const intervalsWithoutQuotes: [number, number][] = [];
    let counter = 0;
    while ((match = re.exec(formula)) != null) {
      if (!counter && match.index > 0)
        intervalsWithoutQuotes.push([0, match.index]);

      if (counter)
        intervalsWithoutQuotes.push([quotesSelection[counter - 1].to, match.index]);

      quotesSelection.push({from: match.index, to: match.index + match[0].length});
      counter++;
    }
    if (counter) {
      if (quotesSelection[counter - 1].to < formula.length)
        intervalsWithoutQuotes.push([quotesSelection[counter - 1].to, formula.length]);
    } else
      intervalsWithoutQuotes.push([0, formula.length]);


    return {quotesSelection, intervalsWithoutQuotes};
  }

  getColumnNamesAndSelections(formula: string): ColumnNamesAndSelections {
    //first check for parts in quotes and collect indexes outside quotes to check them for column names
    const leadingSpaces = formula.length - (DG._isDartium() ? formula.trimLeft().length : formula.trimStart().length);
    const closingSpaces = formula.length - (DG._isDartium() ? formula.trimRight().length : formula.trimEnd().length);
    let isSingleCol = false;
    const columnSelections: {from: number, to: number}[] = [];
    const columnNames: string[] = [];

    const {quotesSelection, intervalsWithoutQuotes} = this.getIntervalsWithinAndOutsideQuotes(formula);

    const openBrackets: number[] = [];
    const closeBrackets: number[] = [];

    const isOpeningBracket = (i: number): string => {
      let bracket = '';
      if (i > 0 && formula[i - 1] === '$')
        bracket = formula[i] === '{' ? '{' : formula[i] === '[' ? '[' : '';
      return bracket;
    };

    const isClosingBracket = (i: number): boolean => {
      return i > 0 && formula[i - 1] !== '\\' && formula[i] === closingBracket;
    };
    const getClosingBracketSym = (sym: string) => {
      return sym === '{' ? '}' : sym === '[' ? ']' : '';
    };

    let openingBracket = '';
    let openingBracketIdx: number | null = null;
    let closingBracket = '';
    for (let i = 0; i < intervalsWithoutQuotes.length; i++) {
      for (let j = intervalsWithoutQuotes[i][0]; j < intervalsWithoutQuotes[i][1]; j++) {
        if (formula[j] === '(') {
          openBrackets.push(j);
          continue;
        }
        if (formula[j] === ')') {
          if (!openBrackets.length)
            closeBrackets.push(j);
          else
            openBrackets.pop();
          continue;
        }

        const bracket = isOpeningBracket(j);
        if (!openingBracket && bracket) {
          openingBracket = bracket;
          openingBracketIdx = j - 1;
          closingBracket = getClosingBracketSym(bracket);
          continue;
        }
        if (openingBracket && isClosingBracket(j)) {
          columnSelections.push({from: openingBracketIdx!, to: j + 1});
          columnNames.push(formula.substring(openingBracketIdx! + 2, j));
          if (openingBracketIdx === 0 + leadingSpaces && j === formula.length - 1 - closingSpaces)
            isSingleCol = true;
          openingBracket = '';
          openingBracketIdx = null;
          closingBracket = '';
        }
      }
    }

    const unmatchedBracketsSelections: {from: number, to: number}[] = [];
    openBrackets.concat(closeBrackets).forEach((it) => {
      unmatchedBracketsSelections.push({from: it, to: it + 1});
    });

    return {quotesSelection, columnSelections, columnNames, isSingleCol, unmatchedBracketsSelections};
  }

  validateFormula(formula: string, columnNames: string[], isSingleCol: boolean): string {
    const unmatchedCols: string[] = [];
    for (const colName of columnNames) {
      const unescapedMatch = grok.functions.handleOuterBracketsInColName(colName, false);
      if (!this.columnNamesLowerCase.includes(unescapedMatch.toLowerCase()) && !unmatchedCols.includes(colName))
        unmatchedCols.push(colName);
    }
    if (unmatchedCols.length) {
      return unmatchedCols.length > 1 ? `Columns ${unmatchedCols.join(',')} are missing` :
        `Column ${unmatchedCols[0]} is missing`;
    }
    //check for single column formula
    if (isSingleCol)
      return '';
    //check syntax errors
    try {
      const funcCall = grok.functions.parse(formula, false);
      if (funcCall.func.outputs.length && (funcCall.func.outputs[0] as DG.Property).type === DG.TYPE.DATA_FRAME)
        this.multipleColsOutput = true;
      const error = this.validateFuncCallTypes(funcCall);
      return this.getErrorMessage(error);
    } catch (e: any) {
      this.multipleColsOutput = false;
      return this.getErrorMessage(e.message ?? e);
    }
  }

  getErrorMessage(message: string) {
    return message?.endsWith(': end of input expected]') ? 'Possible syntax error' : message;
  }

  validateFuncCallTypes(funcCall: DG.FuncCall): string {
    const innerFuncCalls: string[] = [];
    const actualInputParamTypes: { [key: string]: string } = {};
    const actualInputSemTypes: { [key: string]: string } = {};

    if (funcCall.func.name.toLowerCase() === 'if') {
      try {
        this.getIfFuncOutputParam(funcCall);
      } catch (e: any) {
        return e.message ?? e;
      }
    }
    //collect actual input parameter types
    for (const key of Object.keys(funcCall.inputs)) {
      const value = funcCall.inputs[key];
      if (value instanceof DG.FuncCall) {
        //do not allow functions with multiple outputs
        if (value.func.outputs.length > 1)
          return `Function ${value.func} returns multiple values`;
        //treat $COLUMN_NAME as scalar
        if (value.func.name !== COLUMN_FUNCTION_NAME)
          innerFuncCalls.push(key);
        if (value.func.outputs.length) {
          actualInputParamTypes[key] = value.func.outputs[0].propertyType;
          actualInputSemTypes[key] = value.func.outputs[0].semType;
        } else
          actualInputParamTypes[key] = 'dynamic';
      } else {
        actualInputParamTypes[key] = Array.isArray(value) && value.length !== 0 ?
          typeof value[0] : value == null ? 'undefined' : typeof value; // temp for debug
      }
    }

    const mappingMatch = (type: string, actualType: string) => VALIDATION_TYPES_MAPPING[type] &&
      VALIDATION_TYPES_MAPPING[type].includes(actualType);

    //validate types for current function
    for (const property of funcCall.func.inputs) {
      //skip validation of dataframe parameter for vector functions
      if (funcCall.func.options['vectorFunc'] === 'true' && property.propertyType === DG.TYPE.DATA_FRAME)
        continue;
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
          return e.message ?? e;
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
      if (property.propertyType === DG.TYPE.DYNAMIC || (actualInputType === DG.TYPE.DYNAMIC &&
        property.propertyType !== DG.TYPE.COLUMN && property.propertyType !== DG.TYPE.LIST))
        continue;
      //check for optional parameter
      if (property.nullable && actualInputType === 'undefined')
        continue;
      //check for semType match
      if (property.semType && actualSemType && property.semType !== actualSemType)
        // eslint-disable-next-line max-len
        return `Function ${funcCall.func.name} '${property.name}' param should be ${property.semType} type instead of ${actualSemType}`;
      //check column type
      if (property.propertyType === DG.TYPE.COLUMN) {
        if (funcCall.inputs[property.name].func?.name !== COLUMN_FUNCTION_NAME)
          return `Function ${funcCall.func.name} '${property.name}' param should be column type`;
        if (property.propertySubType && property.propertySubType !== actualInputType &&
          !mappingMatch(property.propertySubType, actualInputType))
          // eslint-disable-next-line max-len
          return `Function ${funcCall.func.name} '${property.name}' param should be ${property.propertySubType} column`;
        //check list type
      } else if (property.propertyType === DG.TYPE.LIST) {
        if (!Array.isArray(funcCall.inputs[property.name]))
          return `Function ${funcCall.func.name} '${property.name}' param should be array type`;
        //check type of array if array is not empty
        if (property.propertySubType && property.propertySubType !== actualInputType &&
          !mappingMatch(property.propertySubType, actualInputType) && funcCall.inputs[property.name].length > 0)
          // eslint-disable-next-line max-len
          return `Function ${funcCall.func.name} '${property.name}' param should be array of ${property.propertySubType}`;
        //check for typed lists
      } else if (TYPED_LISTS[property.propertyType]) {
        if (!TYPED_LISTS[property.propertyType].includes(actualInputType))
          // eslint-disable-next-line max-len
          return `Function ${funcCall.func.name} '${property.name}' param type should be ${TYPED_LISTS[property.propertyType].join(' or ')} instead of ${actualInputType}`;
      } else {
        //check for type match
        if (property.propertyType !== actualInputType) {
          //check for type match in mapping
          if (!mappingMatch(property.propertyType, actualInputType))
            // eslint-disable-next-line max-len
            return `Function ${funcCall.func.name} '${property.name}' param should be ${property.propertyType} type instead of ${actualInputType}`;
        }
      }
    }
    //validate inner func calls recursively
    if (innerFuncCalls.length) {
      for (const param of innerFuncCalls) {
        const error = this.validateFuncCallTypes(funcCall.inputs[param]);
        if (error)
          return error;
      }
    }
    return '';
  }

  getFunctionNameAtPosition(view: EditorView, pos: number, side: number,
    packageFunctionsParams: { [key: string]: FuncInfo }, coreFunctionsParams: { [key: string]: FuncInfo },
    withoutSignature?: boolean): { funcName: string, signature?: string, start: number, end: number } | null {
    const {from, to, text} = view.state.doc.lineAt(pos);
    let start = pos; let end = pos;
    while (start > from && /\w|:/.test(text[start - from - 1]))
      start--;
    while (end < to && /\w|:/.test(text[end - from]))
      end++;
    if (start == pos && side < 0 || end == pos && side > 0)
      return null;
    const funcName = text.slice(start - from, end - from);
    if (!packageFunctionsParams[funcName] && !coreFunctionsParams[funcName])
      return {funcName: funcName, start: start, end: end};
    if (withoutSignature)
      return {funcName: funcName, start: start, end: end};
    const funcParams = funcName.includes(':') ? packageFunctionsParams[funcName] : coreFunctionsParams[funcName];
    if (!funcParams)
      return {funcName: funcName, start: start, end: end};
    const funcInputs = funcParams.params.filter((it) => it.propName !== FUNC_OUTPUT_TYPE);
    const funcOutputs = funcParams.params.filter((it) => it.propName === FUNC_OUTPUT_TYPE);
    let funcOutputType = '';
    if (funcOutputs.length)
      funcOutputType = funcOutputs[0].propType;
    return {
      signature: `${funcName}${funcInputs.length ?
        `(${funcInputs.map((it) => `${it.propName}:${it.propType}`).join(', ')})` : ''}: ${funcOutputType}`,
      funcName: funcName,
      start: start,
      end: end,
    };
  }

  hoverTooltipCustom(packageFunctionsParams: { [key: string]: FuncInfo },
    coreFunctionsParams: { [key: string]: FuncInfo }): Extension {
    return hoverTooltip((view: EditorView, pos: number, side: number) => {
      const res = this.getFunctionNameAtPosition(view, pos, side, packageFunctionsParams, coreFunctionsParams);
      if (!res || !res.signature)
        return null;
      return {
        pos: res.start,
        end: res.end,
        above: true,
        create(view) {
          const dom = document.createElement('div');
          dom.textContent = res.signature!;
          return {dom};
        },
      };
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
      ]),
    });
  }


  /** Creates and initializes the Preview Grid. */
  initUiPreview(): HTMLDivElement {
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
    this.widgetColumns = DG.ColumnGrid.popup(this.sourceDf!, {widgetMode: true});

    if (this.widgetColumns!.grid) { //added check for grid property for backward compatibility
      const columnsGrid = this.widgetColumns?.grid!;
      columnsGrid.autoSize(350, 345, undefined, undefined, true);
      columnsGrid.root.classList.add('add-new-column-columns-grid');
      ui.onSizeChanged(this.uiDialog!.root).subscribe(() => {
        const gridEl = this.uiDialog!.root.getElementsByClassName('add-new-column-columns-grid');
        if (gridEl.length) {
          const newHeight = gridEl[0].clientHeight - 5;
          columnsGrid.autoSize(350, newHeight, undefined, undefined, true);
        }
      });
    }

    this.columnsDf = this.widgetColumns!.dfColumns;
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
    });

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
      ignorePackages: PACKAGES_TO_EXCLUDE,
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
    let colName = this.widget ? this.call.getParamValue('name') : this.getResultColumnName().colName;

    if (this.error && this.error !== SYNTAX_ERROR && !this.changedType) {
      //clearing preview df in case of error (use only within dialog)
      if (!this.widget) {
        const rowCount = this.gridPreview!.dataFrame.rowCount;
        this.gridPreview!.dataFrame = DG.DataFrame.fromColumns([DG.Column.string(colName, rowCount)]);
        this.gridPreview!.col(colName)!.backColor = this.newColumnBgColor;
      }
      return;
    }

    this.changedType = false;

    //in case name was changed in nameInput, do not recalculate preview
    if (changeName) {
      if (!this.error)
        ui.empty(this.errorDiv);
      //check if column with the same name already exists
      if (this.sourceDf?.columns.names().some((name) => name.toLowerCase() === colName.toLowerCase()) && !this.error &&
        !(this.edit && this.call.getParamValue('name')?.toLowerCase() === colName.toLowerCase())) {
        const unusedColName = this.getResultColumnName().unusedName;
        this.errorDiv.append(ui.divText(
          `Column named ${colName} already exists. Column will be named ${unusedColName}`, 'cm-warning-div'));
        colName = unusedColName;
      }
      const gridPreviewCol = this.gridPreview!.dataFrame.col(this.currentCalculatedColName);
      if (gridPreviewCol) {
        gridPreviewCol.name = colName;
        this.gridPreview!.invalidate();
      }
      this.currentCalculatedColName = colName;
      return;
    }

    this.currentCalculatedColName = colName;
    // Making the Column List for the Preview dataframe:
    const columnIds = this.findUniqueColumnNamesInExpression(expression);

    const type = this.widget ? this.call.getParamValue('type') : this.getSelectedType()[0];

    // Looking for non-empty rows in columns used in formula
    this.findNonEmptyRowsForPreview(columnIds);
    const potentialColIds: string[] = [];

    //set update indicator only in case we are within dialog
    if (!this.widget)
      ui.setUpdateIndicator(this.gridPreview!.root, true);

    await this.getPreviewResults(colName, type, expression, potentialColIds);

    //do not validate column type in case function returns multiple columns
    if (this.multipleColsOutput)
      this.error = '';

    this.updateError();

    //do not need to create preview grid in case of widget, so return
    if (this.widget)
      return;

    ui.setUpdateIndicator(this.gridPreview!.root, false);

    if (potentialColIds.length === 0)
      potentialColIds[0] = colName;
    columnIds.push(...potentialColIds);
    this.gridPreview!.dataFrame = this.previwDf!.clone(null, columnIds);
    for (const colName of potentialColIds)
      this.gridPreview!.col(colName)!.backColor = this.newColumnBgColor;
    this.resultColumnType = this.previwDf!.col(potentialColIds[0])!.type;

    for (const colName of potentialColIds)
      this.previwDf!.columns.remove(colName);

    //setting format to preview columns
    for (const colName of this.gridPreview!.dataFrame.columns.names()) {
      if (FLOATING_POINT_TYPES.includes(this.gridPreview!.dataFrame.col(colName)!.type))
        this.gridPreview!.dataFrame.col(colName)!.tags[DG.TAGS.FORMAT] = '#.00000';
    }

    this.setAutoType(); // Adding (or removing) the column auto-type caption to "Auto" item in the ChoiceBox.
  }

  async getPreviewResults(colName: string, colType: string, expression: string, potentialColIds: string[]):
    Promise<void> {
    const call = (DG.Func.find({name: 'AddNewColumn'})[0]).prepare({table: this.previwDf!,
      name: colName, expression: expression, type: 'auto'});
    const sub = this.previwDf!.onColumnsAdded.subscribe((args: DG.ColumnsArgs) => {
      potentialColIds[potentialColIds.length] = args.columns[0].name;
      const mappedTypes = VALIDATION_TYPES_MAPPING[colType] ?? [];
      this.error = colType !== 'auto' &&
       args.columns[0].type !== colType && !mappedTypes.includes(args.columns[0].type) ?
        `Result column type (${args.columns[0].type}) doesn't match with current column type (${colType}).
          Change column type ${this.widget ? 'using \'Edit in dialog\'' : ''} or modify formula.` : '';
    });
    await call.call(false, undefined, {processed: true, report: false});
    /*    await this.previwDf!.columns.addNewCalculated(
        colName,
        this.inputExpression!.value,
        ...this.getSelectedType()
    );*/
    sub.unsubscribe();
  }

  findNonEmptyRowsForPreview(columnIds: string[]) {
    let nonEmptyIdx = 0;
    if (columnIds.length) {
      const rowCount = this.sourceDf!.rowCount;
      for (let i = 0; i < rowCount; i++) {
        if (columnIds.every((colName) => !this.sourceDf!.col(colName)!.isNone(i))) {
          nonEmptyIdx = i;
          break;
        }
      }
      const minRowIdx = nonEmptyIdx + this.previwDf!.rowCount <= rowCount ?
        nonEmptyIdx : rowCount - this.previwDf!.rowCount;
      const maxRowIdx = minRowIdx + this.previwDf!.rowCount;
      this.previwDf = this.sourceDf!.clone(DG.BitSet.create(maxRowIdx,
        (idx) => idx >= minRowIdx && idx < maxRowIdx));
    }
  }

  /** Finds all unique column names used in the Expression input field. */
  findUniqueColumnNamesInExpression(expression: string): string[] {
    const columnsAndSelections = this.getColumnNamesAndSelections(expression);

    if (!columnsAndSelections.columnNames.length)
      return [];

    const lcSourceColumnNames = (this.sourceDf!.columns.names() as string[]).map((name) => name.toLowerCase());

    return columnsAndSelections.columnNames
      .filter((v, i, a) => a.indexOf(v) === i) // Selecting only unique names.
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
      const funcName = this.getFunctionNameAtPosition(cm, parenthesesPos, -1,
        this.packageFunctionsParams, this.coreFunctionsParams, true)?.funcName;
      const isAggr = funcName ?
        Object.entries(DG.AGG).map(([key, value]) => value).includes(funcName!.toLocaleLowerCase() as DG.AGG) : false;
      const escapedColName = grok.functions.handleOuterBracketsInColName(x.name, true);
      snippet = isAggr ? `\$[${escapedColName}]` : `\${${escapedColName}}`;
    } else if (this.typeOf(x, DG.Func)) {
      const params = (x as DG.Func).inputs.map((it) => it.semType ?? it.propertyType);
      const colPos = this.findColumnTypeMatchingParam(x);
      if (colPos !== -1) {
        const isAggr = Object.entries(DG.AGG).map(([key, value]) => value)
          .includes((x as DG.Func).name.toLocaleLowerCase() as DG.AGG);
        const escapedColName = grok.functions.handleOuterBracketsInColName(this.selectedColumn!.name, true);
        params[colPos] = isAggr ? `\$[${escapedColName}]` : `\${${escapedColName}}`;
      }
      const paramsStr = params.join(', ');
      const funcName = (x as DG.Func).nqName.startsWith('core:') ? (x as DG.Func).name : (x as DG.Func).nqName;
      snippet = `${funcName}(${paramsStr})`;
    } else
      return;
    const value = cm.state.doc.toString();
    const cursorPos = cm.state.selection.main.head;
    const newValue = value.slice(0, cursorPos) + snippet + value.slice(cursorPos);
    cm.dispatch({
      changes: {
        from: 0,
        to: cm.state.doc.length,
        insert: newValue,
      },
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
  getResultColumnName(): {colName: string, unusedName: string} {
    const input = this.inputName!.input as HTMLInputElement;
    let value: string = input.value;
    if ((value ?? '') == '')
      value = input.placeholder;
    return {colName: value, unusedName: this.sourceDf!.columns.getUnusedName(value)};
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
    selectedColumns.forEach((v) => { types.push(v.type), semTypes.push(v.semType); });
    types = types.filter((v, i, a) => a.indexOf(v) === i);
    semTypes = semTypes.filter((v, i, a) => a.indexOf(v) === i && v != null);
    return [types, semTypes];
  }

  /** Adds a New Column to the source table or edit formula for an existing column. */
  async addNewColumnAction(): Promise<void> {
    if (this.edit) {
      const colToUpdate = this.sourceDf?.col(this.call!.getParamValue('name'));
      if (colToUpdate) {
        const name = this.widget ? this.call!.getParamValue('name') : this.inputName!.value;
        const type = this.widget ? this.call.getParamValue('type') : this.getSelectedType()[0];
        const treatAsString = this.widget ? this.call.getParamValue('treatAsString') : this.getSelectedType()[1];
        await colToUpdate.applyFormula(this.codeMirror!.state.doc.toString().trim(), type, treatAsString);
        if (name !== colToUpdate.name)
          colToUpdate.name = this.sourceDf?.columns.getUnusedName(name) ?? name;
        grok.shell.o = colToUpdate;
      } else
        grok.shell.error(`Column ${this.call!.getParamValue('name')} is missing in the table`);
    } else {
      if (!this.call.getParamValue('table'))
        this.call.setParamValue('table', this.sourceDf);
      this.call.setParamValue('name', this.edit ? this.inputName!.value : this.getResultColumnName().unusedName);
      this.call.setParamValue('expression', this.codeMirror!.state.doc.toString().trim());
      this.call.setParamValue('type', this.getSelectedType()[0]);
      this.call.setParamValue('treatAsString', this.getSelectedType()[1]);
      if (!this.edit)
        this.call.setParamValue('subscribeOnChanges', true);
      await this.call.call(false, undefined, {processed: false});
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

  functionsCompletions(colNames: string[], packageNames: string[], fullPackageFunctionNames: string[],
    coreFunctionsNames: string[], packageFunctionsNames: {[key: string]: string[]},
    packageFunctionsParams: {[key: string]: FuncInfo}, coreFunctionsParams: {[key: string]: FuncInfo}) {
    return (context: CompletionContext) => {
      const word = context.matchBefore(/[\w|:|$|${|$\[|'|"]*/);
      if (!word || word?.from === word?.to && !context.explicit)
        return null;

      //check if word is inside quotes and if yes - do not show autocomplete
      const {quotesSelection} =
        this.getIntervalsWithinAndOutsideQuotes(context.state.doc.toString());
      for (let i = 0; i < quotesSelection.length; i++) {
        if (word.from >= quotesSelection[i].from && word.to <= quotesSelection[i].to)
          return null;
      }

      const options: Completion[] = [];
      let index = word!.from;
      let filter = true;
      const getFuncSignature = (name: string, propInfo: PropInfo[]) => {
        return `${name}(${propInfo
          .filter((it) => it.propName !== FUNC_OUTPUT_TYPE)
          .map((it)=> it.propType).join(',')})`;
      };
      if (word.text.includes(':')) {
        const colonIdx = word.text.indexOf(':');
        const packNameCaseIns = word.text.substring(0, word.text.indexOf(':'));
        const packName = getKeyCaseInsensitive(packageFunctionsNames, packNameCaseIns);
        if (packName) {
          packageFunctionsNames[packName].forEach((name: string, idx: number) => {
            options.push({
              label: name,
              type: 'variable',
              apply: getFuncSignature(name, packageFunctionsParams[`${packName}:${name}`].params),
            });
          });
        }
        index = word!.from + colonIdx + 1;
        filter = !word.text.endsWith(':');
      } else if (word.text.includes('$') || word.text.includes('${') || word.text.includes('$[')) {
        const dollarIdx = word.text.lastIndexOf('$');
        //check if there is a function before dollar sign (only in case index > 2: function name at least one letter and an opening brace) and if it is an aggregation function which requires square braces
        let openingSym = '{';
        if (context.view && index > 1) {
          const funcName = this.getFunctionNameAtPosition(context.view!, index - 2, -1, this.packageFunctionsParams,
            this.coreFunctionsParams, true)?.funcName;
          const isAggr = funcName ?
            Object.entries(DG.AGG).map(([key, value]) => value).includes(funcName!.toLocaleLowerCase() as DG.AGG) : false;
          if (isAggr)
            openingSym = '[';
        }
        const closingSym = openingSym === '{' ? '}' : ']';
        const openingBracketIdx = word.text.indexOf(openingSym) > dollarIdx ? word.text.indexOf(openingSym) : -1;
        const closingBracket = context.state.doc.length > word.text.length ?
          context.state.doc.toString()[word.to] === closingSym : false;
        colNames.forEach((name: string) => options.push({label: name, type: 'variable',
          apply: openingBracketIdx !== -1 ?
            closingBracket ? `${grok.functions.handleOuterBracketsInColName(name, true)}` :
              `${grok.functions.handleOuterBracketsInColName(name, true)}${closingSym}` :
            closingBracket ? `${openingSym}${grok.functions.handleOuterBracketsInColName(name, true)}` :
              `${openingSym}${grok.functions.handleOuterBracketsInColName(name, true)}${closingSym}`}));
        index = word!.from + (openingBracketIdx === -1 ? word.text.lastIndexOf('$') + 1 : openingBracketIdx + 1);
        filter = !word.text.endsWith('$') && !word.text.endsWith(openingSym);
      } else {
        //number of core functions
        const cf = coreFunctionsNames.length;
        //number of core functions + package functions
        const cpf = cf + fullPackageFunctionNames.length;
        coreFunctionsNames.concat(fullPackageFunctionNames).concat(packageNames)
          .forEach((name: string, idx: number) => options.push({
            label: name, type: 'variable',
            apply: idx < cf ? getFuncSignature(name, coreFunctionsParams[name].params) :
              idx < cpf ? getFuncSignature(name, packageFunctionsParams[name].params) :
                `${name}:`,
            detail: idx < cpf ? '' : 'package',
            section: idx < cpf ? '' : 'packages',
          }));
      }
      return {
        from: index,
        options: options,
        filter: filter,
      } as CompletionResult;
    };
  }

  getValueType(value: any): string {
    if (!value)
      return 'null';
    let type: string = typeof value;
    if (typeof value === 'number') {
      const str = value.toString();
      type = str.indexOf('.') === -1 && str.indexOf(',') === -1 ? DG.TYPE.INT : DG.TYPE.FLOAT;
    }
    return type;
  }

  getIfFuncOutputParam(call: DG.FuncCall) {
    let outType = '';

    const values = Object.entries(call.inputParams).map(([key, value]) => value);

    for (const ip of values) {
      if (ip.name == 'ifTrue' || ip.name == 'ifFalse') {
        let pType = '';
        if ((ip as DG.FuncCallParam).value instanceof DG.FuncCall &&
          (ip as DG.FuncCallParam).value.func.name.toLowerCase() === 'if')
          pType = this.getIfFuncOutputParam((ip as DG.FuncCallParam).value);
        else { //in case columns is passed as parameter (GetCurrentRowFieldFunc) - we take column type
          pType = (ip as DG.FuncCallParam).value instanceof DG.FuncCall ?
            (ip as DG.FuncCallParam).value.func.name === 'GetCurrentRowField' ?
              this.sourceDf?.col((ip as DG.FuncCallParam).value.inputs['field'])?.type ?? 'null' :
              Object.keys((ip as DG.FuncCallParam).value.outputParams).length > 0 ?
                // eslint-disable-next-line max-len
                (ip as DG.FuncCallParam).value.outputParams[Object.keys((ip as DG.FuncCallParam).value.outputParams)[0]].property.propertyType : 'null' :
            this.getValueType((ip as DG.FuncCallParam).value);
        }
        if (outType == '') {
          outType = pType;
          continue;
        } else
          outType = this.getOutputParamType(outType, pType);
      }
    }
    return outType;
  }

  getOutputParamType(first: string, second: string) {
    const isNullParam = (param: string) => param == 'null' || param == 'undefined' || param == null;
    if (isNullParam(first))
      return !isNullParam(second) ? second : DG.TYPE.STRING;
    else {
      if (isNullParam(second))
        return first;
      else {
        if (second.toLowerCase() == first.toLowerCase())
          return second;
        else {
          if (first == DG.TYPE.DYNAMIC || second == DG.TYPE.DYNAMIC)
            return DG.TYPE.DYNAMIC;
          if (isNumerical(second) && isNumerical(first)) {
            const resType = numericCast[first][second];
            if (resType == '')
              // eslint-disable-next-line max-len
              throw new Error(`If function params types (${first}, ${second}) do not match and cannot be casted to each other`);
            return resType;
          } else
            // eslint-disable-next-line max-len
            throw new Error(`If function params types (${first}, ${second}) do not match and cannot be casted to each other`);
        }
      }
    }
  }
}

export function prepareAddNewColumnFuncCall(col: DG.Column): DG.FuncCall {
  const expression = col.getTag(DG.Tags.Formula);
  const table = col.dataFrame;
  const f = DG.Func.byName('AddNewColumn');
  const fc = f.prepare({
    'table': table,
    'expression': expression,
    'name': col.name,
    'type': col.type,
  });
  return fc;
}

function getKeyCaseInsensitive(obj: Record<string, any>, keyCaseIns: string): string | undefined {
  let key;
  for (const objKey of Object.keys(obj)) {
    if (objKey.toLowerCase() === keyCaseIns.toLowerCase())
      key = objKey;
  }
  return key;
}
