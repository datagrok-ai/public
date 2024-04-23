import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DimReductionMethods} from './types';
import {DimReductionEditorOptions, DBScanOptions} from '../functionEditors/dimensionality-reduction-editor';
import {IDBScanOptions} from '@datagrok-libraries/math';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import {DIM_RED_DEFAULT_POSTPROCESSING_FUNCTION_META, DIM_RED_POSTPROCESSING_FUNCTION_TAG,
  DIM_RED_PREPROCESSING_FUNCTION_TAG, SUPPORTED_DISTANCE_FUNCTIONS_TAG,
  SUPPORTED_SEMTYPES_TAG, SUPPORTED_TYPES_TAG, SUPPORTED_UNITS_TAG} from '../functionEditors/consts';
import {DistanceAggregationMethods} from '../distance-matrix/types';
import {IDimReductionParam, ITSNEOptions, IUMAPOptions} from './multi-column-dim-reducer';
import {Subject} from 'rxjs';
import '../../css/styles.css';
import {getGPUAdapterDescription} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';

export class UMAPOptions {
  learningRate: IDimReductionParam =
  {uiName: 'Learninig rate', value: 1, tooltip: 'The initial learning rate for the embedding optimization'};
  // nComponents: IDimReductionParam =
  // {uiName: 'Components', value: 2, tooltip: 'The number of components (dimensions) to project the data to'};
  nEpochs: IDimReductionParam =
    {uiName: 'Epochs', value: 0,
      tooltip: 'The number of epochs to optimize embeddings via SGD. Computed automatically if set to 0'};
  nNeighbors: IDimReductionParam =
  {uiName: 'Neighbors', value: 15, tooltip: 'The number of nearest neighbors to construct the fuzzy manifold'};
  spread: IDimReductionParam =
  {uiName: 'Spread', value: 1,
    tooltip:
    `The effective scale of embedded points, used with min distance to control 
    the clumped/dispersed nature of the embedding`};
  minDist: IDimReductionParam =
  {uiName: 'Min distance', value: 0.1,
    tooltip: `The effective minimum distance between embedded points, 
  used with spread to control the clumped/dispersed nature of the embedding`};
  randomSeed: IDimReductionParam<string> = {uiName: 'Random seed', value: null, tooltip: 'Random seed', type: 'string'};
  useWebGPU: IDimReductionParam<boolean> =
    {uiName: 'Use WebGPU', value: false, tooltip: 'Use WebGPU for KNN computations', type: 'boolean'};
  constructor() {
    getGPUAdapterDescription().then((desc) => {
      if (desc) {
        this.useWebGPU.tooltip += ` (${desc})`;
      } else {
        this.useWebGPU.value = false;
        this.useWebGPU.forceRemove = true;
      }
    });
  };
}

export class TSNEOptions {
  epsilon: IDimReductionParam =
  {uiName: 'Epsilon', value: 10, tooltip: 'Epsilon is learning rate'};
  perplexity: IDimReductionParam =
  {uiName: 'Perplexity', value: 30, tooltip: 'Roughly how many neighbors each point influences'};
  // dim: IDimReductionParam = {uiName: 'Dimensionality', value: 2, tooltip: 'Dimensionality of the embedding'};

  constructor() {};
}

export type DimRedSupportedFunctions = {
    func: DG.Func,
    semTypes: string[],
    types: string[],
    units: string[],
    distanceFunctions: string[]
}

export type MultiColDimReductionParams = {
    table: DG.DataFrame,
    columns: DG.ColumnList,
    methodName: DimReductionMethods,
    preprocessingFunctions: DG.Func[],
    distanceMetrics: string[],
    weights: number[],
    options: (IUMAPOptions | ITSNEOptions) & Partial<IDBScanOptions> & {preprocessingFuncArgs: Options[]} & Options,
    plotEmbeddings?: boolean,
    clusterEmbeddings?: boolean,
}

export class MultiColumnDimReductionEditor {
    editorSettings: DimReductionEditorOptions = {};
    tableInput!: DG.InputBase<DG.DataFrame | null>;
    columnsInput!: DG.InputBase<DG.Column[]>;
    columnsInputRoot!: HTMLElement;
    columnOptEditors: DimReductionColumnEditor[] = [];
    columnOptEditorsRoot: HTMLElement = ui.div();
    columnParamsEditorRoot: HTMLElement = ui.div();
    columnParamsEditorAccordion: DG.Accordion;
    weightsEditorRoot: HTMLElement = ui.div();
    columnFunctionsMap: {[key: string]: string[]} = {};
    methodsParams: {[key: string]: UMAPOptions | TSNEOptions} = {
      [DimReductionMethods.UMAP]: new UMAPOptions(),
      [DimReductionMethods.T_SNE]: new TSNEOptions()
    };
    dbScanParams = new DBScanOptions();
    methodSettingsDiv = ui.inputs([]);
    dbScanSettingsDiv = ui.inputs([]);
    supportedFunctions: {[name: string]: DimRedSupportedFunctions} = {};
    methodInput: DG.InputBase<string | null>;
    methodSettingsIcon: HTMLElement;
    dbScanSettingsIcon: HTMLElement;
    plotEmbeddingsInput = ui.boolInput('Plot embeddings', true);
    clusterEmbeddingsInput = ui.boolInput('Cluster embeddings', false);
    postProcessingEditor: PostProcessingFuncEditor;
    aggregationMethodInput = ui.choiceInput('Aggregation', DistanceAggregationMethods.EUCLIDEAN,
      [DistanceAggregationMethods.EUCLIDEAN, DistanceAggregationMethods.MANHATTAN]);

    onColumnsChanged: Subject<void>;
    constructor(editorSettings: DimReductionEditorOptions = {}) {
      this.aggregationMethodInput.setTooltip('Aggregation method for combining distances between columns');
      this.onColumnsChanged = new Subject();
      this.editorSettings = editorSettings;
      this.columnParamsEditorAccordion = ui.accordion();
      const preporcessingFuncs = DG.Func.find({tags: [DIM_RED_PREPROCESSING_FUNCTION_TAG]});
      preporcessingFuncs.forEach((f) => {
        const semTypes: string = f.options.get(SUPPORTED_SEMTYPES_TAG) ?? '';
        const name = f.friendlyName ?? f.name;
        const types: string = f.options.get(SUPPORTED_TYPES_TAG) ?? '';
        const units: string = f.options.get(SUPPORTED_UNITS_TAG) ?? '';
        const distanceFunctions: string = f.options.get(SUPPORTED_DISTANCE_FUNCTIONS_TAG) ?? '';
        if (this.editorSettings.semtype && !semTypes.includes(this.editorSettings.semtype))
          return;
        if (this.editorSettings.type && !types.includes(this.editorSettings.type))
          return;
        if (this.editorSettings.units && !units.includes(this.editorSettings.units))
          return;

        this.supportedFunctions[name] = {
          func: f,
          semTypes: semTypes ? semTypes.split(',') : [],
          types: types ? types.split(',') : [],
          units: units ? units.split(',') : [],
          distanceFunctions: distanceFunctions ? distanceFunctions.split(',') : [],
        };
      });

      this.postProcessingEditor = new PostProcessingFuncEditor();

      this.tableInput = ui.tableInput('Table', grok.shell.tv.dataFrame, grok.shell.tables, () => {
        this.onTableInputChanged();
      });
      this.onTableInputChanged();
      let settingsOpened = false;
      let dbScanSettingsOpened = false;

      this.methodInput = ui.choiceInput('Method', DimReductionMethods.UMAP,
        [DimReductionMethods.UMAP, DimReductionMethods.T_SNE], () => {
          if (settingsOpened)
            this.createAlgorithmSettingsDiv(this.methodSettingsDiv, this.methodsParams[this.methodInput.value!]);
        });
      this.methodSettingsIcon = ui.icons.settings(()=> {
        settingsOpened = !settingsOpened;
        if (!settingsOpened)
          ui.empty(this.methodSettingsDiv);
        else
          this.createAlgorithmSettingsDiv(this.methodSettingsDiv, this.methodsParams[this.methodInput.value!]);
      }, 'Modify methods parameters');
      this.dbScanSettingsIcon = ui.icons.settings(()=> {
        dbScanSettingsOpened = !dbScanSettingsOpened;
        if (!dbScanSettingsOpened)
          ui.empty(this.dbScanSettingsDiv);
        else
          this.createAlgorithmSettingsDiv(this.dbScanSettingsDiv, this.dbScanParams);
      }, 'Modify clustering parameters');
      this.clusterEmbeddingsInput.classList.add('ml-dim-reduction-settings-input');
      this.clusterEmbeddingsInput.root.prepend(this.dbScanSettingsIcon);
      this.methodInput.root.classList.add('ml-dim-reduction-settings-input');
      this.methodInput.root.prepend(this.methodSettingsIcon);
      this.columnParamsEditorAccordion.addPane('Column options', () => this.columnOptEditorsRoot, true, null, false);
      this.columnParamsEditorAccordion.root.style.display = 'none';
      this.columnParamsEditorRoot.appendChild(this.columnParamsEditorAccordion.root);
      this.columnParamsEditorRoot.appendChild(this.weightsEditorRoot);
    }

    onTableInputChanged() {
      const table = this.tableInput.value;
      if (!table)
        return;
      ui.empty(this.columnOptEditorsRoot);
      ui.empty(this.weightsEditorRoot);

      this.columnFunctionsMap = {};
      const columns = table.columns.toList();
      columns.forEach((col) => {
        Object.keys(this.supportedFunctions).forEach((funcName) => {
          const semTypes = this.supportedFunctions[funcName].semTypes;
          const types = this.supportedFunctions[funcName].types;
          const units = this.supportedFunctions[funcName].units;
          const semTypeSupported = !semTypes.length || (col.semType && semTypes.includes(col.semType));
          const typeSuported = !types.length || types.includes(col.type);
          const unitsSupported = !units.length ||
            (col.getTag(DG.TAGS.UNITS) && units.includes(col.getTag(DG.TAGS.UNITS)));
          if (semTypeSupported && typeSuported && unitsSupported) {
            if (!this.columnFunctionsMap[col.name])
              this.columnFunctionsMap[col.name] = [];
            this.columnFunctionsMap[col.name].push(funcName);
          }
        });
      });
      const supportedColNames = Object.keys(this.columnFunctionsMap);
      const columnsInput = ui.columnsInput('Columns', table, () => {
        this.onColumnsChanged.next();
        ui.empty(this.columnOptEditorsRoot);
        ui.empty(this.weightsEditorRoot);
        const cols = columnsInput.value;
        if (!cols || cols?.length < 2)
          this.aggregationMethodInput.root.style.display = 'none';
        else
          this.aggregationMethodInput.root.style.display = 'flex';
        if (!cols || cols.length === 0) {
          this.columnParamsEditorAccordion.root.style.display = 'none';
          return;
        }
        const editorWidths = [20, 30, 30, 20];
        this.columnOptEditors = cols.map((col) => {
          const editorClass = new DimReductionColumnEditor(col, this.columnFunctionsMap[col.name].map((it) =>
            this.supportedFunctions[it]), editorWidths);
          return editorClass;
        });
        const editorTitles = ['Column', 'Encoding function', 'Similarity metric', 'Weight'].map((it, i) =>
          ui.h1(it, {style: {width: `${editorWidths[i]}%`, margin: 0}}));

        ui.tooltip.bind(editorTitles[1], 'Encoding function for the column values');
        ui.tooltip.bind(editorTitles[2], 'Distance/Similarity metric for the encoded column values');
        ui.tooltip.bind(editorTitles[3], 'Weight of the column for combining distances between values');

        const editorTitleRoot = ui.divH(editorTitles, {classes: 'ml-dim-reduction-column-editor-header-root'});
        this.columnOptEditorsRoot.appendChild(editorTitleRoot);
        const editorsRoot = ui.divV([], {style: {maxHeight: '400px', overflow: 'auto'}});
        this.columnOptEditors.forEach((editor) => {
          editorsRoot.appendChild(editor.accordionDiv);
        });
        this.columnOptEditorsRoot.appendChild(editorsRoot);
        if (this.columnOptEditors.length > 0)
          this.columnParamsEditorAccordion.root.style.display = 'flex';
      }, {available: supportedColNames});
      columnsInput.fireChanged();
      if (!this.columnsInputRoot) {
        this.columnsInputRoot = columnsInput.root;
        this.columnsInput = columnsInput;
      } else {
        ui.empty(this.columnsInputRoot);
        this.columnsInput = columnsInput;
        Array.from(this.columnsInput.root.children)
          .forEach((it) => this.columnsInputRoot.appendChild(it));
      }
    }
    private createAlgorithmSettingsDiv(
      paramsForm: HTMLElement, params: UMAPOptions | TSNEOptions | DBScanOptions): HTMLElement {
      ui.empty(paramsForm);
      Object.keys(params).forEach((it: any) => {
        const param: IDimReductionParam | IDimReductionParam<string> | IDimReductionParam<boolean> =
          (params as any)[it];
        if (!param.forceRemove) {
          const input = param.type === 'string' ?
            ui.stringInput(param.uiName, param.value ?? '', () => {
              param.value = (input as DG.InputBase<string>).value;
            }) : param.type === 'boolean' ?
              ui.boolInput(param.uiName, param.value ?? false, () => {
                param.value = (input as DG.InputBase<boolean>).value;
              }) :
              ui.floatInput(param.uiName, param.value as any, () => {
                param.value = input.value;
              });
          ui.tooltip.bind(input.input ?? input.root, param.tooltip);
          paramsForm.append(input.root);
        }
      });
      return paramsForm;
    }

    get algorithmOptions(): IUMAPOptions | ITSNEOptions {
      const algorithmParams: UMAPOptions | TSNEOptions = this.methodsParams[this.methodInput.value!];
      const options: any = {};
      Object.keys(algorithmParams).forEach((key: string) => {
        if ((algorithmParams as any)[key].value != null)
          options[key] = (algorithmParams as any)[key].value;
      });
      return options;
    }

    get dbScanOptions(): IDBScanOptions {
      return {
        dbScanEpsilon: this.dbScanParams.epsilon.value ?? 0.01,
        dbScanMinPts: this.dbScanParams.minPts.value ?? 4
      };
    }

    public getEditor(): HTMLElement {
      const div = ui.div([
        this.tableInput.root,
        this.columnsInputRoot,
        this.columnParamsEditorRoot,
        this.aggregationMethodInput.root,
        this.methodInput.root,
        this.methodSettingsDiv,
        this.plotEmbeddingsInput,
        this.postProcessingEditor.root,
      ], {style: {minWidth: '420px'}, classes: 'ui-form'});
      return div;
    }

    public getParams() {
      return {
        table: this.tableInput.value!,
        columns: this.columnsInput.value!,
        methodName: this.methodInput.value!,
        preprocessingFunctions: this.columnOptEditors.map((it) => it.preProcessingFunction),
        distanceMetrics: this.columnOptEditors.map((it) => it.similarityMetricInput.value!),
        weights: this.columnOptEditors.map((it) => it.weight ?? 1),
        options: {...this.algorithmOptions, ...this.dbScanOptions,
          preprocessingFuncArgs: this.columnOptEditors.map((it) => it.preprocessingFunctionSettings)},
        plotEmbeddings: this.plotEmbeddingsInput.value,
        clusterEmbeddings: false,
        postProcessingFunction: this.postProcessingEditor.postProcessingFunction,
        postProcessingFunctionArgs: this.postProcessingEditor.args,
        aggreaggregationMethod: this.aggregationMethodInput.value
      };
    }
}


class DimReductionColumnEditor {
    preprocessingFuncSettingsDiv = ui.div([]);
    preprocessingFunctionInput: DG.InputBase<string | null>;
    preprocessingFuncSettingsIcon: HTMLElement;
    similarityMetricInput!: DG.InputBase<string | null>;
    similarityMetricInputRoot!: HTMLElement;
    preprocessingFunctionSettings: Options = {};
    accordionDiv: HTMLElement;
    column: DG.Column;
    supportedFunctions: DimRedSupportedFunctions[];
    editorDiv: HTMLElement = ui.div([]);
    hasExtraSettings: boolean = true;
    functionsMap: {[_: string]: DG.Func} = {};
    needsConfiguration: boolean = false;
    weightInput: DG.InputBase<number | null>;
    weight: number = 1;
    constructor(column: DG.Column, supportedFunctions: DimRedSupportedFunctions[], editorWidths: number[]) {
      this.weightInput = ui.floatInput('Weight', 1, () => { this.weight = this.weightInput.value ?? 1; });
      this.column = column;
      // sort by specificity
      this.supportedFunctions = supportedFunctions.sort((a, b) => {
        if ((a.units.length === 0 || b.units.length === 0) && a.units.length !== b.units.length)
          return b.units.length - a.units.length;
        if (a.units.length !== b.units.length)
          return a.units.length - b.units.length;
        if (a.semTypes.length === 0 || b.semTypes.length === 0)
          return b.semTypes.length - a.semTypes.length;
        if (a.semTypes.length !== b.semTypes.length)
          return a.semTypes.length - b.semTypes.length;
        return a.types.length - b.types.length;
      });
      this.supportedFunctions.forEach((f) => {
        this.functionsMap[getFuncName(f.func)] = f.func;
      });
      this.preprocessingFunctionInput = ui.choiceInput('Encoding function',
        getFuncName(this.supportedFunctions[0].func), this.supportedFunctions.map((it) => getFuncName(it.func)),
        () => {
          const val = this.preprocessingFunctionInput.value!;
          const func = this.functionsMap[val];
          this.preprocessingFunctionSettings = {};
          this.hasExtraSettings = func.inputs.length > 2;
          const supF = this.supportedFunctions.find((it) => getFuncName(it.func) === val)!;
          this.getSimilarityMetricInput(supF);
          ui.empty(this.preprocessingFuncSettingsDiv);
          settingsOpened = false;
          if (!this.hasExtraSettings)
            this.preprocessingFuncSettingsIcon.style.display = 'none';
          else
            this.preprocessingFuncSettingsIcon.style.display = 'flex';
        }
      );
      this.preprocessingFunctionInput.root.style.display = 'flex';

      this.createSettingsDiv(this.preprocessingFuncSettingsDiv, this.supportedFunctions[0].func)
        .then(() => { ui.empty(this.preprocessingFuncSettingsDiv); });

      this.getSimilarityMetricInput(this.supportedFunctions[0]);
      this.hasExtraSettings = this.supportedFunctions[0].func.inputs.length > 2;

      let settingsOpened = false;
      this.preprocessingFuncSettingsIcon = ui.icons.settings(async () => {
        settingsOpened = !settingsOpened;
        if (settingsOpened) {
          await this.createSettingsDiv(this.preprocessingFuncSettingsDiv,
            this.functionsMap[this.preprocessingFunctionInput.value!]);
        } else {
          ui.empty(this.preprocessingFuncSettingsDiv);
        }
      }, 'Modify encoding function parameters');

      this.preprocessingFunctionInput.root.classList.add('ml-dim-reduction-settings-input');
      this.preprocessingFunctionInput.root.prepend(this.preprocessingFuncSettingsIcon);

      if (!this.hasExtraSettings)
        this.preprocessingFuncSettingsIcon.style.display = 'none';
      else
        this.preprocessingFuncSettingsIcon.style.display = 'flex';

      this.needsConfiguration = !(supportedFunctions.length < 2 && !this.hasExtraSettings &&
      supportedFunctions[0].distanceFunctions.length < 2);

      const columnTitle = ui.h3(this.column.name, {classes: 'ml-dim-reduction-column-editor-column-title'});
      const colOptEditors = [
        columnTitle, this.preprocessingFunctionInput.root,
        this.similarityMetricInputRoot, this.weightInput.root
      ];
      // assign tooltips
      ui.tooltip.bind(columnTitle, this.column.name);

      colOptEditors.forEach((it, i) => {
        it.style.width = `${editorWidths[i]}%`;
      });

      //add classes
      colOptEditors.forEach((it) => it.classList.add('ml-dim-reduction-column-editor-input-root'));

      const distanceOptionsDiv = ui.divH(
        colOptEditors, {classes: 'ml-dim-reduction-column-editor-root'}
      );

      this.accordionDiv = ui.divV([]);
      this.editorDiv.appendChild(distanceOptionsDiv);
      this.editorDiv.appendChild(this.preprocessingFuncSettingsDiv);
      this.accordionDiv.appendChild(this.editorDiv);
    }

    getSimilarityMetricInput(preFunc: DimRedSupportedFunctions) {
      const input = ui.choiceInput('Similarity metric', preFunc.distanceFunctions[0], preFunc.distanceFunctions);
      if (!this.similarityMetricInputRoot) {
        this.similarityMetricInputRoot = input.root;
        this.similarityMetricInput = input;
      } else {
        ui.empty(this.similarityMetricInputRoot);
        this.similarityMetricInput = input;
        Array.from(this.similarityMetricInput.root.children)
          .forEach((it) => this.similarityMetricInputRoot.appendChild(it));
      }
    }

    get preProcessingFunction() {
      return this.functionsMap[this.preprocessingFunctionInput.value!];
    }
    async createSettingsDiv(paramsForm: HTMLElement, func: DG.Func): Promise<HTMLElement> {
      ui.empty(paramsForm);
      if (func.inputs.length < 3)
        return ui.div();
      const fc = func.prepare();
      const inputs = await fc.buildEditor(ui.div());
      for (let i = 2; i < func.inputs.length; i++) {
        const fInput = func.inputs[i];
        const val = this.preprocessingFunctionSettings[fInput.name] ||
            fc.inputParams[func.inputs[i].name].value || fInput.defaultValue;
        if (val)
          this.preprocessingFunctionSettings[fInput.name] = val;
        const input = inputs.find((inp) => inp.property.name === fInput.name);
        if (!input)
          continue;
        if (this.preprocessingFunctionSettings[fInput.name] !== null &&
            this.preprocessingFunctionSettings[fInput.name] !== undefined)
          input.value = this.preprocessingFunctionSettings[fInput.name];
        input.onChanged(() => { this.preprocessingFunctionSettings[fInput.name] = input.value; });
        paramsForm.append(input.root);
      }
      paramsForm.style.marginBottom = '10px';
      return paramsForm;
    }
}

function getFuncName(func: DG.Func) {
  return func.friendlyName ?? func.name;
}

class PostProcessingFuncEditor {
  postProcessingFunctionsMap: {[key: string]: DG.Func | null} = {};
  postProcessingFunctionInput: DG.ChoiceInput<string | null>;
  private _root: HTMLElement = ui.div([]);
  private _postProcessingArgs: Options = {};
  private _argsElement: HTMLElement = ui.div([]);
  private _settingsIcon: HTMLElement;
  private _settingsOpened: boolean = false;
  constructor() {
    this._settingsIcon = ui.icons.settings(async () => {
      this._settingsOpened = !this._settingsOpened;
      if (this._settingsOpened)
        this._argsElement.style.display = 'block';
      else
        this._argsElement.style.display = 'none';
    });
    this._argsElement.style.display = 'none';

    const postProcessingFuncs = DG.Func.find({tags: [DIM_RED_POSTPROCESSING_FUNCTION_TAG]})
      .filter((f) => f.inputs.length >= 2);
    postProcessingFuncs.forEach((f) => {
      const name = f.friendlyName ?? f.name;
      this.postProcessingFunctionsMap[name] = f;
    });
    this.postProcessingFunctionsMap['None'] = null;

    const defaultPostProcessingFunc = Object.keys(this.postProcessingFunctionsMap).find((it) =>
      !!this.postProcessingFunctionsMap[it]?.options?.[DIM_RED_DEFAULT_POSTPROCESSING_FUNCTION_META]) ?? 'None';
    this.postProcessingFunctionInput =
      ui.choiceInput('Postprocessing', defaultPostProcessingFunc,
        Object.keys(this.postProcessingFunctionsMap), async () => { await this.onFunctionChanged(); },
        {nullable: false});
    this.onFunctionChanged();
    this.postProcessingFunctionInput.nullable = false;
    this.postProcessingFunctionInput.classList.add('ml-dim-reduction-settings-input');
    this.postProcessingFunctionInput.root.prepend(this._settingsIcon);
    this._root.appendChild(this.postProcessingFunctionInput.root);
    this._root.appendChild(this._argsElement);
  }

  get postProcessingFunction() {
    return this.postProcessingFunctionInput.value ?
      this.postProcessingFunctionsMap[this.postProcessingFunctionInput.value] : null;
  }

  async onFunctionChanged() {
    const func = this.postProcessingFunction;
    ui.empty(this._argsElement);
    this._postProcessingArgs = {};
    if (!func || func.inputs.length < 3) {
      this._settingsIcon.style.display = 'none';
      return;
    };
    this._settingsIcon.style.display = 'flex';
    const fc = func.prepare();
    const inputs = await fc.buildEditor(ui.div());
    for (let i = 2; i < func.inputs.length; i++) {
      const fInput = func.inputs[i];
      const val = this._postProcessingArgs[fInput.name] ||
        fc.inputParams[func.inputs[i].name].value || fInput.defaultValue;
      if (val)
        this._postProcessingArgs[fInput.name] = val;
      const input = inputs.find((inp) => inp.property.name === fInput.name);
      if (!input)
        continue;
      input.onChanged(() => { this._postProcessingArgs[fInput.name] = input.value; });
      this._argsElement.append(input.root);
    }
  }

  get root() {
    return this._root;
  }

  get args() {
    return this._postProcessingArgs;
  }
}
