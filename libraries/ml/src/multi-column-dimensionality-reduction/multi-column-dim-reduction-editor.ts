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
import { vectorDistanceMetricsMethods, VectorMetricsNames } from '../typed-metrics';

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
    {uiName: 'Use WebGPU', value: false, tooltip: 'Use WebGPU for Distance and UMAP computations', type: 'boolean',
      disableTooltip: 'WebGPU is not available'};
  constructor() {
    getGPUAdapterDescription().then((desc) => {
      if (desc) {
        this.useWebGPU.tooltip += ` (${desc})`;
        this.useWebGPU.value = true; // enable by default
      } else {
        this.useWebGPU.value = false;
        this.useWebGPU.disable = true;
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
    columnFunctionsMap: {[key: string]: string[]} = {};
    methodsParams: {[key: string]: UMAPOptions | TSNEOptions} = {
      [DimReductionMethods.UMAP]: new UMAPOptions(),
      [DimReductionMethods.T_SNE]: new TSNEOptions()
    };
    dbScanParams = new DBScanOptions();
    methodSettingsDivs: HTMLElement[] = [];
    supportedFunctions: {[name: string]: DimRedSupportedFunctions} = {};
    methodInput: DG.InputBase<string | null>;
    methodSettingsIcon: HTMLElement;
    methodSettingsAnchor: HTMLElement = ui.div();
    plotEmbeddingsInput = ui.input.bool('Plot embeddings', {value: true});
    postProcessingEditor: PostProcessingFuncEditor;
    aggregationMethodInput = ui.input.choice('Aggregation', {value: DistanceAggregationMethods.EUCLIDEAN,
      items: [DistanceAggregationMethods.EUCLIDEAN, DistanceAggregationMethods.MANHATTAN]});
    vectorDistanceInput: DG.InputBase<string | null> = ui.input.choice('Distance metric', {value: VectorMetricsNames.Euclidean, items:[VectorMetricsNames.Euclidean, VectorMetricsNames.Manhattan, VectorMetricsNames.Cosine]});

    onColumnsChanged: Subject<void>;
    constructor(editorSettings: DimReductionEditorOptions = {}) {
      this.aggregationMethodInput.setTooltip('Aggregation method for combining distances between columns');
      this.vectorDistanceInput.root.style.display = 'none';
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

      this.tableInput = ui.input.table('Table', {value: grok.shell.tv.dataFrame, items: grok.shell.tables,
        onValueChanged: () => {
          this.onTableInputChanged();
        }
      });
      this.onTableInputChanged();
      let settingsOpened = false;

      this.methodInput = ui.input.choice('Method', {value: DimReductionMethods.UMAP,
        items: [DimReductionMethods.UMAP, DimReductionMethods.T_SNE], onValueChanged: (value) => {
          if (settingsOpened)
            this.createAlgorithmSettingsDiv(this.methodsParams[value]);
        }});
      this.methodSettingsIcon = ui.icons.settings(()=> {
        settingsOpened = !settingsOpened;
        if (!settingsOpened) {
          this.methodSettingsDivs.forEach((it) => it.remove());
          this.methodSettingsDivs = [];
        } else { this.createAlgorithmSettingsDiv(this.methodsParams[this.methodInput.value!]); }
      }, 'Modify methods parameters');

      this.methodInput.root.classList.add('ml-dim-reduction-settings-input');
      this.methodInput.root.prepend(this.methodSettingsIcon);
      this.columnParamsEditorAccordion.addPane('Column options', () => this.columnOptEditorsRoot, true, null, false);
      this.columnParamsEditorAccordion.root.style.display = 'none';
      this.columnParamsEditorRoot.appendChild(this.columnParamsEditorAccordion.root);
      this.columnParamsEditorRoot.appendChild(this.vectorDistanceInput.root);
    }

    onTableInputChanged() {
      const table = this.tableInput.value;
      if (!table)
        return;
      ui.empty(this.columnOptEditorsRoot);

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
            (col.meta.units && units.includes(col.meta.units));
          if (semTypeSupported && typeSuported && unitsSupported) {
            if (!this.columnFunctionsMap[col.name])
              this.columnFunctionsMap[col.name] = [];
            this.columnFunctionsMap[col.name].push(funcName);
          }
        });
      });
      const supportedColNames = Object.keys(this.columnFunctionsMap);
      const columnsInput = ui.input.columns('Columns', {table: table, onValueChanged: (value) => {
        this.onColumnsChanged.next();
        ui.empty(this.columnOptEditorsRoot);
        if (!value || value?.length < 2)
          this.aggregationMethodInput.root.style.display = 'none';
        else
          this.aggregationMethodInput.root.style.display = 'flex';
        if (!value || value.length === 0) {
          this.columnParamsEditorAccordion.root.style.display = 'none';
          return;
        }

        this.columnOptEditors = value.map((col) => {
          const editorClass = new DimReductionColumnEditor(col, this.columnFunctionsMap[col.name].map((it) =>
            this.supportedFunctions[it]));
          return editorClass;
        });
        const editorsRoot = ui.divV([], {style: {maxHeight: '400px', overflow: 'auto'}});
        this.columnOptEditors.forEach((editor) => {
          editorsRoot.appendChild(editor.accordionDiv);
        });
        //this.columnOptEditorsRoot.appendChild(editorsRoot);
        const doubledColEditors = new Array<any>(this.columnOptEditors.length * 2).fill(null)
          .map((_, i) => i % 2 === 0 ? this.columnOptEditors[i / 2].colOptEditors : []);
        let c = 0;
        const table = ui.table(doubledColEditors, (item) => {
          c++;
          if (item && item.length > 0)
            return item;
          const paramsEditor =
            (this.columnOptEditors[Math.floor((c - 1) / 2)].preprocessingFuncSettingsDiv = ui.div([]));
          return [paramsEditor, ui.div(), ui.div(), ui.div()];
        }, ['Column', 'Encoding function', 'Similarity metric', 'Weight']);
        this.columnOptEditors
          .forEach((it) => {
            it.preprocessingFuncSettingsDiv?.parentElement?.setAttribute('colspan', '4');
            it.preprocessingFuncSettingsDiv?.parentElement?.parentElement?.style?.setProperty('height', 'unset');
          });
        if (this.columnOptEditors.length > 0)
          this.columnParamsEditorAccordion.root.style.display = 'flex';
        table.classList.add('ml-dim-reduction-column-editor-table-root');
        this.columnOptEditorsRoot.appendChild(table);

        if (value.every((it) => it.isNumerical)) {
          this.aggregationMethodInput.root.style.display = 'none';
          this.columnParamsEditorAccordion.root.style.display = 'none';
          this.vectorDistanceInput.root.style.display = 'flex';
        } else {
          this.vectorDistanceInput.root.style.display = 'none';
        }
      }, available: supportedColNames});
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
      params: UMAPOptions | TSNEOptions | DBScanOptions): void {
      this.methodSettingsDivs.forEach((it) => it.remove());
      this.methodSettingsDivs = [];
      const anchor = this.methodSettingsAnchor;
      const parent = anchor.parentElement as HTMLDivElement | null;
      if (!parent)
        return;
      Object.keys(params).forEach((it: any) => {
        const param: IDimReductionParam | IDimReductionParam<string> | IDimReductionParam<boolean> =
          (params as any)[it];

        const input = param.type === 'string' ?
          ui.input.string(param.uiName, {value: param.value ?? '', onValueChanged: (value) => {
            param.value = value;
          }}) : param.type === 'boolean' ?
            ui.input.bool(param.uiName, {value: param.value ?? false, onValueChanged: (value) => {
              param.value = value;
            }}) :
            ui.input.float(param.uiName, {value: param.value as any, onValueChanged: (value) => {
              param.value = value;
            }});
        if (param.disable) {
          input.enabled = false;
          ui.tooltip.bind(input.input ?? input.root, param.disableTooltip ?? '');
        } else { ui.tooltip.bind(input.input ?? input.root, param.tooltip); }
        parent.insertBefore(input.root, anchor);
        this.methodSettingsDivs.push(input.root);
      });
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
        this.methodSettingsAnchor,
        this.plotEmbeddingsInput,
        this.postProcessingEditor.root,
      ], {style: {minWidth: '420px'}, classes: 'dim-reduction-dialog-form'});
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
        aggreaggregationMethod: this.aggregationMethodInput.value,
        vectorDistanceMetric: this.vectorDistanceInput.value as VectorMetricsNames | undefined
      };
    }

    // used by history
    public getInput() {
      return {
        columns: this.columnsInput.value.map((it) => it.name),
        method: this.methodInput.value,
        preprocessingFunctions: this.columnOptEditors.map((it) => it.preprocessingFunctionInput.value),
        distanceMetrics: this.columnOptEditors.map((it) => it.similarityMetricInput.value!),
        weights: this.columnOptEditors.map((it) => it.weight ?? 1),
        options: {...this.algorithmOptions, ...this.dbScanOptions,
          preprocessingFuncArgs: this.columnOptEditors.map((it) => it.preprocessingFunctionSettings)},
        plotEmbeddings: this.plotEmbeddingsInput.value,
        postProcessingFunction: this.postProcessingEditor.postProcessingFunctionInput.value ?? null,
        postProcessingFunctionArgs: this.postProcessingEditor.args,
        aggreaggregationMethod: this.aggregationMethodInput.value,
        vectorDistanceMetric: this.vectorDistanceInput.value
      };
    }

    public getStringInput() {
      return JSON.stringify(this.getInput());
    }

    public async applyStringInput(input: string) {
      try {
        const parsed = JSON.parse(input);
        await this.applyInput(parsed);
      } catch (e) {
        grok.shell.error('Error applying input from history');
        console.error(e);
      }
    }

    public async applyInput(input: ReturnType<typeof this.getInput>) {
      try {
        const cols = input.columns.map((it) => this.tableInput.value!.col(it));
        if (cols.some((it) => !it))
          throw new Error('Some columns are not found');
        this.columnsInput.value = cols as DG.Column[];
        this.columnsInput.fireChanged();
        this.methodInput.value = input.method;
        this.columnOptEditors.forEach((it, i) => {
          it.preprocessingFunctionInput.value = input.preprocessingFunctions[i];
          it.similarityMetricInput.value = input.distanceMetrics[i];
          it.weightInput.value = input.weights[i];
          it.preprocessingFunctionSettings = input.options.preprocessingFuncArgs[i];
        });
        this.plotEmbeddingsInput.value = input.plotEmbeddings;
        this.postProcessingEditor.postProcessingFunctionInput.value = input.postProcessingFunction;
        await this.postProcessingEditor._prevChangePromise;
        this.postProcessingEditor._postProcessingArgs = input.postProcessingFunctionArgs;
        this.aggregationMethodInput.value = input.aggreaggregationMethod;
        const copiedAlgoOptions: Options = {};
        Object.keys(this.methodsParams[input.method!]).forEach((key) => {
          copiedAlgoOptions[key] = input.options[key as keyof typeof input.options];
        });
        Object.keys(this.methodsParams[input.method!]).forEach((key) => {
          ((this.methodsParams[input.method!] as any)[key]).value = copiedAlgoOptions[key];
        });
        const wasSettingsOpened = this.methodSettingsDivs.length > 0;
        this.createAlgorithmSettingsDiv(this.methodsParams[input.method!]);
        if (!wasSettingsOpened) {
          this.methodSettingsDivs.forEach((it) => it.remove());
          this.methodSettingsDivs = [];
        }

        await this.postProcessingEditor.onFunctionChanged(input.postProcessingFunctionArgs);

        this.vectorDistanceInput.value = input.vectorDistanceMetric;
      } catch (e) {
        grok.shell.error('Error applying input from history');
        console.error(e);
      }
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
    colOptEditors: HTMLElement[] = [];
    constructor(column: DG.Column, supportedFunctions: DimRedSupportedFunctions[]) {
      this.weightInput = ui.input.float('Weight',
        {value: 1, onValueChanged: (value) => { this.weight = value ?? 1; }});
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
      this.preprocessingFunctionInput = ui.input.choice('Encoding function',
        {value: getFuncName(this.supportedFunctions[0].func),
          items: this.supportedFunctions.map((it) => getFuncName(it.func)),
          onValueChanged: (value) => {
            const func = this.functionsMap[value];
            this.preprocessingFunctionSettings = {};
            this.hasExtraSettings = func.inputs.length > 2;
            const supF = this.supportedFunctions.find((it) => getFuncName(it.func) === value)!;
            this.getSimilarityMetricInput(supF);
            ui.empty(this.preprocessingFuncSettingsDiv);
            settingsOpened = false;
            if (!this.hasExtraSettings)
              this.preprocessingFuncSettingsIcon.style.display = 'none';
            else
              this.preprocessingFuncSettingsIcon.style.display = 'flex';
          }});
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
      this.colOptEditors = [
        columnTitle, this.preprocessingFunctionInput.root,
        this.similarityMetricInputRoot, this.weightInput.root
      ];
      // assign tooltips
      ui.tooltip.bind(columnTitle, this.column.name);

      // this.colOptEditors.forEach((it, i) => {
      //   it.style.width = `${editorWidths[i]}%`;
      // });

      //add classes
      this.colOptEditors.forEach((it) => it.classList.add('ml-dim-reduction-column-editor-input-root'));

      const distanceOptionsDiv = ui.divH(
        this.colOptEditors, {classes: 'ml-dim-reduction-column-editor-root'}
      );

      this.accordionDiv = ui.divV([]);
      this.editorDiv.appendChild(distanceOptionsDiv);
      this.editorDiv.appendChild(this.preprocessingFuncSettingsDiv);
      this.accordionDiv.appendChild(this.editorDiv);
    }

    getSimilarityMetricInput(preFunc: DimRedSupportedFunctions) {
      const input =
        ui.input.choice('Similarity metric', {value: preFunc.distanceFunctions[0], items: preFunc.distanceFunctions});
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
        input.onChanged.subscribe((value) => { this.preprocessingFunctionSettings[fInput.name] = value; });
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
  public _postProcessingArgs: Options = {};
  private _argsElement: HTMLElement = ui.div([]);
  private _settingsIcon: HTMLElement;
  private _settingsOpened: boolean = false;
  public _prevChangePromise: Promise<void> = Promise.resolve();
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
      ui.input.choice('Postprocessing', {value: defaultPostProcessingFunc,
        items: Object.keys(this.postProcessingFunctionsMap),
        onValueChanged: async () => { await this.onFunctionChanged(); },
        nullable: false});
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

  async onFunctionChanged(presetArgs: Options = {}) {
    await this._prevChangePromise;
    let resolve: () => void;
    this._prevChangePromise = new Promise<void>((res) => resolve = res);
    try {
      const func = this.postProcessingFunction;
      ui.empty(this._argsElement);
      this._postProcessingArgs = {};
      if (!func || func.inputs.length < 3) {
        this._settingsIcon.style.display = 'none';
        return;
      };
      this._settingsIcon.style.display = 'flex';
      const fc = func.prepare(presetArgs);
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
        input.onChanged.subscribe((value) => { this._postProcessingArgs[fInput.name] = value; });
        this._argsElement.append(input.root);
      }
    } catch (e) {
      grok.shell.error('Error applying postprocessing function');
      console.error(e);
    } finally {
      resolve!();
    }
  }

  get root() {
    return this._root;
  }

  get args() {
    return this._postProcessingArgs;
  }
}
