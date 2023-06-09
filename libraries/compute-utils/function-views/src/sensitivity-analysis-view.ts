/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {BehaviorSubject} from 'rxjs';
import {getPropViewers} from './shared/utils';
import {CARD_VIEW_TYPE, FUNCTIONS_VIEW_TYPE, SCRIPTS_VIEW_TYPE, VIEWER_PATH, viewerTypesMapping} from './shared/consts';
import {VarianceBasedSenstivityAnalysis} from './variance-based-analysis/sensitivityAnalysis';

const RUN_NAME_COL_LABEL = 'Run name' as const;
const RUN_ID_COL_LABEL = 'RunId' as const;

enum DISTRIB_TYPE {
  LINEAR = 'Linear',
  NORMAL = 'Normal',
  RANDOM = 'Random',
}
const DISTRIB_TYPES = [
  DISTRIB_TYPE.LINEAR,
  // DISTRIB_TYPE.NORMAL,
  // DISTRIB_TYPE.RANDOM,
];

type SensitivityNumericValues = {
  prop: DG.Property,
  const: number,
  min: number,
  max: number,
  lvl: number,
  distrib: DISTRIB_TYPE,
  isChanging: BehaviorSubject<boolean>,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
};

type SensitivityBoolValues = {
  prop: DG.Property,
  const: boolean,
  isChanging: BehaviorSubject<boolean>,
  type: DG.TYPE.BOOL,
}

type SensitivityConstValues = {
  prop: DG.Property,
  const: DG.InputBase,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL>,
}

type SensitivityValues = SensitivityNumericValues | SensitivityBoolValues | SensitivityConstValues;

type SensitivityNumericInputs = {
  constInput: DG.InputBase,
  minInput: DG.InputBase,
  maxInput: DG.InputBase,
  lvlInput: DG.InputBase,
  distribInput: DG.InputBase,
  isChangingInput: DG.InputBase<boolean | null>,
  prop: DG.Property,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
}

type SensitivityBoolInputs = {
  defaultInput: DG.InputBase,
  isChangingInput: DG.InputBase<boolean | null>,
  prop: DG.Property,
  type: DG.TYPE.BOOL,
}

type SensitivityConstInputs = {
  constInput: DG.InputBase,
  prop: DG.Property,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL>,
}

type SensitvityInputs = SensitivityNumericInputs | SensitivityBoolInputs | SensitivityConstInputs;


export class SensitivityAnalysisView {
  generateFormRows = (func: DG.Func) => {
    const store = (func.inputs.map((inputProp) => {
      switch (inputProp.propertyType) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
      case DG.TYPE.FLOAT:
        return {
          prop: inputProp,
          const: inputProp.defaultValue ?? 0,
          min: inputProp.defaultValue ?? 0,
          max: inputProp.defaultValue ?? 0,
          lvl: 3,
          distrib: DISTRIB_TYPES[0],
          isChanging: new BehaviorSubject(false),
          type: inputProp.propertyType,
        };
      case DG.TYPE.BOOL:
        return {
          prop: inputProp,
          const: inputProp.defaultValue,
          isChanging: new BehaviorSubject(false),
          type: inputProp.propertyType,
        };
      default:
        return {
          prop: inputProp,
          type: inputProp.propertyType,
          const: inputProp.defaultValue,
          isChanging: new BehaviorSubject(false),
        };
      }
    }) as SensitivityValues[]).reduce((acc, values, valIdx) => {
      acc[func.inputs[valIdx].name] = values;
      return acc;
    }, {} as Record<string, SensitivityValues>);

    const inputs = func.inputs.reduce((acc, inputProp) => {
      switch (inputProp.propertyType) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
      case DG.TYPE.FLOAT:
        const numericStore = store[inputProp.name] as SensitivityNumericValues;
        acc[inputProp.name] = {
          constInput: inputProp.propertyType === DG.TYPE.FLOAT ?
            ui.floatInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue, (v: number) => numericStore.const = v):
            ui.intInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue, (v: number) => numericStore.const = v),
          minInput: inputProp.propertyType === DG.TYPE.FLOAT ?
            ui.floatInput(`${inputProp.caption ?? inputProp.name} min`, inputProp.defaultValue, (v: number) => numericStore.min = v):
            ui.intInput(`${inputProp.caption ?? inputProp.name} min`, inputProp.defaultValue, (v: number) => numericStore.min = v),
          maxInput: inputProp.propertyType === DG.TYPE.FLOAT ?
            ui.floatInput(`${inputProp.caption ?? inputProp.name} max`, inputProp.defaultValue, (v: number) => numericStore.max = v):
            ui.intInput(`${inputProp.caption ?? inputProp.name} max`, inputProp.defaultValue, (v: number) => numericStore.max = v),
          lvlInput: ui.intInput('Levels', 3, (v: number) => numericStore.lvl = v),
          distribInput: ui.choiceInput('Distribution', DISTRIB_TYPES[0], DISTRIB_TYPES, (v: DISTRIB_TYPE) => numericStore.distrib = v),
          isChangingInput: (() => {
            const input = ui.switchInput(' ', numericStore.isChanging.value, (v: boolean) => numericStore.isChanging.next(v));            
            $(input.root).css({'min-width': '50px', 'width': '50px'});
            $(input.captionLabel).css({'min-width': '0px', 'max-width': '0px'});
            return input;
          })(),
          type: inputProp.propertyType,
          prop: inputProp,
        };
        const ref = acc[inputProp.name] as SensitivityNumericInputs;
        numericStore.isChanging.subscribe((newVal) => {
          if (newVal) {
            $(ref.constInput.root).hide();
            $(ref.minInput.root).show();
            $(ref.maxInput.root).show();
            if (!this.isVarianceBasedInput.value) {
              $(ref.lvlInput.root).show();
              $(ref.distribInput.root).show();
            } else {
              $(ref.lvlInput.root).hide();
              $(ref.distribInput.root).hide();
            }
          } else {
            $(ref.constInput.root).show();
            $(ref.minInput.root).hide();
            $(ref.maxInput.root).hide();
            $(ref.lvlInput.root).hide();
            $(ref.distribInput.root).hide();
          }
        });
        $([ref.minInput.root, ref.maxInput.root, ref.lvlInput.root, ref.distribInput.root]).css({
          'width': '25%',
        });
        break;
      case DG.TYPE.BOOL:
        const boolStore = store[inputProp.name] as SensitivityBoolValues;
        acc[inputProp.name] = {
          defaultInput: ui.boolInput('Default value', false, (v: boolean) => boolStore.const = v),
          isChangingInput: ui.switchInput('Is changing', boolStore.isChanging.value, (v: boolean) => boolStore.isChanging.next(v)),
          type: inputProp.propertyType,
          prop: inputProp,
        };
        break;
      default:
        const constStore = store[inputProp.name] as SensitivityConstValues;
        acc[inputProp.name] = {
          constInput: ui.input.forProperty(inputProp, undefined, {onValueChanged: (v: any) => constStore.const = v}),
          type: inputProp.propertyType,
          prop: inputProp,
        };
      }

      return acc;
    }, {} as Record<string, SensitvityInputs>);

    return {store, inputs};
  };

  isVarianceBasedInput = ui.switchInput('Use variance based analysis', false);
  samplesInput = ui.intInput('Samples', 100);
  inputConfig = this.generateFormRows(this.func);
  comparisonView!: DG.TableView;

  constructor(
    public func: DG.Func,
    public options: {
      name?: string,
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
      configFunc?: undefined,
    } = {
      name: undefined,
      parentView: undefined,
      parentCall: undefined,
      configFunc: undefined,
    },
  ) {
    const configFunc = options.configFunc ?? func;

    const allParamViewers = [
      ...configFunc.inputs,
      ...configFunc.outputs,
    ]
      .map((prop) => getPropViewers(prop))
      .reduce((acc, config) => {
        if (!acc[config.name])
          acc[config.name] = config.config;
        else
          acc[config.name].push(...config.config);
        return acc;
      }, {} as Record<string, Record<string, string | boolean>[]>);


    const comparisonDf = DG.DataFrame.create(0);

    const addColumnsFromProp = (configProp: DG.Property): DG.Column[] => {
      if (configProp.propertyType === DG.TYPE.DATA_FRAME) {
        const requestedViewersConfigs = allParamViewers[configProp.name];

        const viewerColumns = requestedViewersConfigs.map((config) => {
          let columnName = configProp.caption ?? configProp.name;
          const newColumn = DG.Column.fromType(DG.TYPE.DATA_FRAME, columnName, 0);
          const unusedName = comparisonDf.columns.getUnusedName(newColumn.name);
          newColumn.name = unusedName;
          columnName = unusedName;
          newColumn.temp[VIEWER_PATH] = config;
          comparisonDf.columns.add(newColumn);

          return newColumn;
        });

        return viewerColumns;
      } else {
        let columnName = configProp.caption ?? configProp.name;
        //@ts-ignore
        const newColumn = DG.Column.fromType(configProp.propertyType, columnName, 0);
        const unusedName = comparisonDf.columns.getUnusedName(newColumn.name);
        newColumn.name = unusedName;
        columnName = unusedName;
        comparisonDf.columns.add(newColumn);
        return [newColumn];
      }
    };

    configFunc.inputs.forEach((prop) => addColumnsFromProp(prop));
    configFunc.outputs.forEach((prop) => addColumnsFromProp(prop));

    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-12878
    const cardView = [...grok.shell.views].find((view) => view.type === CARD_VIEW_TYPE || view.type === SCRIPTS_VIEW_TYPE || view.type === FUNCTIONS_VIEW_TYPE);
    if (cardView) grok.shell.v = cardView;

    // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-12879
    this.comparisonView = grok.shell.addTableView(comparisonDf);
    this.comparisonView.temp = {'isComparison': true};

    const form = this.buildFormWithBtn();

    this.comparisonView.dockManager.dock(
      form,
      DG.DOCK_TYPE.LEFT,
      null,
      this.options.name ?? `${this.func.name} - Sensitivity Analysis`,
      0.25,
    );
  }

  private buildFormWithBtn() {
    const adaptStyle = () => {
      if (allPropInputs.some((input) => {
        return $(input).width() < 350 && $(input).width() > 0;
      }))
        $(form.container).addClass('ui-form-condensed');
      else
        $(form.container).removeClass('ui-form-condensed');
    };

    let prevCategory = 'Misc';
    const allPropInputs = [] as HTMLElement[];
    const form = Object.values(this.inputConfig.inputs)
      .reduce(({container, switches}, propConfig) => {
        let propInputs = [] as HTMLElement[];

        switch (propConfig.type) {
        case DG.TYPE.INT:
        case DG.TYPE.BIG_INT:
        case DG.TYPE.FLOAT:
          propConfig = propConfig as SensitivityNumericInputs;
          propInputs = [
            propConfig.isChangingInput.root,
            propConfig.constInput.root,
            propConfig.minInput.root,
            propConfig.maxInput.root,
            propConfig.lvlInput.root,
            propConfig.distribInput.root,
          ];
          propConfig.isChangingInput.onChanged(adaptStyle);
          allPropInputs.push(...propInputs.filter((_, idx) => idx > 0));
          switches.push(propConfig.isChangingInput);
          break;
        case DG.TYPE.BOOL:
          propConfig = propConfig as SensitivityBoolInputs;
          propInputs = [
            propConfig.isChangingInput.root,
            propConfig.defaultInput.root,
          ];
          propConfig.isChangingInput.onChanged(adaptStyle);
          allPropInputs.push(...propInputs.filter((_, idx) => idx > 0));
          switches.push(propConfig.isChangingInput);
          break;
        default:
          propInputs = [
            propConfig.constInput.root,
          ];
          allPropInputs.push(...propInputs);
        }

        const prop = propConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.h2(prop.category, {style: {'margin-top': '10px'}}));
          prevCategory = prop.category;
        }
        container.appendChild(ui.divH(propInputs, {style: {'max-width': '100%'}}));

        return {container, switches};
      }, {
        container: ui.divV([ui.divH([
          this.isVarianceBasedInput.root,
          this.samplesInput.root,
        ])], 'ui-form ui-form-wide ui-form-left'),
        switches: [] as DG.InputBase[],
      });

    const buttons = ui.buttonsInput([ui.bigButton('Run sensitivity analysis', async () => {
      if (!this.isVarianceBasedInput.value)
        this.runSimpleAnalysis();
      else
        this.runVarianceAnalysis();
    })]);
    $(buttons).css({
      'align-items': 'start',
    });
    form.container.appendChild(
      buttons,
    );

    form.switches.forEach((switchEl) => {
      $(switchEl.root).css({
        'width': '50px',
        'max-width': '50px',
        'min-width': '50px',
        'padding-top': '0px',
      });
    });

    $(form.container).css({
      'width': '100%',
      'max-width': '100%',
      'padding-top': '0px',
      'overflow-y': 'scroll',
    });

    ui.tools.handleResize(form.container, adaptStyle);

    return form.container;
  }

  private async runVarianceAnalysis() {
    const options = {
      func: this.func,
      fixedInputs: Object.keys(this.inputConfig.store).filter((propName) => {
        switch (this.inputConfig.store[propName].type) {
        case DG.TYPE.INT:
        case DG.TYPE.BIG_INT:
        case DG.TYPE.FLOAT:
          const numPropConfig = this.inputConfig.store[propName] as SensitivityNumericValues;
          return numPropConfig.isChanging.value === false;
        default:
          return true;
        }
      }).map((propName) => ({
        name: propName,
        value: this.inputConfig.store[propName].const,
      })),
      variedInputs: Object.keys(this.inputConfig.store).filter((propName) => {
        switch (this.inputConfig.store[propName].type) {
        case DG.TYPE.INT:
        case DG.TYPE.BIG_INT:
        case DG.TYPE.FLOAT:
          const numPropConfig = this.inputConfig.store[propName] as SensitivityNumericValues;
          return numPropConfig.isChanging.value === true;
        default:
          return false;
        }
      }).map((propName) => {
        const propConfig = this.inputConfig.store[propName] as SensitivityNumericValues;

        return {
          prop: propConfig.prop,
          min: propConfig.min,
          max: propConfig.max,
        };
      }),
      samplesCount: this.samplesInput.value || 1,
    };
    console.log(options);
    const analysis = new VarianceBasedSenstivityAnalysis(options.func, options.fixedInputs, options.variedInputs, options.samplesCount);
    const analysisResults = await analysis.perform();

    const funcEvalResults = analysisResults.funcEvalResults;
    const firstOrderIndeces = analysisResults.firstOrderSobolIndeces;
    const totalOrderIndeces = analysisResults.totalOrderSobolIndeces;

    const outoutNames = firstOrderIndeces.columns.names();
    const evalDataframeNames = funcEvalResults.columns.names();    

    this.comparisonView.dataFrame = funcEvalResults;

    // add correlation plot
    this.comparisonView.addViewer(DG.Viewer.correlationPlot(funcEvalResults));

    // add scatterplot    
    this.comparisonView.addViewer(DG.Viewer.scatterPlot(funcEvalResults, 
      {x: evalDataframeNames[0], 
       y: outoutNames[1]
      }
    ));

    // add barchart with 1-st order Sobol' indeces    
    this.comparisonView.addViewer(DG.Viewer.barChart(firstOrderIndeces, 
      { title: firstOrderIndeces.name, 
        split: outoutNames[0], 
        value: outoutNames[1], 
        valueAggrType: 'avg'
      }
    ));

    // add barchart with total order Sobol' indeces    
    this.comparisonView.addViewer(DG.Viewer.barChart(totalOrderIndeces, 
      { title: totalOrderIndeces.name, 
        split: outoutNames[0], 
        value: outoutNames[1], 
        valueAggrType: 'avg'
      }
    ));        
  }

  private async runSimpleAnalysis() {
    const paramValues = Object.keys(this.inputConfig.store).reduce((acc, propName) => {
      switch (this.inputConfig.store[propName].type) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
        const numPropConfig = this.inputConfig.store[propName] as SensitivityNumericValues;
        const intStep = (numPropConfig.max - numPropConfig.min) / (numPropConfig.lvl - 1);
        acc[propName] = numPropConfig.isChanging.value ?
          Array.from({length: numPropConfig.lvl}, (_, i) => Math.round(numPropConfig.min + i*intStep)) :
          [numPropConfig.const];
        break;
      case DG.TYPE.FLOAT:
        const floatPropConfig = this.inputConfig.store[propName] as SensitivityNumericValues;
        const floatStep = (floatPropConfig.max - floatPropConfig.min) / (floatPropConfig.lvl - 1);
        acc[propName] = floatPropConfig.isChanging.value ?
          Array.from({length: floatPropConfig.lvl}, (_, i) => floatPropConfig.min + i*floatStep) :
          [floatPropConfig.const];
        break;
      case DG.TYPE.BOOL:
        const boolPropConfig = this.inputConfig.store[propName] as SensitivityBoolValues;
        acc[propName] = boolPropConfig.isChanging.value ?
          [boolPropConfig.const, !boolPropConfig.const]:
          [boolPropConfig.const];
        break;
      default:
        const constPropConfig = this.inputConfig.store[propName] as SensitivityConstValues;
        acc[propName] = [constPropConfig.const];
      }

      return acc;
    }, {} as Record<string, any[]>);

    let runParams = Object.values(paramValues)[0].map((item) => [item]) as any[][];
    for (let i = 1; i < Object.values(paramValues).length; i++) {
      const values = Object.values(paramValues)[i];

      const newRunParams = [] as any[][];
      for (const accVal of runParams) {
        for (const val of values)
          newRunParams.push([...accVal, val]);
      }

      runParams = newRunParams;
    }

    const funccalls = runParams.map((runParams) => this.func.prepare(
      this.func.inputs
        .map((input, idx) => ({name: input.name, idx}))
        .reduce((acc, {name, idx}) => {
          acc[name] = runParams[idx];
          return acc;
        }, {} as Record<string, any>),
    ));

    const pi = DG.TaskBarProgressIndicator.create(`Running ${funccalls.length} function calls...`);
    const calledFuncCalls = await Promise.all(funccalls.map(async (funccall) => await funccall.call()));
    pi.close();

    const configFunc = calledFuncCalls[0].func;

    const allParamViewers = [
      ...configFunc.inputs,
      ...configFunc.outputs,
    ]
      .map((prop) => getPropViewers(prop))
      .reduce((acc, config) => {
        if (!acc[config.name])
          acc[config.name] = config.config;
        else
          acc[config.name].push(...config.config);
        return acc;
      }, {} as Record<string, Record<string, string | boolean>[]>);

    const addColumnsFromProp = (configProp: DG.Property): DG.Column[] => {
      if (configProp.propertyType === DG.TYPE.DATA_FRAME) {
        const requestedViewersConfigs = allParamViewers[configProp.name];

        const viewerColumns = requestedViewersConfigs.map((config) => {
          let columnName = configProp.caption ?? configProp.name;
          const newColumn = DG.Column.fromType(DG.TYPE.DATA_FRAME, columnName, calledFuncCalls.length);
          newColumn.init(
            (idx: number) => calledFuncCalls[idx].inputs[configProp.name] ?? calledFuncCalls[idx].outputs[configProp.name],
          );
          const unusedName = comparisonDf.columns.getUnusedName(newColumn.name);
          newColumn.name = unusedName;
          columnName = unusedName;
          newColumn.temp[VIEWER_PATH] = config;
          comparisonDf.columns.add(newColumn);

          return newColumn;
        });

        return viewerColumns;
      } else {
        let columnName = configProp.caption ?? configProp.name;
        //@ts-ignore
        const newColumn = DG.Column.fromType(configProp.propertyType, columnName, calledFuncCalls.length);
        newColumn.init(
          (idx: number) => calledFuncCalls[idx].inputs[configProp.name] ?? calledFuncCalls[idx].outputs[configProp.name],
        );
        const unusedName = comparisonDf.columns.getUnusedName(newColumn.name);
        newColumn.name = unusedName;
        columnName = unusedName;
        comparisonDf.columns.add(newColumn);
        return [newColumn];
      }
    };

    const comparisonDf = DG.DataFrame.create(calledFuncCalls.length);
    const uniqueRunNames = [] as string[];
    calledFuncCalls.forEach((run) => {
      let defaultRunName = run.options['title'] ?? `${run.func.name} - ${new Date(run.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})}`;
      let idx = 2;
      while (uniqueRunNames.includes(defaultRunName)) {
        defaultRunName = `${run.func.name} - ${new Date(run.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})} - ${idx}`;
        idx++;
      }
      uniqueRunNames.push(defaultRunName);
    });

    comparisonDf.columns.add(DG.Column.fromStrings(
      RUN_NAME_COL_LABEL,
      uniqueRunNames,
    ));
    comparisonDf.name = this.options.parentCall?.func.name ? `${this.options.parentCall?.func.name} - comparison` : `${calledFuncCalls[0].func.name} - comparison`;

    configFunc.inputs.forEach((prop) => addColumnsFromProp(prop));
    configFunc.outputs.forEach((prop) => addColumnsFromProp(prop));

    setTimeout(async () => {
      this.defaultCustomize();
    }, 0);

    this.comparisonView.dataFrame = comparisonDf;
  }

  private defaultCustomize(): void {
    this.comparisonView.grid.props.rowHeight = 180;
    this.comparisonView.grid.props.showAddNewRowIcon = false;
    this.comparisonView.grid.props.allowEdit = false;

    this.comparisonView.grid.columns.byName(RUN_NAME_COL_LABEL)!.width = 70;

    for (let i = 0; i < this.comparisonView.grid.columns.length; i++) {
      const gridCol = this.comparisonView.grid.columns.byIndex(i)!;
      if (gridCol.column?.temp[VIEWER_PATH]) {
        gridCol.width = 350;
        gridCol.cellType = 'html';
      };
    }

    this.comparisonView.grid.onCellPrepare((gc) => {
      if (gc.isColHeader || gc.isRowHeader) return;

      if (gc.tableColumn!.name === RUN_NAME_COL_LABEL)
        gc.style.textVertical = true;

      if (gc.tableColumn && gc.tableColumn.temp[VIEWER_PATH]) {
        const initialValue = gc.cell.value;

        const viewerConfig = gc.tableColumn.temp[VIEWER_PATH];
        const viewerType = viewerConfig['type'] as string;
        gc.element =
          ui.waitBox(async () => {
            const viewer = Object.values(viewerTypesMapping).includes(viewerType) ?
              DG.Viewer.fromType(viewerType, initialValue) :
              await initialValue.plot.fromType(viewerType) as DG.Viewer;

            // Workaround required since getOptions and setOptions are not symmetrical
            if (!gc.tableColumn!.temp[VIEWER_PATH]['look']) {
              viewer.setOptions(viewerConfig);
              gc.tableColumn!.temp[VIEWER_PATH] = viewer.getOptions();
            } else
              viewer.setOptions(gc.tableColumn!.temp[VIEWER_PATH]['look']);

            return viewer.root;
          });

        gc.element.style.width = '100%';
        gc.element.style.height = '100%';
      }
    });

    // Catching events to render context panel
    grok.events.onCurrentObjectChanged.subscribe(({sender}) => {
      if (
        sender instanceof DG.Column &&
        sender.type === DG.TYPE.DATA_FRAME &&
        grok.shell.tv &&
        grok.shell.tv.temp['isComparison'] &&
        [DG.VIEWER.LINE_CHART, DG.VIEWER.SCATTER_PLOT, DG.VIEWER.HISTOGRAM, DG.VIEWER.BOX_PLOT].includes(sender.temp[VIEWER_PATH]['type'])
      ) {
        grok.shell.windows.showProperties = true;

        const getAppendedDfs = (column: DG.Column) => {
          const appendedDf = column.get(0).clone() as DG.DataFrame;
          appendedDf.columns.addNew(RUN_ID_COL_LABEL, DG.TYPE.STRING).init(column.dataFrame.get(RUN_NAME_COL_LABEL, 0));

          for (let i = 1; i < column.length; i++) {
            const newRunDf = column.get(i).clone() as DG.DataFrame;
            newRunDf.columns.addNew(RUN_ID_COL_LABEL, DG.TYPE.STRING).init(column.dataFrame.get(RUN_NAME_COL_LABEL, i));

            // If one of the columns is parsed as int, it could be converted into double for proper append
            const convertibleTypes = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT] as DG.ColumnType[];
            for (let j = 0; j < newRunDf.columns.length; j++) {
              const newDfColumn = newRunDf.columns.byIndex(j);
              const appendDfColumn = appendedDf.columns.byIndex(j);

              if (
                newDfColumn.type !== appendDfColumn.type &&
                convertibleTypes.includes(newDfColumn.type) &&
                convertibleTypes.includes(appendDfColumn.type)
              ) {
                if (newDfColumn.type !== DG.COLUMN_TYPE.FLOAT) newRunDf.columns.replace(newDfColumn, newDfColumn.convertTo(DG.COLUMN_TYPE.FLOAT));
                if (appendDfColumn.type !== DG.COLUMN_TYPE.FLOAT) appendedDf.columns.replace(appendDfColumn, appendDfColumn.convertTo(DG.COLUMN_TYPE.FLOAT));
              }
            }

            appendedDf.append(newRunDf, true);
          }

          return appendedDf;
        };
        const unitedDf = sender.temp['unitedDf'] as DG.DataFrame ?? getAppendedDfs(sender);
        sender.temp['unitedDf'] = unitedDf;

        const config = sender.temp[VIEWER_PATH]['look'];

        // Avoiding cycling event emission
        setTimeout(() => {
          switch (sender.temp[VIEWER_PATH]['type']) {
          case DG.VIEWER.LINE_CHART:
            grok.shell.o = ui.box(unitedDf.plot.line({...config, 'split': RUN_ID_COL_LABEL}).root);
            break;
          case DG.VIEWER.SCATTER_PLOT:
            grok.shell.o = ui.box(unitedDf.plot.scatter({...config, 'color': RUN_ID_COL_LABEL}).root);
            break;
          case DG.VIEWER.HISTOGRAM:
            grok.shell.o = ui.box(unitedDf.plot.histogram({...config, 'split': RUN_ID_COL_LABEL}).root);
            break;
          case DG.VIEWER.BOX_PLOT:
            grok.shell.o = ui.box(unitedDf.plot.box({...config, 'category': RUN_ID_COL_LABEL}).root);
            break;
          }
        });
      }
    });
  }
}
