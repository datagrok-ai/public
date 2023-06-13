/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {BehaviorSubject} from 'rxjs';
import {getDfFromRuns, getPropViewers} from './shared/utils';
import {CARD_VIEW_TYPE, FUNCTIONS_VIEW_TYPE, SCRIPTS_VIEW_TYPE, VIEWER_PATH, viewerTypesMapping} from './shared/consts';
import {VarianceBasedSenstivityAnalysis} from './variance-based-analysis/sensitivityAnalysis';
import {RunComparisonView} from './run-comparison-view';

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

enum ANALYSIS_TYPE {
  SIMPLE_ANALYSIS = 'Simple',
  VARIANCE_ANALYSIS = 'Variance-based',
}

type AnalysisProps = {
  analysisType: InputWithValue<ANALYSIS_TYPE>,
  samplesCount: InputWithValue,
}

type InputWithValue<T = number> = {input: DG.InputBase, value: T};

type InputValues = {
  isChanging: InputWithValue<BehaviorSubject<boolean>>,
  const: InputWithValue<boolean | number>,
  min: InputWithValue,
  max: InputWithValue,
  lvl: InputWithValue,
  distrib: InputWithValue<DISTRIB_TYPE>,
}

type SensitivityNumericStore = {
  prop: DG.Property,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
} & InputValues;

type SensitivityBoolStore = {
  prop: DG.Property,
  type: DG.TYPE.BOOL,
} & InputValues;

type SensitivityConstStore = {
  prop: DG.Property,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL>,
} & InputValues

type SensitivityStore = SensitivityNumericStore | SensitivityBoolStore | SensitivityConstStore;

export class SensitivityAnalysisView {
  generateForm = (func: DG.Func) => {
    const analysisInputs = {
      analysisType: {
        input: ui.choiceInput(
          'Analysis type', ANALYSIS_TYPE.SIMPLE_ANALYSIS, [ANALYSIS_TYPE.SIMPLE_ANALYSIS, ANALYSIS_TYPE.VARIANCE_ANALYSIS], (v: ANALYSIS_TYPE) => {
            analysisInputs.analysisType.value = v;

            if (v === ANALYSIS_TYPE.VARIANCE_ANALYSIS)
              $(analysisInputs.samplesCount.input.root).show();
            else
              $(analysisInputs.samplesCount.input.root).hide();
          }),
        value: ANALYSIS_TYPE.SIMPLE_ANALYSIS,
      },
      samplesCount: {
        input: ui.intInput('Samples', 100, (v: number) => analysisInputs.samplesCount.value = v),
        value: 100,
      },
    } as AnalysisProps;

    $(analysisInputs.samplesCount.input.root).hide();

    const inputs = func.inputs.reduce((acc, inputProp) => {
      switch (inputProp.propertyType) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
      case DG.TYPE.FLOAT:
        acc[inputProp.name] = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input:
              inputProp.propertyType === DG.TYPE.FLOAT ?
                ui.floatInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue, (v: number) => acc[inputProp.name].const.value = v):
                ui.intInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue, (v: number) => acc[inputProp.name].const.value = v),
            value: inputProp.defaultValue,
          },
          min: {
            input:
              inputProp.propertyType === DG.TYPE.FLOAT ?
                ui.floatInput(`${inputProp.caption ?? inputProp.name} min`, inputProp.defaultValue, (v: number) => acc[inputProp.name].min.value = v):
                ui.intInput(`${inputProp.caption ?? inputProp.name} min`, inputProp.defaultValue, (v: number) => acc[inputProp.name].min.value = v),
            value: inputProp.defaultValue,
          },
          max: {
            input:
              inputProp.propertyType === DG.TYPE.FLOAT ?
                ui.floatInput(`${inputProp.caption ?? inputProp.name} max`, inputProp.defaultValue, (v: number) => acc[inputProp.name].max.value = v):
                ui.intInput(`${inputProp.caption ?? inputProp.name} max`, inputProp.defaultValue, (v: number) => acc[inputProp.name].max.value = v),
            value: inputProp.defaultValue,
          },
          lvl: {
            input: ui.intInput('Levels', 3, (v: number) => acc[inputProp.name].lvl.value = v),
            value: 3,
          },
          distrib: {
            input: ui.choiceInput('Distribution', DISTRIB_TYPES[0], DISTRIB_TYPES, (v: DISTRIB_TYPE) => acc[inputProp.name].distrib.value = v),
            value: inputProp.defaultValue,
          },
          isChanging: {
            input: (() => {
              const input = ui.switchInput(' ', false, (v: boolean) => acc[inputProp.name].isChanging.value.next(v));
              $(input.root).css({'min-width': '50px', 'width': '50px'});
              $(input.captionLabel).css({'min-width': '0px', 'max-width': '0px'});
              return input;
            })(),
            value: new BehaviorSubject<boolean>(false),
          },
        } as SensitivityNumericStore;
        const ref = acc[inputProp.name];
        ref.isChanging.value.subscribe((newVal) => {
          if (newVal) {
            $(ref.const.input.root).hide();
            $(ref.min.input.root).show();
            $(ref.max.input.root).show();
            if (analysisInputs.analysisType.value === ANALYSIS_TYPE.SIMPLE_ANALYSIS) {
              $(ref.lvl.input.root).show();
              $(ref.distrib.input.root).show();
            } else {
              $(ref.lvl.input.root).hide();
              $(ref.distrib.input.root).hide();
            }
          } else {
            $(ref.const.input.root).show();
            $(ref.min.input.root).hide();
            $(ref.max.input.root).hide();
            $(ref.lvl.input.root).hide();
            $(ref.distrib.input.root).hide();
          }
        });
        $([ref.min.input.root, ref.max.input.root, ref.lvl.input.root, ref.distrib.input.root]).css({
          'width': '25%',
        });
        break;
      case DG.TYPE.BOOL:
        acc[inputProp.name] = {
          const: {
            input: ui.boolInput('Default value', false, (v: boolean) => acc[inputProp.name].const.value = v),
            value: false,
          } as InputWithValue<boolean>,
          isChanging: {
            input: ui.switchInput('Is changing', false, (v: boolean) => acc[inputProp.name].isChanging.value.next(v)),
            value: new BehaviorSubject<boolean>(false),
          },
          type: inputProp.propertyType,
          prop: inputProp,
        } as SensitivityBoolStore;
        break;
      default:
        acc[inputProp.name] = {
          const: {
            input: ui.input.forProperty(inputProp, undefined, {onValueChanged: (v: any) => acc[inputProp.name].const.value = v}),
            value: inputProp.defaultValue,
          },
          type: inputProp.propertyType,
          prop: inputProp,
        } as SensitivityConstStore;
      }

      return acc;
    }, {} as Record<string, SensitivityStore>);

    return {analysisInputs, inputs};
  };

  store = this.generateForm(this.func);
  comparisonView: DG.TableView;

  static async fromEmpty(
    func: DG.Func,
    options: {
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
    } = {
      parentView: undefined,
      parentCall: undefined,
    },
  ) {
    const cardView = [...grok.shell.views].find((view) => view.type === CARD_VIEW_TYPE);

    const v = await RunComparisonView.fromComparedRuns([], func,
      {
        parentView: cardView,
        parentCall: options.parentCall,
      });
    grok.shell.addView(v);

    new this(
      func,
      v,
      options,
    );
  }

  constructor(
    public func: DG.Func,
    baseView: DG.TableView,
    public options: {
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
      configFunc?: undefined,
    } = {
      parentView: undefined,
      parentCall: undefined,
      configFunc: undefined,
    },
  ) {
    const form = this.buildFormWithBtn();
    this.comparisonView = baseView;

    this.comparisonView.dockManager.dock(
      form,
      DG.DOCK_TYPE.LEFT,
      null,
      `${this.func.name} - Sensitivity Analysis`,
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
    const form = Object.values(this.store.inputs)
      .reduce(({container, switches}, inputConfig) => {
        let propInputs = [] as HTMLElement[];

        switch (inputConfig.type) {
        case DG.TYPE.INT:
        case DG.TYPE.BIG_INT:
        case DG.TYPE.FLOAT:
          inputConfig = inputConfig as SensitivityNumericStore;
          propInputs = [
            inputConfig.isChanging.input.root,
            inputConfig.const.input.root,
            inputConfig.min.input.root,
            inputConfig.max.input.root,
            inputConfig.lvl.input.root,
            inputConfig.distrib.input.root,
          ];
          inputConfig.isChanging.input.onChanged(adaptStyle);
          allPropInputs.push(...propInputs.filter((_, idx) => idx > 0));
          switches.push(inputConfig.isChanging.input);
          break;
        case DG.TYPE.BOOL:
          inputConfig = inputConfig as SensitivityBoolStore;
          propInputs = [
            inputConfig.isChanging.input.root,
            inputConfig.const.input.root,
          ];
          inputConfig.isChanging.input.onChanged(adaptStyle);
          allPropInputs.push(...propInputs.filter((_, idx) => idx > 0));
          switches.push(inputConfig.isChanging.input);
          break;
        default:
          propInputs = [
            inputConfig.const.input.root,
          ];
          allPropInputs.push(...propInputs);
        }

        const prop = inputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.h2(prop.category, {style: {'margin-top': '10px'}}));
          prevCategory = prop.category;
        }
        container.appendChild(ui.divH(propInputs, {style: {'max-width': '100%'}}));

        return {container, switches};
      }, {
        container: ui.divV([ui.divH([
          this.store.analysisInputs.analysisType.input.root,
          this.store.analysisInputs.samplesCount.input.root,
        ])], 'ui-form ui-form-wide ui-form-left'),
        switches: [] as DG.InputBase[],
      });

    const buttons = ui.buttonsInput([ui.bigButton('Run sensitivity analysis', async () => {
      if (this.store.analysisInputs.analysisType.value === ANALYSIS_TYPE.SIMPLE_ANALYSIS)
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
      fixedInputs: Object.keys(this.store.inputs).filter((propName) => {
        switch (this.store.inputs[propName].type) {
        case DG.TYPE.INT:
        case DG.TYPE.BIG_INT:
        case DG.TYPE.FLOAT:
          const numPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
          return numPropConfig.isChanging.value.value === false;
        default:
          return true;
        }
      }).map((propName) => ({
        name: propName,
        value: this.store.inputs[propName].const.value,
      })),
      variedInputs: Object.keys(this.store.inputs).filter((propName) => {
        switch (this.store.inputs[propName].type) {
        case DG.TYPE.INT:
        case DG.TYPE.BIG_INT:
        case DG.TYPE.FLOAT:
          const numPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
          return numPropConfig.isChanging.value.value === true;
        default:
          return false;
        }
      }).map((propName) => {
        const propConfig = this.store.inputs[propName] as SensitivityNumericStore;

        return {
          prop: propConfig.prop,
          min: propConfig.min.value,
          max: propConfig.max.value,
        };
      }),
      samplesCount: this.store.analysisInputs.samplesCount.value || 1,
    };

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
        y: outoutNames[1],
      },
    ));

    // add barchart with 1-st order Sobol' indeces
    this.comparisonView.addViewer(DG.Viewer.barChart(firstOrderIndeces,
      {title: firstOrderIndeces.name,
        split: outoutNames[0],
        value: outoutNames[1],
        valueAggrType: 'avg',
      },
    ));

    // add barchart with total order Sobol' indeces
    this.comparisonView.addViewer(DG.Viewer.barChart(totalOrderIndeces,
      {title: totalOrderIndeces.name,
        split: outoutNames[0],
        value: outoutNames[1],
        valueAggrType: 'avg',
      },
    ));
  }

  private async runSimpleAnalysis() {
    const paramValues = Object.keys(this.store.inputs).reduce((acc, propName) => {
      switch (this.store.inputs[propName].type) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
        const numPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
        const intStep = (numPropConfig.max.value - numPropConfig.min.value) / (numPropConfig.lvl.value - 1);
        acc[propName] = numPropConfig.isChanging.value.value ?
          Array.from({length: numPropConfig.lvl.value}, (_, i) => Math.round(numPropConfig.min.value + i*intStep)) :
          [numPropConfig.const.value];
        break;
      case DG.TYPE.FLOAT:
        const floatPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
        const floatStep = (floatPropConfig.max.value - floatPropConfig.min.value) / (floatPropConfig.lvl.value - 1);
        acc[propName] = floatPropConfig.isChanging.value.value ?
          Array.from({length: floatPropConfig.lvl.value}, (_, i) => floatPropConfig.min.value + i*floatStep) :
          [floatPropConfig.const.value];
        break;
      case DG.TYPE.BOOL:
        const boolPropConfig = this.store.inputs[propName] as SensitivityBoolStore;
        acc[propName] = boolPropConfig.isChanging.value.value ?
          [boolPropConfig.const.value, !boolPropConfig.const.value]:
          [boolPropConfig.const.value];
        break;
      default:
        const constPropConfig = this.store.inputs[propName] as SensitivityConstStore;
        acc[propName] = [constPropConfig.const.value];
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

    const comparisonDf = getDfFromRuns(
      calledFuncCalls,
      this.func,
      this.options,
    );

    this.comparisonView.dataFrame = comparisonDf;

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
  }
}
