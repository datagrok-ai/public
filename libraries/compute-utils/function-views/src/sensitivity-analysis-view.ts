/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {BehaviorSubject} from 'rxjs';
import {getDfFromRuns} from './shared/utils';
import {CARD_VIEW_TYPE, VIEWER_PATH, viewerTypesMapping} from './shared/consts';
import {SobolAnalysis} from './variance-based-analysis/sobol-sensitivity-analysis';
import {RandomAnalysis} from './variance-based-analysis/random-sensitivity-analysis';
import {RunComparisonView} from './run-comparison-view';
import {combineLatest} from 'rxjs';

const RUN_NAME_COL_LABEL = 'Run name' as const;

enum DISTRIB_TYPE {
  UNIFORM = 'Uniform',
  NORMAL = 'Normal',
  RANDOM = 'Random',
}
const DISTRIB_TYPES = [
  DISTRIB_TYPE.UNIFORM,
  // DISTRIB_TYPE.NORMAL,
  // DISTRIB_TYPE.RANDOM,
];

enum ANALYSIS_TYPE {
  SIMPLE_ANALYSIS = 'Elementary',
  RANDOM_ANALYSIS = 'Monte Carlo',
  SOBOL_ANALYSIS = 'Sobol method',
}

enum DF_OPTIONS {
  LAST_ROW = 'Last row',
  FIRST_ROW = 'First row',
  ALL_COLUMNS = 'All'
}

type AnalysisProps = {
  analysisType: InputWithValue<BehaviorSubject<ANALYSIS_TYPE>>,
  samplesCount: InputWithValue,
}

type InputWithValue<T = number> = {input: DG.InputBase, value: T};

type InputValues = {
  isChanging: InputWithValue<BehaviorSubject<boolean>>,
  const: InputWithValue<boolean | number>,
  constForm: HTMLElement,
  saForm: HTMLElement,
}

type SensitivityNumericStore = {
  prop: DG.Property,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
  min: InputWithValue,
  max: InputWithValue,
  lvl: InputWithValue,
  distrib: InputWithValue<DISTRIB_TYPE>,
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
  generateInputFields = (func: DG.Func) => {
    const analysisInputs = {
      analysisType: {
        input: ui.choiceInput(
          'Analysis type', ANALYSIS_TYPE.SIMPLE_ANALYSIS, [ANALYSIS_TYPE.SIMPLE_ANALYSIS, ANALYSIS_TYPE.RANDOM_ANALYSIS, ANALYSIS_TYPE.SOBOL_ANALYSIS],
          (v: ANALYSIS_TYPE) => {
            analysisInputs.analysisType.value.next(v);
            this.updateRunButtonText();
          }),
        value: new BehaviorSubject(ANALYSIS_TYPE.SIMPLE_ANALYSIS),
      },
      samplesCount: {
        input: ui.intInput('Samples', 100, (v: number) => {
          analysisInputs.samplesCount.value = v;
          this.updateRunButtonText();
        }),
        value: 100,
      },
    } as AnalysisProps;

    const getInputValue = (input: DG.Property, key: string) => (
      input.options[key] === undefined ? input.defaultValue : Number(input.options[key])
    );

    const inputs = func.inputs.reduce((acc, inputProp) => {
      switch (inputProp.propertyType) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
      case DG.TYPE.FLOAT:
        const temp = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input:
              inputProp.propertyType === DG.TYPE.FLOAT ?
                ui.floatInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue, (v: number) => ref.const.value = v):
                ui.intInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue, (v: number) => ref.const.value = v),
            value: inputProp.defaultValue,
          },
          min: {
            input:
              (() => {
                const inp = inputProp.propertyType === DG.TYPE.FLOAT ?
                  ui.floatInput(`${inputProp.caption ?? inputProp.name}`, getInputValue(inputProp, 'min'), (v: number) => ref.min.value = v):
                  ui.intInput(`${inputProp.caption ?? inputProp.name}`, getInputValue(inputProp, 'min'), (v: number) => ref.min.value = v);

                (inp.input as HTMLInputElement).placeholder = 'Min';
                $(inp.input).css({'padding-right': '0px', 'margin-right': '4px'});
                return inp;
              })(),
            value: getInputValue(inputProp, 'min'),
          },
          max: {
            input:
              (() => {
                const inp = inputProp.propertyType === DG.TYPE.FLOAT ?
                  ui.floatInput(` `, getInputValue(inputProp, 'max'), (v: number) => ref.max.value = v):
                  ui.intInput(` `, getInputValue(inputProp, 'max'), (v: number) => ref.max.value = v);

                (inp.input as HTMLInputElement).placeholder = 'Max';
                $(inp.input).css({'padding-left': '0px', 'margin-left': '4px'});
                return inp;
              })(),
            value: getInputValue(inputProp, 'max'),
          },
          lvl: {
            input: ui.intInput('Levels', 3, (v: number) => {
              ref.lvl.value = v;
              this.updateRunButtonText();
            }),
            value: 3,
          },
          distrib: {
            input: ui.choiceInput('Grid', DISTRIB_TYPES[0], DISTRIB_TYPES, (v: DISTRIB_TYPE) => ref.distrib.value = v),
            value: inputProp.defaultValue,
          },
          isChanging: {
            input: (() => {
              const input = ui.switchInput(' ', false, (v: boolean) => {
                ref.isChanging.value.next(v);
                this.updateRunButtonText();
              });
              $(input.root).css({'min-width': '50px', 'width': '50px'});
              $(input.captionLabel).css({'min-width': '0px', 'max-width': '0px'});
              return input;
            })(),
            value: new BehaviorSubject<boolean>(false),
          },
        };

        temp.min.input.root.append(temp.max.input.input);
        const simpleSa = [temp.lvl.input, temp.distrib.input];
        acc[inputProp.name] = {
          ...temp,
          constForm: ui.form([temp.const.input], {style: {flexGrow: '1'}}),
          saForm: ui.form([
            temp.min.input,
            ...simpleSa,
          ], {style: {flexGrow: '1'}}),
        } as SensitivityNumericStore;

        const ref = acc[inputProp.name] as SensitivityNumericStore;
        combineLatest([
          temp.isChanging.value, analysisInputs.analysisType.value,
        ]).subscribe(([isChanging, analysisType]) => {
          if (isChanging) {
            $(ref.constForm).hide();
            if (analysisType === ANALYSIS_TYPE.SIMPLE_ANALYSIS) {
              $(ref.saForm).show();
              simpleSa.forEach((input) => $(input.root).show());
            } else {
              $(ref.saForm).show();
              simpleSa.forEach((input) => $(input.root).hide());
            }
          } else {
            $(ref.constForm).show();
            $(ref.saForm).hide();
          }
        });
        break;
      case DG.TYPE.BOOL:
        const tempBool = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input: ui.boolInput(' ', false, (v: boolean) => ref.const.value = v),
            value: false,
          } as InputWithValue<boolean>,
          isChanging: {
            input: ui.switchInput('Is changing', false, (v: boolean) => ref.isChanging.value.next(v)),
            value: new BehaviorSubject<boolean>(false),
          },
        };

        acc[inputProp.name] = {
          ...tempBool,
          constForm: ui.form([tempBool.const.input], {style: {flexGrow: '1'}}),
          saForm: ui.form([tempBool.const.input], {style: {flexGrow: '1'}}),
        } as SensitivityBoolStore;
        break;
      default:
        acc[inputProp.name] = {
          const: {
            input: ui.input.forProperty(inputProp, undefined, {onValueChanged: (v: any) => ref.const.value = v}),
            value: inputProp.defaultValue,
          },
          type: inputProp.propertyType,
          prop: inputProp,
        } as SensitivityConstStore;
      }

      return acc;
    }, {} as Record<string, SensitivityStore>);

    const outputs = func.outputs.reduce((acc, outputProp) => {
      const temp = {
        prop: outputProp,
        input:
          (() => {
            const input = outputProp.propertyType === DG.TYPE.DATA_FRAME ?
              ui.choiceInput(outputProp.caption ?? outputProp.name, DF_OPTIONS.LAST_ROW, [DF_OPTIONS.LAST_ROW, DF_OPTIONS.FIRST_ROW], (v: DF_OPTIONS) => {
                temp.value.returning = v;
              }):
              ui.input.forProperty(outputProp);

          // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13004
          input.captionLabel.firstChild!.replaceWith(ui.span([outputProp.caption ?? outputProp.name]));

          return input;
          })(),
        analysisInputs:
          outputProp.propertyType === DG.TYPE.DATA_FRAME ? [(() => {
            const input = ui.stringInput('Column', DF_OPTIONS.ALL_COLUMNS, (v: string) => {
              temp.value.colName = v;
            });

            return input;
          })()]: [],
        value: {
          returning: DF_OPTIONS.LAST_ROW,
          colName: DF_OPTIONS.ALL_COLUMNS as string,
        },
        isInterest: {
          input: (() => {
            const input = ui.switchInput(' ', false, (v: boolean) => temp.isInterest.value.next(v));
            $(input.root).css({'min-width': '50px', 'width': '50px'});
            $(input.captionLabel).css({'min-width': '0px', 'max-width': '0px'});
            return input;
          })(),
          value: new BehaviorSubject<boolean>(false),
        },
      };
      $(temp.input.input).hide();

      if (temp.prop.propertyType === DG.TYPE.DATA_FRAME) {
        temp.isInterest.value.subscribe((isInterest) => {
          if (isInterest) {
            temp.analysisInputs.forEach((input) => $(input.root).show());
            $(temp.input.input).show();
          } else {
            temp.analysisInputs.forEach((input) => $(input.root).hide());
            $(temp.input.input).hide();
          }
        });
      }

      switch (outputProp.propertyType) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
      case DG.TYPE.FLOAT:
      case DG.TYPE.BOOL:
      case DG.TYPE.DATA_FRAME:
        break;
      default:
        $(temp.isInterest.input.root).css({'visibility': 'hidden'});
      }

      acc[outputProp.name] = temp;

      return acc;
    }, {} as Record<string, {
      prop: DG.Property,
      input: DG.InputBase,
      analysisInputs: DG.InputBase[],
      value: {
        returning: DF_OPTIONS,
        colName: string | null
        //identifier: string | null
      }
      isInterest: {input: DG.InputBase, value: BehaviorSubject<boolean>}
    }>);

    return {analysisInputs, inputs, outputs};
  };

  private openedViewers = [] as DG.Viewer[];
  private runButton: HTMLButtonElement;

  store = this.generateInputFields(this.func);
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
    this.runButton = this.buildRunButton();
    const form = this.buildFormWithBtn();
    this.addTooltips();
    this.comparisonView = baseView;

    this.comparisonView.dockManager.dock(
      form,
      DG.DOCK_TYPE.LEFT,
      null,
      `${this.func.name} - Sensitivity Analysis`,
      0.25,
    );

    this.comparisonView.grid.columns.byName(RUN_NAME_COL_LABEL)!.visible = false;
  }

  private closeOpenedViewers() {
    for (const v of this.openedViewers)
      v.close();

    this.openedViewers.splice(0);
  }

  private getFuncCallCount(analysisInputs: AnalysisProps, inputs: Record<string, SensitivityStore>): number {
    let variedInputsCount = 0;

    switch (analysisInputs.analysisType.value.value) {
    case ANALYSIS_TYPE.SIMPLE_ANALYSIS:
      let product = 1;

      for (const name of Object.keys(inputs)) {
        if (inputs[name].isChanging.value.value) {
          product *= (inputs[name] as SensitivityNumericStore).lvl.value;
          ++variedInputsCount;
        }
      }

      if (variedInputsCount === 0)
        return 0;

      return product;

    case ANALYSIS_TYPE.RANDOM_ANALYSIS:
      return analysisInputs.samplesCount.value;

    case ANALYSIS_TYPE.SOBOL_ANALYSIS:
      for (const name of Object.keys(inputs)) {
        if (inputs[name].isChanging.value.value)
          ++variedInputsCount;
      }

      if (variedInputsCount === 0)
        return 0;

      return (variedInputsCount + 2) * analysisInputs.samplesCount.value;

    default:
      return 0;
    }
  }

  private getFuncCallCountAsString(): string {
    const funcCallCount = this.getFuncCallCount(this.store.analysisInputs, this.store.inputs);

    if (funcCallCount < 1000)
      return String(funcCallCount);

    if (funcCallCount < 1000000)
      return String(Math.ceil(funcCallCount / 10) / 100) + 'k';

    if (funcCallCount < 1000000000)
      return String(Math.ceil(funcCallCount / 10000) / 100) + 'm';

    return String(Math.ceil(funcCallCount / 10000000) / 100) + 'b';
  }

  private updateRunButtonText(): void {
    this.runButton.textContent = `Run (${this.getFuncCallCountAsString()})`;
  }

  private buildFormWithBtn() {
    // const adaptStyle = () => {
    //   if (allPropInputs.some((input) => {
    //     return $(input).width() < 350 && $(input).width() > 0;
    //   }))
    //     $(form.container).addClass('ui-form-condensed');
    //   else
    //     $(form.container).removeClass('ui-form-condensed');
    // };

    let prevCategory = 'Misc';
    const form = Object.values(this.store.inputs)
      .reduce(({container, switches}, inputConfig) => {
        const prop = inputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.h2(prop.category));
          prevCategory = prop.category;
        }

        const inputRow = ui.divH([
          ui.form([inputConfig.isChanging.input]),
          inputConfig.constForm,
          inputConfig.saForm,
        ], 'ui-form ui-form-wide');
        $(inputRow).removeClass('ui-div');
        container.appendChild(inputRow);

        return {container, switches};
      }, {
        container: ui.form([
          this.store.analysisInputs.analysisType.input,
          this.store.analysisInputs.samplesCount.input,
        ], 'ui-form-wide'),
        switches: [] as DG.InputBase[],
      });

    const outputsTitle = ui.h2('Outputs');
    form.container.appendChild(outputsTitle);
    prevCategory = 'Misc';

    const outputForm = Object.values(this.store.outputs)
      .reduce(({container, switches}, outputConfig) => {
        const prop = outputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.h2(prop.category));
          prevCategory = prop.category;
        }

        const inputRow = ui.divH([
          ui.form([outputConfig.isInterest.input]),
          ui.form([outputConfig.input, ...outputConfig.analysisInputs], {style: {flexGrow: '1'}}),
        ], 'ui-form ui-form-wide');
        $(inputRow).removeClass('ui-div');
        container.appendChild(inputRow);

        return {container, switches};
      }, {
        container: ui.form([]),
        switches: [] as DG.InputBase[],
      });

    this.store.analysisInputs.analysisType.value.subscribe((analysisType) => {
      if (analysisType === ANALYSIS_TYPE.SIMPLE_ANALYSIS) {
        $(outputsTitle).hide();
        $(outputForm.container).hide();
        $(this.store.analysisInputs.samplesCount.input.root).hide();
      } else {
        $(outputsTitle).show();
        $(outputForm.container).show();
        $(this.store.analysisInputs.samplesCount.input.root).show();
      }
    });

    this.updateRunButtonText();

    const buttons = ui.buttonsInput([this.runButton, /*
      ui.bigButton('Run sensitivity analysis', async () => {
      this.closeOpenedViewers();

      if (this.store.analysisInputs.analysisType.value.value === ANALYSIS_TYPE.SIMPLE_ANALYSIS)
        this.runSimpleAnalysis();
      else if (this.store.analysisInputs.analysisType.value.value === ANALYSIS_TYPE.VARIANCE_ANALYSIS)
        this.runVarianceAnalysis();
      else
        this.runRandomAnalysis();
    })*/]);
    $(buttons).css({
      'align-items': 'start',
    });

    form.container.appendChild(outputForm.container);
    form.container.appendChild(
      buttons,
    );

    $(form.container).css({
      'overflow-y': 'scroll',
    });

    // add tooltips


    return form.container;
  }

  private addTooltips(): void {
    // type of analysis
    ui.tooltip.bind(this.store.analysisInputs.analysisType.input.root, () => {
      switch (this.store.analysisInputs.analysisType.value.value) {
      case ANALYSIS_TYPE.SIMPLE_ANALYSIS:
        return 'Elementary analysis: the function is evaluated with respect to the selected inputs varying within the specified ranges';
      case ANALYSIS_TYPE.RANDOM_ANALYSIS:
        return 'Monte Carlo simulation: the function is evaluated with respect to random variation of the selected inputs within the specified ranges';
      case ANALYSIS_TYPE.SOBOL_ANALYSIS:
        return 'Variance-based sensitivity analysis: the Sobol\' indices are computed';
      default:
        return 'Unknown method!';
      }
    });

    // run button
    ui.tooltip.bind(this.runButton, () => {
      return `Run sensitivity analysis: the function is evaluated ${this.getFuncCallCount(this.store.analysisInputs, this.store.inputs)} times`;
    });

    // samples count
    ui.tooltip.bind(this.store.analysisInputs.samplesCount.input.root, () => {
      switch (this.store.analysisInputs.analysisType.value.value) {
      case ANALYSIS_TYPE.RANDOM_ANALYSIS:
        return 'Input parameters sets count';
      case ANALYSIS_TYPE.SOBOL_ANALYSIS:
        return 'Sample size for the Sobol\' indices computation';
      default:
        return 'Unknown method!';
      }
    });

    // switchInputs for inputs
    for (const propName of Object.keys(this.store.inputs)) {
      const propConfig = this.store.inputs[propName];

      ui.tooltip.bind(propConfig.isChanging.input.root, () => {
        if (propConfig.isChanging.value.value === false)
          return 'Switch to mark input as mutable';
        return 'Switch to mark input as immutable';
      });

      ui.tooltip.bind(propConfig.const.input.root, 'Fixed input value');
      ui.tooltip.bind((propConfig as SensitivityNumericStore).min.input.root, 'Varied input range');
      ui.tooltip.bind((propConfig as SensitivityNumericStore).max.input.root, 'Varied input range');
      ui.tooltip.bind((propConfig as SensitivityNumericStore).lvl.input.root, 'Range values count');
      ui.tooltip.bind((propConfig as SensitivityNumericStore).distrib.input.root, 'Range values grid type');
    }

    // switchInputs for outputs
    for (const propName of Object.keys(this.store.outputs)) {
      const propConfig = this.store.outputs[propName];

      ui.tooltip.bind(propConfig.isInterest.input.root, () => {
        if (propConfig.isInterest.value.value === false)
          return 'Switch to mark output as requiring analysis';
        return 'Switch to mark output as not requiring analysis';
      });

      switch (propConfig.prop.propertyType) {
      case DG.TYPE.DATA_FRAME:
        ui.tooltip.bind(propConfig.input.root, () => {
          if (propConfig.isInterest.value.value === false)
            return 'Dataframe';
          return 'Specify dataframe part that requires analysis';
        });
        break;
      case DG.TYPE.INT:
      case DG.TYPE.FLOAT:
      case DG.TYPE.BIG_INT:
        ui.tooltip.bind(propConfig.input.root, 'Scalar');
        break;
      default:
        break;
      }
    }
  }

  private isAnyInputSelected(): boolean {
    for (const propName of Object.keys(this.store.inputs)) {
      if (this.store.inputs[propName].isChanging.value.value === true)
        return true;
    }
    return false;
  }

  private isAnyOutputSelected(): boolean {
    for (const propName of Object.keys(this.store.outputs)) {
      if (this.store.outputs[propName].isInterest.value.value === true)
        return true;
    }
    return false;
  }

  private buildRunButton(): HTMLButtonElement {
    return ui.bigButton('Run', async () => {
      this.closeOpenedViewers();

      if (!this.isAnyInputSelected()) {
        grok.shell.error('None of inputs is marked as mutable.');
        return;
      }

      switch (this.store.analysisInputs.analysisType.value.value) {
      case ANALYSIS_TYPE.SIMPLE_ANALYSIS:
        this.runSimpleAnalysis();
        break;
      case ANALYSIS_TYPE.RANDOM_ANALYSIS:
        if (!this.isAnyInputSelected()) {
          grok.shell.error('None of inputs is marked as mutable.');
          return;
        }
        this.runRandomAnalysis();
        break;
      case ANALYSIS_TYPE.SOBOL_ANALYSIS:
        if (!this.isAnyInputSelected()) {
          grok.shell.error('None of inputs is marked as mutable.');
          return;
        }
        this.runSobolAnalysis();
        break;
      default:
        break;
      }
    });
  }

  private getFixedInputColumns(rowCount: number): DG.Column [] {
    return Object.keys(this.store.inputs).filter((propName) => {
      switch (this.store.inputs[propName].type) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
      case DG.TYPE.FLOAT:
        const numPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
        return numPropConfig.isChanging.value.value === false;
      default:
        return true;
      }
    }).map((propName) => (DG.Column.fromList(
      this.store.inputs[propName].type as unknown as DG.COLUMN_TYPE,
      this.store.inputs[propName].prop.caption ?? propName,
      Array(rowCount).fill(this.store.inputs[propName].const.value),
    )));
  }

  private async runSobolAnalysis() {
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

    const outputsOfInterest = [];

    for (const outputName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[outputName];

      if (output.isInterest.value.value) {
        outputsOfInterest.push({
          prop: output.prop,
          value: {
            row: output.value.returning === DF_OPTIONS.LAST_ROW ? -1 : 1,
            columns: output.value.colName,
          }});
      }
    }

    const analysis = new SobolAnalysis(options.func, options.fixedInputs, options.variedInputs, outputsOfInterest, options.samplesCount);
    const analysisResults = await analysis.perform();

    const funcEvalResults = analysisResults.funcEvalResults;
    const firstOrderIndeces = analysisResults.firstOrderSobolIndices;
    const totalOrderIndeces = analysisResults.totalOrderSobolIndices;
    const outputNames = firstOrderIndeces.columns.names();
    const evalDataframeNames = funcEvalResults.columns.names();
    this.comparisonView.dataFrame = funcEvalResults;
    const colNamesToShow = funcEvalResults.columns.names();
    const fixedInputs = this.getFixedInputColumns(funcEvalResults.rowCount);

    // add columns with fixed inputs
    for (const col of fixedInputs)
      funcEvalResults.columns.add(col);

    // hide columns with fixed inputs
    grok.shell.tableView(funcEvalResults.name).grid.columns.setVisible(colNamesToShow);

    // add correlation plot
    const corPlot = this.comparisonView.addViewer(DG.Viewer.correlationPlot(funcEvalResults));

    this.comparisonView.dockManager.dock(corPlot, 'right', undefined, '', 0.25);

    // add scatterplot
    const scatPlot = this.comparisonView.addViewer(DG.Viewer.scatterPlot(funcEvalResults,
      {x: evalDataframeNames[0],
        y: outputNames[1],
      },
    ));

    // add barchart with 1-st order Sobol' indices
    const bChartSobol1 = this.comparisonView.addViewer(DG.Viewer.barChart(firstOrderIndeces,
      {title: firstOrderIndeces.name,
        split: outputNames[0],
        value: outputNames[1],
        valueAggrType: 'avg',
      },
    ));

    this.comparisonView.dockManager.dock(bChartSobol1, 'right', undefined, '', 0.2);

    // add barchart with total order Sobol' indices
    const bChartSobolT = this.comparisonView.addViewer(DG.Viewer.barChart(totalOrderIndeces,
      {title: totalOrderIndeces.name,
        split: outputNames[0],
        value: outputNames[1],
        valueAggrType: 'avg',
      },
    ));

    this.openedViewers = this.openedViewers.concat([corPlot, scatPlot, bChartSobol1, bChartSobolT]);
  }

  private async runRandomAnalysis() {
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

    const outputsOfInterest = [];

    for (const outputName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[outputName];

      if (output.isInterest.value.value) {
        outputsOfInterest.push({
          prop: output.prop,
          value: {
            row: output.value.returning === DF_OPTIONS.LAST_ROW ? -1 : 1,
            columns: output.value.colName,
          }});
      }
    }

    const analysis = new RandomAnalysis(options.func, options.fixedInputs, options.variedInputs, outputsOfInterest, options.samplesCount);
    const funcEvalResults = await analysis.perform();
    this.comparisonView.dataFrame = funcEvalResults;
    const colNamesToShow = funcEvalResults.columns.names();
    const fixedInputs = this.getFixedInputColumns(funcEvalResults.rowCount);

    // add columns with fixed inputs
    for (const col of fixedInputs)
      funcEvalResults.columns.add(col);

    // hide columns with fixed inputs
    grok.shell.tableView(funcEvalResults.name).grid.columns.setVisible(colNamesToShow);

    // add correlation plot
    const corPlot = this.comparisonView.addViewer(DG.Viewer.correlationPlot(funcEvalResults));

    this.comparisonView.dockManager.dock(corPlot, 'right', undefined, '', 0.4);

    const names = funcEvalResults.columns.names();

    // add scatterplot
    const scatPlot = this.comparisonView.addViewer(DG.Viewer.scatterPlot(funcEvalResults,
      {x: names[0],
        y: names[options.variedInputs.length],
      }));

    this.openedViewers = this.openedViewers.concat([corPlot, scatPlot]);
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
    this.comparisonView.grid.columns.byName(RUN_NAME_COL_LABEL)!.visible = false;

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
