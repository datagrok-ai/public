/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {BehaviorSubject} from 'rxjs';
import {getDfFromRuns, getPropViewers} from './shared/utils';
import {SobolAnalysis} from './variance-based-analysis/sobol-sensitivity-analysis';
import {RandomAnalysis} from './variance-based-analysis/random-sensitivity-analysis';
import {getOutput} from './variance-based-analysis/sa-outputs-routine';
import {getCalledFuncCalls} from './variance-based-analysis/utils';
import {RunComparisonView} from './run-comparison-view';
import {combineLatest} from 'rxjs';
import '../css/sens-analysis.css';
import {CARD_VIEW_TYPE, VIEWER_PATH, viewerTypesMapping} from '../../shared-utils/consts';

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
  GRID_ANALYSIS = 'Grid',
  RANDOM_ANALYSIS = 'Monte Carlo',
  SOBOL_ANALYSIS = 'Sobol',
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
  isChanging: BehaviorSubject<boolean>,
  const: InputWithValue<boolean | number | string>,
  constForm: HTMLElement,
  saForm: DG.InputBase[],
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

type SensitivityStrStore = {
  prop: DG.Property,
  type: DG.TYPE.STRING,
} & InputValues;

type SensitivityConstStore = {
  prop: DG.Property,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL | DG.TYPE.STRING>,
} & InputValues;

type SensitivityStore = SensitivityNumericStore | SensitivityBoolStore | SensitivityConstStore | SensitivityStrStore;

const getWorkaroundElement = () => ui.div([], 'sa-switch-input');

export class SensitivityAnalysisView {
  generateInputFields = (func: DG.Func) => {
    const analysisInputs = {
      analysisType: {
        input: ui.choiceInput(
          'Method', ANALYSIS_TYPE.GRID_ANALYSIS, [ANALYSIS_TYPE.GRID_ANALYSIS, ANALYSIS_TYPE.RANDOM_ANALYSIS, ANALYSIS_TYPE.SOBOL_ANALYSIS],
          (v: ANALYSIS_TYPE) => {
            analysisInputs.analysisType.value.next(v);
            this.updateRunButtonText();
          }),
        value: new BehaviorSubject(ANALYSIS_TYPE.GRID_ANALYSIS),
      },
      samplesCount: {
        input: ui.intInput('Samples', 100, (v: number) => {
          analysisInputs.samplesCount.value = v;
          this.updateRunButtonText();
        }),
        value: 100,
      },
    } as AnalysisProps;

    analysisInputs.analysisType.input.root.insertBefore(getWorkaroundElement(), analysisInputs.analysisType.input.captionLabel);
    analysisInputs.samplesCount.input.root.insertBefore(getWorkaroundElement(), analysisInputs.samplesCount.input.captionLabel);

    const getInputValue = (input: DG.Property, key: string) => (
      input.options[key] === undefined ? input.defaultValue : Number(input.options[key])
    );

    const inputs = func.inputs.reduce((acc, inputProp) => {
      switch (inputProp.propertyType) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
      case DG.TYPE.FLOAT:
        const isChangingInputMin = (() => {
          const input = ui.switchInput(' ', false, (v: boolean) => {
            ref.isChanging.next(v);
            this.updateRunButtonText();
          });
          $(input.root).addClass('sa-switch-input');
          $(input.captionLabel).hide();
          return input;
        })();

        const isChangingInputConst = (() => {
          const input = ui.switchInput(' ', false, (v: boolean) => {
            ref.isChanging.next(v);
            this.updateRunButtonText();
          });
          $(input.root).addClass('sa-switch-input');
          $(input.captionLabel).hide();
          return input;
        })();

        const temp = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input:
            (() => {
              const inp = ui.intInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue, (v: number) => ref.const.value = v);

              inp.root.insertBefore(isChangingInputConst.root, inp.captionLabel);
              return inp;
            })(),
            value: inputProp.defaultValue,
          },
          min: {
            input:
              (() => {
                const inp = ui.floatInput(`${inputProp.caption ?? inputProp.name} min`, getInputValue(inputProp, 'min'), (v: number) => (ref as SensitivityNumericStore).min.value = v);

                (inp.input as HTMLInputElement).placeholder = 'Min';

                inp.root.insertBefore(isChangingInputMin.root, inp.captionLabel);
                return inp;
              })(),
            value: getInputValue(inputProp, 'min'),
          },
          max: {
            input: ui.floatInput(`${inputProp.caption ?? inputProp.name} max`, getInputValue(inputProp, 'max'), (v: number) => (ref as SensitivityNumericStore).max.value = v),
            value: getInputValue(inputProp, 'max'),
          },
          lvl: {
            input: ui.intInput('Samples', 3, (v: number) => {
              (ref as SensitivityNumericStore).lvl.value = v;
              this.updateRunButtonText();
            }),
            value: 3,
          },
          distrib: {
            input: ui.choiceInput('Grid', DISTRIB_TYPES[0], DISTRIB_TYPES, (v: DISTRIB_TYPE) => (ref as SensitivityNumericStore).distrib.value = v),
            value: inputProp.defaultValue,
          },
          isChanging: new BehaviorSubject<boolean>(false),
        };

        [temp.max.input, temp.lvl.input, temp.distrib.input].forEach((input) => {
          input.root.insertBefore(getWorkaroundElement(), input.captionLabel);
        });

        const simpleSa = [temp.lvl.input, temp.distrib.input];
        acc[inputProp.name] = {
          ...temp,
          constForm: temp.const.input.root,
          saForm: [
            temp.min.input,
            temp.max.input,
            ...simpleSa,
          ],
        } as SensitivityNumericStore;

        const ref = acc[inputProp.name] as SensitivityNumericStore;
        ref.isChanging.subscribe((val) => {
          isChangingInputMin.notify = false;
          isChangingInputMin.value = val;
          isChangingInputMin.notify = true;
        });
        ref.isChanging.subscribe((val) => {
          isChangingInputConst.notify = false;
          isChangingInputConst.value = val;
          isChangingInputConst.notify = true;
        });
        combineLatest([
          temp.isChanging, analysisInputs.analysisType.value,
        ]).subscribe(([isChanging, analysisType]) => {
          if (isChanging) {
            $(ref.constForm).hide();
            if (analysisType === ANALYSIS_TYPE.GRID_ANALYSIS) {
              ref.saForm.forEach((input) => $(input.root).show());
              simpleSa.forEach((input) => $(input.root).show());
            } else {
              ref.saForm.forEach((input) => $(input.root).show());
              simpleSa.forEach((input) => $(input.root).hide());
            }
          } else {
            $(ref.constForm).show();
            ref.saForm.forEach((input) => $(input.root).hide());
          }
        });
        break;

      case DG.TYPE.BOOL:
        const isChangingInputBoolConst = (() => {
          const input = ui.switchInput(' ', false, (v: boolean) => {
            boolRef.isChanging.next(v);
            this.updateRunButtonText();
          });
          $(input.root).addClass('sa-switch-input');
          $(input.captionLabel).hide();
          return input;
        })();

        const tempBool = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input: (() => {
              const temp = ui.boolInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue ?? false, (v: boolean) => boolRef.const.value = v);
              temp.root.insertBefore(isChangingInputBoolConst.root, temp.captionLabel);

              ui.tooltip.bind(isChangingInputBoolConst.root, () => {
                if (boolRef.isChanging.value === false)
                  return 'Switch to mark input as mutable';
                return 'Switch to mark input as immutable';
              });

              return temp;
            })(),
            value: false,
          } as InputWithValue<boolean>,
          changing: {
            input: ui.boolInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue ?? false, (v: boolean) => boolRef.const.value = v),
            value: false,
          } as InputWithValue<boolean>,
          isChanging: new BehaviorSubject<boolean>(false),
        };

        acc[inputProp.name] = {
          ...tempBool,
          constForm: tempBool.const.input.root,
          saForm: [tempBool.changing.input],
        } as SensitivityBoolStore;

        const boolRef = acc[inputProp.name] as SensitivityBoolStore;
        combineLatest([
          tempBool.isChanging, analysisInputs.analysisType.value,
        ]).subscribe(([isChanging, analysisType]) => {
          if (isChanging) {
            $(boolRef.constForm).hide();
            if (analysisType === ANALYSIS_TYPE.GRID_ANALYSIS) {
              boolRef.saForm.forEach((input) => $(input.root).css('visibility', 'visible'));
              simpleSa.forEach((input) => $(input.root).css('visibility', 'visible'));
            } else {
              boolRef.saForm.forEach((input) => $(input.root).css('visibility', 'visible'));
              simpleSa.forEach((input) => $(input.root).css('visibility', 'hidden'));
            }
          } else {
            $(boolRef.constForm).show();
            boolRef.saForm.forEach((input) => $(input.root).css('visibility', 'hidden'));
          }
        });
        boolRef.isChanging.subscribe((val) => {
          isChangingInputBoolConst.notify = false;
          isChangingInputBoolConst.value = val;
          isChangingInputBoolConst.notify = true;
        });
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

            input.addCaption(outputProp.caption ?? outputProp.name);

            const isInterestInput = (() => {
              const input = ui.switchInput(' ', outputProp.propertyType !== DG.TYPE.DATA_FRAME, (v: boolean) => {
                temp.isInterest.next(v);
                this.updateRunButtonText();
              });
              $(input.root).addClass('sa-switch-input');
              $(input.captionLabel).hide();
              return input;
            })();

            ui.tooltip.bind(isInterestInput.root, () => {
              if (temp.isInterest.value === false)
                return 'Switch to mark output as requiring analysis';
              return 'Switch to mark output as not requiring analysis';
            });
            // TODO: Replace by proper check
            switch (outputProp.propertyType) {
            case DG.TYPE.INT:
            case DG.TYPE.BIG_INT:
            case DG.TYPE.FLOAT:
            case DG.TYPE.BOOL:
            case DG.TYPE.DATA_FRAME:
              break;
            default:
              $(isInterestInput.root).css({'visibility': 'hidden'});
            }
            input.root.insertBefore(isInterestInput.root, input.captionLabel);

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
        isInterest: new BehaviorSubject<boolean>(outputProp.propertyType !== DG.TYPE.DATA_FRAME),
      };
      $(temp.input.input).css('visibility', 'hidden');

      temp.analysisInputs.forEach((input) => {
        input.root.insertBefore(getWorkaroundElement(), input.captionLabel);
      });

      if (temp.prop.propertyType === DG.TYPE.DATA_FRAME) {
        temp.isInterest.subscribe((isInterest) => {
          if (isInterest) {
            temp.analysisInputs.forEach((input) => $(input.root).show());
            $(temp.input.input).css('visibility', 'visible');
          } else {
            temp.analysisInputs.forEach((input) => $(input.root).hide());
            $(temp.input.input).css('visibility', 'hidden');
          }
        });
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
      }
      isInterest: BehaviorSubject<boolean>
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

    const saDock = this.comparisonView.dockManager.dock(
      form,
      DG.DOCK_TYPE.LEFT,
      null,
      `${this.func.name} - Sensitivity Analysis`,
      0.25,
    );
    saDock.container.containerElement.style.minWidth = '220px';
    saDock.container.containerElement.style.maxWidth = '390px';

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
    case ANALYSIS_TYPE.GRID_ANALYSIS:
      let product = 1;

      for (const name of Object.keys(inputs)) {
        if (inputs[name].isChanging.value) {
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
        if (inputs[name].isChanging.value)
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
    if (!this.canEvaluationBeRun())
      return '0';

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
    let prevCategory = 'Misc';
    const form = Object.values(this.store.inputs)
      .reduce(({container, switches}, inputConfig) => {
        const prop = inputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.h2(prop.category));
          prevCategory = prop.category;
        }

        container.append(
          inputConfig.constForm,
          ...inputConfig.saForm.map((input) => input.root),
        );

        return {container, switches};
      }, {
        container: ui.form([
          this.store.analysisInputs.analysisType.input,
          this.store.analysisInputs.samplesCount.input,
        ], {style: {'overflow-y': 'scroll', 'padding-right': '4px'}}),
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

        container.append(
          outputConfig.input.root,
          ...outputConfig.analysisInputs.map((input) => input.root),
        );

        return {container, switches};
      }, {
        container: form.container,
        switches: [] as DG.InputBase[],
      });

    this.store.analysisInputs.analysisType.value.subscribe((analysisType) => {
      if (analysisType === ANALYSIS_TYPE.GRID_ANALYSIS)
        $(this.store.analysisInputs.samplesCount.input.root).hide();
      else {
        $(outputsTitle).show();
        $(outputForm.container).show();
        $(this.store.analysisInputs.samplesCount.input.root).show();
      }
    });

    // make at least one output of interest
    let isAnyOutputSelectedAsOfInterest = false;

    for (const name of Object.keys(this.store.outputs)) {
      if (this.store.outputs[name].isInterest.value === true) {
        isAnyOutputSelectedAsOfInterest = true;
        break;
      }
    }

    if (!isAnyOutputSelectedAsOfInterest) {
      const firstOutput = this.store.outputs[Object.keys(this.store.outputs)[0]];
      firstOutput.isInterest.next(true);
      // firstOutput.isInterest.input.value = true;
    }

    this.updateRunButtonText();

    const buttons = ui.buttonsInput([this.runButton]);

    form.container.appendChild(
      buttons,
    );

    ui.tools.handleResize(form.container, (w: number) => {
      if (w < 320)
        $(form.container).addClass('ui-form-condensed');
      else
        $(form.container).removeClass('ui-form-condensed');
    });

    return form.container;
  }

  private addTooltips(): void {
    // type of analysis
    ui.tooltip.bind(this.store.analysisInputs.analysisType.input.root, () => {
      switch (this.store.analysisInputs.analysisType.value.value) {
      case ANALYSIS_TYPE.GRID_ANALYSIS:
        return 'Grid analysis: the function is evaluated with respect to the selected inputs varying within the specified ranges';
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
      if (!this.isAnyInputSelected())
        return 'Select mutable input(s) to run sensitivity analysis';

      if (!this.isAnyOutputSelected())
        return 'Select output(s) requiring analysis';

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
      const inpType = this.store.inputs[propName].prop.propertyType;
      if (inpType === DG.TYPE.BOOL)
        continue;

      if (inpType === DG.TYPE.STRING)
        continue;

      const propConfig = this.store.inputs[propName];

      const name = propConfig.prop.caption ?? propConfig.prop.name;

      ui.tooltip.bind(propConfig.const.input.root, 'Input value');
      ui.tooltip.bind((propConfig as SensitivityNumericStore).min.input.root, `Min & Max values of ${name}`);
      ui.tooltip.bind((propConfig as SensitivityNumericStore).max.input.root, `Min & Max values of ${name}`);
      ui.tooltip.bind((propConfig as SensitivityNumericStore).lvl.input.root, `Number of samples along the axis ${name}`);
      ui.tooltip.bind((propConfig as SensitivityNumericStore).distrib.input.root, 'Type of grid');
    }

    // switchInputs for outputs
    for (const propName of Object.keys(this.store.outputs)) {
      const propConfig = this.store.outputs[propName];

      switch (propConfig.prop.propertyType) {
      case DG.TYPE.DATA_FRAME:
        ui.tooltip.bind(propConfig.input.root, () => {
          if (propConfig.isInterest.value === false)
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
      if (this.store.inputs[propName].isChanging.value === true)
        return true;
    }
    return false;
  }

  private isAnyOutputSelected(): boolean {
    for (const propName of Object.keys(this.store.outputs)) {
      if (this.store.outputs[propName].isInterest.value === true)
        return true;
    }
    return false;
  }

  private canEvaluationBeRun(): boolean {
    return this.isAnyInputSelected() && this.isAnyOutputSelected();
  }

  private buildRunButton(): HTMLButtonElement {
    return ui.bigButton('Run', async () => {
      if (!this.canEvaluationBeRun())
        return;

      switch (this.store.analysisInputs.analysisType.value.value) {
      case ANALYSIS_TYPE.GRID_ANALYSIS:
        this.runElementaryAnalysis();
        break;
      case ANALYSIS_TYPE.RANDOM_ANALYSIS:
        this.runRandomAnalysis();
        break;
      case ANALYSIS_TYPE.SOBOL_ANALYSIS:
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
        return numPropConfig.isChanging.value === false;
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
          return numPropConfig.isChanging.value === false;
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
          return numPropConfig.isChanging.value === true;
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

      if (output.isInterest.value) {
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
    this.closeOpenedViewers();
    const funcEvalResults = analysisResults.funcEvalResults;
    const calledFuncCalls = analysisResults.funcCalls;
    const firstOrderIndeces = analysisResults.firstOrderSobolIndices;
    const totalOrderIndeces = analysisResults.totalOrderSobolIndices;
    const outputNames = firstOrderIndeces.columns.names();
    this.comparisonView.dataFrame = funcEvalResults;
    const colNamesToShow = funcEvalResults.columns.names();
    const fixedInputs = this.getFixedInputColumns(funcEvalResults.rowCount);

    // add columns with fixed inputs
    for (const col of fixedInputs)
      funcEvalResults.columns.add(col);

    const ID_COLUMN_NAME = 'ID';
    funcEvalResults.columns.add(DG.Column.fromStrings(ID_COLUMN_NAME, calledFuncCalls.map((call) => call.id)));

    // hide columns with fixed inputs
    this.comparisonView.grid.columns.setVisible([colNamesToShow[0]]); // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13450
    this.comparisonView.grid.columns.setVisible(colNamesToShow);

    // add correlation plot
    const corPlot = this.comparisonView.addViewer(DG.Viewer.correlationPlot(funcEvalResults));
    this.comparisonView.dockManager.dock(corPlot, 'right', undefined, '', 0.4);
    this.openedViewers.push(corPlot);

    const nameOfNonFixedOutput = this.getOutputNameForScatterPlot(colNamesToShow, funcEvalResults, options.variedInputs.length);

    // add other vizualizations depending on the varied inputs dimension
    if (options.variedInputs.length === 1) {
      const lineChart = this.comparisonView.addViewer(
        DG.Viewer.lineChart(DG.DataFrame.fromColumns(funcEvalResults.columns.byNames(colNamesToShow)), {
          x: colNamesToShow[0],
          markerSize: 1,
          markerType: 'gradient',
          sharex: true,
          multiAxis: true,
          multiAxisLegendPosition: 'RightCenter',
        }));
      this.openedViewers.push(lineChart);
    } else {
      const scatterPlot = this.comparisonView.addViewer(DG.Viewer.scatterPlot( funcEvalResults, {
        x: colNamesToShow[0],
        y: colNamesToShow[1],
        color: nameOfNonFixedOutput,
        size: nameOfNonFixedOutput,
        markerMaxSize: 12,
        jitterSize: 5,
      }));
      this.openedViewers.push(scatterPlot);
    }

    // add barchart with 1-st order Sobol' indices
    const bChartSobol1 = this.comparisonView.addViewer(DG.Viewer.barChart(firstOrderIndeces,
      {title: firstOrderIndeces.name,
        split: outputNames[0],
        value: nameOfNonFixedOutput, //outputNames[1],
        valueAggrType: 'avg',
      },
    ));

    this.comparisonView.dockManager.dock(bChartSobol1, 'right', undefined, '', 0.2);

    // add barchart with total order Sobol' indices
    const bChartSobolT = this.comparisonView.addViewer(DG.Viewer.barChart(totalOrderIndeces,
      {title: totalOrderIndeces.name,
        split: outputNames[0],
        value: nameOfNonFixedOutput, //outputNames[1],
        valueAggrType: 'avg',
      },
    ));

    this.openedViewers = this.openedViewers.concat([bChartSobol1, bChartSobolT]);

    this.comparisonView.grid.onCellClick.subscribe((cell: DG.GridCell) => {
      const selectedRunId = cell.tableRow?.get(ID_COLUMN_NAME);
      const selectedRun = calledFuncCalls.find((call) => call.id === selectedRunId);

      if (!selectedRun) return;

      const scalarParams = ([...selectedRun.outputParams.values()] as DG.FuncCallParam[])
        .filter((param) => DG.TYPES_SCALAR.has(param.property.propertyType));
      const scalarTable = DG.HtmlTable.create(
        scalarParams,
        (scalarVal: DG.FuncCallParam) =>
          [scalarVal.property.caption ?? scalarVal.property.name, selectedRun.outputs[scalarVal.property.name], scalarVal.property.options['units']],
      ).root;

      const dfParams = ([...selectedRun.outputParams.values()] as DG.FuncCallParam[])
        .filter((param) => param.property.propertyType === DG.TYPE.DATA_FRAME);
      const dfPanes = dfParams.reduce((acc, param) => {
        const configs = getPropViewers(param.property).config;

        const dfValue = selectedRun.outputs[param.name];
        const paneName = param.property.caption ?? param.property.name;
        configs.map((config) => {
          const viewerType = config['type'] as string;
          const viewer = DG.Viewer.fromType(viewerType, dfValue);
          viewer.setOptions(config);
          $(viewer.root).css({'width': '100%'});
          if (acc[paneName])
            acc[paneName].push(viewer.root);
          else acc[paneName] = [viewer.root];
        });

        return acc;
      }, {} as {[name: string]: HTMLElement[]});

      const overviewPanelConfig = {
        'Output scalars': [scalarTable],
        ...dfPanes,
      };
      const overviewPanel = ui.accordion();
      $(overviewPanel.root).css({'width': '100%'});
      Object.entries(overviewPanelConfig).map((e) => {
        overviewPanel.addPane(e[0], () => ui.divV(e[1]));
      });

      this.comparisonView.grid.props.rowHeight = 25;

      grok.shell.o = overviewPanel.root;
    });
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
          return numPropConfig.isChanging.value === false;
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
          return numPropConfig.isChanging.value === true;
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

      if (output.isInterest.value) {
        outputsOfInterest.push({
          prop: output.prop,
          value: {
            row: output.value.returning === DF_OPTIONS.LAST_ROW ? -1 : 1,
            columns: output.value.colName,
          }});
      }
    }

    const analysis = new RandomAnalysis(options.func, options.fixedInputs, options.variedInputs, outputsOfInterest, options.samplesCount);
    const analysiResults = await analysis.perform();
    const funcEvalResults = analysiResults.funcEvalResults;
    const calledFuncCalls = analysiResults.funcCalls;

    this.closeOpenedViewers();
    this.comparisonView.dataFrame = funcEvalResults;
    //this.comparisonView.grid.col(ID_COLUMN_NAME)!.visible = false;
    const colNamesToShow = funcEvalResults.columns.names();
    const fixedInputs = this.getFixedInputColumns(funcEvalResults.rowCount);

    // add columns with fixed inputs
    for (const col of fixedInputs)
      funcEvalResults.columns.add(col);

    const ID_COLUMN_NAME = 'ID';
    funcEvalResults.columns.add(DG.Column.fromStrings(ID_COLUMN_NAME, calledFuncCalls.map((call) => call.id)));

    // hide columns with fixed inputs
    this.comparisonView.grid.columns.setVisible([colNamesToShow[0]]); // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13450
    this.comparisonView.grid.columns.setVisible(colNamesToShow);

    // add correlation plot
    const corPlot = this.comparisonView.addViewer(DG.Viewer.correlationPlot(funcEvalResults));
    this.comparisonView.dockManager.dock(corPlot, 'right', undefined, '', 0.4);
    this.openedViewers.push(corPlot);
    this.comparisonView.grid.props.rowHeight = 25;

    const nameOfNonFixedOutput = this.getOutputNameForScatterPlot(colNamesToShow, funcEvalResults, options.variedInputs.length);

    // add other vizualizations depending on the varied inputs dimension
    if (options.variedInputs.length === 1) {
      const lineChart = this.comparisonView.addViewer(
        DG.Viewer.lineChart(DG.DataFrame.fromColumns(funcEvalResults.columns.byNames(colNamesToShow)), {
          x: colNamesToShow[0],
          markerSize: 1,
          markerType: 'gradient',
          sharex: true,
          multiAxis: true,
          multiAxisLegendPosition: 'RightCenter',
        }));
      this.openedViewers.push(lineChart);
    } else {
      const scatterPlot = this.comparisonView.addViewer(DG.Viewer.scatterPlot( funcEvalResults, {
        x: colNamesToShow[0],
        y: colNamesToShow[1],
        color: nameOfNonFixedOutput,
        size: nameOfNonFixedOutput,
        markerMaxSize: 12,
        jitterSize: 5,
      }));
      this.openedViewers.push(scatterPlot);
    }

    this.comparisonView.grid.onCellClick.subscribe((cell: DG.GridCell) => {
      const selectedRunId = cell.tableRow?.get(ID_COLUMN_NAME);
      const selectedRun = calledFuncCalls.find((call) => call.id === selectedRunId);

      if (!selectedRun) return;

      const scalarParams = ([...selectedRun.outputParams.values()] as DG.FuncCallParam[])
        .filter((param) => DG.TYPES_SCALAR.has(param.property.propertyType));
      const scalarTable = DG.HtmlTable.create(
        scalarParams,
        (scalarVal: DG.FuncCallParam) =>
          [scalarVal.property.caption ?? scalarVal.property.name, selectedRun.outputs[scalarVal.property.name], scalarVal.property.options['units']],
      ).root;

      const dfParams = ([...selectedRun.outputParams.values()] as DG.FuncCallParam[])
        .filter((param) => param.property.propertyType === DG.TYPE.DATA_FRAME);
      const dfPanes = dfParams.reduce((acc, param) => {
        const configs = getPropViewers(param.property).config;

        const dfValue = selectedRun.outputs[param.name];
        const paneName = param.property.caption ?? param.property.name;
        configs.map((config) => {
          const viewerType = config['type'] as string;
          const viewer = DG.Viewer.fromType(viewerType, dfValue);
          viewer.setOptions(config);
          $(viewer.root).css({'width': '100%'});
          if (acc[paneName])
            acc[paneName].push(viewer.root);
          else acc[paneName] = [viewer.root];
        });

        return acc;
      }, {} as {[name: string]: HTMLElement[]});

      const overviewPanelConfig = {
        'Output scalars': [scalarTable],
        ...dfPanes,
      };
      const overviewPanel = ui.accordion();
      $(overviewPanel.root).css({'width': '100%'});
      Object.entries(overviewPanelConfig).map((e) => {
        overviewPanel.addPane(e[0], () => ui.divV(e[1]));
      });

      grok.shell.o = overviewPanel.root;
    });
  }

  private async runSimpleAnalysis() {
    const paramValues = Object.keys(this.store.inputs).reduce((acc, propName) => {
      switch (this.store.inputs[propName].type) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
        const numPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
        const intStep = (numPropConfig.max.value - numPropConfig.min.value) / (numPropConfig.lvl.value - 1);
        acc[propName] = numPropConfig.isChanging.value ?
          Array.from({length: numPropConfig.lvl.value}, (_, i) => Math.round(numPropConfig.min.value + i*intStep)) :
          [numPropConfig.const.value];
        break;
      case DG.TYPE.FLOAT:
        const floatPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
        const floatStep = (floatPropConfig.max.value - floatPropConfig.min.value) / (floatPropConfig.lvl.value - 1);
        acc[propName] = floatPropConfig.isChanging.value ?
          Array.from({length: floatPropConfig.lvl.value}, (_, i) => floatPropConfig.min.value + i*floatStep) :
          [floatPropConfig.const.value];
        break;
      case DG.TYPE.BOOL:
        const boolPropConfig = this.store.inputs[propName] as SensitivityBoolStore;
        acc[propName] = boolPropConfig.isChanging.value ?
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

    this.comparisonView.grid.props.rowHeight = 25;
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

    this.comparisonView.grid.props.rowHeight = 25;
  }

  private async runElementaryAnalysis() {
    const paramValues = Object.keys(this.store.inputs).reduce((acc, propName) => {
      switch (this.store.inputs[propName].type) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
        const numPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
        const intStep = (numPropConfig.max.value - numPropConfig.min.value) / (numPropConfig.lvl.value - 1);
        acc[propName] = numPropConfig.isChanging.value ?
          Array.from({length: numPropConfig.lvl.value}, (_, i) => Math.round(numPropConfig.min.value + i*intStep)) :
          [numPropConfig.const.value];
        break;
      case DG.TYPE.FLOAT:
        const floatPropConfig = this.store.inputs[propName] as SensitivityNumericStore;
        const floatStep = (floatPropConfig.max.value - floatPropConfig.min.value) / (floatPropConfig.lvl.value - 1);
        acc[propName] = floatPropConfig.isChanging.value ?
          Array.from({length: floatPropConfig.lvl.value}, (_, i) => floatPropConfig.min.value + i*floatStep) :
          [floatPropConfig.const.value];
        break;
      case DG.TYPE.BOOL:
        const boolPropConfig = this.store.inputs[propName] as SensitivityBoolStore;
        acc[propName] = boolPropConfig.isChanging.value ?
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

    const calledFuncCalls = await getCalledFuncCalls(funccalls);

    this.closeOpenedViewers();

    const variedInputsColumns = [] as DG.Column[];
    const rowCount = calledFuncCalls.length;
    const fixedInputsColumns = this.getFixedInputColumns(rowCount);
    const colNamesToShow = [] as string[];

    for (const inputName of Object.keys(this.store.inputs)) {
      const input = this.store.inputs[inputName];
      const prop = input.prop;

      if (input.isChanging.value) {
        colNamesToShow.push(prop.caption ?? prop.name);
        variedInputsColumns.push(DG.Column.fromType(
          prop.propertyType as unknown as DG.COLUMN_TYPE,
          prop.caption ?? prop.name,
          rowCount,
        ));
      }
    }

    const ID_COLUMN_NAME = 'ID';
    const funcEvalResults = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(ID_COLUMN_NAME, calledFuncCalls.map((call) => call.id)),
      ...variedInputsColumns,
    ]);

    for (let row = 0; row < rowCount; ++row) {
      for (const inputName of Object.keys(this.store.inputs)) {
        const input = this.store.inputs[inputName];
        const prop = input.prop;

        if (input.isChanging.value)
          funcEvalResults.set(prop.caption ?? prop.name, row, calledFuncCalls[row].inputs[inputName]);
      }
    }

    const outputsOfInterest = [];

    for (const outputName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[outputName];

      if (output.isInterest.value) {
        outputsOfInterest.push({
          prop: output.prop,
          elements: [],
          row: output.value.returning === DF_OPTIONS.LAST_ROW ? -1 : 1,
        });
      }
    }

    const outputsOfInterestColumns = getOutput(calledFuncCalls, outputsOfInterest).columns;

    for (const col of outputsOfInterestColumns) {
      funcEvalResults.columns.add(col);
      colNamesToShow.push(col.name);
    }

    for (const col of fixedInputsColumns)
      funcEvalResults.columns.add(col);

    this.comparisonView.dataFrame = funcEvalResults;
    this.comparisonView.grid.col(ID_COLUMN_NAME)!.visible = false;
    this.comparisonView.grid.onCellClick.subscribe((cell: DG.GridCell) => {
      const selectedRunId = cell.tableRow?.get(ID_COLUMN_NAME);
      const selectedRun = calledFuncCalls.find((call) => call.id === selectedRunId);

      if (!selectedRun) return;

      const scalarParams = ([...selectedRun.outputParams.values()] as DG.FuncCallParam[])
        .filter((param) => DG.TYPES_SCALAR.has(param.property.propertyType));
      const scalarTable = DG.HtmlTable.create(
        scalarParams,
        (scalarVal: DG.FuncCallParam) =>
          [scalarVal.property.caption ?? scalarVal.property.name, selectedRun.outputs[scalarVal.property.name], scalarVal.property.options['units']],
      ).root;

      const dfParams = ([...selectedRun.outputParams.values()] as DG.FuncCallParam[])
        .filter((param) => param.property.propertyType === DG.TYPE.DATA_FRAME);
      const dfPanes = dfParams.reduce((acc, param) => {
        const configs = getPropViewers(param.property).config;

        const dfValue = selectedRun.outputs[param.name];
        const paneName = param.property.caption ?? param.property.name;
        configs.map((config) => {
          const viewerType = config['type'] as string;
          const viewer = DG.Viewer.fromType(viewerType, dfValue);
          viewer.setOptions(config);
          $(viewer.root).css({'width': '100%'});
          if (acc[paneName])
            acc[paneName].push(viewer.root);
          else acc[paneName] = [viewer.root];
        });

        return acc;
      }, {} as {[name: string]: HTMLElement[]});

      const overviewPanelConfig = {
        'Output scalars': [scalarTable],
        ...dfPanes,
      };
      const overviewPanel = ui.accordion();
      $(overviewPanel.root).css({'width': '100%'});
      Object.entries(overviewPanelConfig).map((e) => {
        overviewPanel.addPane(e[0], () => ui.divV(e[1]));
      });

      grok.shell.o = overviewPanel.root;
    });

    // hide columns with fixed inputs
    this.comparisonView.grid.columns.setVisible([colNamesToShow[0]]); // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13450
    this.comparisonView.grid.columns.setVisible(colNamesToShow);

    // add correlation plot
    const corPlot = this.comparisonView.addViewer(DG.Viewer.correlationPlot(funcEvalResults));
    this.comparisonView.dockManager.dock(corPlot, 'right', undefined, '', 0.4);
    this.openedViewers.push(corPlot);
    this.comparisonView.grid.props.rowHeight = 25;

    const nameOfNonFixedOutput = this.getOutputNameForScatterPlot(colNamesToShow, funcEvalResults, variedInputsColumns.length);

    console.log(nameOfNonFixedOutput);

    // add other vizualizations depending on the varied inputs dimension
    switch (variedInputsColumns.length) {
    case 1:
      const lineChart = this.comparisonView.addViewer(
        DG.Viewer.lineChart(DG.DataFrame.fromColumns(funcEvalResults.columns.byNames(colNamesToShow)), {
          x: colNamesToShow[0],
          markerSize: 1,
          markerType: 'gradient',
          sharex: true, multiAxis: true,
          multiAxisLegendPosition: 'RightCenter',
        }));
      this.openedViewers.push(lineChart);
      break;

    case 2:
      const surfacePlot = this.comparisonView.addViewer(DG.VIEWER.SURFACE_PLOT, {
        X: colNamesToShow[0], // here, captials are used due to features of surface plot
        Y: colNamesToShow[1],
        Z: nameOfNonFixedOutput,
      });
      this.openedViewers.push(surfacePlot);
      break;

    default:
      const scatterPlot = this.comparisonView.addViewer(DG.Viewer.scatterPlot(
        funcEvalResults, {
          x: colNamesToShow[0],
          y: colNamesToShow[1],
          color: nameOfNonFixedOutput,
          size: nameOfNonFixedOutput,
          markerMaxSize: 12,
          jitterSize: 5,
        }));
      this.openedViewers.push(scatterPlot);
      break;
    }
  }

  private getOutputNameForScatterPlot(names: string[], table: DG.DataFrame, start: number): string {
    for (let i = start; i < names.length; ++i) {
      const min = table.col(names[i])?.min ?? 0;
      const max = table.col(names[i])?.max ?? 0;

      if (min < max)
        return names[i];
    }

    return names[start];
  }
}
