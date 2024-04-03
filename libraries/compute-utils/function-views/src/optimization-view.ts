/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {BehaviorSubject} from 'rxjs';
import {getPropViewers} from './shared/utils';
import {RandomAnalysis} from './variance-based-analysis/random-sensitivity-analysis';
import {getOutput} from './variance-based-analysis/sa-outputs-routine';
import {getCalledFuncCalls} from './variance-based-analysis/utils';
import {RunComparisonView} from './run-comparison-view';
import {combineLatest} from 'rxjs';
import '../css/sens-analysis.css';
import {CARD_VIEW_TYPE} from '../../shared-utils/consts';
import {DOCK_RATIO, ROW_HEIGHT, STARTING_HELP} from './optimization/constants';

const RUN_NAME_COL_LABEL = 'Run name' as const;
const supportedInputTypes = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.BOOL, DG.TYPE.DATA_FRAME];
const supportedOutputTypes = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.BOOL, DG.TYPE.DATA_FRAME];

enum OPTIMIZATION_TYPE {
  MIN = 'min',
  MAX = 'max',
}

enum DF_OPTIONS {
  LAST_ROW = 'Last row',
  FIRST_ROW = 'First row',
  ALL_COLUMNS = '',
  BY_COL_VAL = 'By value in column',
}

type AnalysisProps = {
  analysisType: InputWithValue<BehaviorSubject<OPTIMIZATION_TYPE>>,
  samplesCount: InputWithValue,
}

type InputWithValue<T = number> = {input: DG.InputBase, value: T};

type InputValues = {
  isChanging: BehaviorSubject<boolean>,
  const: InputWithValue<boolean | number | string>,
  constForm: DG.InputBase[],
  saForm: DG.InputBase[],
}

type SensitivityNumericStore = {
  prop: DG.Property,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
  min: InputWithValue,
  max: InputWithValue,
  lvl: InputWithValue,
} & InputValues;

type SensitivityBoolStore = {
  prop: DG.Property,
  type: DG.TYPE.BOOL,
  lvl: number,
} & InputValues;

type SensitivityConstStore = {
  prop: DG.Property,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL | DG.TYPE.STRING>,
  lvl: 1,
} & InputValues;

type SensitivityStore = SensitivityNumericStore | SensitivityBoolStore | SensitivityConstStore;

const getSwitchMock = () => ui.div([], 'sa-switch-input');

const isNumericProp = (prop: DG.Property) => ((prop.propertyType === DG.TYPE.INT) || (prop.propertyType === DG.TYPE.FLOAT));

export class OptimizationView {
  generateInputFields = (func: DG.Func) => {
    const analysisInputs = {
      analysisType: {
        input: ui.choiceInput(
          'Goal', OPTIMIZATION_TYPE.MAX, [OPTIMIZATION_TYPE.MAX, OPTIMIZATION_TYPE.MIN],
          (v: OPTIMIZATION_TYPE) => {
            analysisInputs.analysisType.value.next(v);
            this.updateRunWidgetsState();
            this.setAnalysisInputTooltip();
            this.store.analysisInputs.samplesCount.input.setTooltip('Number of the function evaluations');
          }),
        value: new BehaviorSubject(OPTIMIZATION_TYPE.MAX),
      },
      samplesCount: {
        input: ui.intInput('Evaluations', 10, (v: number) => {
          analysisInputs.samplesCount.value = v;
          this.updateRunWidgetsState();
        }),
        value: 10,
      },
    } as AnalysisProps;

    analysisInputs.analysisType.input.root.insertBefore(getSwitchMock(), analysisInputs.analysisType.input.captionLabel);
    analysisInputs.samplesCount.input.root.insertBefore(getSwitchMock(), analysisInputs.samplesCount.input.captionLabel);

    const getInputValue = (input: DG.Property, key: string) => (
      input.options[key] === undefined ? input.defaultValue : Number(input.options[key])
    );

    const getSwitchElement = (defaultValue: boolean, f: (v: boolean) => any, isInput = true) => {
      const input = ui.switchInput(' ', defaultValue, f);
      $(input.root).addClass('sa-switch-input');
      $(input.captionLabel).hide();

      ui.tooltip.bind(input.root, () => {
        if (isInput) {
          return (input.value) ?
            `Switch to mark as immutable` :
            `Switch to mark as mutable`;
        } else {
          return !input.value ?
            `Switch to mark as a target` :
            `Switch to mark as not a target`;
        }
      });

      return input;
    };

    const inputs = func.inputs.reduce((acc, inputProp) => {
      switch (inputProp.propertyType) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
      case DG.TYPE.FLOAT:
        const isChangingInputMin = getSwitchElement(false, (v: boolean) => {
          ref.isChanging.next(v);
          this.updateRunWidgetsState();
        });

        const isChangingInputConst = getSwitchElement(false, (v: boolean) => {
          ref.isChanging.next(v);
          this.updateRunWidgetsState();
        });

        const temp = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input:
            (() => {
              const inp = ui.intInput(inputProp.caption ?? inputProp.name, inputProp.defaultValue, (v: number) => ref.const.value = v);
              inp.root.insertBefore(isChangingInputConst.root, inp.captionLabel);
              inp.addPostfix(inputProp.options['units']);
              return inp;
            })(),
            value: inputProp.defaultValue,
          },
          min: {
            input:
              (() => {
                const inp = ui.floatInput(`${inputProp.caption ?? inputProp.name} (min)`, getInputValue(inputProp, 'min'), (v: number) => (ref as SensitivityNumericStore).min.value = v);
                inp.root.insertBefore(isChangingInputMin.root, inp.captionLabel);
                inp.addPostfix(inputProp.options['units']);
                return inp;
              })(),
            value: getInputValue(inputProp, 'min'),
          },
          max: {
            input: (() => {
              const inp = ui.floatInput(`${inputProp.caption ?? inputProp.name} (max)`, getInputValue(inputProp, 'max'), (v: number) => (ref as SensitivityNumericStore).max.value = v);
              inp.addPostfix(inputProp.options['units']);
              return inp;
            })(),
            value: getInputValue(inputProp, 'max'),
          },
          lvl: {
            input: ui.intInput('Samples', 3, (v: number) => {
              (ref as SensitivityNumericStore).lvl.value = v;
              this.updateRunWidgetsState();
            }),
            value: 3,
          },
          isChanging: new BehaviorSubject<boolean>(false),
        };

        [temp.max.input, temp.lvl.input].forEach((input) => {
          input.root.insertBefore(getSwitchMock(), input.captionLabel);
          $(input.root).removeProp('display');
        });

        const simpleSa = [temp.lvl.input];
        acc[inputProp.name] = {
          ...temp,
          constForm: [temp.const.input],
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

          isChangingInputConst.notify = false;
          isChangingInputConst.value = val;
          isChangingInputConst.notify = true;
        });
        combineLatest([
          temp.isChanging, analysisInputs.analysisType.value,
        ]).subscribe(([isChanging, analysisType]) => {
          if (isChanging) {
            ref.constForm.forEach((input) => $(input.root).hide());
            ref.saForm.forEach((input) => $(input.root).css('display', 'flex'));
            simpleSa.forEach((input) => analysisType === OPTIMIZATION_TYPE.MIN ? $(input.root).css('display', 'flex'): $(input.root).hide());
          } else {
            ref.constForm.forEach((input) => $(input.root).css('display', 'flex'));
            ref.saForm.forEach((input) => $(input.root).hide());
          }
        });
        break;

      case DG.TYPE.BOOL:
        const isChangingInputBoolConst = getSwitchElement(false, (v: boolean) => {
          boolRef.isChanging.next(v);
          this.updateRunWidgetsState();
        });

        const tempBool = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input: (() => {
              const temp = ui.boolInput(`${inputProp.caption ?? inputProp.name}`, inputProp.defaultValue ?? false, (v: boolean) => boolRef.const.value = v);
              temp.root.insertBefore(isChangingInputBoolConst.root, temp.captionLabel);

              return temp;
            })(),
            value: false,
          } as InputWithValue<boolean>,
          isChanging: new BehaviorSubject<boolean>(false),
          lvl: 1,
        };

        acc[inputProp.name] = {
          ...tempBool,
          constForm: [tempBool.const.input],
          saForm: [],
        } as SensitivityBoolStore;
        const boolRef = acc[inputProp.name] as SensitivityBoolStore;
        boolRef.isChanging.subscribe((v) => {
          boolRef.lvl = v ? 2: 1;
        });
        break;
      default:
        const switchMock = getSwitchMock();

        const tempDefault = {
          input: (() => {
            const temp = ui.input.forProperty(inputProp, undefined, {onValueChanged: (v: DG.InputBase) => tempDefault.value = v.value});
            temp.root.insertBefore(switchMock, temp.captionLabel);

            temp.addPostfix(inputProp.options['units']);

            return temp;
          })(),
          value: inputProp.defaultValue,
        };
        acc[inputProp.name] = {
          const: tempDefault,
          constForm: [tempDefault.input],
          saForm: [] as DG.InputBase[],
          type: inputProp.propertyType,
          prop: inputProp,
          isChanging: new BehaviorSubject(false),
        } as SensitivityConstStore;
      }

      return acc;
    }, {} as Record<string, SensitivityStore>);

    const outputs = func.outputs.filter((outputProp) => isNumericProp(outputProp))
      .reduce((acc, outputProp) => {
        const temp = {
          prop: outputProp,
          input:
          (() => {
            const caption = outputProp.caption ?? outputProp.name;
            const input = ui.input.forProperty(outputProp);
            input.setTooltip('Output scalar');
            input.value = 0;

            input.addCaption(caption);

            const isInterestInput = supportedOutputTypes.includes(outputProp.propertyType) ?
              getSwitchElement(
                this.toSetSwitched,
                (v: boolean) => {
                  this.targetCount += v ? 1 : -1;
                  console.log('============================================');
                  this.store.analysisInputs.analysisType.input.root.hidden = (this.targetCount !== 1);
                  temp.isInterest.next(v);
                  temp.analysisInputs.forEach((inp) => {
                    inp.root.hidden = (temp.value.returning !== DF_OPTIONS.BY_COL_VAL);
                  });
                  this.updateRunWidgetsState();
                },
                false,
              ).root: getSwitchMock();
            input.root.insertBefore(isInterestInput, input.captionLabel);
            this.toSetSwitched = false;

            return input;
          })(),
          analysisInputs:
          outputProp.propertyType === DG.TYPE.DATA_FRAME ? [(() => {
            const input = ui.stringInput('Column', DF_OPTIONS.ALL_COLUMNS, (v: string) => {
              temp.value.colName = v;
            });
            input.root.insertBefore(getSwitchMock(), input.captionLabel);
            input.root.hidden = true;
            input.setTooltip(`Name of column of the '${outputProp.caption ?? outputProp.name}' dataframe`);
            return input;
          })(),
          (() => {
            const input = ui.floatInput('Value', 0, (v: number) => {temp.value.colValue = v;});
            input.root.insertBefore(getSwitchMock(), input.captionLabel);
            input.root.hidden = true;
            input.setTooltip(`Value specifying the '${outputProp.caption ?? outputProp.name}' dataframe row`);
            return input;
          })(),
          ]: [],

          value: {
            returning: DF_OPTIONS.LAST_ROW,
            colName: DF_OPTIONS.ALL_COLUMNS as string,
            colValue: 0,
          },
          isInterest: new BehaviorSubject<boolean>(this.targetCount === 1),
        };
        $(temp.input.input).css('visibility', 'hidden');

        acc[outputProp.name] = temp;

        return acc;
      }, {} as Record<string, {
      prop: DG.Property,
      input: DG.InputBase,
      analysisInputs: DG.InputBase[],
      value: {
        returning: DF_OPTIONS,
        colName: string,
        colValue: number,
      }
      isInterest: BehaviorSubject<boolean>
    }>);

    return {analysisInputs, inputs, outputs};
  };

  private openedViewers = [] as DG.Viewer[];
  private runButton = ui.bigButton('Run', async () => await this.runOptimization(), 'Run optimization');
  private runIcon = ui.iconFA('play', async () => await this.runOptimization(), 'Run optimization');
  private helpIcon = ui.iconFA('question', () => {
    window.open('https://datagrok.ai/help/compute.md#optimization', '_blank');
  }, 'Open help in a new tab');
  private tableDockNode: DG.DockNode | undefined;
  private helpMdNode: DG.DockNode | undefined;
  private gridSubscription: any = null;
  private targetCount = 1;
  private toSetSwitched = true;

  store = this.generateInputFields(this.func);
  comparisonView!: DG.TableView;

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
    if (!this.isOptimizationApplicable(func)) {
      grok.shell.warning('Optimization is not applicable: the function has no scalar outputs.');
      baseView.close();
      return;
    }

    this.runIcon = ui.iconFA('play', async () => await this.runOptimization(), 'Run optimization');
    this.runIcon.style.color = 'var(--green-2)';
    this.runIcon.classList.add('fas');

    this.helpIcon = ui.iconFA('question', () => {
      window.open('https://datagrok.ai/help/compute.md#optimization', '_blank');
    }, 'Open help in a new tab');

    const form = this.buildFormWithBtn();
    this.runButton.disabled = !this.canEvaluationBeRun();
    this.runIcon.hidden = this.runButton.disabled;
    this.addTooltips();
    this.comparisonView = baseView;

    const opt = this.comparisonView.dockManager.dock(
      form,
      DG.DOCK_TYPE.LEFT,
      null,
      `${this.func.name} - Sensitivity Analysis`,
      0.25,
    );
    /*saDock.container.containerElement.style.minWidth = '220px';
    saDock.container.containerElement.style.maxWidth = '390px';*/

    this.comparisonView.grid.columns.byName(RUN_NAME_COL_LABEL)!.visible = false;

    const rbnPanels = this.comparisonView.getRibbonPanels();
    rbnPanels.push([this.helpIcon, this.runIcon]);
    this.comparisonView.setRibbonPanels(rbnPanels);

    this.comparisonView.name = this.comparisonView.name.replace('comparison', 'optimization');
    this.comparisonView.helpUrl = 'https://datagrok.ai/help/compute.md#optimization';
    this.tableDockNode = this.comparisonView.dockManager.findNode(this.comparisonView.grid.root);
    const helpMD = ui.markdown(STARTING_HELP);
    helpMD.style.padding = '10px';
    helpMD.style.overflow = 'auto';
    this.helpMdNode = this.comparisonView.dockManager.dock(helpMD, DG.DOCK_TYPE.FILL, this.tableDockNode, 'About');
  }

  private closeOpenedViewers() {
    for (const v of this.openedViewers)
      v.close();

    this.openedViewers.splice(0);

    if (this.helpMdNode) {
      this.helpMdNode.detachFromParent();
      this.helpMdNode.container.destroy();
      this.helpMdNode = undefined;
    }

    if (this.gridSubscription) {
      this.gridSubscription.unsubscribe();
      this.gridSubscription = null;
    }
  }

  private isOptimizationApplicable(func: DG.Func): boolean {
    for (const output of func.outputs) {
      if (isNumericProp(output))
        return true;
    }

    return false;
  }

  private updateRunWidgetsState(): void {
    this.runButton.disabled = !this.canEvaluationBeRun();
    this.runIcon.hidden = this.runButton.disabled;
  }

  private buildFormWithBtn() {
    let prevCategory = 'Misc';

    const optSettingsDiv = ui.div([this.store.analysisInputs.samplesCount.input.root]);
    optSettingsDiv.hidden = true;

    const optSettingsIcon = ui.iconFA('cog', () => optSettingsDiv.hidden = !optSettingsDiv.hidden, 'Modify optimization settings');

    const form = Object.values(this.store.outputs)
      .reduce((container, outputConfig) => {
        const prop = outputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.p(prop.category));
          prevCategory = prop.category;
        }

        container.append(
          outputConfig.input.root,
          ...outputConfig.analysisInputs.map((input) => input.root),
        );

        return container;
      }, ui.div([
        ui.h2('Optimize'),
        ui.div([optSettingsIcon], {style: {'margin-top': '-22px', 'text-align': 'right'}}),
      ], {style: {'overflow-y': 'scroll', 'width': '100%'}}));

    form.appendChild(this.store.analysisInputs.analysisType.input.root);
    form.appendChild(optSettingsDiv);
    form.appendChild(ui.h2('with respect to'));
    prevCategory = 'Misc';

    Object.values(this.store.inputs)
      .reduce((container, inputConfig) => {
        const prop = inputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.p(prop.category));
          prevCategory = prop.category;
        }

        container.append(
          ...inputConfig.constForm.map((input) => input.root),
          ...inputConfig.saForm.map((input) => input.root),
        );

        return container;
      }, form);

    $(form).addClass('ui-form');

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

    this.updateRunWidgetsState();

    const buttons = ui.buttonsInput([this.runButton]);

    form.appendChild(
      buttons,
    );

    $(form).css({
      'padding-left': '12px',
      'overflow-y': 'scroll',
      'padding-right': '4px',
    });
    return form;
  }

  private setAnalysisInputTooltip(): void {
    let msg: string;

    switch (this.store.analysisInputs.analysisType.value.value) {
    case OPTIMIZATION_TYPE.MIN:
      msg = 'Minimize target';
      break;
    case OPTIMIZATION_TYPE.MAX:
      msg = 'Maximize target';
      break;
    default:
      msg = 'Unknown method!';
      break;
    }

    this.store.analysisInputs.analysisType.input.setTooltip(msg);
  }

  private addTooltips(): void {
    // type of analysis
    this.setAnalysisInputTooltip();

    // samples count
    this.store.analysisInputs.samplesCount.input.setTooltip('Number of the function evaluations');

    // switchInputs for inputs
    for (const propName of Object.keys(this.store.inputs)) {
      const prop = this.store.inputs[propName].prop;
      if (!isNumericProp(prop))
        continue;

      const propConfig = this.store.inputs[propName];
      const name = propConfig.prop.caption ?? propConfig.prop.name;
      propConfig.const.input.setTooltip(`'${name}' value`);
      (propConfig as SensitivityNumericStore).min.input.setTooltip(`Min value of '${name}'`);
      (propConfig as SensitivityNumericStore).max.input.setTooltip(`Max value of '${name}'`);
      (propConfig as SensitivityNumericStore).lvl.input.setTooltip(`Number of samples along the axis '${name}'`);
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

  private async runOptimization(): Promise<void> {
    if (!this.canEvaluationBeRun())
      return;

    if (this.targetCount > 1)
      await this.runMonteCarloOptimization();
    else
      await this.runMonteCarloOptimization(); // TODO: to replace it by Leonid's optimizer
  }

  private getFixedInputColumns(rowCount: number): DG.Column[] {
    return Object.values(this.store.inputs).filter((input) => {
      if (input.type === DG.TYPE.DATA_FRAME)
        return false;

      if (!supportedInputTypes.includes(input.type))
        return true;
      return !input.isChanging.value;
    }).map((input) => {
      return DG.Column.fromList(
        input.type as unknown as DG.COLUMN_TYPE,
        input.prop.caption ?? input.prop.name,
        Array(rowCount).fill(input.const.value),
      );
    });
  }

  private getFixedInputs() {
    return Object.keys(this.store.inputs).filter((propName) => {
      if (supportedInputTypes.includes(this.store.inputs[propName].type))
        return !this.store.inputs[propName].isChanging.value;

      return true;
    });
  }

  private getVariedInputs() {
    return Object.keys(this.store.inputs).filter((propName) => {
      if (supportedInputTypes.includes(this.store.inputs[propName].type))
        return this.store.inputs[propName].isChanging.value;

      return false;
    });
  }

  private async runMonteCarloOptimization() {
    const options = {
      func: this.func,
      fixedInputs: this.getFixedInputs().map((propName) => ({
        name: propName,
        value: this.store.inputs[propName].const.value,
      })),
      variedInputs: this.getVariedInputs().map((propName) => {
        const propConfig = this.store.inputs[propName] as SensitivityNumericStore;

        return {
          prop: propConfig.prop,
          min: propConfig.min.value,
          max: propConfig.max.value,
        };
      }),
      samplesCount: this.store.analysisInputs.samplesCount.value || 1,
    };

    const outputsOfInterest = this.getOutputsOfInterest();
    const analysis = new RandomAnalysis(options.func, options.fixedInputs, options.variedInputs, outputsOfInterest, options.samplesCount);
    const analysiResults = await analysis.perform();
    const funcEvalResults = analysiResults.funcEvalResults;
    const calledFuncCalls = analysiResults.funcCalls;

    this.closeOpenedViewers();
    this.comparisonView.dataFrame = funcEvalResults;
    const colNamesToShow = funcEvalResults.columns.names();
    const fixedInputs = this.getFixedInputColumns(funcEvalResults.rowCount);

    // add columns with fixed inputs & mark them as fixed
    for (const col of fixedInputs) {
      col.name = funcEvalResults.columns.getUnusedName(`${col.name} (fixed)`);
      funcEvalResults.columns.add(col);
    }

    // hide columns with fixed inputs
    this.comparisonView.grid.columns.setVisible([colNamesToShow[0]]); // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13450
    this.comparisonView.grid.columns.setVisible(colNamesToShow);
    this.comparisonView.grid.props.rowHeight = ROW_HEIGHT;

    // add PC plot
    const pcPlot = DG.Viewer.pcPlot(funcEvalResults, {columnNames: colNamesToShow});
    this.comparisonView.dockManager.dock(pcPlot, DG.DOCK_TYPE.LEFT, this.tableDockNode, '', DOCK_RATIO.PC_PLOT);
    this.openedViewers.push(pcPlot);

    const nameOfNonFixedOutput = this.getOutputNameForScatterPlot(colNamesToShow, funcEvalResults, options.variedInputs.length);

    // other vizualizations depending on the varied inputs dimension
    const graphViewer = (options.variedInputs.length === 1) ?
      DG.Viewer.lineChart(funcEvalResults, this.getLineChartOpt(colNamesToShow)) :
      DG.Viewer.scatterPlot(funcEvalResults, this.getScatterOpt(colNamesToShow, nameOfNonFixedOutput));

    this.openedViewers.push(graphViewer);
    this.comparisonView.dockManager.dock(graphViewer, DG.DOCK_TYPE.DOWN, this.tableDockNode, '', DOCK_RATIO.GRAPH);

    this.gridSubscription = this.comparisonView.grid.onCellClick.subscribe((cell: DG.GridCell) => {
      const selectedRun = calledFuncCalls[cell.tableRowIndex ?? 0];

      const scalarParams = ([...selectedRun.outputParams.values()])
        .filter((param) => DG.TYPES_SCALAR.has(param.property.propertyType));
      const scalarTable = DG.HtmlTable.create(
        scalarParams,
        (scalarVal: DG.FuncCallParam) =>
          [scalarVal.property.caption ?? scalarVal.property.name, selectedRun.outputs[scalarVal.property.name], scalarVal.property.options['units']],
      ).root;

      const dfParams = ([...selectedRun.outputParams.values()])
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

      let overviewPanelConfig: Object;
      let paneToExpandIdx: number;

      if (scalarParams.length > 0) {
        paneToExpandIdx = 1;
        overviewPanelConfig = {
          'Output scalars': [scalarTable],
          ...dfPanes,
        };
      } else {
        paneToExpandIdx = 0;
        overviewPanelConfig = {
          ...dfPanes,
        };
      }

      const overviewPanel = ui.accordion();
      $(overviewPanel.root).css({'width': '100%'});
      Object.entries(overviewPanelConfig).map((e) => {
        overviewPanel.addPane(e[0], () => ui.divV(e[1]));
      });

      overviewPanel.panes[paneToExpandIdx].expanded = true;

      grok.shell.o = overviewPanel.root;
    });
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

  private getOutputsOfInterest() {
    const outputsOfInterest = [];

    for (const outputName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[outputName];

      if (output.isInterest.value) {
        let rowVal: number | null;

        switch (output.value.returning) {
        case DF_OPTIONS.LAST_ROW:
          rowVal = -1;
          break;
        case DF_OPTIONS.FIRST_ROW:
          rowVal = 0;
          break;

        default:
          rowVal = null;
        }

        outputsOfInterest.push({
          prop: output.prop,
          value: {
            row: rowVal,
            colName: output.value.colName,
            colValue: output.value.colValue,
          }});
      }
    }

    return outputsOfInterest;
  }

  private getScatterOpt(colNamesToShow: string[], nameOfNonFixedOutput: string): Object {
    return {
      xColumnName: colNamesToShow[0],
      yColumnName: colNamesToShow[1],
      color: nameOfNonFixedOutput,
      size: nameOfNonFixedOutput,
      markerMaxSize: 12,
      jitterSize: 5,
    };
  }

  private getLineChartOpt(colNamesToShow: string[]): Object {
    return {
      xColumnName: colNamesToShow[0],
      yColumnNames: colNamesToShow.slice(1, Math.min(colNamesToShow.length, 8)),
      markerSize: 1,
      markerType: DG.MARKER_TYPE.GRADIENT,
      sharex: true,
      multiAxis: true,
      multiAxisLegendPosition: 'RightCenter',
    };
  }
}
