/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {BehaviorSubject, Subject} from 'rxjs';
import {RunComparisonView} from './run-comparison-view';
import {combineLatest} from 'rxjs';
import {take, filter, debounceTime} from 'rxjs/operators';
import '../css/sens-analysis.css';
import {CARD_VIEW_TYPE} from '../../shared-utils/consts';
import {getDefaultValue} from './shared/utils';
import {STARTING_HELP, TITLE, GRID_SIZE, METHOD, methodTooltip, LOSS, lossTooltip, FITTING_UI,
  INDICES, NAME, LOSS_FUNC_CHART_OPTS, SIZE, TIMEOUT} from './fitting/constants';

import {nelderMeadSettingsOpts} from './fitting/optimizer-nelder-mead';
import {getCategoryWidget, getShowInfoWidget, getLossFuncDf, lightenRGB,
  getScalarsGoodnessOfFitViewer, getHelpIcon, getRadarTooltip, toUseRadar, getRandomSeedSettings,
  getEarlyStoppingInputs,
  makeGetCalledFuncCall} from './fitting/fitting-utils';
import {OptimizationResult, Extremum, TargetTableOutput, ValueBoundsData, OutputTargetItem} from './fitting/optimizer-misc';
import {getLookupChoiceInput} from './shared/lookup-tools';

import {IVP, IVP2WebWorker, PipelineCreator} from 'diff-grok';
import {getFittedParams} from './fitting/diff-studio/nelder-mead';
import {getNonSimilar} from './fitting/similarity-utils';
import {ScalarsFitRadar} from './fitting/scalars-fit-radar';
import {getPropViewers} from '../../shared-utils/utils';
import {runFormula} from './fitting/formulas-resolver';
import {runOptimizer} from './fitting/optimizer-api';

const colors = DG.Color.categoricalPalette;
const colorsCount = colors.length;

const miscName = 'Misc';
const RUN_NAME_COL_LABEL = 'Run name' as const;
const SIDE_INPUT_CLASS = 'side-input';
const SIDE_ICON_CLASS = 'side-icon';
const FORM_TITLE_CLASS = 'form-title';
const FORM_SECTION_CLASS = 'form-section';

const supportedOutputTypes = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.DATA_FRAME];
type OutputTarget = number | DG.DataFrame | null;

type InputWithValue<T = number> = {input: DG.InputBase, value: T};

type InputRangeBounds = {
  input: DG.InputBase<number>,
  formula: DG.InputBase<string>,
  useFormula: DG.InputBase<boolean>,
  value: number
};

type InputValues = {
  isChanging: BehaviorSubject<boolean>,
  const: InputWithValue<boolean | number | string | DG.DataFrame>,
  isChangingInput?: DG.InputBase<boolean>,
  constForm: DG.InputBase[],
  saForm: DG.InputBase[],
}

type FittingNumericStore = {
  prop: DG.Property,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
  min: InputRangeBounds,
  max: InputRangeBounds,
} & InputValues;

type FittingBoolStore = {
  prop: DG.Property,
  type: DG.TYPE.BOOL,
} & InputValues;

type FittingConstStore = {
  prop: DG.Property,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL | DG.TYPE.STRING>,
} & InputValues;

/** Goodness of fit (GoF) viewer */
type GoFtable = {
  caption: string,
  table: DG.DataFrame,
  chart: DG.VIEWER,
  opts: Partial<DG.IBarChartSettings | DG.IScatterPlotSettings>,
};

export type RangeDescription = {
  default?: number;
  min?: number | string;
  max?: number | string;
  step?: number;
  isPrimary?: boolean;
  enabled?: boolean;
}

export type TargetDescription = {
  default?: any;
  argumentCol?: string;
  isPrimary?: boolean;
  enabled?: boolean;
}

type FittingInputsStore = FittingNumericStore | FittingBoolStore | FittingConstStore;

type FittingOutputsStore = {
  prop: DG.Property,
  input: DG.InputBase,
  isInterest: BehaviorSubject<boolean>,
  showInfoWidget: HTMLElement,
  isInterestInput: DG.InputBase<boolean>,
  target: OutputTarget,
  argName: string,
  argColInput: DG.ChoiceInput<string | null>,
  funcColsInput: DG.InputBase<DG.Column[]>,
};

type ValidationInfo = {
  isValid: boolean,
  msg: string,
};

/** Check validness of outputs */
function getValidation(store: FittingOutputsStore): ValidationInfo {
  if (store.prop.propertyType !== DG.TYPE.DATA_FRAME)
    return {isValid: true, msg: ''};

  const argName = store.argName;
  const funcColsVals = store.funcColsInput.value;

  if ((funcColsVals === null) || (funcColsVals === null)) {
    return {
      isValid: false,
      msg: `Incomplete the "${store.input.caption}" target`,
    };
  }

  if (funcColsVals.map((col) => col.name).includes(argName)) {
    return {
      isValid: false,
      msg: `Invalid ${store.input.caption}: functions should not contain argument ('${argName}')`,
    };
  }

  return {isValid: true, msg: ''};
}

export type DiffGrok = {
  ivp: IVP,
  ivpWW: IVP2WebWorker,
  pipelineCreator: PipelineCreator,
};

const isValidForFitting = (prop: DG.Property) => ((prop.propertyType === DG.TYPE.INT) || (prop.propertyType === DG.TYPE.FLOAT) || (prop.propertyType === DG.TYPE.DATA_FRAME));

export class FittingView {
  private validationCheckRequests$ = new Subject<true>();

  generateInputFields = (func: DG.Func) => {
    const getInputValue = (inputProp: DG.Property, key: keyof RangeDescription) => {
      const range = this.options.ranges?.[inputProp.name];
      // treat strings as formulas for valid fitting target types, otherwise use string as a value
      if (range?.[key] != undefined && (typeof range?.[key] !== 'string' || !isValidForFitting(inputProp)))
        return range[key];
      switch (key) {
      case 'min':
        return inputProp.min;
      case 'max':
        return inputProp.max;
      default:
        return getDefaultValue(inputProp);
      }
    };

    const getRangeFormula = (inputProp: DG.Property, key: keyof RangeDescription) => {
      const range = this.options.ranges?.[inputProp.name];
      if (range?.[key] != undefined && typeof range?.[key] === 'string')
        return range[key];

      if (key === 'min' && inputProp.options.minFormula)
        return inputProp.options.minFormula;
      if (key === 'max' && inputProp.options.maxFormula)
        return inputProp.options.maxFormula;
    };

    const getFormulaInput = (inputProp: DG.Property, key: keyof RangeDescription) => {
      const formula = getRangeFormula(inputProp, key) ?? '';
      const caption = `${inputProp.caption ?? inputProp.name} (${key})`;
      const it = this;
      const inp = ui.input.textArea(caption, {
        value: formula,
        onValueChanged() {
          it.updateApplicabilityState();
        },
      });

      inp.setTooltip(`Formula for '${caption}', variable: ${inputProp.name}`);
      return inp;
    };

    const getFormulaToggleInput = (inputProp: DG.Property, key: 'min' | 'max', inputNumber: DG.InputBase, inputFormula: DG.InputBase) => {
      const formula = getRangeFormula(inputProp, key);
      const it = this;
      const toggleInputs = (val: boolean) => {
        if (!this.allowFormulas) {
          $(boolInput.root).hide();
          $(inputFormula.root).hide();
        } else if (val) {
          $(inputNumber.root).hide();
          $(inputFormula.root).show();
        } else {
          $(inputNumber.root).show();
          $(inputFormula.root).hide();
        }
      };
      const boolInput = ui.input.bool('Use formula', {
        value: !!formula,
        onValueChanged(val) {
          toggleInputs(val);
          it.updateApplicabilityState();
        },
      });
      toggleInputs(!!formula);

      const icon = ui.iconFA('question', () => {
        const avaliableVars = this.getAvaliableVars(inputProp.name).flatMap(([, { prop }]) => [ui.divText(prop.name), ui.divText(prop.caption)]);
        const varsGrid = ui.div([ui.h3('Variable'), ui.h3('Caption'), ...avaliableVars], {style: {display: 'grid', 'gridTemplateColumns': '50% 50%'}});
        const popupHeader = `Set formula for ${key} bound for ${inputProp.caption ?? inputProp.name}\n` +
          `Avaliable variables:\n`;
        ui.showPopup(ui.div([ui.divText(popupHeader), varsGrid], {
          style: {maxWidth: '500px', maxHeight: '800px', overflowY: 'auto', userSelect: 'text', padding: '10px'}
        }), icon);
      });
      boolInput.addOptions(icon)
      return boolInput;
    };

    const getSwitchElement = (defaultValue: boolean, f: (v: boolean) => any, isInput: boolean = true) => {
      const input = ui.input.toggle(' ', {value: defaultValue, onValueChanged: (value) => f(value)});
      $(input.root).addClass('sa-switch-input');
      $(input.captionLabel).hide();

      ui.tooltip.bind(input.root, () => {
        if (isInput) {
          return (input.value) ?
            `Switch OFF fitting the input` :
            `Switch ON fitting the input`;
        } else {
          return !input.value ?
            `Mark the output as a target` :
            `Mark the output as not a target`;
        }
      });

      return input;
    };

    const inputs = func.inputs.reduce((acc, inputProp) => {
      const defaultValue = getInputValue(inputProp, 'default');

      if (inputProp.propertyType === DG.TYPE.FLOAT) {
        const isChangingInput = getSwitchElement(false, (v: boolean) => {
          ref.isChanging.next(v);
          this.updateApplicabilityState();
        });

        isChangingInput.classList.add(SIDE_INPUT_CLASS);

        const caption = inputProp.caption ?? inputProp.name;

        const constInput = (() => {
          const inp = ui.input.float(caption, {
            value: defaultValue, onValueChanged: (value) => {
              ref.const.value = value;
              this.updateApplicabilityState();
            },
          });
          inp.addPostfix(inputProp.options['units']);
          inp.setTooltip(`Value of '${caption}', variable: ${inputProp.name}`);
          inp.nullable = false;
          return inp;
        })();

        const minInput = (() => {
          const inp = ui.input.float(`${caption} (min)`, {
            value: getInputValue(inputProp, 'min'), onValueChanged: (value) => {
              (ref).min.value = value;
              this.updateApplicabilityState();
            },
          });
          inp.addValidator((s: string) => (Number(s) > ref.max.value) ? 'Greater than max' : null);
          inp.addPostfix(inputProp.options['units']);
          inp.setTooltip(`Min value '${caption}', variable: ${inputProp.name}`);
          inp.nullable = false;
          return inp as DG.InputBase<number>;
        })();

        const minFormulaInput = getFormulaInput(inputProp, 'min');
        const minUseFormulaInput = getFormulaToggleInput(inputProp, 'min', minInput, minFormulaInput);

        const maxInput = (() => {
          const inp = ui.input.float(`${caption} (max)`, {
            value: getInputValue(inputProp, 'max'), onValueChanged: (value) => {
              (ref).max.value = value;
              this.updateApplicabilityState();
            },
          });
          inp.addValidator((s: string) => (Number(s) < ref.min.value) ? 'Smaller than min' : null);
          inp.addPostfix(inputProp.options['units']);
          inp.setTooltip(`Max value of '${caption}', variable: ${inputProp.name}`);
          inp.nullable = false;
          return inp as DG.InputBase<number>;
        })();

        const maxFormulaInput = getFormulaInput(inputProp, 'max');
        const maxUseFormulaInput = getFormulaToggleInput(inputProp, 'max', maxInput, maxFormulaInput);

        const ref = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input: constInput,
            value: defaultValue,
          },
          min: {
            input: minInput,
            value: getInputValue(inputProp, 'min'),
            formula: minFormulaInput,
            useFormula: minUseFormulaInput,
          },
          max: {
            input: maxInput,
            value: getInputValue(inputProp, 'max'),
            formula: maxFormulaInput,
            useFormula: maxUseFormulaInput,
          },
          isChanging: new BehaviorSubject<boolean>(false),
          isChangingInput,
          constForm: [constInput],
          saForm: [
            minInput,
            minFormulaInput,
            minUseFormulaInput,
            maxInput,
            maxFormulaInput,
            maxUseFormulaInput,
          ],
        };

        acc[inputProp.name] = ref;

        ref.isChanging.subscribe((val) => {
          isChangingInput.notify = false;
          isChangingInput.value = val;
          isChangingInput.notify = true;
        });

        combineLatest([
          ref.isChanging,
        ]).subscribe(([isChanging]) => {
          if (isChanging) {
            ref.constForm.forEach((input) => $(input.root).hide());
            ref.saForm.forEach((input) => $(input.root).css('display', 'flex'));
            minUseFormulaInput.fireChanged();
            maxUseFormulaInput.fireChanged();
          } else {
            ref.constForm.forEach((input) => $(input.root).css('display', 'flex'));
            ref.saForm.forEach((input) => $(input.root).hide());
          }
        });
      } else {
        const tempDefault = {
          input: (() => {
            const temp = ui.input.forProperty(inputProp);
            temp.caption = inputProp.caption ?? inputProp.name;
            temp.onInput.subscribe(() => {
              tempDefault.value = temp.value;
              this.updateApplicabilityState();
            });
            temp.addPostfix(inputProp.options['units']);
            temp.nullable = false;

            return temp;
          })(),
          value: defaultValue,
        };
        tempDefault.input.value = defaultValue;
        acc[inputProp.name] = {
          const: tempDefault,
          constForm: [tempDefault.input],
          saForm: [] as DG.InputBase[],
          type: inputProp.propertyType as any,
          prop: inputProp,
          isChanging: new BehaviorSubject(false),
        };
      }

      return acc;
    }, {} as Record<string, FittingInputsStore>);

    const outputs = func.outputs.filter((prop) => supportedOutputTypes.includes(prop.propertyType))
      .reduce((acc, outputProp) => {
        const defaultValue = (this.options.targets?.[outputProp.name]?.default != null) ?
          this.options.targets?.[outputProp.name]?.default :
          (outputProp.propertyType === DG.TYPE.DATA_FRAME) ? null : getInputValue(outputProp, 'default');

        let defaultArgColName: string | null = null;
        let defaultColNames: Array<string | null> = [null];
        let defaultFuncCols: DG.Column[] | undefined = undefined;
        let defaultTable: DG.DataFrame | undefined = undefined;

        if ((defaultValue != null) && (defaultValue instanceof DG.DataFrame)) {
          defaultTable = defaultValue;
          const numericalCols = defaultTable.columns.toList().filter((col) => col.isNumerical);
          defaultColNames = numericalCols.map((col) => col.name);
          defaultArgColName = this.options.targets?.[outputProp.name]?.argumentCol ?? defaultColNames[0];
          defaultFuncCols =numericalCols.filter((col) => col.name !== defaultArgColName);
        }

        const validator = (_: string) => {
          const validation = getValidation(temp);

          if (validation.isValid) {
            temp.argColInput.input.classList.remove('d4-invalid');
            temp.funcColsInput.input.classList.remove('d4-invalid');
            ui.tooltip.bind(temp.argColInput.input, 'Independent variable');
            ui.tooltip.bind(temp.funcColsInput.input, 'Target dependent variables');

            return null;
          }

          temp.argColInput.input.classList.add('d4-invalid');
          temp.funcColsInput.input.classList.add('d4-invalid');
          ui.tooltip.bind(temp.argColInput.input, validation.msg);
          ui.tooltip.bind(temp.funcColsInput.input, validation.msg);

          return validation.msg;
        };

        const isInterest = new BehaviorSubject(false);
        const isInterestInput = getSwitchElement(
          isInterest.value,
          (v: boolean) => {
            temp.isInterest.next(v);
            this.updateApplicabilityState();
          },
          false,
        );
        isInterestInput.classList.add(SIDE_INPUT_CLASS);


        isInterest.subscribe((val) => {
          isInterestInput.notify = false;
          try {
            isInterestInput.value = val;
          } finally {
            isInterestInput.notify = true;
          }
        });

        const showInfoWidget = getShowInfoWidget(
          isInterestInput.root,
          outputProp.caption ?? outputProp.name,
        );
        showInfoWidget.classList.add(SIDE_ICON_CLASS);

        const temp: FittingOutputsStore = {
          prop: outputProp,
          input:
          (() => {
            const caption = outputProp.caption ?? outputProp.name;
            const input = ui.input.forProperty(outputProp);
            input.addCaption(caption);
            input.value = defaultValue;

            /** dataframe input icons*/
            let dfInputIcons: HTMLElement | null = null;

            if (outputProp.propertyType === DG.TYPE.DATA_FRAME) {
              input.setTooltip('Target dataframe');
              dfInputIcons = input.root.querySelector('div[class = "ui-input-options"]');
            } else
              input.setTooltip('Target scalar');

            isInterest.subscribe((val) => {
              input.input.hidden = !val;
              if (input.root.lastChild)
                (input.root.lastChild as HTMLElement).hidden = !val;
              showInfoWidget.hidden = !val;

              if (dfInputIcons != null)
                dfInputIcons.hidden = !val;
            });

            input.nullable = false;
            ui.tooltip.bind(input.captionLabel, (outputProp.propertyType === DG.TYPE.DATA_FRAME) ? 'Dataframe output of the model' : 'Scalar output of the model');

            if (this.options.targets?.[outputProp.name]?.default != null)
              setTimeout(() => input.value = this.options.targets?.[outputProp.name]?.default, 0);

            input.onChanged.subscribe((value) => {
              temp.target = input.value; // fixing the bug https://reddata.atlassian.net/browse/GROK-16642

              const clearInputs = () => {
                temp.argColInput.items = [null];
                temp.argColInput.value = null;

                ui.input.setColumnsInputTable(temp.funcColsInput, grok.data.demo.randomWalk(2, 2), (col: DG.Column) => false);
              };

              if (input.value === null)
                clearInputs();
              else if (outputProp.propertyType === DG.TYPE.DATA_FRAME) {
                if (temp.target) {
                  const colNames: string[] = [];
                  const numericalCols: DG.Column[] = [];

                  for (const col of (input.value as DG.DataFrame).columns) {
                    if (col.isNumerical) {
                      colNames.push(col.name);
                      numericalCols.push(col);
                    }
                  }

                  if (colNames.length > 1) {
                    temp.argColInput.items = colNames;
                    temp.argColInput.value = colNames[0];

                    ui.input.setColumnsInputTable(temp.funcColsInput, input.value as DG.DataFrame, (col: DG.Column) => col.isNumerical);
                    temp.funcColsInput.value = numericalCols.slice(1);
                  } else {
                    grok.shell.warning(`Not enough of numerical columns: ${colNames.length}, expected: > 1`);
                    clearInputs();
                  }
                } else
                  clearInputs();
              }

              this.updateApplicabilityState();
            });

            return input;
          })(),
          isInterestInput,
          showInfoWidget,
          argName: defaultArgColName ?? '_',
          argColInput: (() => {
            const input = ui.input.choice<string | null>('argument', {
              value: defaultArgColName,
              items: defaultColNames,
              tooltipText: 'Independent variable',
              onValueChanged: (value) => {
                if (value !== null)
                  temp.argName = value;
                else
                  temp.argName = '_';

                this.updateApplicabilityState();
              },
            });

            ui.tooltip.bind(input.captionLabel, 'Column with values of the independent variable');

            const infoIcon = ui.icons.info(() => alert('Hello!'));
            infoIcon.classList.add('sa-switch-input');
            input.nullable = false;
            input.addValidator(validator);
            isInterest.subscribe((val) => input.root.hidden = !val || (outputProp.propertyType !== DG.TYPE.DATA_FRAME));

            return input;
          })(),
          isInterest,
          target: defaultValue,
          funcColsInput: (() => {
            const input = ui.input.columns('functions', {
              table: defaultTable,
              value: defaultFuncCols,
              nullable: false,
              tooltipText: 'Target dependent variables',
              onValueChanged: (cols) => this.updateApplicabilityState(),
            });
            isInterest.subscribe((val) => input.root.hidden = !val || (outputProp.propertyType !== DG.TYPE.DATA_FRAME));
            ui.tooltip.bind(input.captionLabel, 'Columns with values of target dependent variables:');
            input.addValidator(validator);

            return input;
          })(),
        };

        acc[outputProp.name] = temp;

        return acc;
      }, {} as Record<string, FittingOutputsStore>);

    return {inputs, outputs};
  };

  private readyToRun = false;
  private isFittingRunning = false;

  private runIcon = ui.iconFA('play', async () => {
    if (this.readyToRun) {
      this.isFittingRunning = true;
      this.updateApplicabilityState();
      this.updateRunIconDisabledTooltip('In progress...');

      await this.runOptimization();

      this.isFittingRunning = false;
      this.updateApplicabilityState();
    }
  });

  private acceptIcon = ui.iconFA('ballot-check', async () => {
    const choiceItems = Array.from({length: this.currentFuncCalls.length}, (_, i) => i + 1);
    if (choiceItems.length === 0) {
      grok.shell.warning('No fittings');
      return;
    }
    let chosenItem = 1;
    const input = ui.input.choice('Select fitting', {items: choiceItems, value: chosenItem, onValueChanged: (x) => chosenItem = x});
    const confirmed = await new Promise((resolve, _reject) => {
      ui.dialog({title: 'Accept fitting'})
        .add(ui.div([input]))
        .onOK(() => resolve(true))
        .onCancel(() => resolve(false))
        .show({modal: true, fullScreen: true, width: 600, height: 200, center: true});
    });
    if (!confirmed || chosenItem < 0)
      return;

    this.acceptIcon.remove();
    const chosenCall = this.currentFuncCalls[chosenItem-1];
    this.isFittingAccepted = true;
    this.acceptedFitting$.next(chosenCall);
    this.baseView.close();
  });

  private helpIcon = getHelpIcon();

  private randInputs: ReturnType<typeof getRandomSeedSettings> = getRandomSeedSettings();
  private earlyStoppingInputs: ReturnType<typeof getEarlyStoppingInputs> = getEarlyStoppingInputs();

  private showFailsBtn = ui.button('Issues', () => {
    switch (this.method) {
    case METHOD.NELDER_MEAD:
      this.showNelderMeadFails();
      break;

    default:
      throw new Error(`Not implemented the '${this.method}' method`);
    }
  }, 'Show issues');

  private gridClickSubscription: any = null;
  private gridCellChangeSubscription: any = null;

  private failsDF: DG.DataFrame | null = null;

  private samplesCount = FITTING_UI.SAMPLES;
  private samplesCountInput = ui.input.forProperty(DG.Property.fromOptions({
    name: TITLE.SAMPLES,
    inputType: 'Int',
    defaultValue: this.samplesCount,
    min: 1,
    max: 1000,
  }));
  private loss = LOSS.RMSE;
  private lossInput = ui.input.choice(TITLE.LOSS_LOW, {value: this.loss, items: [LOSS.MAD, LOSS.RMSE], onValueChanged: (value) => {
    this.loss = value;
    this.lossInput.setTooltip(lossTooltip.get(this.loss)!);
  }});
  private similarity = FITTING_UI.SIMILARITY_DEFAULT;
  private similarityInput = ui.input.forProperty(DG.Property.fromOptions({
    name: TITLE.SIMILARITY,
    inputType: 'Float',
    defaultValue: this.similarity,
    min: FITTING_UI.SIMILARITY_MIN,
    max: FITTING_UI.SIMILARITY_MAX,
    units: '%',
  }));

  // Enable formulas, unless annotation says otherwise
  private allowFormulas = true;

  // Auxiliary dock nodes with results
  private helpDN: DG.DockNode | undefined = undefined;

  private method = METHOD.NELDER_MEAD;
  private methodInput = ui.input.choice(TITLE.METHOD, {value: this.method, items: [METHOD.NELDER_MEAD], onValueChanged: (value) => {
    this.method = value;
    this.showHideSettingInputs();
    this.methodInput.setTooltip(methodTooltip.get(this.method)!);
  }});

  // The Nelder-Mead method settings
  private nelderMeadSettings = new Map<string, number>();

  // Gradient descent settings
  private gradDescentSettings = {
    iterCount: 10,
    learningRate: 0.0001,
  };

  private defaultsOverrides: Record<string, any> = {};

  private settingsInputs = new Map<METHOD, DG.InputBase[]>();

  store = this.generateInputFields(this.func);
  comparisonView!: DG.TableView;

  private fittingSettingsIcon = ui.iconFA('cog', () => {
    const prevState = this.fittingSettingsDiv.hidden;
    this.fittingSettingsDiv.hidden = !prevState;
  }, 'Modify fitting settings');

  private fittingSettingsDiv = ui.divV([]);

  private diffGrok: DiffGrok | undefined = undefined;

  private currentFuncCalls: DG.FuncCall[] = [];
  private isFittingAccepted = false;
  public acceptedFitting$ = new Subject<DG.FuncCall | null>();

  static async fromEmpty(
    func: DG.Func,
    options: {
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
      inputsLookup?: string,
      ranges?: Record<string, RangeDescription>,
      targets?: Record<string, TargetDescription>,
      acceptMode?: boolean,
      diffGrok?: DiffGrok,
    } = {
      parentView: undefined,
      parentCall: undefined,
      inputsLookup: undefined,
      ranges: undefined,
      targets: undefined,
      acceptMode: false,
      diffGrok: undefined,
    },
  ) {
    const cardView = [...grok.shell.views].find((view) => view.type === CARD_VIEW_TYPE);

    const v = await RunComparisonView.fromComparedRuns([], func,
      {
        parentView: cardView,
        parentCall: options.parentCall,
      });
    grok.shell.addView(v);

    return new this(
      func,
      v,
      options,
    );
  }

  private constructor(
    public func: DG.Func,
    private baseView: RunComparisonView,
    public options: {
      parentView?: DG.View,
      parentCall?: DG.FuncCall,
      configFunc?: undefined,
      inputsLookup?: string,
      ranges?: Record<string, RangeDescription>,
      targets?: Record<string, TargetDescription>,
      acceptMode?: boolean,
      diffGrok?: DiffGrok,
    } = {
      parentView: undefined,
      parentCall: undefined,
      configFunc: undefined,
      inputsLookup: undefined,
      ranges: undefined,
      targets: undefined,
      acceptMode: false,
      diffGrok: undefined,
    },
  ) {
    if (!this.isOptimizationApplicable(func)) {
      grok.shell.warning('Fitting is not applicable: the function has no scalar outputs.');
      baseView.close();
      return;
    }

    nelderMeadSettingsOpts.forEach((vals, key) => this.nelderMeadSettings.set(key, vals.default));
    this.parseDefaultsOverrides();

    this.randInputs = getRandomSeedSettings(this.defaultsOverrides);
    this.earlyStoppingInputs = getEarlyStoppingInputs(this.defaultsOverrides);
    if (this.defaultsOverrides.samplesCount) {
      this.samplesCount = this.defaultsOverrides.samplesCount;
      this.samplesCountInput.value = this.samplesCount;
    }
    if (this.defaultsOverrides.similarity) {
      this.similarity = this.defaultsOverrides.similarity;
      this.similarityInput.value = this.similarity;
    }
    if (this.defaultsOverrides.allowFormulas != null)
      this.allowFormulas = this.defaultsOverrides.allowFormulas;


    grok.events.onViewRemoved.pipe(filter((v) => v.id === baseView.id), take(1)).subscribe(() => {
      if (options.acceptMode && !this.isFittingAccepted) {
        this.acceptedFitting$.next(null);
        this.isFittingAccepted = true;
      }
    });

    this.validationCheckRequests$.pipe(
      debounceTime(0),
    ).subscribe(() => {
      this.readyToRun = this.canFittingBeRun() && (!this.isFittingRunning);
      this.updateRunIconStyle();
    });

    this.buildForm(options.inputsLookup).then((form) => {
      this.comparisonView = baseView;

      this.comparisonView.dockManager.dock(
        form,
        DG.DOCK_TYPE.LEFT,
        null,
        `${this.func.name} - Fitting`,
        0.25,
      );

      this.comparisonView.grid.columns.byName(RUN_NAME_COL_LABEL)!.visible = false;

      const rbnPanels = [[this.helpIcon, this.runIcon, ...(this.options.acceptMode ? [this.acceptIcon] : [])]];
      this.comparisonView.setRibbonPanels(rbnPanels);
      this.fittingSettingsDiv.hidden = true;

      this.comparisonView.name = this.comparisonView.name.replace('comparison', 'fitting');
      this.comparisonView.helpUrl = '/help/compute/function-analysis#parameter-optimization';
      const helpMD = ui.markdown(STARTING_HELP);
      helpMD.style.padding = '10px';
      helpMD.style.overflow = 'auto';
      this.helpDN = this.comparisonView.dockManager.dock(
        helpMD,
        DG.DOCK_TYPE.FILL,
        this.comparisonView.dockManager.findNode(this.comparisonView.grid.root),
        'About',
      );
      this.methodInput.setTooltip(methodTooltip.get(this.method)!);
      ui.tooltip.bind(this.methodInput.captionLabel, 'Method for finding optimal model parameters');
      this.lossInput.setTooltip(lossTooltip.get(this.loss)!);
      ui.tooltip.bind(this.lossInput.captionLabel, 'Method for measuring the difference between model outputs and target values');

      this.showFailsBtn.style.padding = '0px';

      this.samplesCountInput.addCaption(TITLE.SAMPLES);
      this.samplesCountInput.onChanged.subscribe((value) => {
        this.samplesCount = value;
        this.updateApplicabilityState();
      });
      this.samplesCountInput.setTooltip('Number of initial points used to start the optimization process. Higher value = higher chance of finding more distinct points');

      this.similarityInput.onChanged.subscribe((value) => {
        this.similarity = value;
        this.updateApplicabilityState();
      });
      this.similarityInput.addCaption(TITLE.SIMILARITY);
      this.similarityInput.setTooltip('Maximum relative difference between points to consider them identical. Lower = more unique results');
      ui.tooltip.bind(this.similarityInput.captionLabel, `Maximum relative deviation (%) between similar fitted points`);

      this.updateRunIconStyle();
      this.updateRunIconDisabledTooltip('Select inputs for fitting');
      this.runIcon.classList.add('fas');

      this.updateApplicabilityState();
      //this.processOptionTargets();
    });

    this.diffGrok = options.diffGrok;
  } // constructor

  private parseDefaultsOverrides() {
    this.defaultsOverrides = JSON.parse(this.func.options['fittingSettings'] || '{}');
  }

  /** Check fiiting applicability to the function */
  private isOptimizationApplicable(func: DG.Func): boolean {
    for (const output of func.outputs) {
      if (isValidForFitting(output))
        return true;
    }

    return false;
  } // isOptimizationApplicable

  /** Generate UI for the Nelder-Mead method settings */
  private generateNelderMeadSettingsInputs(): void {
    const inputs: DG.InputBase[] = [];

    let needsInitalUpdate = false;

    nelderMeadSettingsOpts.forEach((opts, key) => {
      const inp = ui.input.forProperty(DG.Property.fromOptions({
        name: opts.caption,
        inputType: opts.inputType,
        defaultValue: this.defaultsOverrides[key] ?? opts.default,
        min: opts.min,
        max: opts.max,
      }));

      inp.setTooltip(opts.tooltipText);

      inp.addCaption(opts.caption);
      inp.nullable = false;

      if (this.defaultsOverrides[key] != null) {
        this.nelderMeadSettings.set(key, this.defaultsOverrides[key]);
        needsInitalUpdate = true;
      }

      inp.onChanged.subscribe((value) => {
        this.nelderMeadSettings.set(key, value);
        this.updateApplicabilityState();
      });

      inputs.push(inp);
    });

    this.settingsInputs.set(METHOD.NELDER_MEAD, inputs);
    if (needsInitalUpdate)
      this.updateApplicabilityState();
  } // generateNelderMeadSettingsInputs

  /** Check correctness of the Nelder-Mead settings */
  private areNelderMeadSettingsCorrect(): boolean {
    for (const [key, val] of this.nelderMeadSettings) {
      if ((val === null) || (val === undefined) || (val > nelderMeadSettingsOpts.get(key)!.max) || (val < nelderMeadSettingsOpts.get(key)!.min)) {
        this.updateRunIconDisabledTooltip(`Invalid "${key}": check method settings`);
        return false;
      }
    }

    return true;
  }

  /** Generate UI for the gradient descent method */
  private generateGradDescentSettingsInputs(): void {
    const iterInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'iterations',
      inputType: 'Int',
      defaultValue: this.defaultsOverrides.iterCount ?? this.gradDescentSettings.iterCount,
      min: 1,
      max: 10000,
    }));
    iterInp.onChanged.subscribe((value) => this.gradDescentSettings.iterCount = value);

    const learningRateInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'Learning rate',
      inputType: 'Float',
      defaultValue: this.defaultsOverrides.learningRate ?? this.gradDescentSettings.learningRate,
      min: 1e-6,
      max: 1000,
    }));
    learningRateInp.onChanged.subscribe((value) => this.gradDescentSettings.learningRate = value);

    this.settingsInputs.set(METHOD.GRAD_DESC, [iterInp, learningRateInp]);
  } // generateGradDescentSettingsInputs

  /** Check correctness of the gradient descent settings */
  private areGradDescentSettingsCorrect(): boolean {
    return false;
  }

  /** Check correctness of the method settings */
  private areMethodSettingsCorrect(): boolean {
    switch (this.method) {
    case METHOD.NELDER_MEAD:
      return this.areNelderMeadSettingsCorrect();

    case METHOD.GRAD_DESC:
      return this.areGradDescentSettingsCorrect();

    default:
      return true;
    }
  } // areMethodSettingsCorrect

  /** Create UI for each method */
  private generateSettingInputs(): void {
    this.generateNelderMeadSettingsInputs();
    this.generateGradDescentSettingsInputs();
  }

  /** Show settings UI of the current method */
  private showHideSettingInputs(): void {
    this.settingsInputs.forEach((inputsArray, method) => inputsArray.forEach((input) => input.root.hidden = method !== this.method));
  }

  /** Return value lookup widget */
  private async getLookupElement(inputsLookup?: string) {
    if (inputsLookup === undefined)
      return null;

    const constIputs = new Map<string, DG.InputBase>();
    Object.keys(this.store.inputs).forEach((name) => constIputs.set(name, this.store.inputs[name].constForm[0]));

    const lookupElement = await getLookupChoiceInput(inputsLookup, constIputs);

    return lookupElement;
  }

  private addFormInputs(
    form: HTMLElement,
    inputsByCategories: Map<string, HTMLElement[]>,
    topCategory: string | null,
    expandHandler: (r: HTMLElement, isExpanded: boolean, category: string) => void,
  ) {
    const newCategoriesElements = new Map<string, HTMLElement>();
    if (inputsByCategories.size > 1) {
      if (topCategory !== null) {
        const roots = inputsByCategories.get(topCategory);
        const catEl = getCategoryWidget(topCategory, roots!, expandHandler);
        catEl.classList.add(FORM_SECTION_CLASS);
        newCategoriesElements.set(topCategory, catEl);
        form.append(catEl);
        form.append(...roots!);
      }

      inputsByCategories.forEach((roots, category) => {
        if ((category !== miscName) && (category !== topCategory)) {
          const catEl = getCategoryWidget(category, roots, expandHandler);
          catEl.classList.add(FORM_SECTION_CLASS);
          newCategoriesElements.set(category, catEl);
          form.append(catEl);
          form.append(...roots);
        }
      });

      if (topCategory !== miscName) {
        const miscRoots = inputsByCategories.get(miscName);

        if (miscRoots!.length > 0) {
          const roots = inputsByCategories.get(miscName)!;
          const catEl = getCategoryWidget(miscName, roots, expandHandler);
          catEl.classList.add(FORM_SECTION_CLASS);
          newCategoriesElements.set(miscName, catEl);
          form.append(catEl);
          form.append(...roots);
        }
      }
    } else
      form.append(...inputsByCategories.get(miscName)!);
    return newCategoriesElements;
  }

  /** Build form with inputs */
  private async buildForm(inputsLookup?: string) {
    // the main form
    const fitHeader = ui.h1(TITLE.FIT);
    fitHeader.classList.add(FORM_TITLE_CLASS);

    ui.tooltip.bind(fitHeader, 'Select inputs to be fitted');
    const form = ui.form([], {style: {padding: '5px'}});

    form.append(fitHeader);

    //0. Handling primary/secondary inputs outputs, but only if at
    //least one is set as primray
    const primaryInputRoots = new Set<HTMLElement>();

    // inputs/outputs grouped by categories
    const inputsByCategories = new Map<string, HTMLElement[]>([[miscName, []]]);
    const outputsByCategories = new Map<string, HTMLElement[]>([[miscName, []]]);
    const collapsedInputCategories = new Set<string>();
    const collapsedOutputCategories = new Set<string>();

    const primaryToggle = ui.input.toggle('Show only primary', {
      value: false, onValueChanged: (showPrimaryOnly) => {
        if (primaryInputRoots.size === 0)
          return;
        const iterPayload = [
          [inputsByCategories, inputCategoriesByName, collapsedInputCategories],
          [outputsByCategories, outputCategoriesByName, collapsedOutputCategories],
        ] as const;
        for (const [rootsMap, catMap, collapsedCats] of iterPayload) {
          for (const [catName, roots] of rootsMap.entries()) {
            const hideCat = showPrimaryOnly && !roots.some((root) => primaryInputRoots.has(root));
            const catEl = catMap.get(catName);
            if (catEl)
              catEl.hidden = hideCat;
            for (const root of roots)
              root.hidden = collapsedCats.has(catName) || (showPrimaryOnly && !primaryInputRoots.has(root));
          }
        }
      },
    });

    primaryToggle.root.hidden = true;
    form.append(primaryToggle.root);

    //1. Inputs of the function

    // group inputs by categories
    Object.values(this.store.inputs).forEach((inputConfig) => {
      const category = inputConfig.prop.category;
      const propName = inputConfig.prop.name;
      const roots = [
        ...(inputConfig.isChangingInput ? [inputConfig.isChangingInput.root] : []),
        ...inputConfig.constForm.map((input) => input.root),
        ...inputConfig.saForm.map((input) => input.root),
      ];
      const isPrimary = this.options.ranges?.[propName]?.isPrimary;
      if (isPrimary)
        roots.forEach((item) => primaryInputRoots.add(item));

      if (inputsByCategories.has(category))
        inputsByCategories.get(category)!.push(...roots);
      else
        inputsByCategories.set(category, roots);
    });

    const lookupElement = await this.getLookupElement(inputsLookup);
    let topCategory: string | null = null;

    if (lookupElement !== null) {
      const inputs = inputsByCategories.get(lookupElement.category);
      topCategory = lookupElement.category;

      if (inputs !== undefined)
        inputsByCategories.set(topCategory, [lookupElement.input.root].concat(inputs));
      else
        inputsByCategories.set(topCategory, [lookupElement.input.root]);
    }

    // add inputs to the main form (grouped by categories)
    const inputCategoriesByName = this.addFormInputs(form, inputsByCategories, topCategory, (el, isExpanded, catName) => {
      isExpanded ? collapsedInputCategories.delete(catName) : collapsedInputCategories.add(catName);
      el.hidden = !isExpanded || (primaryToggle.value && !primaryInputRoots.has(el));
    });

    // Check if isChanging is enabled
    for (const name of Object.keys(this.store.inputs)) {
      if (this.options.ranges?.[name]?.enabled)
        this.store.inputs[name].isChanging.next(true);
    }

    // 2. Outputs of the function
    const toGetHeader = ui.h1(TITLE.TARGET);
    toGetHeader.classList.add(FORM_TITLE_CLASS);

    ui.tooltip.bind(toGetHeader, 'Select target outputs');
    form.appendChild(toGetHeader);

    // group outputs by categories
    Object.values(this.store.outputs).forEach((outputConfig) => {
      const category = outputConfig.prop.category;
      const roots = [outputConfig.isInterestInput.root, outputConfig.input.root];
      if (outputConfig.prop.type === DG.TYPE.DATA_FRAME) {
        roots.push(
          outputConfig.showInfoWidget,
          outputConfig.argColInput.root,
          outputConfig.funcColsInput.root,
        );
      }
      const propName = outputConfig.prop.name;
      const isPrimary = this.options.targets?.[propName]?.isPrimary;
      if (isPrimary)
        roots.forEach((item) => primaryInputRoots.add(item));

      if (outputsByCategories.has(category))
        outputsByCategories.get(category)!.push(...roots);
      else
        outputsByCategories.set(category, roots);
    });

    // add outputs to the main form (grouped by categories)
    const outputCategoriesByName = this.addFormInputs(form, outputsByCategories, null, (el, isExpanded, catName) => {
      isExpanded ? collapsedOutputCategories.delete(catName) : collapsedOutputCategories.add(catName);
      el.hidden = !isExpanded || (primaryToggle.value && !primaryInputRoots.has(el));
    });

    if (primaryInputRoots.size > 0) {
      primaryToggle.root.hidden = false;
      primaryToggle.value = true;
    }

    // 3. Make one output of interest
    let isAnyOutputSelectedAsOfInterest = false;

    for (const name of Object.keys(this.store.outputs)) {
      if (this.options.targets?.[name]?.enabled)
        this.store.outputs[name].isInterest.next(true);

      if (this.store.outputs[name].isInterest.value === true)
        isAnyOutputSelectedAsOfInterest = true;
    }

    if (!isAnyOutputSelectedAsOfInterest) {
      const firstOutput = this.store.outputs[Object.keys(this.store.outputs)[0]];
      firstOutput.isInterest.next(true);
    }

    // 4. Method & settings block
    const usingHeader = ui.h1('Using');
    usingHeader.classList.add(FORM_TITLE_CLASS);
    ui.tooltip.bind(usingHeader, 'Specify fitting method');
    form.appendChild(usingHeader);
    this.fittingSettingsIcon.style.minWidth = '50px';
    this.fittingSettingsIcon.classList.add(SIDE_ICON_CLASS);

    form.appendChild(this.fittingSettingsIcon);
    form.appendChild(this.methodInput.root);

    // Add general settings
    [this.lossInput, this.samplesCountInput, this.similarityInput].forEach((inp) => {
      inp.nullable = false;
      this.fittingSettingsDiv.appendChild(inp.root);
    });

    // Add random generator settings
    this.fittingSettingsDiv.appendChild(this.randInputs.reproducibility.root);

    this.fittingSettingsDiv.appendChild(this.randInputs.seed.root);

    // Add filtering results inputs
    this.fittingSettingsDiv.appendChild(this.earlyStoppingInputs.earlyStopping.root);

    this.fittingSettingsDiv.appendChild(this.earlyStoppingInputs.threshold.root);

    this.fittingSettingsDiv.appendChild(this.earlyStoppingInputs.stopAtFirst.root);

    // Add input related to the methods
    this.generateSettingInputs();
    this.showHideSettingInputs();
    form.appendChild(this.fittingSettingsDiv);
    this.settingsInputs.forEach((array) => array.forEach((input) => {
      this.fittingSettingsDiv.append(input.root);
    }));

    // deal with side inputs
    for (const item of Array.from(form.children) as HTMLElement[]) {
      if (item.classList.contains(SIDE_INPUT_CLASS)) {
        item.style.height = '0';
        item.style.top = '22px';
        item.style.position = 'relative';
      } else if (item.classList.contains(SIDE_ICON_CLASS)) {
        item.style.height = '0';
        item.style.top = '14px';
        item.style.left = '14px';
        item.style.position = 'relative';
      } else if (item.classList.contains(FORM_TITLE_CLASS) || item.classList.contains(FORM_SECTION_CLASS))
        continue;
      else
        item.style.marginLeft = '50px';
    }

    this.updateApplicabilityState();

    return form;
  } // buildForm

  /** Update run icon tooltip: disabled case */
  private updateRunIconDisabledTooltip(msg: string): void {
    ui.tooltip.bind(this.runIcon, () => {
      const label = ui.label(msg);
      label.style.color = '#FF0000';
      return label;
    });
  } // updateRunIconDisabledTooltip

  /** Check applicability of fitting */
  private updateApplicabilityState(): void {
    this.validationCheckRequests$.next(true);
  } // updateApplicabilityState

  private updateRunIconStyle(): void {
    if (this.readyToRun) {
      this.runIcon.style.color = 'var(--green-2)';
      ui.tooltip.bind(this.runIcon, 'Run fitting');
    } else
      this.runIcon.style.color = 'var(--grey-3)';
  } // updateRunIconStyle

  private getInputValue(input: DG.InputBase) {
    return input.inputType === 'Choice' ? input.stringValue : input.value;
  }

  // make inputs with formulas bounds last
  private splitInputs() {
    return Object.entries(this.store.inputs).reduce((acc, item) => {
      const [, val] = item;
      if (!val.isChanging.value ||
        (!(val as FittingNumericStore).min.useFormula.value && !(val as FittingNumericStore).max.useFormula.value))
        acc.nonFormulaInputs.push(item);
      else
        acc.formulaInputs.push(item);
      return acc;
    }, {nonFormulaInputs: [], formulaInputs: []} as { nonFormulaInputs: [string, FittingInputsStore][], formulaInputs: [string, FittingInputsStore][] });
  } // splitInputs


  private getStoreInputsAndBounds(fittingStore: FittingNumericStore, minContext: Record<string, any>, maxContext: Record<string, any>) {
    const minInp = fittingStore.min.useFormula.value ? fittingStore.min.formula : fittingStore.min.input;
    const maxInp = fittingStore.max.useFormula.value ? fittingStore.max.formula : fittingStore.max.input;
    const minFormula = fittingStore.min.useFormula.value ? fittingStore.min.formula.value : undefined;
    const maxFormula = fittingStore.max.useFormula.value ? fittingStore.max.formula.value : undefined;
    const min = fittingStore.min.useFormula.value ? runFormula(minFormula!, minContext) : fittingStore.min.input.value;
    const max = fittingStore.max.useFormula.value ? runFormula(maxFormula!, maxContext) : fittingStore.max.input.value;
    return {minInp, maxInp, min, max, minFormula, maxFormula};
  } // getStoreInputsAndBounds

  private getAvaliableVars(name: string) {
    const {nonFormulaInputs, formulaInputs} = this.splitInputs();
    const idx = formulaInputs.findIndex(([n]) => n === name);
    if (idx >= 0)
      return [...nonFormulaInputs, ...formulaInputs.slice(0, idx)];
    return nonFormulaInputs;
  }

  /** Check inputs */
  private areInputsReady(): boolean {
    let isAnySelected = false;

    const errors: string[] = [];

    // make inputs with formulas bounds last
    const {nonFormulaInputs, formulaInputs} = this.splitInputs();

    const minContext: Record<string, any> = {};
    const maxContext: Record<string, any> = {};

    // deal with non-formula inputs first, same approach as in sampleParamsWithFormulaBounds,
    // only checking top and bottom points
    for (const [name, input] of [...nonFormulaInputs, ...formulaInputs]) {
      const caption = input.prop.caption ?? input.prop.name;

      if (input.isChanging.value === true) {
        isAnySelected = true;
        const fittingStore = (input as FittingNumericStore);
        const {minInp, maxInp, min, max} = this.getStoreInputsAndBounds(fittingStore, minContext, maxContext);

        minContext[name] = min;
        maxContext[name] = max;

        if (min != null && max != null) {
          if (min > max) {
            errors.push(`Invalid min & max of "${caption}"`);
            minInp.input.classList.add('d4-invalid');
            maxInp.input.classList.add('d4-invalid');
          } else {
            minInp.input.classList.remove('d4-invalid');
            maxInp.input.classList.remove('d4-invalid');
          }
        } else
          errors.push(`Incomplete "${caption}"`);
      } else {
        const val = this.getInputValue(input.const.input);
        minContext[name] = val;
        maxContext[name] = val;

        if (val == null && !JSON.parse(input.prop.options.optional ?? 'false'))
          errors.push(`Incomplete "${caption}"`);
      }
    }

    if (!isAnySelected)
      errors.push(`No parameters for fitting selected`);

    if (errors.length > 0)
      this.updateRunIconDisabledTooltip(errors.join('\n'));

    return errors.length === 0;
  } // areInputsReady

  /** Check outputs */
  private areOutputsReady(): boolean {
    let isAnySelected = false;
    let areSelectedFilled = true;
    let cur: boolean;

    for (const propName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[propName];

      if (output.isInterest.value === true) {
        isAnySelected = true;
        const val = output.input.value;
        cur = (val !== null) && (val !== undefined);

        if (output.prop.propertyType === DG.TYPE.DATA_FRAME) {
          const validation = getValidation(output);

          if (!validation.isValid) {
            this.updateRunIconDisabledTooltip(validation.msg);
            return false;
          }

          cur &&= (output.funcColsInput.value !== null);

          if (cur)
            cur &&= (output.argColInput.value !== null) && (output.argName !== '') && (output.funcColsInput.value.length > 0);
        }

        if (!cur)
          this.updateRunIconDisabledTooltip(`Incomplete the "${output.input.caption}" target`);

        areSelectedFilled = areSelectedFilled && cur;
      }
    }

    if (!isAnySelected)
      this.updateRunIconDisabledTooltip(`No targets selected`);

    return isAnySelected && areSelectedFilled;
  } // areOutputsReady

  /** Check samples count */
  private isSamplesCountValid(): boolean {
    if (this.samplesCount > 0)
      return true;

    this.updateRunIconDisabledTooltip(`Invalid "${TITLE.SAMPLES}"`);
    return false;
  } // isSamplesCountValid

  /** Check similarity */
  private isSimilarityValid(): boolean {
    if ((this.similarity >= FITTING_UI.SIMILARITY_MIN) && (this.similarity <= FITTING_UI.SIMILARITY_MAX))
      return true;

    this.updateRunIconDisabledTooltip(`Invalid "${TITLE.SIMILARITY}"`);
    return false;
  } // isSamplesCountValid

  /** Check inputs/outputs/settings */
  private canFittingBeRun(): boolean {
    return this.areInputsReady() && this.areOutputsReady() &&
      this.isSamplesCountValid() && this.areMethodSettingsCorrect() && this.isSimilarityValid();
  }

  /** Return names of the fixed inputs */
  private getFixedInputs() {
    return Object.keys(this.store.inputs).filter((propName) => !this.store.inputs[propName].isChanging.value);
  }

  /** Return names of the fitted inputs */
  private getFittedInputs() {
    return Object.keys(this.store.inputs)
      .filter((propName) => (this.store.inputs[propName].type === DG.TYPE.FLOAT) && this.store.inputs[propName].isChanging.value);
  }

  private makeOptimizerInputsConfig() {
    // make input bounds payload
    const {nonFormulaInputs, formulaInputs} = this.splitInputs();
    const inputsBounds: Record<string, ValueBoundsData> = {};
    const minContext: Record<string, any> = {};
    const maxContext: Record<string, any> = {};
    for (const [name, input] of [...nonFormulaInputs, ...formulaInputs]) {
      if (input.isChanging.value === true) {
        const fittingStore = (input as FittingNumericStore);
        const {min, max, minFormula, maxFormula} = this.getStoreInputsAndBounds(fittingStore, minContext, maxContext);
        minContext[name] = min;
        maxContext[name] = max;

        const bottom = minFormula ? {
          name,
          type: 'formula' as const,
          formula: minFormula,
        } : {
          name,
          type: 'value' as const,
          value: min!,
        };

        const top = maxFormula ? {
          name,
          type: 'formula' as const,
          formula: maxFormula,
        } : {
          name,
          type: 'value' as const,
          value: max!,
        };

        inputsBounds[name] = {
          type: 'changing',
          top,
          bottom,
        };
      } else {
        const value = this.getInputValue(input.const.input);
        inputsBounds[name] = {
          type: 'const',
          value,
        };
        minContext[name] = value;
        maxContext[name] = value;
      }
    }

    // use original order
    const inputsBoundsOrdered: Record<string, ValueBoundsData> = {};
    const variedInputs = this.getFittedInputs();
    const fixedInputs = this.getFixedInputs();

    for (const name of [...fixedInputs, ...variedInputs])
      inputsBoundsOrdered[name] = inputsBounds[name];


    return inputsBoundsOrdered;
  }

  private makeOptimizerOutputsConfig(): OutputTargetItem[] {
    // make outputs payload
    const outputsOfInterest = this.getOutputsOfInterest();
    const outputTargets: OutputTargetItem[] = outputsOfInterest.map((item) => {
      if (item.prop.propertyType !== DG.TYPE.DATA_FRAME) {
        return {
          type: item.prop.propertyType as DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
          propName: item.prop.name,
          target: item.target as number,
        };
      } else {
        return {
          type: item.prop.propertyType,
          propName: item.prop.name,
          target: item.target as DG.DataFrame,
          argName: item.argName,
          cols: item.funcColsInput.value,
        };
      }
    });
    return outputTargets;
  }

  /** Perform optimization */
  private async runOptimization(): Promise<void> {
    try {
    // check applicability
      if (!this.canFittingBeRun())
        return;

      if (this.samplesCount < 1)
        return;

      if ((this.similarity < FITTING_UI.SIMILARITY_MIN) || (this.similarity > FITTING_UI.SIMILARITY_MAX))
        return;

      this.failsDF = null;

      // inputs of the source function
      const inputs: any = {};

      // add fixed inputs
      const fixedInputs = this.getFixedInputs();
      fixedInputs.forEach((name) => inputs[name] = this.store.inputs[name].const.value);

      // get varied inputs, optimization is performed with respect to them
      const variedInputs = this.getFittedInputs();
      const dim = variedInputs.length;

      // varied inputs specification
      const variedInputNames: string[] = [];
      const variedInputsCaptions = new Array<string>(dim);

      // set varied inputs specification
      variedInputs.forEach((name, idx) => {
        const propConfig = this.store.inputs[name] as FittingNumericStore;
        variedInputNames.push(name);
        variedInputsCaptions[idx] = propConfig.prop.caption ?? propConfig.prop.name;
      });

      // get selected output
      const outputsOfInterest = this.getOutputsOfInterest();
      const outputsCount = outputsOfInterest.length;
      if (outputsCount < 1) {
        grok.shell.error('No output is selected for optimization.');
        return;
      }

      const toShowTableName = outputsOfInterest.map((item) => item.prop.propertyType).filter((type) => type === DG.TYPE.DATA_FRAME).length > 1;

      /** Get call funcCall with the specified inputs */
      const getCalledFuncCall = makeGetCalledFuncCall(this.func, inputs, variedInputNames, true);

      // get optimizer inputs/outputs config
      const inputBounds = this.makeOptimizerInputsConfig();
      const outputTargets = this.makeOptimizerOutputsConfig();

      const costTooltip = this.loss === LOSS.MAD ? 'scaled maximum absolute deviation' : 'scaled root mean square error';

      let optResult: OptimizationResult;

      if (this.method !== METHOD.NELDER_MEAD)
        throw new Error(`Not implemented the '${this.method}' method`);

      // Perform optimization
      if (this.diffGrok !== undefined) {
        try {
          const index = INDICES.DIFF_STUDIO_OUTPUT;
          optResult = await getFittedParams(
            {
              loss: this.loss,
              ivp: this.diffGrok.ivp,
              ivp2ww: this.diffGrok.ivpWW,
              pipelineCreator: this.diffGrok.pipelineCreator,
              settings: this.nelderMeadSettings,
              variedInputNames,
              bounds: inputBounds,
              fixedInputs: inputs,
              argColName: outputsOfInterest[index].argName,
              funcCols: outputsOfInterest[index].funcColsInput.value,
              target: outputsOfInterest[index].target as DG.DataFrame,
              samplesCount: this.samplesCount,
              reproSettings: this.randInputs.settings,
              earlyStoppingSettings: this.earlyStoppingInputs.settings,
            },
          );
        } catch (err) { // run fitting in the main thread if in-webworker run failed
          console.error(err);
          [optResult] = await runOptimizer({
            lossType: this.loss,
            func: this.func,
            inputBounds,
            outputTargets,
            samplesCount: this.samplesCount,
            similarity: this.similarity,
            settings: this.nelderMeadSettings,
            reproSettings: this.randInputs.settings,
            earlyStoppingSettings: this.earlyStoppingInputs.settings,
          });
        }
      } else {
        [optResult] = await runOptimizer({
          lossType: this.loss,
          func: this.func,
          inputBounds,
          outputTargets,
          samplesCount: this.samplesCount,
          similarity: this.similarity,
          settings: this.nelderMeadSettings,
          reproSettings: this.randInputs.settings,
          earlyStoppingSettings: this.earlyStoppingInputs.settings,
        });
      }

      const extrema = optResult.extremums;
      const allExtrCount = extrema.length;

      // Process fails
      if (optResult.fails != null) {
        this.failsDF = optResult.fails;
        const cols = this.failsDF.columns;

        variedInputsCaptions.forEach((cap, idx) => cols.byIndex(idx).name = cols.getUnusedName(cap));

        fixedInputs.forEach((name) => {
          const inp = this.store.inputs[name];

          if (inp.prop.propertyType === DG.TYPE.DATA_FRAME)
            return;

          cols.add(DG.Column.fromList(
            inp.prop.propertyType as unknown as DG.COLUMN_TYPE,
            cols.getUnusedName(inp.prop.caption ?? inp.prop.name),
            new Array(this.failsDF?.rowCount).fill(inp.const.value),
          ));
        });
      }

      // Sort all extrema with respect to the loss function
      extrema.sort((a: Extremum, b: Extremum) => a.cost - b.cost);

      // Extract target dataframes
      const targetDfs: TargetTableOutput[] = outputsOfInterest
        .filter((output) => output.prop.propertyType === DG.TYPE.DATA_FRAME)
        .map((output) => {
          return {name: output.prop.name, target: output.target as DG.DataFrame, argColName: output.argName};
        });

      // Get non-similar points
      const nonSimilarExtrema = await getNonSimilar(extrema, this.similarity, getCalledFuncCall, targetDfs);
      const rowCount = nonSimilarExtrema.length;

      this.clearPrev();

      // Show info/warning reporting results
      if (optResult.fails != null) {
        if (allExtrCount < 1) {
          grok.shell.warning(ui.divV([
            ui.label(`Failed to find ${this.samplesCount} point${this.samplesCount > 1 ? 's' : ''}`),
            ui.div([this.showFailsBtn]),
          ]));

          return;
        } else {
          grok.shell.info(ui.divV([
            ui.label(`Found ${rowCount} point${rowCount > 1 ? 's' : ''}`),
            ui.div([this.showFailsBtn]),
          ]));
        }
      } else
        grok.shell.info(`Found ${rowCount} point${rowCount > 1 ? 's' : ''}`);

      if (allExtrCount < 1) {
        this.currentFuncCalls = [];
        this.comparisonView.dataFrame = DG.DataFrame.fromColumns([DG.Column.fromFloat64Array(TITLE.LOSS, [] as any)]);
        return;
      }

      const lossVals = new Float64Array(rowCount);
      const grid = this.comparisonView.grid;
      const tooltips = new Map([[TITLE.LOSS as string, `The final loss obtained: ${costTooltip}`]]);

      nonSimilarExtrema.forEach((extr, idx) => lossVals[idx] = extr.cost);

      // Add fitting results to the table: iteration & loss
      const reportTable = DG.DataFrame.fromColumns([DG.Column.fromFloat64Array(TITLE.LOSS, lossVals)]);
      this.comparisonView.dataFrame = reportTable;
      const reportColumns = reportTable.columns;

      // Add fitting results to the table: fitted parameters
      variedInputsCaptions.forEach((cap, idx, arr) => {
        cap = reportColumns.getUnusedName(cap);
        arr[idx] = cap;
        const raw = new Float64Array(rowCount);

        for (let j = 0; j < rowCount; ++j)
          raw[j] = nonSimilarExtrema[j].point[idx];

        reportColumns.add(DG.Column.fromFloat64Array(cap, raw));
        tooltips.set(cap, `Obtained values of '${cap}'`);
      });

      // Set properties of the grid
      grid.props.rowHeight = GRID_SIZE.ROW_HEIGHT;
      grid.props.showAddNewRowIcon = false;
      grid.props.allowEdit = false;

      // Add goodness of fit viewers
      const gofTables = new Array<Map<string, GoFtable>>(rowCount);

      this.currentFuncCalls = [];
      for (let idx = 0; idx < rowCount; ++idx) {
        const gofDfs = new Map<string, GoFtable>();
        const scalarGoFs: GoFtable[] = [];

        const calledFuncCall = await getCalledFuncCall(nonSimilarExtrema[idx].point);
        this.currentFuncCalls.push(calledFuncCall);

        outputsOfInterest.forEach((output) => {
          const gofs = this.getOutputGofTable(output.prop, output.target, calledFuncCall, output.argName, toShowTableName);
          gofs.forEach((item) => {
            if (item.chart !== DG.VIEWER.BAR_CHART)
              gofDfs.set(item.caption, item);
            else
              scalarGoFs.push(item);
          });
        });

        if (scalarGoFs.length > 0) {
          const gof = this.getScalarsGofTable(scalarGoFs);
          const name = gof.table.columns.length > 2 ? NAME.SCALARS : gof.table.columns.byIndex(0).name;
          gofDfs.set(name, gof);
        }

        gofTables[idx] = gofDfs;
      }

      // Add goodness of fit columns
      gofTables[0].forEach((_, name) => {
        reportColumns.addNew(name, DG.COLUMN_TYPE.DATA_FRAME).init((row: number) => gofTables[row].get(name)!.table);
        tooltips.set(name, 'Goodness of fit');
        const lossFuncGraphGridCol = grid.columns.byName(name);
        lossFuncGraphGridCol!.cellType = 'html';
        lossFuncGraphGridCol!.width = GRID_SIZE.LOSS_GRAPH_WIDTH;
      });

      // Add dataframe with loss function values
      const lossGraphColName = reportColumns.getUnusedName(`${this.loss} by iterations`);
      tooltips.set(lossGraphColName, `Minimizing ${costTooltip}`);
      reportColumns.addNew(lossGraphColName, DG.COLUMN_TYPE.DATA_FRAME).init((row: number) => getLossFuncDf(nonSimilarExtrema[row]));
      const lossFuncGraphGridCol = grid.columns.byName(lossGraphColName);
      lossFuncGraphGridCol!.cellType = 'html';
      lossFuncGraphGridCol!.width = GRID_SIZE.LOSS_GRAPH_WIDTH;

      // Check whether to use the radar viewer
      let toAddRadars = false;
      gofTables[0].forEach((gof) => toAddRadars ||= toUseRadar(gof.table));

      // Add viewers to the grid
      let toReorderCols = true;

      grid.onCellPrepare(async (gc: DG.GridCell) => {
        if (toReorderCols) {
          grid.columns.setOrder(reportColumns.names().filter((name) => name !== lossGraphColName).concat([lossGraphColName]));
          toReorderCols = false;
        }

        if (gc.isColHeader || gc.isRowHeader) return;

        if (gc.isTableCell && gc.gridColumn.name === lossGraphColName && gc.cell.value !== null)
          gc.style.element = gc.cell.value.plot.line(LOSS_FUNC_CHART_OPTS).root;

        const row = gc.gridRow;

        if (gc.isTableCell && (row < rowCount)) {
          const gof = gofTables[row].get(gc.gridColumn.name);

          if (gof !== undefined) {
            if (gof.chart === DG.VIEWER.SCATTER_PLOT)
              gc.style.element = gc.cell.value.plot.scatter(gof.opts).root;
            else {
              if (toAddRadars) {
                const container = ui.divV([], {style: {height: `${GRID_SIZE.ROW_HEIGHT}px`}});
                ui.tooltip.bind(container, () => getRadarTooltip());
                gc.style.element = container;
                setTimeout(async () => {
                  const rViewer = new ScalarsFitRadar(gc.cell.value as DG.DataFrame);
                  const root = await rViewer.getRoot();
                  root.style.marginTop = '0px';
                  root.style.marginLeft = '0px';
                  container.append(root);
                }, TIMEOUT.RADAR);
              } else
                gc.style.element = getScalarsGoodnessOfFitViewer(gc.cell.value as DG.DataFrame);
            }
          }
        }
      });

      // Set loss-column format & sort by loss-vals
      grid.columns.byName(TITLE.LOSS)!.format = 'scientific';
      grid.columns.byName(TITLE.LOSS)!.width = GRID_SIZE.LOSS_COL_WIDTH;

      // Add tooltips
      grid.onCellTooltip(function(cell, x, y) {
        if (cell.isColHeader) {
          const cellCol = cell.tableColumn;
          if (cellCol) {
            const name = cell.tableColumn.name;

            if ((name === NAME.SCALARS) && toAddRadars)
              ui.tooltip.show(getRadarTooltip(), x, y);
            else {
              const msg = tooltips.get(name);
              ui.tooltip.show(msg ?? '', x, y);
            }

            return true;
          }
        }
      });

      // Set current row with min loss
      reportTable.currentCell = reportTable.cell(0, TITLE.LOSS);

      // Add grid cell effects: show funccall outputs in the context panel
      const cellEffect = async (cell: DG.GridCell) => {
        if ((cell === null) || (cell === undefined))
          return;

        const row = cell.tableRowIndex ?? 0;
        //const selectedRun = calledFuncCalls[row];
        const selectedRun = await getCalledFuncCall(nonSimilarExtrema[row].point);

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

        if (scalarParams.length > 0) {
          overviewPanelConfig = {
            'Output scalars': [scalarTable],
            ...dfPanes,
          };
        } else {
          overviewPanelConfig = {
            ...dfPanes,
          };
        }

        const overviewPanel = ui.accordion();
        $(overviewPanel.root).css({'width': '100%'});
        Object.entries(overviewPanelConfig).map((e) => {
          overviewPanel.addPane(e[0], () => ui.divV(e[1]));
        });

        const paneToExpandIdx = overviewPanel.panes.length - 1;

        if (paneToExpandIdx < 0)
          throw new Error('Function has no outputs');

        overviewPanel.panes[overviewPanel.panes.length - 1].expanded = true;

        grok.shell.o = overviewPanel.root;
      };

      this.gridClickSubscription = grid.onCellClick.subscribe(cellEffect);
      this.gridCellChangeSubscription = grid.onCurrentCellChanged.subscribe(cellEffect);
    } catch (error) {
      grok.shell.error(error instanceof Error ? error.message : 'The platform issue');
    }
  } // runOptimization

  /** Get outputs selected as target */
  private getOutputsOfInterest() {
    const outputsOfInterest = [];

    for (const outputName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[outputName];

      if (output.isInterest.value)
        outputsOfInterest.push(output);
    }

    return outputsOfInterest;
  } // getOutputsOfInterest

  /** Return a table with scalar output results */
  private getScalarsGofTable(gofs: GoFtable[]): GoFtable {
    const cols: DG.Column[] = [];
    gofs.forEach((gof) => cols.push(gof.table.columns.byIndex(INDICES.SCALAR_VALS_COL)));
    cols.push(DG.Column.fromStrings(NAME.CATEGORY, [NAME.SIMULATION, NAME.TARGET]));

    return {
      caption: NAME.SCALARS,
      table: DG.DataFrame.fromColumns(cols),
      chart: cols.length > 2 ? DG.VIEWER.RADAR_VIEWER : DG.VIEWER.BAR_CHART,
      opts: {},
    };
  } // getScalarsGofTable

  /** Return output goodness of fit (GOF) table */
  private getOutputGofTable(prop: DG.Property, target: OutputTarget, call: DG.FuncCall, argColName: string, toShowDfCaption: boolean): GoFtable[] {
    const type = prop.propertyType;
    const name = prop.name;
    const caption = prop.caption ?? name;

    // Provide an appropriate viewer for each output
    switch (type) {
    // bar chart for each scalar
    case DG.TYPE.FLOAT:
    case DG.TYPE.INT:
    case DG.TYPE.BIG_INT:
      return [{
        caption: caption,
        chart: DG.VIEWER.BAR_CHART,
        table: DG.DataFrame.fromColumns([
          DG.Column.fromStrings(NAME.SCALARS, [NAME.SIMULATION, NAME.TARGET]),
          DG.Column.fromList(type as unknown as DG.COLUMN_TYPE, caption, [call.getParamValue(name), target]),
        ]),
        opts: {},
      }];

    // a set of linecharts comparing columns of output dataframe to their targets
    case DG.TYPE.DATA_FRAME:
      const result: GoFtable[] = [];
      const simDf = call.getParamValue(name) as DG.DataFrame;
      const expDf = target as DG.DataFrame;
      const simArgCol = simDf.col(argColName);
      const expArgCol = expDf.col(argColName);

      if ((simArgCol === null) || (expArgCol === null))
        throw new Error('Creating viewer fails: incorrect argument column name');

      // Create a table with results: experiment vs. simulation
      const expCols = expDf.columns;
      const expVsSimDf = simDf.clone(null, expCols.names());

      // Add columns from a table with experiment
      expVsSimDf.append(DG.DataFrame.fromColumns(expCols.toList().map((col) => {
        // To avoid columns of different type
        const simCol = expVsSimDf.col(col.name);

        if (simCol === null)
          throw new Error(`Inconsistent dataframes: no column "${col.name}" in output of of the function`);

        return col.convertTo(simCol.type);
      })), true);

      // Add category column (for splitting data)
      const funcColNames = expVsSimDf.columns.names().filter((name) => name !== argColName);

      // Add style cols
      const categories = new Array<string>(simDf.rowCount).fill(NAME.SIMULATION).concat(new Array<string>(expDf.rowCount).fill(NAME.TARGET));
      const rowsCount = categories.length;

      funcColNames.forEach((name) => {
        expVsSimDf.columns.add(DG.Column.fromStrings(`${NAME.CATEGORY}_${name}`, categories));
      });

      const sizes = new Int32Array(rowsCount).fill(SIZE.SIMULATION, 0, simDf.rowCount).fill(SIZE.TARGET, simDf.rowCount);
      expVsSimDf.columns.add(DG.Column.fromInt32Array(NAME.SIZE, sizes));

      // Coloring scheme
      const simFuncColNames = simDf.columns.names().filter((name) => name !== argColName);

      // Create linecharts
      funcColNames.forEach((name) => {
        const catColName = `${NAME.CATEGORY}_${name}`;
        const catCol = expVsSimDf.col(catColName);
        const colorIdx = simFuncColNames.indexOf(name) % colorsCount;
        const color = colors[colorIdx];
        const rgb = DG.Color.toRgb(color);

        catCol!.colors.setCategorical({
          'Simulation': color,
          'Target': simFuncColNames.length > 1 ? lightenRGB(rgb, SIZE.LIGHTER_PERC) : colors[colorIdx + 1],
        });

        result.push({
          caption: toShowDfCaption ? `${caption}: [${name}]` : `[${name}]`,
          table: expVsSimDf,
          chart: DG.VIEWER.SCATTER_PLOT,
          opts: {
            xColumnName: argColName,
            yColumnName: name,
            showColorSelector: false,
            showSizeSelector: false,
            legendVisibility: 'Always',
            legendPosition: 'Top',
            linesWidth: SIZE.LINE_CHART_LINE_WIDTH,
            linesOrderColumnName: argColName,
            markerType: DG.MARKER_TYPE.CIRCLE,
            markerMinSize: SIZE.MIN_MARKER,
            markerMaxSize: SIZE.MAX_MARKER,
            colorColumnName: catColName,
            sizeColumnName: NAME.SIZE,
          } as DG.IScatterPlotSettings,
        });
      });

      return result;

    default:
      throw new Error('Unsupported output type');
    }
  } // getOutputGofTable

  /** Clear previous results: close dock nodes and unsubscribe from events */
  private clearPrev(): void {
    if (this.helpDN !== undefined) {
      this.comparisonView.dockManager.close(this.helpDN);
      this.helpDN = undefined;
    }

    if (this.gridClickSubscription) {
      this.gridClickSubscription.unsubscribe();
      this.gridClickSubscription = null;
    }

    if (this.gridCellChangeSubscription) {
      this.gridCellChangeSubscription.unsubscribe();
      this.gridCellChangeSubscription = null;
    }
  } // clearPrev

  /** Show fails of the Nelder-Mead fails */
  private showNelderMeadFails(): void {
    if (this.failsDF === null)
      return;

    // show fails info: inputs, issue message
    const view = grok.shell.addTableView(this.failsDF);

    // add method's settings
    view.addViewer(DG.Viewer.form(DG.DataFrame.fromColumns(
      Object.entries(this.nelderMeadSettings).map((e) => DG.Column.fromFloat64Array(nelderMeadSettingsOpts.get(e[0])?.caption ?? e[0], new Float64Array([e[1]]))),
    )), {description: 'The Nelder-Mead method settings', showNavigation: false});

    // create tooltips
    const count = this.getFittedInputs().length;
    const tooltips = new Map<string, string>();
    this.failsDF.columns.names().forEach((name, idx) => {
      if (idx < count)
        tooltips.set(name, `Initial value of "${name}" used in the Nelder-Mead method`);
      else if (idx > count)
        tooltips.set(name, `Value of "${name}" (fixed input)`);
      else
        tooltips.set(name, 'Error message raised');
    });

    // add tooltips
    view.grid.onCellTooltip(function(cell, x, y) {
      if (cell.isColHeader) {
        const cellCol = cell.tableColumn;
        if (cellCol) {
          const name = cell.tableColumn.name;
          const msg = tooltips.get(name);
          ui.tooltip.show(msg ?? '', x, y);
          return true;
        }
      }
    });
  } // showNelderMeadFails
}
