/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {BehaviorSubject} from 'rxjs';
import {RunComparisonView} from './run-comparison-view';
import {combineLatest} from 'rxjs';
import '../css/sens-analysis.css';
import {CARD_VIEW_TYPE} from '../../shared-utils/consts';
import {STARTING_HELP} from './fitting/constants';
import {optimize} from './fitting/optimizer';

import {NELDER_MEAD_DEFAULTS} from './fitting/optimizer-nelder-mead';

const RUN_NAME_COL_LABEL = 'Run name' as const;
const supportedOutputTypes = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.DATA_FRAME];

enum METHOD {
  NELDER_MEAD = 'Nelder-Mead',
  GRAD_DESC = 'Gradient descent',
};

type InputWithValue<T = number> = {input: DG.InputBase, value: T};

type InputValues = {
  isChanging: BehaviorSubject<boolean>,
  const: InputWithValue<boolean | number | string | DG.DataFrame>,
  constForm: DG.InputBase[],
  saForm: DG.InputBase[],
}

type SensitivityNumericStore = {
  prop: DG.Property,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
  min: InputWithValue,
  max: InputWithValue,
} & InputValues;

type SensitivityBoolStore = {
  prop: DG.Property,
  type: DG.TYPE.BOOL,
} & InputValues;

type SensitivityConstStore = {
  prop: DG.Property,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL | DG.TYPE.STRING>,
} & InputValues;

type FittingInputsStore = SensitivityNumericStore | SensitivityBoolStore | SensitivityConstStore;

const getSwitchMock = () => ui.div([], 'sa-switch-input');

const isNumericProp = (prop: DG.Property) => ((prop.propertyType === DG.TYPE.INT) || (prop.propertyType === DG.TYPE.FLOAT));

export class FittingView {
  generateInputFields = (func: DG.Func) => {
    const getInputValue = (input: DG.Property, key: string) => (
      input.options[key] === undefined ? input.defaultValue : Number(input.options[key])
    );

    const getSwitchElement = (defaultValue: boolean, f: (v: boolean) => any, isInput: boolean = true) => {
      const input = ui.switchInput(' ', defaultValue, f);
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
      if (inputProp.propertyType === DG.TYPE.FLOAT) {
        const isChangingInputMin = getSwitchElement(false, (v: boolean) => {
          ref.isChanging.next(v);
          this.updateRunWidgetsState();
        });

        const isChangingInputConst = getSwitchElement(false, (v: boolean) => {
          ref.isChanging.next(v);
          this.updateRunWidgetsState();
        });

        const caption = inputProp.caption ?? inputProp.name;

        const temp = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input:
            (() => {
              const inp = ui.intInput(caption, inputProp.defaultValue, (v: number) => ref.const.value = v);
              inp.root.insertBefore(isChangingInputConst.root, inp.captionLabel);
              inp.addPostfix(inputProp.options['units']);
              inp.setTooltip(`Value of the '${caption}' input`);
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
                inp.setTooltip(`Min value of the '${caption}' input`);
                return inp;
              })(),
            value: getInputValue(inputProp, 'min'),
          },
          max: {
            input: (() => {
              const inp = ui.floatInput(`${inputProp.caption ?? inputProp.name} (max)`, getInputValue(inputProp, 'max'), (v: number) => (ref as SensitivityNumericStore).max.value = v);
              inp.addPostfix(inputProp.options['units']);
              inp.setTooltip(`Max value of the '${caption}' input`);
              return inp;
            })(),
            value: getInputValue(inputProp, 'max'),
          },
          isChanging: new BehaviorSubject<boolean>(false),
        };

        [temp.max.input].forEach((input) => {
          input.root.insertBefore(getSwitchMock(), input.captionLabel);
          $(input.root).removeProp('display');
        });

        acc[inputProp.name] = {
          ...temp,
          constForm: [temp.const.input],
          saForm: [
            temp.min.input,
            temp.max.input,
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
          temp.isChanging,
        ]).subscribe(([isChanging]) => {
          if (isChanging) {
            ref.constForm.forEach((input) => $(input.root).hide());
            ref.saForm.forEach((input) => $(input.root).css('display', 'flex'));
          } else {
            ref.constForm.forEach((input) => $(input.root).css('display', 'flex'));
            ref.saForm.forEach((input) => $(input.root).hide());
          }
        });
      } else {
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
    }, {} as Record<string, FittingInputsStore>);

    const outputs = func.outputs.filter((prop) => supportedOutputTypes.includes(prop.propertyType))
      .reduce((acc, outputProp) => {
        const temp = {
          prop: outputProp,
          input:
          (() => {
            const caption = outputProp.caption ?? outputProp.name;
            const input = ui.input.forProperty(outputProp);
            input.addCaption(caption);
            input.setTooltip(this.toSetSwitched ?
              'Target value' :
              (outputProp.propertyType === DG.TYPE.DATA_FRAME) ? 'Output dataframe' : 'Output scalar');
            input.input.hidden = !this.toSetSwitched;
            input.nullable = true;

            input.onChanged(() => {
              temp.target = input.value;

              if (outputProp.propertyType === DG.TYPE.DATA_FRAME) {
                const colNames: string[] = [];

                for (const col of (input.value as DG.DataFrame).columns) {
                  if (col.type === DG.COLUMN_TYPE.FLOAT || col.type === DG.COLUMN_TYPE.INT)
                    colNames.push(col.name);
                }

                temp.colNameInput.items = colNames;
                temp.colNameInput.value = colNames[0];
              }
            });

            if (outputProp.propertyType === DG.TYPE.DATA_FRAME)
              (input.root.lastElementChild as HTMLDivElement).hidden = !this.toSetSwitched;

            const isInterestInput = getSwitchElement(
              this.toSetSwitched,
              (v: boolean) => {
                temp.isInterest.next(v);
                this.updateRunWidgetsState();
                input.input.hidden = !v;
                input.setTooltip(v ? 'Target value' :
                  (outputProp.propertyType === DG.TYPE.DATA_FRAME) ? 'Output dataframe' : 'Output scalar');
                if (outputProp.propertyType === DG.TYPE.DATA_FRAME) {
                  (input.root.lastElementChild as HTMLDivElement).hidden = !v;
                  temp.colNameInput.root.hidden = !v;
                }
              },
              false,
            ).root;
            input.root.insertBefore(isInterestInput, input.captionLabel);

            return input;
          })(),
          colNameInput: (() => {
            const input = ui.choiceInput('argument', '', [''], (v: string) => temp.colName = v);
            input.setTooltip('Column with argument values');
            input.root.insertBefore(getSwitchMock(), input.captionLabel);
            input.root.hidden = outputProp.propertyType !== DG.TYPE.DATA_FRAME || !this.toSetSwitched;
            this.toSetSwitched = false;

            return input;
          })(),
          isInterest: new BehaviorSubject<boolean>(this.toSetSwitched),
          target: null,
          colName: '',
        };

        acc[outputProp.name] = temp;

        return acc;
      }, {} as Record<string, {
      prop: DG.Property,
      input: DG.InputBase,
      isInterest: BehaviorSubject<boolean>,
      target: number | DG.DataFrame | null,
      colName: string,
      colNameInput: DG.InputBase,
    }>);

    return {inputs, outputs};
  };

  private runButton = ui.bigButton('Run', async () => await this.runOptimization(), 'Run fitting');
  private runIcon = ui.iconFA('play', async () => await this.runOptimization(), 'Run fitting');
  private helpIcon = ui.iconFA('question', () => {
    window.open('https://datagrok.ai/help/compute/#input-parameter-optimization', '_blank');
  }, 'Open help in a new tab');
  private tableDockNode: DG.DockNode | undefined;
  private helpMdNode: DG.DockNode | undefined;
  private gridSubscription: any = null;
  private toSetSwitched = true;

  private method = METHOD.NELDER_MEAD;
  private methodInput = ui.choiceInput('method', this.method, [METHOD.NELDER_MEAD, METHOD.GRAD_DESC], () => {
    this.method = this.methodInput.value!;
    this.showHideSettingInputs();
  });

  // The Nelder-Mead method settings
  private nelderMeadSettings = {
    tolerance: NELDER_MEAD_DEFAULTS.TOLERANCE,
    maxIter: NELDER_MEAD_DEFAULTS.MAX_ITER,
    nonZeroParam: NELDER_MEAD_DEFAULTS.NON_ZERO_PARAM,
    initialScale: NELDER_MEAD_DEFAULTS.INITIAL_SCALE,
    scaleReflaction: NELDER_MEAD_DEFAULTS.SCALE_REFLECTION,
    scaleExpansion: NELDER_MEAD_DEFAULTS.SCALE_EXPANSION,
    scaleContraction: NELDER_MEAD_DEFAULTS.SCALE_CONTRACTION,
  };

  // Foo settings
  private gradDescentSettings = {
    iterCount: 10,
    learningRate: 0.0001,
  };

  private settingsInputs = new Map<METHOD, DG.InputBase[]>();

  store = this.generateInputFields(this.func);
  comparisonView!: DG.TableView;

  private fittingSettingsIcon = ui.iconFA('cog', () => {
    const prevState = this.fittingSettingsDiv.hidden;
    this.fittingSettingsDiv.hidden = !prevState;
  }, 'Modify fitting settings');

  private fittingSettingsDiv = ui.divV([]);

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
      grok.shell.warning('Fitting is not applicable: the function has no scalar outputs.');
      baseView.close();
      return;
    }

    this.runIcon = ui.iconFA('play', async () => await this.runOptimization(), 'Run fitting');
    this.runIcon.style.color = 'var(--green-2)';
    this.runIcon.classList.add('fas');

    this.fittingSettingsDiv.hidden = true;

    const form = this.buildFormWithBtn();
    this.runButton.disabled = !this.canEvaluationBeRun();
    this.runIcon.hidden = this.runButton.disabled;
    this.comparisonView = baseView;

    this.comparisonView.dockManager.dock(
      form,
      DG.DOCK_TYPE.LEFT,
      null,
      `${this.func.name} - Fitting`,
      0.25,
    );

    this.comparisonView.grid.columns.byName(RUN_NAME_COL_LABEL)!.visible = false;

    const rbnPanels = this.comparisonView.getRibbonPanels();
    rbnPanels.push([this.helpIcon, this.runIcon]);
    this.comparisonView.setRibbonPanels(rbnPanels);

    this.comparisonView.name = this.comparisonView.name.replace('comparison', 'fitting');
    this.comparisonView.helpUrl = 'https://datagrok.ai/help/compute/#input-parameter-optimization';
    this.tableDockNode = this.comparisonView.dockManager.findNode(this.comparisonView.grid.root);
    const helpMD = ui.markdown(STARTING_HELP);
    helpMD.style.padding = '10px';
    helpMD.style.overflow = 'auto';
    this.helpMdNode = this.comparisonView.dockManager.dock(helpMD, DG.DOCK_TYPE.FILL, this.tableDockNode, 'About');
    this.methodInput.setTooltip('Numerical method for minimizing objective function');
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

  private generateNelderMeadSettingsInputs(): void {
    const tolInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'tolerance',
      inputType: 'Float',
      defaultValue: NELDER_MEAD_DEFAULTS.TOLERANCE,
      min: 1e-20,
      max: 1e-1,
    }));
    tolInp.onChanged(() => this.nelderMeadSettings.tolerance = tolInp.value);

    const maxIterInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'Max iterations',
      inputType: 'Int',
      defaultValue: NELDER_MEAD_DEFAULTS.MAX_ITER,
      min: 1,
      max: 10000,
    }));
    maxIterInp.onChanged(() => this.nelderMeadSettings.maxIter = maxIterInp.value);

    const nonZeroParamInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'Non-zero param',
      inputType: 'Float',
      defaultValue: NELDER_MEAD_DEFAULTS.NON_ZERO_PARAM,
      min: 1e-20,
      max: 1e-1,
    }));
    nonZeroParamInp.onChanged(() => this.nelderMeadSettings.nonZeroParam = nonZeroParamInp.value);

    const initialScaleInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'Initial scale',
      inputType: 'Float',
      defaultValue: NELDER_MEAD_DEFAULTS.INITIAL_SCALE,
      min: 1e-20,
      max: 1e-1,
    }));
    initialScaleInp.onChanged(() => this.nelderMeadSettings.initialScale = initialScaleInp.value);

    const scaleReflactionInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'Scale reflection',
      inputType: 'Float',
      defaultValue: NELDER_MEAD_DEFAULTS.SCALE_REFLECTION,
    }));
    scaleReflactionInp.onChanged(() => this.nelderMeadSettings.scaleReflaction = scaleReflactionInp.value);

    const scaleExpansionInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'Scale expansion',
      inputType: 'Float',
      defaultValue: NELDER_MEAD_DEFAULTS.SCALE_EXPANSION,
    }));
    scaleExpansionInp.onChanged(() => this.nelderMeadSettings.scaleExpansion = scaleExpansionInp.value);

    const scaleContractionInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'Scale contraction',
      inputType: 'Float',
      defaultValue: NELDER_MEAD_DEFAULTS.SCALE_CONTRACTION,
    }));
    scaleContractionInp.onChanged(() => this.nelderMeadSettings.scaleContraction = scaleContractionInp.value);

    this.settingsInputs.set(METHOD.NELDER_MEAD, [
      tolInp,
      maxIterInp,
      nonZeroParamInp,
      initialScaleInp,
      scaleReflactionInp,
      scaleExpansionInp,
      scaleContractionInp,
    ]);
  }

  private generateGradDescentSettingsInputs(): void {
    const iterInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'iterations',
      inputType: 'Int',
      defaultValue: this.gradDescentSettings.iterCount,
      min: 1,
      max: 10000,
    }));
    iterInp.onChanged(() => this.gradDescentSettings.iterCount = iterInp.value);

    const learningRateInp = ui.input.forProperty(DG.Property.fromOptions({
      name: 'Learning rate',
      inputType: 'Float',
      defaultValue: this.gradDescentSettings.learningRate,
      min: 1e-6,
      max: 1000,
    }));
    learningRateInp.onChanged(() => this.gradDescentSettings.learningRate = learningRateInp.value);

    this.settingsInputs.set(METHOD.GRAD_DESC, [iterInp, learningRateInp]);
  }

  private generateSettingInputs(): void {
    this.generateNelderMeadSettingsInputs();
    this.generateGradDescentSettingsInputs();
  }

  private showHideSettingInputs(): void {
    this.settingsInputs.forEach((inputsArray, method) => inputsArray.forEach((input) => input.root.hidden = method !== this.method));
  }

  private buildFormWithBtn() {
    let prevCategory = 'Misc';
    const fitHeader = ui.h2('Fit');
    ui.tooltip.bind(fitHeader, 'Select inputs to be fitted');

    const form = Object.values(this.store.inputs)
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
      }, ui.div([
        fitHeader,
      ], {style: {'overflow-y': 'scroll', 'width': '100%'}}));

    const toGetHeader = ui.h2('Target');
    ui.tooltip.bind(toGetHeader, 'Select target outputs');
    form.appendChild(toGetHeader);
    prevCategory = 'Misc';

    Object.values(this.store.outputs)
      .reduce((container, outputConfig) => {
        const prop = outputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.p(prop.category));
          prevCategory = prop.category;
        }

        container.append(outputConfig.input.root);
        container.append(outputConfig.colNameInput.root);

        return container;
      }, form);

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

    const usingHeader = ui.h2('Using');
    ui.tooltip.bind(usingHeader, 'Specify fitting method');
    form.appendChild(usingHeader);
    this.methodInput.root.insertBefore(this.fittingSettingsIcon, this.methodInput.captionLabel);
    this.fittingSettingsIcon.style.minWidth = '50px';
    form.appendChild(this.methodInput.root);

    const settingsHeader = ui.h3('with settings');
    ui.tooltip.bind(settingsHeader, () => `Settings of the ${this.method} method`);
    this.fittingSettingsDiv.appendChild(settingsHeader);
    this.generateSettingInputs();
    this.showHideSettingInputs();
    form.appendChild(this.fittingSettingsDiv);
    this.settingsInputs.forEach((array) => array.forEach((input) => {
      input.root.insertBefore(getSwitchMock(), input.captionLabel);
      this.fittingSettingsDiv.append(input.root);
    }));

    $(form).addClass('ui-form');

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
    console.log(this.store);

    if (!this.canEvaluationBeRun())
      return;
    //await this.runLocalMinimumOptimization();
  }

  private getFixedInputs() {
    return Object.keys(this.store.inputs)
      .filter((propName) => (this.store.inputs[propName].type === DG.TYPE.FLOAT) && !this.store.inputs[propName].isChanging.value);
  }

  private getVariedInputs() {
    return Object.keys(this.store.inputs)
      .filter((propName) => (this.store.inputs[propName].type === DG.TYPE.FLOAT) && this.store.inputs[propName].isChanging.value);
  }

  /** Perform Nelder-Mead method */
  private async runLocalMinimumOptimization() {
    // inputs of the source function
    const inputs: any = {};

    // add fixed inputs
    this.getFixedInputs().forEach((name) => inputs[name] = this.store.inputs[name].const.value);

    // get varied inputs, optimization is performed with respect to them
    const variedInputs = this.getVariedInputs();
    const dim = variedInputs.length;

    // varied inputs specification
    const variedInputNames = [] as string[];
    const minVals = new Float32Array(dim);
    const maxVals = new Float32Array(dim);

    // set varied inputs specification
    variedInputs.forEach((name, idx) => {
      const propConfig = this.store.inputs[name] as SensitivityNumericStore;
      minVals[idx] = propConfig.min.value;
      maxVals[idx] = propConfig.max.value;
      variedInputNames.push(name);
    });

    // get selected output
    const outputsOfInterest = this.getOutputsOfInterest();
    if (outputsOfInterest.length !== 1) {
      grok.shell.error('No output is selected for optimization.');
      return;
    }
    const outputName = outputsOfInterest[0].name;

    /** Cost function to be optimized */
    const costFunc = async (x: Float32Array): Promise<number> => {
      x.forEach((val, idx) => inputs[variedInputNames[idx]] = val);
      const funcCall = this.func.prepare(inputs);
      const calledFuncCall = await funcCall.call();

      return calledFuncCall.getParamValue(outputName);
    };

    const extr = await optimize(costFunc,
      minVals,
      maxVals,
      100,
    );

    console.log(extr);
  } // runNelderMeadMethod

  private getOutputsOfInterest() {
    const outputsOfInterest = [];

    for (const outputName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[outputName];

      if (output.isInterest.value)
        outputsOfInterest.push(output.prop);
    }

    return outputsOfInterest;
  }
}
