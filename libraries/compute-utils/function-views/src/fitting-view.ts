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
import {STARTING_HELP, TITLE, DOCK_RATIO, REPORT_DF_TOOLTIP, COL_WIDTH} from './fitting/constants';
import {performNelderMeadOptimization} from './fitting/optimizer';

import {NELDER_MEAD_DEFAULTS, NelderMeadSettings} from './fitting/optimizer-nelder-mead';
import {getErrors} from './fitting/fitting-utils';
import {Extremum} from './fitting/optimizer-misc';

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
            const temp = ui.input.forProperty(inputProp);
            temp.onChanged(() => tempDefault.value = temp.value);
            temp.root.insertBefore(switchMock, temp.captionLabel);
            temp.addPostfix(inputProp.options['units']);

            return temp;
          })(),
          value: inputProp.defaultValue,
        };
        tempDefault.input.value = null;
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
            input.value = 0;
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
          target: (outputProp.propertyType !== DG.TYPE.DATA_FRAME) ? 0 : null,
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

  private openedViewers: DG.Viewer[] = [];

  private runButton = ui.bigButton('Run', async () => await this.runOptimization(), 'Run fitting');
  private runIcon = ui.iconFA('play', async () => await this.runOptimization(), 'Run fitting');
  private helpIcon = ui.iconFA('question', () => {
    window.open('https://datagrok.ai/help/compute/#input-parameter-optimization', '_blank');
  }, 'Open help in a new tab');
  private gridClickSubscription: any = null;
  private gridCellChangeSubscription: any = null;
  private toSetSwitched = true;
  private samplesCount = 10;
  private samplesCountInput = ui.input.forProperty(DG.Property.fromOptions({
    name: TITLE.SAMPLES,
    inputType: 'Int',
    defaultValue: this.samplesCount,
    min: 1,
    max: 1000,
  }));

  // Dock Nodes with results
  private tableDockNode: DG.DockNode | undefined = undefined;
  private helpMdNode: DG.DockNode | undefined = undefined;
  private goodnessOfFitNode: DG.DockNode | undefined = undefined;
  private pcPlotNode: DG.DockNode | undefined = undefined;
  private lossFuncLineChartsNode: DG.DockNode | undefined = undefined;
  private tempDockNodes: DG.DockNode[] = [];

  private method = METHOD.NELDER_MEAD;
  private methodInput = ui.choiceInput(TITLE.METHOD, this.method, [METHOD.NELDER_MEAD], () => {
    this.method = this.methodInput.value!;
    this.showHideSettingInputs();
  });

  // The Nelder-Mead method settings
  private nelderMeadSettings: NelderMeadSettings = {
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
    this.tempDockNodes.push(this.comparisonView.dockManager.dock(helpMD, DG.DOCK_TYPE.FILL, this.tableDockNode, 'About'));
    this.methodInput.setTooltip('Numerical method for minimizing loss function');

    this.samplesCountInput.addCaption(TITLE.SAMPLES);
    this.samplesCountInput.onChanged(() => {
      this.samplesCount = this.samplesCountInput.value;
      this.updateRunWidgetsState();
    });
    this.samplesCountInput.setTooltip('Number fitted inputs sets');
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
    const fitHeader = ui.h2(TITLE.FIT);
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

    const toGetHeader = ui.h2(TITLE.TARGET);
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

    this.samplesCountInput.root.insertBefore(getSwitchMock(), this.samplesCountInput.captionLabel);
    form.appendChild(this.samplesCountInput.root);

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
    return this.isAnyInputSelected() && this.isAnyOutputSelected() && (this.samplesCount > 0);
  }

  private getFixedInputs() {
    return Object.keys(this.store.inputs).filter((propName) => !this.store.inputs[propName].isChanging.value);
  }

  private getVariedInputs() {
    return Object.keys(this.store.inputs)
      .filter((propName) => (this.store.inputs[propName].type === DG.TYPE.FLOAT) && this.store.inputs[propName].isChanging.value);
  }

  /** Perform the Nelder-Mead method */
  private async runOptimization(): Promise<void> {
    try {
    // check applicability
      if (!this.canEvaluationBeRun())
        return;

      if (this.samplesCount < 1)
        return;

      // inputs of the source function
      const inputs: any = {};

      // add fixed inputs
      this.getFixedInputs().forEach((name) => inputs[name] = this.store.inputs[name].const.value);
      //console.log(inputs);

      // get varied inputs, optimization is performed with respect to them
      const variedInputs = this.getVariedInputs();
      const dim = variedInputs.length;

      // varied inputs specification
      const variedInputNames: string[] = [];
      const minVals = new Float32Array(dim);
      const maxVals = new Float32Array(dim);

      const variedInputsCaptions = new Array<string>(dim);

      // set varied inputs specification
      variedInputs.forEach((name, idx) => {
        const propConfig = this.store.inputs[name] as SensitivityNumericStore;
        minVals[idx] = propConfig.min.value;
        maxVals[idx] = propConfig.max.value;
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

      /** Get call funcCall with the specified inputs */
      const getCalledFuncCall = async (x: Float32Array): Promise<DG.FuncCall> => {
        x.forEach((val, idx) => inputs[variedInputNames[idx]] = val);
        const funcCall = this.func.prepare(inputs);
        return await funcCall.call();
      };

      /** Cost function to be optimized */
      const costFunc = async (x: Float32Array): Promise<number> => {
        x.forEach((val, idx) => inputs[variedInputNames[idx]] = val);
        const funcCall = this.func.prepare(inputs);
        const calledFuncCall = await funcCall.call();

        //console.log(x);

        let cost = 0;

        outputsOfInterest.forEach((output) => {
          if (output.prop.propertyType !== DG.TYPE.DATA_FRAME)
            cost += (output.target as number - calledFuncCall.getParamValue(output.prop.name)) ** 2;
          //console.log(output.target as number - calledFuncCall.getParamValue(output.prop.name));
          else {
          //const errs = getErrors(output.colName, output.target as DG.DataFrame, calledFuncCall.getParamValue(output.prop.name));
          //console.log(errs);
            getErrors(output.colName, output.target as DG.DataFrame, calledFuncCall.getParamValue(output.prop.name)).forEach((err) => cost += err ** 2);
          }
        //console.log(output.target);
        //console.log('---------------------------------------------------------------');
        //console.log(calledFuncCall.getParamValue(output.prop.name));
        });

        //console.log(cost);

        return cost;
      };

      //console.log(await costFunc(new Float32Array([1, 5])));

      let extremums: Extremum[];

      // Perform optimization
      if (this.method === METHOD.NELDER_MEAD)
        extremums = await performNelderMeadOptimization(costFunc, minVals, maxVals, this.nelderMeadSettings, this.samplesCount);
      else
        throw new Error(`The '${this.method}' method has not been implemented.`);

      this.closeDockNodes();

      const rowCount = this.samplesCount;
      const iterCounts = new Int32Array(rowCount);
      const lossVals = new Float32Array(rowCount);
      const grid = this.comparisonView.grid;

      const lossFuncLineChartsRoots = new Map<number, HTMLElement>();
      const goodnessOfFitDivs = new Map<number, HTMLDivElement>();

      extremums.forEach((extr, idx) => {
        iterCounts[idx] = extr.iterCount;
        lossVals[idx] = extr.cost;
        lossFuncLineChartsRoots.set(idx, DG.Viewer.lineChart(
          DG.DataFrame.fromColumns([
            DG.Column.fromList(DG.COLUMN_TYPE.INT, TITLE.ITER, [...Array(extr.iterCount).keys()].map((i) => i + 1)),
            DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, TITLE.LOSS, extr.iterCosts.slice(0, extr.iterCount))]),
          {
            description: 'Loss: sum of squared errors',
            showYAxis: true,
            showXAxis: true,
          }).root);

        goodnessOfFitDivs.set(idx, ui.divV([ui.label(`Goodness of fit: ${idx}`)]));
      });

      // Add fitting results to the table: iteration & loss
      const reportTable = DG.DataFrame.fromColumns([
        DG.Column.fromInt32Array(TITLE.ITERATIONS, iterCounts),
        DG.Column.fromFloat32Array(TITLE.LOSS, lossVals),
      ]);
      this.comparisonView.dataFrame = reportTable;
      const reportColumns = reportTable.columns;

      // Add fitting results to the table: fitted parameters
      variedInputsCaptions.forEach((cap, idx, arr) => {
        arr[idx] = reportColumns.getUnusedName(cap);
        const raw = new Float32Array(rowCount);

        for (let j = 0; j < rowCount; ++j)
          raw[j] = extremums[j].point[idx];

        reportColumns.add(DG.Column.fromFloat32Array(variedInputsCaptions[idx], raw));
      });

      // Add id's of obtained points
      const idName = reportColumns.getUnusedName(TITLE.ID);
      const idVals = new Int32Array(rowCount);
      idVals.forEach((_, idx, arr) => arr[idx] = idx);
      reportColumns.add(DG.Column.fromInt32Array(idName, idVals));

      // Hide id-column
      grid.columns.setVisible([idName]); // DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13450
      grid.columns.setVisible(reportColumns.names().filter((name) => name !== idName));

      // Set loss-column format & sort by loss-vals
      grid.columns.byName(TITLE.LOSS)!.format = 'scientific';
      grid.sort([TITLE.LOSS]);
      grid.columns.byName(TITLE.LOSS)!.width = COL_WIDTH;

      // Add tooltips
      grid.onCellTooltip(function(cell, x, y) {
        if (cell.isColHeader) {
          const cellCol = cell.tableColumn;
          if (cellCol) {
            const name = cell.tableColumn.name;
            const msg = REPORT_DF_TOOLTIP.get(name);
            if (msg)
              ui.tooltip.show(msg, x, y);
            else
              ui.tooltip.show(`Obtained values of '${name}'`, x, y);
            return true;
          }
        }
      });

      // Add PC plot
      const pcPlot = DG.Viewer.pcPlot(reportTable, {
        columnNames: variedInputsCaptions.concat([TITLE.LOSS]),
        description: 'Variation',
      });
      const pcPlotNode = this.comparisonView.dockManager.dock(pcPlot, DG.DOCK_TYPE.LEFT, this.tableDockNode, '', DOCK_RATIO.PC_PLOT);

      // Add loss linecharts
      const lossFuncLineChartsDiv = ui.divV([]);
      lossFuncLineChartsRoots.forEach((r) => lossFuncLineChartsDiv.appendChild(r));
      const lossFuncLineChartsNode = this.comparisonView.dockManager.dock(lossFuncLineChartsDiv, DG.DOCK_TYPE.DOWN, pcPlotNode, '', DOCK_RATIO.LOSS_PLOT);

      // Add goodness of fit divs
      const gofDiv = ui.divV([]);
      goodnessOfFitDivs.forEach((d) => gofDiv.appendChild(d));
      const goodnessOfFitNode = this.comparisonView.dockManager.dock(gofDiv, DG.DOCK_TYPE.DOWN, this.tableDockNode, '', DOCK_RATIO.FIT_DIV);

      // Update nodes store & hide titles
      this.tempDockNodes.push(pcPlotNode, lossFuncLineChartsNode, goodnessOfFitNode);
      this.tempDockNodes.forEach((node) => {
        const title = node.container.dart.elementTitle;
        if (title)
          title.hidden = true;
      });

      // Find item with min loss
      let minLossExtrIdx = 0;
      let minLoss = extremums[0].cost;
      for (let i = 1; i < rowCount; ++i) {
        if (extremums[i].cost < minLoss) {
          minLoss = extremums[i].cost;
          minLossExtrIdx = i;
        }
      }
      // Set current row with min loss
      reportTable.currentCell = reportTable.cell(minLossExtrIdx, TITLE.LOSS);

      // Set visibilities
      lossFuncLineChartsRoots.forEach((r, i) => {
        r.hidden = i !== minLossExtrIdx;
        r.style.width = '100%';
        r.style.height = '100%';
        goodnessOfFitDivs.get(i)!.hidden = i !== minLossExtrIdx;
      });

      // Add grid cell effects
      const cellEffect = async (cell: DG.GridCell) => {
        const row = cell.tableRowIndex ?? 0;
        lossFuncLineChartsRoots.forEach((r, i) => {
          r.hidden = i !== row;
          goodnessOfFitDivs.get(i)!.hidden = i !== row;
        });
      };

      this.gridClickSubscription = grid.onCellClick.subscribe(cellEffect);
      this.gridCellChangeSubscription = grid.onCurrentCellChanged.subscribe(cellEffect);
    } catch (error) {
      grok.shell.error(error instanceof Error ? error.message : 'The platform issue');
    }
  } // runOptimization

  /** Perform the Nelder-Mead method */
  private async runOptimizationOld(): Promise<void> {
    try {
    // check applicability
      if (!this.canEvaluationBeRun())
        return;

      // inputs of the source function
      const inputs: any = {};

      // add fixed inputs
      this.getFixedInputs().forEach((name) => inputs[name] = this.store.inputs[name].const.value);
      //console.log(inputs);

      // get varied inputs, optimization is performed with respect to them
      const variedInputs = this.getVariedInputs();
      const dim = variedInputs.length;

      // varied inputs specification
      const variedInputNames: string[] = [];
      const minVals = new Float32Array(dim);
      const maxVals = new Float32Array(dim);

      const variedInputsCaptions = new Array<string>(dim);

      // set varied inputs specification
      variedInputs.forEach((name, idx) => {
        const propConfig = this.store.inputs[name] as SensitivityNumericStore;
        minVals[idx] = propConfig.min.value;
        maxVals[idx] = propConfig.max.value;
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

      /** Get call funcCall with the specified inputs */
      const getCalledFuncCall = async (x: Float32Array): Promise<DG.FuncCall> => {
        x.forEach((val, idx) => inputs[variedInputNames[idx]] = val);
        const funcCall = this.func.prepare(inputs);
        return await funcCall.call();
      };

      /** Cost function to be optimized */
      const costFunc = async (x: Float32Array): Promise<number> => {
        x.forEach((val, idx) => inputs[variedInputNames[idx]] = val);
        const funcCall = this.func.prepare(inputs);
        const calledFuncCall = await funcCall.call();

        //console.log(x);

        let cost = 0;

        outputsOfInterest.forEach((output) => {
          if (output.prop.propertyType !== DG.TYPE.DATA_FRAME)
            cost += (output.target as number - calledFuncCall.getParamValue(output.prop.name)) ** 2;
          //console.log(output.target as number - calledFuncCall.getParamValue(output.prop.name));
          else {
          //const errs = getErrors(output.colName, output.target as DG.DataFrame, calledFuncCall.getParamValue(output.prop.name));
          //console.log(errs);
            getErrors(output.colName, output.target as DG.DataFrame, calledFuncCall.getParamValue(output.prop.name)).forEach((err) => cost += err ** 2);
          }
        //console.log(output.target);
        //console.log('---------------------------------------------------------------');
        //console.log(calledFuncCall.getParamValue(output.prop.name));
        });

        //console.log(cost);

        return cost;
      };

      //console.log(await costFunc(new Float32Array([1, 5])));

      let extremums: Extremum[];

      // Perform optimization
      if (this.method === METHOD.NELDER_MEAD)
        extremums = await performNelderMeadOptimization(costFunc, minVals, maxVals, this.nelderMeadSettings);
      else
        throw new Error(`The '${this.method}' method has not been implemented.`);

      this.closeDockNodes();

      const extr = extremums[0];
      const rowCount = extr.iterCount;

      // Add fitting results to the table: iteration & loss
      const reportTable = DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.COLUMN_TYPE.INT, TITLE.ITERATIONS, [...Array(rowCount).keys()].map((i) => i + 1)),
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, TITLE.LOSS, extr.iterCosts.slice(0, rowCount)),
      ]);
      this.comparisonView.dataFrame = reportTable;
      const reportColumns = reportTable.columns;
      const fittedRaw: (Float32Array | Int32Array | Float64Array | Uint32Array)[] = [];

      // Add fitting results to the table: fitted parameters
      extr.iterPoints.forEach((vals, idx) => {
        variedInputsCaptions[idx] = reportColumns.getUnusedName(variedInputsCaptions[idx]);
        const col = DG.Column.fromFloat32Array(variedInputsCaptions[idx], new Float32Array(vals.buffer, 0, rowCount));
        fittedRaw.push(col.getRawData());
        reportColumns.add(col);
      });

      // Add PC plot
      const pcPlot = DG.Viewer.pcPlot(reportTable, {
        columnNames: reportTable.columns.names().filter((name) => name !== TITLE.ITERATIONS),
        description: 'Variation',
      });
      const pcPlotDockNode = this.comparisonView.dockManager.dock(pcPlot, DG.DOCK_TYPE.LEFT, this.tableDockNode, '', DOCK_RATIO.PC_PLOT);
      if (pcPlotDockNode.container.dart.elementTitle)
        pcPlotDockNode.container.dart.elementTitle.hidden = true;
      this.openedViewers.push(pcPlot);

      // Add loss function linechart
      const lossLinechart = DG.Viewer.lineChart(reportTable, {
        xColumnName: TITLE.ITERATIONS,
        yColumnNames: [TITLE.LOSS],
        showXSelector: true,
        showYSelectors: false,
        showXAxis: true,
        showYAxis: true,
        description: `${TITLE.LOSS}: sum of squared errors`,
      });
      const lossLinechartDockNode = this.comparisonView.dockManager.dock(lossLinechart, DG.DOCK_TYPE.TOP, pcPlotDockNode, '', DOCK_RATIO.LOSS_PLOT);
      if (lossLinechartDockNode.container.dart.elementTitle)
        lossLinechartDockNode.container.dart.elementTitle.hidden = true;
      this.openedViewers.push(lossLinechart);

      // Goodness of fit (GOF)
      const outputCaptions: string[] = [];
      const outputDf = new Map<string, DG.DataFrame>();
      const outputViewer = new Map<string, DG.Viewer>();

      // Initialize output dataframes & viewers
      outputsOfInterest.forEach((item, idx) => {
        const name = item.prop.name;
        const caption = item.prop.caption ?? name;
        outputCaptions.push(caption);
        const type = item.prop.propertyType;

        switch (type) {
        case DG.TYPE.INT:
        case DG.TYPE.BIG_INT:
        case DG.TYPE.FLOAT:
          const df = DG.DataFrame.fromColumns([
            DG.Column.fromStrings(caption, ['obtained', 'target']),
            DG.Column.fromList(type, TITLE.VALUE, [idx + 1, item.target]),
          ]);
          df.name = caption;
          outputDf.set(name, df);
          outputViewer.set(name, DG.Viewer.barChart(df, {
            valueColumnName: TITLE.VALUE,
            splitColumnName: caption,
            valueAggrType: DG.AGG.AVG,
            showValueSelector: false,
            showCategorySelector: false,
            showStackSelector: false,
          }));
          break;

        case DG.TYPE.DATA_FRAME:
          break;

        default:
          throw new Error(`Unsupported output type: ${item.prop.propertyType}`);
        }
      });

      const x = new Float32Array(dim);

      /** Fit output dataframes with respect to the row of the report table */
      const fillOutputDFs = async (rowIdx: number) => {
        fittedRaw.forEach((arr, idx) => x[idx] = arr[rowIdx]);
        const calledFuncCall = await getCalledFuncCall(x);

        outputsOfInterest.forEach((output) => {
          const name = output.prop.name;
          if (output.prop.propertyType !== DG.TYPE.DATA_FRAME)
            outputDf.get(name)!.set(TITLE.VALUE, 0, calledFuncCall.getParamValue(output.prop.name));
          else {}
        });
      };

      // Fill using data from the last raw
      fillOutputDFs(rowCount - 1);

      // Build UI
      let currentOutput = outputCaptions[0];
      const outputInput = ui.choiceInput(TITLE.OBJECTIVE, currentOutput, outputCaptions, () => {
        currentOutput = outputInput.value!;
        outputViewer.forEach((v, caption) => v.root.hidden = caption !== currentOutput);
      }, {nullable: false});
      outputInput.setTooltip('Output name');
      outputInput.root.hidden = outputsCount < 2;
      const header = ui.h2(currentOutput);
      ui.tooltip.bind(header, 'Output of the function');
      header.style.marginRight = '6px';
      header.hidden = outputsCount > 1;
      const gofItems = [outputInput.root, header];
      outputViewer.forEach((v) => gofItems.push(v.root));
      const gofDiv = ui.divV(gofItems);
      this.goodnessOfFitNode = this.comparisonView.dockManager.dock(gofDiv, DG.DOCK_TYPE.TOP, this.tableDockNode, undefined, DOCK_RATIO.FIT_DIV);
      if (this.goodnessOfFitNode.container.dart.elementTitle)
        this.goodnessOfFitNode.container.dart.elementTitle.hidden = true;

      // Add grid cell effects
      const cellEffect = async (cell: DG.GridCell) => {
        const selectedRow = cell.tableRowIndex ?? 0;
        await fillOutputDFs(selectedRow);
      };

      this.gridClickSubscription = this.comparisonView.grid.onCellClick.subscribe(cellEffect);
      this.gridCellChangeSubscription = this.comparisonView.grid.onCurrentCellChanged.subscribe(cellEffect);
    } catch (error) {
      grok.shell.error(error instanceof Error ? error.message : 'The platform issue');
    }
  } // runOptimizationOld

  /** Get outputs selected as target */
  private getOutputsOfInterest() {
    const outputsOfInterest = [];

    for (const outputName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[outputName];

      if (output.isInterest.value)
        outputsOfInterest.push(output);
    }

    return outputsOfInterest;
  }

  /** Close dock nodes with viewers (visulalizing prev results) */
  private closeDockNodes() {
    this.tempDockNodes.forEach((node) => this.comparisonView.dockManager.close(node));
    this.tempDockNodes.slice(0);
    console.log(this.tempDockNodes);

    /*
    if (this.helpMdNode) {
      this.comparisonView.dockManager.close(this.helpMdNode);
      this.helpMdNode = undefined;
    }

    if (this.goodnessOfFitNode) {
      this.comparisonView.dockManager.close(this.goodnessOfFitNode);
      this.goodnessOfFitNode = undefined;
    }

    if (this.pcPlotNode) {
      this.comparisonView.dockManager.close(this.pcPlotNode);
      this.pcPlotNode = undefined;
    }

    if (this.lossFuncLineChartsNode) {
      this.comparisonView.dockManager.close(this.lossFuncLineChartsNode);
      this.lossFuncLineChartsNode = undefined;
    }

    if (this.gridClickSubscription) {
      this.gridClickSubscription.unsubscribe();
      this.gridClickSubscription = null;
    }

    if (this.gridCellChangeSubscription) {
      this.gridCellChangeSubscription.unsubscribe();
      this.gridCellChangeSubscription = null;
    }*/
  }
}
