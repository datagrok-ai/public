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

const RUN_NAME_COL_LABEL = 'Run name' as const;
const supportedOutputTypes = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.DATA_FRAME];

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
            `Switch OFF fitting this input` :
            `Switch ON fitting this input`;
        } else {
          return !input.value ?
            `Mark this output as a target` :
            `Mark this output as not a target`;
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

            input.onChanged(() => temp.target = input.value);

            if (outputProp.propertyType === DG.TYPE.DATA_FRAME)
              (input.root.lastElementChild as HTMLDivElement).hidden = !this.toSetSwitched;

            const isInterestInput = supportedOutputTypes.includes(outputProp.propertyType) ?
              getSwitchElement(
                this.toSetSwitched,
                (v: boolean) => {
                  temp.isInterest.next(v);
                  this.updateRunWidgetsState();
                  input.input.hidden = !v;
                  input.setTooltip(v ? 'Target value' :
                    (outputProp.propertyType === DG.TYPE.DATA_FRAME) ? 'Output dataframe' : 'Output scalar');
                  if (outputProp.propertyType === DG.TYPE.DATA_FRAME)
                    (input.root.lastElementChild as HTMLDivElement).hidden = !v;
                },
                false,
              ).root: getSwitchMock();
            input.root.insertBefore(isInterestInput, input.captionLabel);
            this.toSetSwitched = false;

            return input;
          })(),
          isInterest: new BehaviorSubject<boolean>(this.toSetSwitched),
          target: null,
        };

        acc[outputProp.name] = temp;

        return acc;
      }, {} as Record<string, {
      prop: DG.Property,
      input: DG.InputBase,
      isInterest: BehaviorSubject<boolean>,
      target: number | DG.DataFrame | null,
    }>);

    return {inputs, outputs};
  };

  private runButton = ui.bigButton('Run', async () => await this.runOptimization(), 'Run fitting');
  private runIcon = ui.iconFA('play', async () => await this.runOptimization(), 'Run fitting');
  private helpIcon = ui.iconFA('question', () => {
    window.open('https://datagrok.ai/help/compute.md#fitting', '_blank');
  }, 'Open help in a new tab');
  private tableDockNode: DG.DockNode | undefined;
  private helpMdNode: DG.DockNode | undefined;
  private gridSubscription: any = null;
  private toSetSwitched = true;

  store = this.generateInputFields(this.func);
  comparisonView!: DG.TableView;

  // Optimization settings: TO ADD METHOD's SETTINGS HERE


  private methodSettingsDiv = ui.divV([ui.label('TO BE ADDED')]); // HERE, TO ADD UI for modifying optimizer settings

  private optSettingsIcon = ui.iconFA('cog', () => {
    const prevState = this.optSettingsDiv.hidden;
    this.optSettingsDiv.hidden = !prevState;
  }, 'Modify optimization settings');

  private optSettingsDiv = ui.divV([
    ui.h2('Settings'),
    this.methodSettingsDiv,
  ]);

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

    this.runIcon = ui.iconFA('play', async () => await this.runOptimization(), 'Run fitting');
    this.runIcon.style.color = 'var(--green-2)';
    this.runIcon.classList.add('fas');

    this.helpIcon = ui.iconFA('question', () => {
      window.open('https://datagrok.ai/help/compute.md#fitting', '_blank');
    }, 'Open help in a new tab');

    this.optSettingsDiv.hidden = true;

    const form = this.buildFormWithBtn();
    this.runButton.disabled = !this.canEvaluationBeRun();
    this.runIcon.hidden = this.runButton.disabled;
    this.comparisonView = baseView;

    this.comparisonView.dockManager.dock(
      form,
      DG.DOCK_TYPE.LEFT,
      null,
      `${this.func.name} - Sensitivity Analysis`,
      0.25,
    );

    this.comparisonView.grid.columns.byName(RUN_NAME_COL_LABEL)!.visible = false;

    const rbnPanels = this.comparisonView.getRibbonPanels();
    rbnPanels.push([this.helpIcon, this.runIcon]);
    this.comparisonView.setRibbonPanels(rbnPanels);

    this.comparisonView.name = this.comparisonView.name.replace('comparison', 'optimization');
    this.comparisonView.helpUrl = 'https://datagrok.ai/help/compute.md#fitting';
    this.tableDockNode = this.comparisonView.dockManager.findNode(this.comparisonView.grid.root);
    const helpMD = ui.markdown(STARTING_HELP);
    helpMD.style.padding = '10px';
    helpMD.style.overflow = 'auto';
    this.helpMdNode = this.comparisonView.dockManager.dock(helpMD, DG.DOCK_TYPE.FILL, this.tableDockNode, 'About');
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
        ui.h2('Fit'),
      ], {style: {'overflow-y': 'scroll', 'width': '100%'}}));

    form.appendChild(this.optSettingsDiv);
    form.appendChild(ui.h2('To get'));
    prevCategory = 'Misc';

    Object.values(this.store.outputs)
      .reduce((container, outputConfig) => {
        const prop = outputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.p(prop.category));
          prevCategory = prop.category;
        }

        container.append(outputConfig.input.root);

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

    form.appendChild(ui.h2('Using'));

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
    return true; //this.targetCount > 0;
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
