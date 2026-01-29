/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {BehaviorSubject, Subject} from 'rxjs';
import {RunComparisonView} from './run-comparison-view';
import {combineLatest} from 'rxjs';
import {take, filter} from 'rxjs/operators';
import '../css/sens-analysis.css';
import {CARD_VIEW_TYPE} from '../../shared-utils/consts';
import {getDefaultValue} from './shared/utils';

import {getCategoryWidget} from './fitting/fitting-utils';
import {getLookupChoiceInput} from './shared/lookup-tools';

import {DiffGrok} from './fitting-view';
import {getOptTypeInput, HELP_URL, STARTING_HELP, TITLE, METHOD, METHOD_HINT,
  getHelpIcon,
  getColorScaleDiv} from './multi-objective-optimization/ui-tools';
import {MoeadManager} from './multi-objective-optimization/moead-manager';
import {OptimizeManager} from './multi-objective-optimization/optimize-manager';
import {getFloatArrays} from './multi-objective-optimization/utils';
import {OPT_TYPE, TooltipInfo} from './multi-objective-optimization/defs';

const RUN_NAME_COL_LABEL = 'Run name' as const;
const supportedOutputTypes = [DG.TYPE.INT, DG.TYPE.BIG_INT, DG.TYPE.FLOAT, DG.TYPE.DATA_FRAME];
type OutputTarget = number | DG.DataFrame | null;

type InputWithValue<T = number> = {input: DG.InputBase, value: T};

type InputValues = {
  isChanging: BehaviorSubject<boolean>,
  const: InputWithValue<boolean | number | string | DG.DataFrame>,
  constForm: DG.InputBase[],
  saForm: DG.InputBase[],
}

type FittingNumericStore = {
  prop: DG.Property,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
  min: InputWithValue,
  max: InputWithValue,
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
  min?: number;
  max?: number;
  step?: number;
}

export type TargetDescription = {
  default?: any;
  argumentCol?: string;
}

type FittingInputsStore = FittingNumericStore | FittingBoolStore | FittingConstStore;

type FittingOutputsStore = {
  prop: DG.Property,
  input: DG.InputBase,
  isInterest: BehaviorSubject<boolean>,
  target: OutputTarget,
};

const getSwitchMock = () => ui.div([], 'sa-switch-input');

export class OptimizationView {
  generateInputFields = (func: DG.Func) => {
    const getInputValue = (input: DG.Property, key: keyof RangeDescription) => {
      const range = this.options.ranges?.[input.name];
      if (range?.[key] != undefined)
        return range[key];
      return input.options[key] === undefined ? getDefaultValue(input) : Number(input.options[key]);
    };

    const getSwitchElement = (defaultValue: boolean, f: (v: boolean) => any, isInput: boolean = true) => {
      const input = ui.input.toggle(' ', {value: defaultValue, onValueChanged: (value) => f(value)});
      $(input.root).addClass('sa-switch-input');
      $(input.captionLabel).hide();

      ui.tooltip.bind(input.root, () => {
        if (isInput) {
          return (input.value) ?
            `Set the input as constant` :
            `Set the input as variable`;
        } else {
          return !input.value ?
            `Mark the output as an objective` :
            `Mark the output as not an objective`;
        }
      });

      return input;
    };

    const inputs = func.inputs.reduce((acc, inputProp) => {
      const defaultValue = getInputValue(inputProp, 'default');

      if (inputProp.propertyType === DG.TYPE.FLOAT) {
        const isChangingInputMin = getSwitchElement(false, (v: boolean) => {
          ref.isChanging.next(v);
          this.updateApplicabilityState();
        });

        const isChangingInputConst = getSwitchElement(false, (v: boolean) => {
          ref.isChanging.next(v);
          this.updateApplicabilityState();
        });

        const caption = inputProp.caption ?? inputProp.name;

        const temp = {
          type: inputProp.propertyType,
          prop: inputProp,
          const: {
            input:
            (() => {
              const inp = ui.input.float(caption, {value: defaultValue, onValueChanged: (value) => {
                ref.const.value = value;
                this.updateApplicabilityState();
              }});
              inp.root.insertBefore(isChangingInputConst.root, inp.captionLabel);
              inp.addPostfix(inputProp.options['units']);
              inp.setTooltip(`Value of '${caption}'`);
              inp.nullable = false;
              return inp;
            })(),
            value: defaultValue,
          },
          min: {
            input:
              (() => {
                const inp = ui.input.float(`${caption} (min)`, {value: getInputValue(inputProp, 'min'), onValueChanged: (value) => {
                  (ref as FittingNumericStore).min.value = value;
                  this.updateApplicabilityState();
                }});
                inp.addValidator((s:string) => (Number(s) > temp.max.value) ? 'Greater than max': null);
                inp.root.insertBefore(isChangingInputMin.root, inp.captionLabel);
                inp.addPostfix(inputProp.options['units']);
                inp.setTooltip(`Min value '${caption}'`);
                inp.nullable = false;
                return inp;
              })(),
            value: getInputValue(inputProp, 'min'),
          },
          max: {
            input: (() => {
              const inp = ui.input.float(`${caption} (max)`, {value: getInputValue(inputProp, 'max'), onValueChanged: (value) => {
                (ref as FittingNumericStore).max.value = value;
                this.updateApplicabilityState();
              }});
              inp.addValidator((s:string) => (Number(s) < temp.min.value) ? 'Smaller than min': null);
              inp.addPostfix(inputProp.options['units']);
              inp.setTooltip(`Max value of '${caption}'`);
              inp.nullable = false;
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
        } as FittingNumericStore;

        const ref = acc[inputProp.name] as FittingNumericStore;
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
            temp.caption = inputProp.caption ?? inputProp.name;
            temp.onInput.subscribe(() => {
              tempDefault.value = temp.value;
              this.updateApplicabilityState();
            });
            temp.root.insertBefore(switchMock, temp.captionLabel);
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
          type: inputProp.propertyType,
          prop: inputProp,
          isChanging: new BehaviorSubject(false),
        } as FittingConstStore;
      }

      return acc;
    }, {} as Record<string, FittingInputsStore>);

    const outputs = func.outputs.filter((prop) => supportedOutputTypes.includes(prop.propertyType))
      .reduce((acc, outputProp) => {
        const temp: FittingOutputsStore = {
          prop: outputProp,
          input:
          (() => {
            const input = ui.input.int(outputProp.caption ?? outputProp.name);
            input.setTooltip((outputProp.propertyType === DG.TYPE.DATA_FRAME) ? 'Target dataframe' : 'Target scalar');
            ui.tooltip.bind(input.captionLabel, (outputProp.propertyType === DG.TYPE.DATA_FRAME) ? 'Dataframe' : 'Scalar');

            if (this.options.targets?.[outputProp.name]?.default != null)
              setTimeout(() => input.value = this.options.targets?.[outputProp.name]?.default, 0);

            const isInterestInput = getSwitchElement(
              true,
              (v: boolean) => {
                temp.isInterest.next(v);
                this.updateApplicabilityState();
              },
              false,
            ).root;
            input.root.insertBefore(isInterestInput, input.captionLabel);
            input.input.hidden = true;

            return input;
          })(),
          isInterest: new BehaviorSubject<boolean>(true),
          target: (outputProp.propertyType !== DG.TYPE.DATA_FRAME) ? 0 : null,
        };

        acc[outputProp.name] = temp;

        return acc;
      }, {} as Record<string, FittingOutputsStore>);

    return {inputs, outputs};
  };

  private openedViewers: DG.Viewer[] = [];
  private readyToRun = false;
  private isOptimizationRunning = false;

  private optTypeInput = getOptTypeInput();

  private progressIndicator: DG.TaskBarProgressIndicator | undefined = undefined;

  private runIcon = ui.iconFA('play', async () => {
    if (this.readyToRun) {
      this.isOptimizationRunning = true;
      this.updateApplicabilityState();
      this.updateRunIconDisabledTooltip('In progress...');
      this.progressIndicator = DG.TaskBarProgressIndicator.create('Optimization...', {cancelable: true});

      await this.runOptimization();

      this.progressIndicator.close();
      this.isOptimizationRunning = false;
      this.updateApplicabilityState();
    }
  });

  private acceptIcon = ui.iconFA('ballot-check', async () => {
    const choiceItems = Array.from({length: this.currentFuncCalls.length}, (_, i) => i + 1);
    let chosenItem = -1;
    const input = ui.input.choice('Select optimization', {items: choiceItems, onValueChanged: (x) => chosenItem = x});
    const confirmed = await new Promise((resolve, _reject) => {
      ui.dialog({title: 'Accept optimization'})
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
  });

  private helpIcon = getHelpIcon();

  private gridClickSubscription: any = null;
  private gridCellChangeSubscription: any = null;

  // Auxiliary dock nodes with results
  private helpDN: DG.DockNode | undefined = undefined;


  store = this.generateInputFields(this.func);
  comparisonView!: DG.TableView;


  private diffGrok: DiffGrok | undefined = undefined;

  private currentFuncCalls: DG.FuncCall[] = [];
  private isFittingAccepted = false;
  public acceptedFitting$ = new Subject<DG.FuncCall | null>();

  private method: METHOD = METHOD.MOEAD;
  private methodInput = ui.input.choice<METHOD>(TITLE.METHOD, {
    value: METHOD.MOEAD,
    nullable: false,
    items: [METHOD.MOEAD],
    tooltipText: METHOD_HINT.get(METHOD.MOEAD),
    onValueChanged: (val) => {
      this.method = val;
      this.methodInput.setTooltip(METHOD_HINT.get(val) ?? '');
      this.updateMetSetInputs();
      this.updateApplicabilityState();
    },
  });

  private optManagers = new Map<METHOD, OptimizeManager>([
    [METHOD.MOEAD, new MoeadManager(this)],
  ]);

  private optSetInps = new Map<METHOD, HTMLElement[]>();

  private settingsHeader = ui.h1(TITLE.SET);
  private areSettingsInputsShown = false;
  private fittingSettingsIcon = ui.iconFA('cog', () => {
    this.areSettingsInputsShown = !this.areSettingsInputsShown;
    this.updateMetSetInputs();
  });

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
    baseView: RunComparisonView,
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
    grok.events.onViewRemoved.pipe(filter((v) => v.id === baseView.id), take(1)).subscribe(() => {
      if (options.acceptMode && !this.isFittingAccepted) {
        this.acceptedFitting$.next(null);
        this.isFittingAccepted = true;
      }
    });

    this.buildForm(options.inputsLookup).then((form) => {
      this.comparisonView = baseView;

      this.comparisonView.dockManager.dock(
        form,
        DG.DOCK_TYPE.LEFT,
        null,
        `${this.func.name} - Optimization`,
        0.25,
      );

      this.comparisonView.grid.columns.byName(RUN_NAME_COL_LABEL)!.visible = false;

      const rbnPanels = [[this.helpIcon, this.runIcon, ...(this.options.acceptMode ? [this.acceptIcon] : [])]];
      this.comparisonView.setRibbonPanels(rbnPanels);

      this.comparisonView.name = this.comparisonView.name.replace('comparison', 'optimization');
      this.comparisonView.helpUrl = HELP_URL;
      const helpMD = ui.markdown(STARTING_HELP);
      helpMD.style.padding = '10px';
      helpMD.style.overflow = 'auto';
      this.helpDN = this.comparisonView.dockManager.dock(
        helpMD,
        DG.DOCK_TYPE.FILL,
        this.comparisonView.dockManager.findNode(this.comparisonView.grid.root),
        'About',
      );

      this.updateRunIconStyle();
      this.runIcon.classList.add('fas');

      this.updateApplicabilityState();
    });

    this.diffGrok = options.diffGrok;
  } // constructor

  /** Check correctness of the method settings */
  private areMethodSettingsCorrect(): boolean {
    return true;
  } // areMethodSettingsCorrect

  /** Return value lookup widget */
  private async getLookupElement(inputsLookup?: string) {
    if (inputsLookup === undefined)
      return null;

    const constIputs = new Map<string, DG.InputBase>();
    Object.keys(this.store.inputs).forEach((name) => constIputs.set(name, this.store.inputs[name].constForm[0]));

    const lookupElement = await getLookupChoiceInput(inputsLookup, constIputs);

    if (lookupElement !== null)
      lookupElement.input.root.insertBefore(getSwitchMock(), lookupElement.input.captionLabel);

    return lookupElement;
  }

  /** Build form with inputs */
  private async buildForm(inputsLookup?: string) {
    // 0. The main form
    const form = ui.div([this.optTypeInput], {style: {overflowY: 'scroll', width: '100%'}});

    // 1. Outputs of the function
    let prevCategory = 'Misc';

    Object.values(this.store.outputs)
      .reduce((container, outputConfig) => {
        const prop = outputConfig.prop;
        if (prop.category !== prevCategory) {
          container.append(ui.h3(prop.category));
          prevCategory = prop.category;
        }

        container.append(outputConfig.input.root);

        return container;
      }, form);

    // 2. Inputs of the function
    const fitHeader = ui.h1(TITLE.CONSTR);
    ui.tooltip.bind(fitHeader, 'Set inputs constraints');

    // inputs grouped by categories
    const inputsByCategories = new Map<string, HTMLElement[]>([['Misc', []]]);

    // group inputs by categories
    Object.values(this.store.inputs).forEach((inputConfig) => {
      const category = inputConfig.prop.category;
      const roots = [...inputConfig.constForm.map((input) => input.root), ...inputConfig.saForm.map((input) => input.root)];

      if (inputsByCategories.has(category))
        inputsByCategories.get(category)!.push(...roots);
      else
        inputsByCategories.set(category, roots);
    });

    form.append(fitHeader);

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
    if (inputsByCategories.size > 1) {
      if (topCategory !== null) {
        const roots = inputsByCategories.get(topCategory);
        form.append(getCategoryWidget(topCategory, roots!));
        form.append(...roots!);
      }

      inputsByCategories.forEach((roots, category) => {
        if ((category !== 'Misc') && (category !== topCategory)) {
          form.append(getCategoryWidget(category, roots));
          form.append(...roots);
        }
      });

      if (topCategory !== 'Misc') {
        const miscRoots = inputsByCategories.get('Misc');

        if (miscRoots!.length > 0) {
          const roots = inputsByCategories.get('Misc')!;
          form.append(getCategoryWidget('Misc', roots));
          form.append(...roots);
        }
      }
    } else
      form.append(...inputsByCategories.get('Misc')!);

    // 3. Make one output of interest
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
    }

    // 4. Method
    form.append(ui.h1(TITLE.USING));
    this.fittingSettingsIcon.style.minWidth = '50px';
    this.methodInput.root.insertBefore(this.fittingSettingsIcon, this.methodInput.captionLabel);
    form.append(this.methodInput.root);

    // 5. Method settings
    form.append(this.settingsHeader);

    this.optManagers.forEach((manager, name) => {
      const inputs = manager.getInputs();

      inputs.forEach((inp) => {
        inp.root.insertBefore(getSwitchMock(), inp.captionLabel);
        form.append(inp.root);
      });

      this.optSetInps.set(name, inputs.map((inp) => inp.root));
    });

    this.updateMetSetInputs();

    $(form).addClass('ui-form');

    this.updateApplicabilityState();

    $(form).css({
      'padding-left': '12px',
      'overflow-y': 'scroll',
      'padding-right': '4px',
      'position': 'relative',
      'z-index': '2',
      'background-color': 'white',
    });
    return form;
  } // buildForm

  /** Update inputs for methods settings */
  private updateMetSetInputs(): void {
    this.updateSettingsIconsHint();
    this.settingsHeader.hidden = !this.areSettingsInputsShown;

    this.optSetInps.forEach((inputs, method) => {
      inputs.forEach((root) => root.hidden = !((method === this.method) && this.areSettingsInputsShown));
    });
  }

  /** Update settings icons tooltip */
  private updateSettingsIconsHint(): void {
    ui.tooltip.bind(
      this.fittingSettingsIcon,
      this.areSettingsInputsShown ?
        'Hide settings' :
        'Modify settings of the method',
    );
  }

  /** Update run icon tooltip: disabled case */
  private updateRunIconDisabledTooltip(msg: string): void {
    ui.tooltip.bind(this.runIcon, () => {
      const label = ui.label(msg);
      label.style.color = '#FF0000';
      return label;
    });
  } // updateRunIconDisabledTooltip

  /** Check applicability of fitting */
  public updateApplicabilityState(): void {
    this.readyToRun = this.canOptimizationBeRun() && (!this.isOptimizationRunning);
    this.updateRunIconStyle();
  } // updateApplicabilityState

  private updateRunIconStyle(): void {
    if (this.readyToRun) {
      this.runIcon.style.color = 'var(--green-2)';
      ui.tooltip.bind(this.runIcon, 'Run optimization');
    } else
      this.runIcon.style.color = 'var(--grey-3)';
  } // updateRunIconStyle

  /** Check inputs */
  private areInputsReady(): boolean {
    let isAnySelected = false;
    let areSelectedFilled = true;
    let cur: boolean;

    for (const propName of Object.keys(this.store.inputs)) {
      const input = this.store.inputs[propName];
      const caption = input.prop.caption ?? input.prop.name;

      if (input.isChanging.value === true) {
        isAnySelected = true;
        const minInp = (input as FittingNumericStore).min.input;
        const maxInp = (input as FittingNumericStore).max.input;
        const min = minInp.value;
        const max = maxInp.value;

        cur = (min !== null) && (min !== undefined) && (max !== undefined) && (max !== undefined);

        if (cur) {
          if (min > max) {
            cur = false;
            this.updateRunIconDisabledTooltip(`Invalid min & max of "${caption}"`);
            minInp.input.classList.add('d4-invalid');
            maxInp.input.classList.add('d4-invalid');
          } else {
            minInp.input.classList.remove('d4-invalid');
            maxInp.input.classList.remove('d4-invalid');
          }
        } else
          this.updateRunIconDisabledTooltip(`Incomplete "${caption}"`);

        areSelectedFilled = areSelectedFilled && cur;
      } else {
        const val = input.const.input.value;
        cur = (val !== null) && (val !== undefined);

        if (!cur)
          this.updateRunIconDisabledTooltip(`Incomplete "${caption}"`);

        areSelectedFilled = areSelectedFilled && cur;
      }
    }

    if (!isAnySelected)
      this.updateRunIconDisabledTooltip(`No inputs selected as variable`);

    return isAnySelected && areSelectedFilled;
  } // areInputsReady

  /** Check outputs */
  private areOutputsSelected(): boolean {
    for (const propName of Object.keys(this.store.outputs)) {
      const output = this.store.outputs[propName];

      if (output.isInterest.value === true)
        return true;
    }

    this.updateRunIconDisabledTooltip(`No objectives selected to be optimized`);

    return false;
  } // areOutputsReady

  /** Check method's settings */
  private areSettingsOfMethodValid(): boolean {
    const optManager = this.optManagers.get(this.method);

    if (optManager === undefined) {
      this.updateRunIconDisabledTooltip('Invalid method');
      return true;
    }

    const validRes = optManager.areSettingsValid();

    if (validRes.res)
      return true;

    this.updateRunIconDisabledTooltip(validRes.msg);
    return false;
  }

  /** Check inputs/outputs/settings */
  private canOptimizationBeRun(): boolean {
    return this.areSettingsOfMethodValid() && this.areInputsReady() &&
      this.areOutputsSelected() && this.areMethodSettingsCorrect();
  }

  /** Return names of the fixed inputs */
  private getFixedInputs() {
    return Object.keys(this.store.inputs).filter((propName) => !this.store.inputs[propName].isChanging.value);
  }

  /** Return names of the varied inputs */
  private getVariadInputs() {
    return Object.keys(this.store.inputs)
      .filter((propName) => (this.store.inputs[propName].type === DG.TYPE.FLOAT) && this.store.inputs[propName].isChanging.value);
  }

  /** Perform optimization */
  private async runOptimization(): Promise<void> {
    try {
      // check applicability
      if (!this.canOptimizationBeRun())
        return;

      // inputs of the source function
      const inputs: any = {};

      // add fixed inputs
      const fixedInputs = this.getFixedInputs();
      fixedInputs.forEach((name) => inputs[name] = this.store.inputs[name].const.value);

      // get varied inputs, optimization is performed with respect to them
      const variedInputs = this.getVariadInputs();
      const inputDim = variedInputs.length;

      // varied inputs specification
      const variedInputNames: string[] = [];
      const minVals = new Float32Array(inputDim);
      const maxVals = new Float32Array(inputDim);
      const variedInputsCaptions = new Array<string>(inputDim);
      const variedInputsTypes = new Array<string>(inputDim);

      // set varied inputs specification
      variedInputs.forEach((name, idx) => {
        const propConfig = this.store.inputs[name] as FittingNumericStore;
        minVals[idx] = propConfig.min.value ?? 0;
        maxVals[idx] = propConfig.max.value ?? 0;
        variedInputNames.push(name);
        variedInputsCaptions[idx] = propConfig.prop.caption ?? propConfig.prop.name;
        variedInputsTypes[idx] = propConfig.prop.propertyType;
      });

      // get selected output
      const outputsOfInterest = this.getOutputsOfInterest();
      const outputsCount = outputsOfInterest.length;
      if (outputsCount < 1) {
        grok.shell.error('No output is selected for optimization.');
        return;
      }

      let outputDim = 0;
      const outputInfo: {caption: string, hint: string}[] = [];

      // Single call for output specification
      minVals.forEach((val, idx) => inputs[variedInputNames[idx]] = val);
      const funcCall = this.func.prepare(inputs);
      const calledFuncCall = await funcCall.call(undefined, undefined, {processed: true, report: false});

      outputsOfInterest.forEach((output) => {
        const caption = output.prop.caption ?? output.prop.name;

        if (output.prop.propertyType !== DG.TYPE.DATA_FRAME) {
          outputInfo.push({
            caption: caption,
            hint: `**Objective**\n\nValues of the output scalar **"${caption}"**`,
          });
          ++outputDim;
        } else {
          const df = calledFuncCall.getParamValue(output.prop.name) as DG.DataFrame;
          const rowCount = df.rowCount;

          for (const col of df.columns) {
            if (!col.isNumerical)
              continue;

            for (let i = 0; i < rowCount; ++i) {
              outputInfo.push({
                caption: `${caption}["${col.name}", ${i + 1}]`,
                hint: `**Objective**\n\nValues of the output dataframe **"${caption}"**:\n\n - column: "${col.name}"\n\n - row: ${i + 1}`,
              });
              ++outputDim;
            }
          }
        }
      });

      const optType: OPT_TYPE = this.optTypeInput.value!;
      const sign = (optType === OPT_TYPE.MIN) ? 1 : -1;

      /** The optimization objective */
      const objective = async (x: Float32Array): Promise<Float32Array> => {
        // return new Float32Array([
        //   -(x[0]**2 + x[1]**2) * sign,
        //   -((x[0] - 2)**2 + (x[1] - 1)**2) * sign,
        //   -((x[0] - 1)**2 + (x[1] - 2)**2) * sign,
        //   -((x[0] - 0.3)**2 + (x[1] - 0.7)**2) * sign,
        // ]);

        x.forEach((val, idx) => inputs[variedInputNames[idx]] = val);
        const funcCall = this.func.prepare(inputs);
        const calledFuncCall = await funcCall.call(undefined, undefined, {processed: true, report: false});

        const result = new Float32Array(outputDim);
        let curOutIdx = 0;

        outputsOfInterest.forEach((output) => {
          if (output.prop.propertyType !== DG.TYPE.DATA_FRAME) {
            result[curOutIdx] = sign * (calledFuncCall.getParamValue(output.prop.name) as number);
            ++curOutIdx;
          } else {
            const df = calledFuncCall.getParamValue(output.prop.name) as DG.DataFrame;
            const rowCount = df.rowCount;

            for (const col of df.columns) {
              if (!col.isNumerical)
                continue;

              const raw = col.getRawData();

              for (let i = 0; i < rowCount; ++i) {
                result[curOutIdx] = sign * raw[i];
                ++curOutIdx;
              }
            }
          }
        });

        if (outputDim !== curOutIdx)
          throw new Error('Inconsistent outputs: dataframes must be of the same size');

        return new Float32Array(result);
      };

      // Perform optimization
      const optManager = this.optManagers.get(this.method);
      if (optManager === undefined)
        throw new Error('Non-supported method');

      const solution = await optManager.perform(
        objective,
        {
          dim: inputDim,
          mins: minVals,
          maxs: maxVals,
        },
        outputDim,
        this.progressIndicator,
      );

      this.clearPrev();

      // Prepare resulting dataframe
      const solutionsCount = solution.length;
      const inpRaw = getFloatArrays(inputDim, solutionsCount);
      const outRaw = getFloatArrays(outputDim, solutionsCount);

      for (let rowIdx = 0; rowIdx < solutionsCount; ++rowIdx) {
        const input = solution[rowIdx];

        for (let colIdx = 0; colIdx < inputDim; ++colIdx)
          inpRaw[colIdx][rowIdx] = input[colIdx];

        const output = await objective(input);

        for (let colIdx = 0; colIdx < outputDim; ++colIdx)
          outRaw[colIdx][rowIdx] = output[colIdx] * sign;
      }

      const resulDataframe = DG.DataFrame.fromColumns(
        inpRaw.map((raw, idx) => DG.Column.fromFloat32Array(`inp ${idx}`, raw, solutionsCount))
          .concat(outRaw.map((raw, idx) => DG.Column.fromFloat32Array(`out ${idx}`, raw, solutionsCount))),
      );

      // Set columns' names & tooltips
      const colTooltips = new Map<string, TooltipInfo>();

      const resCols = resulDataframe.columns;
      variedInputsCaptions.forEach((name, idx) => {
        const unusedName = resCols.getUnusedName(name);
        resCols.byName(`inp ${idx}`).name = unusedName;
        colTooltips.set(unusedName, {
          msg: `**Decision variable**\n\nValues of the input scalar **"${name}"**`,
          isInput: true,
        });
      });

      outputInfo.forEach((info, idx) => {
        const col = resCols.byName(`out ${idx}`);
        const unusedName = resCols.getUnusedName(info.caption);
        col.name = unusedName;
        colTooltips.set(unusedName, {
          msg: info.hint,
          isInput: false,
        });
      });

      this.comparisonView.dataFrame = resulDataframe;
      const grid = this.comparisonView.grid;

      // Add tooltips
      grid.onCellTooltip(function(cell, x, y) {
        if (cell.isColHeader) {
          const cellCol = cell.tableColumn;
          if (cellCol) {
            const info = colTooltips.get(cell.tableColumn.name) ?? {msg: '', isInput: true};
            const elems = [ui.markdown(info.msg)];

            if (!info.isInput)
              elems.push(getColorScaleDiv(optType));

            ui.tooltip.show(ui.divV(elems), x, y);
            return true;
          }
        }
      });

      // Visualize results
      this.openedViewers = optManager.visualize(this.comparisonView, inputDim, outputDim, optType);
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

  /** Clear previous results: close dock nodes and unsubscribe from events */
  private clearPrev(): void {
    for (const v of this.openedViewers)
      v.close();

    this.openedViewers.splice(0);

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

    this.comparisonView.grid.sort([], []);
  } // clearPrev
}
