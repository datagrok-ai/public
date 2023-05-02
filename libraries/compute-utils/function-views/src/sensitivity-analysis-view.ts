/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {RunComparisonView} from './run-comparison-view';
import {BehaviorSubject} from 'rxjs';

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
  min: number, max: number, lvl: number, distrib: DISTRIB_TYPE, isChanging: BehaviorSubject<boolean>,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
};

type SensitivityBoolValues = {
  default: boolean, isChanging: BehaviorSubject<boolean>,
  type: DG.TYPE.BOOL,
}

type SensitivityConstValues = {
  const: DG.InputBase,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL>,
}

type SensitivityValues = SensitivityNumericValues | SensitivityBoolValues | SensitivityConstValues;

type SensitivityNumericInputs = {
  minInput: DG.InputBase,
  maxInput: DG.InputBase,
  lvlInput: DG.InputBase,
  distribInput: DG.InputBase,
  isChangingInput: DG.InputBase<boolean | null>,
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
}

type SensitivityBoolInputs = {
  defaultInput: DG.InputBase,
  isChangingInput: DG.InputBase<boolean | null>,
  type: DG.TYPE.BOOL,
}

type SensitivityConstInputs = {
  constInput: DG.InputBase,
  type: Exclude<DG.TYPE, DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT | DG.TYPE.BOOL>,
}

type SensitvityInputs = SensitivityNumericInputs | SensitivityBoolInputs | SensitivityConstInputs;

const generateFormRows = (func: DG.Func) => {
  const store = (func.inputs.map((inputProp) => {
    switch (inputProp.propertyType) {
    case DG.TYPE.INT:
    case DG.TYPE.BIG_INT:
    case DG.TYPE.FLOAT:
      return {
        min: inputProp.defaultValue ?? 0,
        max: inputProp.defaultValue ?? 0,
        lvl: 3,
        distrib: DISTRIB_TYPES[0],
        isChanging: new BehaviorSubject(false),
        type: inputProp.propertyType,
      };
    case DG.TYPE.BOOL:
      return {
        default: inputProp.defaultValue,
        isChanging: new BehaviorSubject(false),
        type: inputProp.propertyType,
      };
    default:
      return {
        constValue: inputProp.defaultValue,
        isChanging: new BehaviorSubject(false),
        type: inputProp.propertyType,
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
        minInput: inputProp.propertyType === DG.TYPE.FLOAT ?
          ui.floatInput(`${inputProp.caption ?? inputProp.name} min`, inputProp.defaultValue, (v: number) => numericStore.min = v):
          ui.intInput(`${inputProp.caption ?? inputProp.name} min`, inputProp.defaultValue, (v: number) => numericStore.min = v),
        maxInput: inputProp.propertyType === DG.TYPE.FLOAT ?
          ui.floatInput(`${inputProp.caption ?? inputProp.name} max`, inputProp.defaultValue, (v: number) => numericStore.max = v):
          ui.intInput(`${inputProp.caption ?? inputProp.name} max`, inputProp.defaultValue, (v: number) => numericStore.max = v),
        lvlInput: ui.intInput('Levels', 3, (v: number) => numericStore.lvl = v),
        distribInput: ui.choiceInput('Distribution', DISTRIB_TYPES[0], DISTRIB_TYPES, (v: DISTRIB_TYPE) => numericStore.distrib = v),
        isChangingInput: ui.boolInput('Is changing', numericStore.isChanging.value, (v: boolean) => numericStore.isChanging.next(v)),
        type: inputProp.propertyType,
      };
      const ref = acc[inputProp.name] as SensitivityNumericInputs;
      numericStore.isChanging.subscribe((newVal) => {
        if (newVal) {
          $(ref.maxInput.root).show();
          $(ref.lvlInput.root).show();
          $(ref.distribInput.root).show();
        } else {
          $(ref.maxInput.root).hide();
          $(ref.lvlInput.root).hide();
          $(ref.distribInput.root).hide();
        }
      });
      break;
    case DG.TYPE.BOOL:
      const boolStore = store[inputProp.name] as SensitivityBoolValues;
      acc[inputProp.name] = {
        defaultInput: ui.boolInput('Default value', false, (v: boolean) => boolStore.default = v),
        isChangingInput: ui.boolInput('Is changing', boolStore.isChanging.value, (v: boolean) => boolStore.isChanging.next(v)),
        type: inputProp.propertyType,
      };
      break;
    default:
      const constStore = store[inputProp.name] as SensitivityConstValues;
      acc[inputProp.name] = {
        constInput: ui.input.forProperty(inputProp, undefined, {onValueChanged: (v: any) => constStore.const = v}),
        type: inputProp.propertyType,
      };
    }

    return acc;
  }, {} as Record<string, SensitvityInputs>);

  return {store, inputs};
};

export class SensitivityAnalysisView extends DG.ViewBase {
  static openInDockMode(
    func: DG.Func,
    options: {
      name?: string,
    } = {name: undefined},
  ) {
    const saView = new this(func, options);
    grok.shell.dockManager.dock(saView.root, DG.DOCK_TYPE.LEFT, null, this.name, 0.25);
  }

  inputConfig = generateFormRows(this.func);
  leftSide = ui.div();

  constructor(
    public func: DG.Func,
    public options: {
      name?: string,
    } = {name: undefined},
  ) {
    super();
    this.name = this.options.name ?? `${this.func.name} - Sensitivity Analysis`;
    this.box = true;

    this.leftSide = this.buildLeftSide();
    this.leftSide.style.maxWidth = '100%';
    this.root.appendChild(this.leftSide);
  }

  private buildLeftSide() {
    return ui.divV([
      Object.values(this.inputConfig.inputs).reduce((acc, propConfig) => {
        let propInputs = [];

        switch (propConfig.type) {
        case DG.TYPE.INT:
        case DG.TYPE.BIG_INT:
        case DG.TYPE.FLOAT:
          propConfig = propConfig as SensitivityNumericInputs;
          propInputs = [
            propConfig.isChangingInput.root,
            propConfig.minInput.root,
            propConfig.maxInput.root,
            propConfig.lvlInput.root,
            propConfig.distribInput.root,
          ];
          break;
        case DG.TYPE.BOOL:
          propConfig = propConfig as SensitivityBoolInputs;
          propInputs = [
            propConfig.isChangingInput.root,
            propConfig.defaultInput.root,
          ];
          break;
        default:
          propInputs = [
            propConfig.constInput.root,
          ];
        }

        acc.appendChild(ui.divH(propInputs));
        return acc;
      }, ui.divV([], 'ui-form ui-form-condensed')),
      ui.bigButton('Run sensitivity analysis', async () => this.run()),
    ]);
  }

  private async run() {
    const paramValues = Object.keys(this.inputConfig.store).reduce((acc, propName) => {
      switch (this.inputConfig.store[propName].type) {
      case DG.TYPE.INT:
      case DG.TYPE.BIG_INT:
        const numPropConfig = this.inputConfig.store[propName] as SensitivityNumericValues;
        const intStep = (numPropConfig.max - numPropConfig.min) / (numPropConfig.lvl - 1);
        acc[propName] = numPropConfig.isChanging.value ?
          Array.from({length: numPropConfig.lvl}, (_, i) => Math.round(numPropConfig.min + i*intStep)) :
          [numPropConfig.min];
        break;
      case DG.TYPE.FLOAT:
        const floatPropConfig = this.inputConfig.store[propName] as SensitivityNumericValues;
        const floatStep = (floatPropConfig.max - floatPropConfig.min) / (floatPropConfig.lvl - 1);
        acc[propName] = floatPropConfig.isChanging.value ?
          Array.from({length: floatPropConfig.lvl}, (_, i) => floatPropConfig.min + i*floatStep) :
          [floatPropConfig.min];
        break;
      case DG.TYPE.BOOL:
        const boolPropConfig = this.inputConfig.store[propName] as SensitivityBoolValues;
        acc[propName] = boolPropConfig.isChanging.value ?
          [boolPropConfig.default, !boolPropConfig.default]:
          [boolPropConfig.default];
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

    const v = await RunComparisonView.fromComparedRuns(calledFuncCalls);
    grok.shell.addView(v);
  }
}
