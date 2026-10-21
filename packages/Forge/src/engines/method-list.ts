import * as DG from 'datagrok-api/dg';
import {boolColumn, textColumn} from '../columns';
import {defaultHyperparameters, Engine, hasTrainingRows, hyperparametersOf, isComplete, isServerEngine, rolesOf}
  from './engine';
import {EngineRegistry} from './engine-registry';

/** One setting (a `train` input) of a method, as `Forge:listMethods` describes it; absent keys are not declared. */
export interface SettingDetails {
  name: string;
  type: string;
  default?: number | string | boolean;
  min?: number;
  max?: number;
  choices?: string[];
  description?: string;
}

/** The app's five method columns, one row per method: **Method**, **Package**, **Method type**, **Roles**,
 * **Hyperparameters**. */
export function methodColumns(engines: Engine[]): DG.Column[] {
  return [
    textColumn('Method', engines.map((e) => e.name)),
    textColumn('Package', engines.map((e) => e.namespace)),
    textColumn('Method type', engines.map((e) => e.kind)),
    textColumn('Roles', engines.map((e) => rolesOf(e).join(', '))),
    textColumn('Hyperparameters', engines.map((e) => hyperparametersOf(e).map((p) => p.name).join(', '))),
  ];
}

/** One row per method with the {@link methodColumns}. */
export function methodsFrame(engines: Engine[]): DG.DataFrame {
  return DG.DataFrame.fromColumns(methodColumns(engines));
}

function methodSettings(engine: Engine): SettingDetails[] {
  const defaults = defaultHyperparameters(engine);
  return hyperparametersOf(engine).map((p) => {
    const setting: SettingDetails = {name: p.name, type: p.propertyType};
    const value = defaults[p.name];
    if (value !== undefined)
      setting.default = value;
    if (typeof p.min === 'number')
      setting.min = p.min;
    if (typeof p.max === 'number')
      setting.max = p.max;
    if (Array.isArray(p.choices) && p.choices.length > 0)
      setting.choices = p.choices.map(String);
    if (typeof p.description === 'string' && p.description !== '')
      setting.description = p.description;
    return setting;
  });
}

/** `Forge:listMethods` without a table: [engines] (by default the complete methods in discovery order) with the app's
 * five columns, **Runs on**, **Keeps training rows**, **Retrains live** and **Settings** (json text). */
export function methodList(engines: Engine[] = EngineRegistry.discover().filter(isComplete)): DG.DataFrame {
  const [method, pkg, kind, roles, hyperparameters] = methodColumns(engines);
  const frame = DG.DataFrame.fromColumns([method, pkg, kind,
    textColumn('Runs on', engines.map((e) => isServerEngine(e) ? 'server' : 'browser')),
    boolColumn('Keeps training rows', engines.map(hasTrainingRows)),
    boolColumn('Retrains live', engines.map((e) => e.isLiveUpdate)),
    roles, hyperparameters,
    textColumn('Settings', engines.map((e) => JSON.stringify(methodSettings(e)))),
  ]);
  frame.name = 'Methods';
  return frame;
}
