import * as DG from 'datagrok-api/dg';
import {InputToggle} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';

export type ToggledInputs = {
  hidden: Set<string>,
  counts: Map<string, number>,
};

const isLinked = (consistency?: ConsistencyInfo) =>
  consistency?.restriction === 'disabled' || consistency?.restriction === 'restricted';

const hasIssue = (validation?: ValidationResult, consistency?: ConsistencyInfo) =>
  !!validation?.errors?.length || !!validation?.warnings?.length || !!consistency?.inconsistent;

const matches = (toggle: InputToggle, prop: DG.Property, consistency?: ConsistencyInfo) => {
  if (toggle.linked && !isLinked(consistency))
    return false;
  if (!toggle.inputs && !toggle.categories)
    return true;
  return !!toggle.inputs?.includes(prop.name) || !!toggle.categories?.includes(prop.category);
};

export function getToggledInputs(
  props: DG.Property[],
  toggles: InputToggle[],
  shown: Set<string>,
  meta: Record<string, any>,
  validations?: Record<string, ValidationResult>,
  consistencies?: Record<string, ConsistencyInfo>,
): ToggledInputs {
  const hidden = new Set<string>();
  const counts = new Map<string, number>();
  for (const prop of props) {
    const {name} = prop;
    if (meta[name]?.hidden || meta[name]?.keepVisible || hasIssue(validations?.[name], consistencies?.[name]))
      continue;
    const toggle = toggles.find((t) => matches(t, prop, consistencies?.[name]));
    if (!toggle)
      continue;
    counts.set(toggle.caption, (counts.get(toggle.caption) ?? 0) + 1);
    if (!shown.has(toggle.caption))
      hidden.add(name);
  }
  return {hidden, counts};
}
