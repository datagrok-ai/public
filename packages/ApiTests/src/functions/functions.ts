import { category, test } from '@datagrok-libraries/utils/src/test';
import { findFuncsByTags, Indexable, DEFAULT_PARAM_VALUES, SPECIAL_CASES, TAGS } from './utils';


/** The objective is to test standard functions generically. These tests check that the functions can be
 * called with proper inputs via the JS API, while specific unit tests for them are in the core codebase. */
category('Functions', () => {
  const functions = findFuncsByTags(TAGS);
  for (const f of functions) {
    let hasUnsupportedType = false;

    let params = f.inputs.reduce((obj, p) => {
      obj[p.name] = p.defaultValue ?? DEFAULT_PARAM_VALUES[p.propertyType];
      if (obj[p.name] == null)
        hasUnsupportedType = true;
      return obj;
    }, <Indexable>{});

    if (hasUnsupportedType)
      continue;

    if (f.name in SPECIAL_CASES)
      params = SPECIAL_CASES[f.name];

    test(f.name, () => f.apply(params));
  }
});
