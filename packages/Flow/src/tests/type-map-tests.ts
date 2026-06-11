import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {areTypesCompatible, dgTypeToSlotType, getSlotColor, DG_TYPE_MAP} from '../types/type-map';

category('Flow: type-map', () => {
  test('identical types are compatible', async () => {
    expect(areTypesCompatible('dataframe', 'dataframe'), true);
    expect(areTypesCompatible('string', 'string'), true);
    expect(areTypesCompatible('column', 'column'), true);
  });

  test('dynamic and object are wildcards', async () => {
    expect(areTypesCompatible('dynamic', 'dataframe'), true);
    expect(areTypesCompatible('dataframe', 'dynamic'), true);
    expect(areTypesCompatible('object', 'int'), true);
    expect(areTypesCompatible('int', 'object'), true);
  });

  test('numeric widening is symmetric', async () => {
    expect(areTypesCompatible('int', 'double'), true);
    expect(areTypesCompatible('double', 'int'), true);
    expect(areTypesCompatible('num', 'int'), true);
    expect(areTypesCompatible('int', 'num'), true);
  });

  test('list and string_list interconvert', async () => {
    expect(areTypesCompatible('list', 'string_list'), true);
    expect(areTypesCompatible('string_list', 'list'), true);
  });

  test('incompatible types are rejected', async () => {
    expect(areTypesCompatible('dataframe', 'string'), false);
    expect(areTypesCompatible('column', 'dataframe'), false);
    expect(areTypesCompatible('bool', 'datetime'), false);
    expect(areTypesCompatible('string', 'int'), false);
  });

  test('dgTypeToSlotType maps known and unknown types', async () => {
    expect(dgTypeToSlotType('dataframe'), 'dataframe');
    expect(dgTypeToSlotType('blob'), 'byte_array'); // blob slot type is byte_array
    expect(dgTypeToSlotType('totally_unknown'), 'totally_unknown'); // passthrough
  });

  test('every mapped type has a non-empty color', async () => {
    for (const [dgType, def] of Object.entries(DG_TYPE_MAP)) {
      expect(typeof def.color === 'string' && def.color.length > 0, true, `color for ${dgType}`);
      expect(getSlotColor(dgType).length > 0, true, `getSlotColor for ${dgType}`);
    }
  });
});
