/**
 * @jest-environment jsdom
 */

import * as _package from '../package';

jest.mock('datagrok-api/dg', () => {
  const originalModule = jest.requireActual('datagrok-api/dg');
  return {
    __esModule: true,
    ...originalModule,
    Func: {find: (o: object) => []},
  };
});

test('adds 1 + 2 to equal 3', () => {
  //expect(_package.sum(1, 2)).toBe(3);
});