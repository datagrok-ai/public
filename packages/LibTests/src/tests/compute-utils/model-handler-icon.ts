import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {ModelHandler} from '@datagrok-libraries/compute-utils/model-catalog';

category('ComputeUtils: ModelHandler', async () => {
  test('renderIcon for a function without a package', async () => {
    const func = DG.Func.find({name: 'Abs'}).find((f) => f.package == null)!;
    const icon = new ModelHandler().renderIcon(func);
    expect(icon instanceof HTMLElement, true);
  });
});
