import * as DG from 'datagrok-api/dg';
import {category, expect, expectFloat, test} from '@datagrok-libraries/test/src/test';

category('AI: Qnum: static helpers', () => {
  test('exact / less / greater set qualifier', async () => {
    for (var v of [0, 1.5, -3.25, 1e-6, 1234567.89]) {
      const ex = DG.Qnum.exact(v);
      const lt = DG.Qnum.less(v);
      const gt = DG.Qnum.greater(v);
      expect(DG.Qnum.getQ(ex), DG.QNUM_EXACT);
      expect(DG.Qnum.getQ(lt), DG.QNUM_LESS);
      expect(DG.Qnum.getQ(gt), DG.QNUM_GREATER);
      expectFloat(DG.Qnum.getValue(ex), v, 0.01);
      expectFloat(DG.Qnum.getValue(lt), v, 0.01);
      expectFloat(DG.Qnum.getValue(gt), v, 0.01);
    }
  });

  test('exact / less / greater equal create with constant', async () => {
    for (var v of [0, 1.5, -3.25, 100]) {
      expect(DG.Qnum.exact(v), DG.Qnum.create(v, DG.QNUM_EXACT));
      expect(DG.Qnum.less(v), DG.Qnum.create(v, DG.QNUM_LESS));
      expect(DG.Qnum.greater(v), DG.Qnum.create(v, DG.QNUM_GREATER));
    }
  });

  test('parse compact forms', async () => {
    const e = DG.Qnum.parse('10');
    const l = DG.Qnum.parse('<10');
    const g = DG.Qnum.parse('>10');
    expect(DG.Qnum.getQ(e), DG.QNUM_EXACT);
    expect(DG.Qnum.getQ(l), DG.QNUM_LESS);
    expect(DG.Qnum.getQ(g), DG.QNUM_GREATER);
    expectFloat(DG.Qnum.getValue(e), 10, 0.01);
    expectFloat(DG.Qnum.getValue(l), 10, 0.01);
    expectFloat(DG.Qnum.getValue(g), 10, 0.01);
  });

  test('parse tolerates whitespace', async () => {
    const a = DG.Qnum.parse(' < 10');
    const b = DG.Qnum.parse('>  10');
    expect(DG.Qnum.getQ(a), DG.QNUM_LESS);
    expect(DG.Qnum.getQ(b), DG.QNUM_GREATER);
    expectFloat(DG.Qnum.getValue(a), 10, 0.01);
    expectFloat(DG.Qnum.getValue(b), 10, 0.01);
  });

  test('parse rejects invalid input', async () => {
    expect(DG.Qnum.parse('abc') == null, true);
    expect(DG.Qnum.parse('') == null, true);
  });

  test('qualifier returns =, <, >', async () => {
    expect(DG.Qnum.qualifier(DG.Qnum.exact(1.5)), '=');
    expect(DG.Qnum.qualifier(DG.Qnum.less(1.5)), '<');
    expect(DG.Qnum.qualifier(DG.Qnum.greater(1.5)), '>');
  });

  test('toString shape', async () => {
    const sExact = DG.Qnum.toString(DG.Qnum.exact(1.5));
    const sLess = DG.Qnum.toString(DG.Qnum.less(1.5));
    const sGreater = DG.Qnum.toString(DG.Qnum.greater(1.5));
    expect(sExact.startsWith('<') || sExact.startsWith('>'), false);
    expect(sLess.startsWith('<'), true);
    expect(sGreater.startsWith('>'), true);
    expectFloat(parseFloat(sExact), 1.5, 0.01);
    expectFloat(parseFloat(sLess.replace(/^[<>=]/, '')), 1.5, 0.01);
    expectFloat(parseFloat(sGreater.replace(/^[<>=]/, '')), 1.5, 0.01);
  });
});
