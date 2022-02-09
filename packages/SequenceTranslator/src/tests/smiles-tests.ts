import {category, expect, testExpectFinish} from '@datagrok-libraries/utils/src/test';
import {sequenceToSmiles} from '../package';

category('sequence-translator', () => {
  testExpectFinish('basic sequence', async () => {
    const expected = 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O';
    expect(sequenceToSmiles('AGGTCCTCTTGACTTAGGCC'), expected);
  });

  testExpectFinish('invabasic/galnac1', async () => {
    const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
    'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
    expect(sequenceToSmiles('(invabasic)sgg(invabasic)(GalNAc-2-JNJ)'), expected);
  });

  testExpectFinish('invabasic/galnac2', async () => {
    const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
    'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
    expect(sequenceToSmiles('(invabasic)sgsg(invabasic)(GalNAc-2-JNJ)'), expected);
  });

  testExpectFinish('invabasic/galnac3', async () => {
    const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
    'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
    'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
    expect(sequenceToSmiles('(invabasic)sggs(invabasic)(GalNAc-2-JNJ)'), expected);
  });

  testExpectFinish('invabasic/galnac4', async () => {
    const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
    'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
    '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
    expect(sequenceToSmiles('(invabasic)sgg(invabasic)s(GalNAc-2-JNJ)'), expected);
  });
});

