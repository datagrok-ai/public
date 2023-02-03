import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {sequenceToSmiles} from '../utils/structures-works/from-monomers';
import {SYNTHESIZERS} from '../hardcode-to-be-eliminated/map';

category('sequence-translator', () => {
  // test('AGGTCCTCTTGACTTAGGCC', async () => {
  //   const expected = 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=C(C)C(=O)NC2(=O))C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)C[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1OP(=O)(O)OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))C[C@@H]1O';

  //   expect(sequenceToSmiles('AGGTCCTCTTGACTTAGGCC', false, SYNTHESIZERS.RAW_NUCLEOTIDES), expected);
  // });

  // test('invabasic/galnac1', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)sgg(invabasic)(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('invabasic/galnac2', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)sgsg(invabasic)(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('invabasic/galnac3', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)sggs(invabasic)(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('invabasic/galnac4', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)sgg(invabasic)s(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  test('usCfCfUfGfAf', async () => {
    const expected = 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1O';
    expect(sequenceToSmiles('usCfCfUfGfAf', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('usAfsusgsgsg', async () => {
    const expected = 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('usAfsusgsgsg', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('UfUfUfsCfsuacg', async () => {
    const expected = 'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('UfUfUfsCfsuacg', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('susususauasu', async () => {
    const expected = 'OP(=O)(S)OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('susususauasu', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('CfGfCfsGfsCf', async () => {
    const expected = 'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1O';
    expect(sequenceToSmiles('CfGfCfsGfsCf', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('acacacsacsac', async () => {
    const expected = 'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('acacacsacsac', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('cccgggusug', async () => {
    const expected =
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('cccgggusug', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('UfAfCfGfGfCfAfUf', async () => {
    const expected =
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O';
    expect(sequenceToSmiles('UfAfCfGfGfCfAfUf', false, SYNTHESIZERS.AXOLABS), expected);
  });

  // test('(invabasic)sucuCfuUf', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O';
  //   expect(sequenceToSmiles('(invabasic)sucuCfuUf', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invabasic)sAfgcugUf', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O';
  //   expect(sequenceToSmiles('(invabasic)sAfgcugUf', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  test('(invabasic)cuCfuUfsc', async () => {
    const expected =
    'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('(invabasic)cuCfuUfsc', false, SYNTHESIZERS.AXOLABS), expected);
  });

  // test('(invabasic)scususu(GalNAc-2-JNJ)', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)scususu(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  test('(invabasic)usAfsucuCfuUfAfgcugUfgCfacususu', async () => {
    const expected =
    'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('(invabasic)usAfsucuCfuUfAfgcugUfgCfacususu', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('(invabasic)asacgGfuGfCfAfacucuauuca', async () => {
    const expected =
    'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('(invabasic)asacgGfuGfCfAfacucuauuca', false, SYNTHESIZERS.AXOLABS), expected);
  });

  // test('(invabasic)scsgguGfcAfAfCfucuauucuga', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O';
  //   expect(sequenceToSmiles('(invabasic)scsgguGfcAfAfCfucuauucuga', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invabasic)scsaacUfcUfAfUfucuggacuua', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O';
  //   expect(sequenceToSmiles('(invabasic)scsaacUfcUfAfUfucuggacuua', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invabasic)sasacuCfuAfUfUfcuggacuuua', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1O';
  //   expect(sequenceToSmiles('(invabasic)sasacuCfuAfUfUfcuggacuuua', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invabasic)usAfscug(invabasic)(GalNAc-2-JNJ)', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)usAfscug(invabasic)(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invabasic)sasuaaCfcUf(GalNAc-2-JNJ)', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)sasuaaCfcUf(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invabasic)sasuaaCfcUfCfUfuguaguuaua(GalNAc-2-JNJ)', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)sasuaaCfcUfCfUfuguaguuaua(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS),
  //     expected);
  // });

  // test('(invabasic)scsaucacguUfGfCfagccgucuua(invabasic)(GalNAc-2-JNJ)', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)scsaucacguUfGfCfagccgucuua(invabasic)(GalNAc-2-JNJ)',
  //     false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invabasic)susguuUfgCfCfUfacaucuacua(GalNAc-2-JNJ)', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invabasic)susguuUfgCfCfUfacaucuacua(GalNAc-2-JNJ)', false, SYNTHESIZERS.AXOLABS),
  //     expected);
  // });
});
