import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {sequenceToSmiles} from '../model/sequence-to-molfile-utils/from-monomers';
import {SYNTHESIZERS} from '../model/const';

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

  // test('invAb/galnac1', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invAb)sgg(invAb)(GalNAc2)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('invAb/galnac2', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invAb)sgsg(invAb)(GalNAc2)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('invAb/galnac3', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invAb)sggs(invAb)(GalNAc2)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('invAb/galnac4', async () => {
  //   const expected = 'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invAb)sgg(invAb)s(GalNAc2)', false, SYNTHESIZERS.AXOLABS), expected);
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

  // test('(invAb)sucuCfuUf', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O';
  //   expect(sequenceToSmiles('(invAb)sucuCfuUf', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invAb)sAfgcugUf', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C3N=CN=C(N)C=3N=C2)[C@H](F)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C3N=C(N)NC(=O)C=3N=C2)[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1O';
  //   expect(sequenceToSmiles('(invAb)sAfgcugUf', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  test('(invAb)cuCfuUfsc', async () => {
    const expected =
    'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](F)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
    'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](F)[C@@H]1OP(=O)(S)' +
    'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1O';
    expect(sequenceToSmiles('(invAb)cuCfuUfsc', false, SYNTHESIZERS.AXOLABS), expected);
  });

  // test('(invAb)scususu(GalNAc2)', async () => {
  //   const expected =
  //   'O[C@@H]1C[C@@H]O[C@H]1COP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(N)=NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(S)' +
  //   'OC[C@H]1O[C@@H](N2C=CC(=O)NC2(=O))[C@H](OC)[C@@H]1OP(=O)(O)' +
  //   'OCC(O)CNC(=O)CCCC(=O)NC(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)' +
  //   '(COCCC(=O)NCCCNC(=O)CCCCOC2OC(CO)C(O)C(O)C2NC(=O)C)';
  //   expect(sequenceToSmiles('(invAb)scususu(GalNAc2)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  test('(invAb)usAfsucuCfuUfAfgcugUfgCfacususu', async () => {
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
    expect(sequenceToSmiles('(invAb)usAfsucuCfuUfAfgcugUfgCfacususu', false, SYNTHESIZERS.AXOLABS), expected);
  });

  test('(invAb)asacgGfuGfCfAfacucuauuca', async () => {
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
    expect(sequenceToSmiles('(invAb)asacgGfuGfCfAfacucuauuca', false, SYNTHESIZERS.AXOLABS), expected);
  });

  // test('(invAb)scsgguGfcAfAfCfucuauucuga', async () => {
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
  //   expect(sequenceToSmiles('(invAb)scsgguGfcAfAfCfucuauucuga', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invAb)scsaacUfcUfAfUfucuggacuua', async () => {
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
  //   expect(sequenceToSmiles('(invAb)scsaacUfcUfAfUfucuggacuua', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invAb)sasacuCfuAfUfUfcuggacuuua', async () => {
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
  //   expect(sequenceToSmiles('(invAb)sasacuCfuAfUfUfcuggacuuua', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invAb)usAfscug(invAb)(GalNAc2)', async () => {
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
  //   expect(sequenceToSmiles('(invAb)usAfscug(invAb)(GalNAc2)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invAb)sasuaaCfcUf(GalNAc2)', async () => {
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
  //   expect(sequenceToSmiles('(invAb)sasuaaCfcUf(GalNAc2)', false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invAb)sasuaaCfcUfCfUfuguaguuaua(GalNAc2)', async () => {
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
  //   expect(sequenceToSmiles('(invAb)sasuaaCfcUfCfUfuguaguuaua(GalNAc2)', false, SYNTHESIZERS.AXOLABS),
  //     expected);
  // });

  // test('(invAb)scsaucacguUfGfCfagccgucuua(invAb)(GalNAc2)', async () => {
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
  //   expect(sequenceToSmiles('(invAb)scsaucacguUfGfCfagccgucuua(invAb)(GalNAc2)',
  //     false, SYNTHESIZERS.AXOLABS), expected);
  // });

  // test('(invAb)susguuUfgCfCfUfacaucuacua(GalNAc2)', async () => {
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
  //   expect(sequenceToSmiles('(invAb)susguuUfgCfCfUfacaucuacua(GalNAc2)', false, SYNTHESIZERS.AXOLABS),
  //     expected);
  // });
});
