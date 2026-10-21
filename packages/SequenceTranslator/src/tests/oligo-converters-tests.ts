/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {OLIGO_SEM_TYPE} from '../oligo-renderer/types';

import {_package} from '../package-test';

category('OligoRenderer: converters', () => {
  test('Combine Sense + Antisense on sirna-demo.csv', async () => {
    const df = DG.DataFrame.fromCsv(await _package.files.readAsText('samples/sirna-demo.csv'));
    await grok.data.detectSemanticTypes(df);
    const col: DG.Column<string> = await grok.functions.call(`${_package.name}:combineSenseAntisenseToOligoNucleotide`,
      {table: df, senseCol: df.getCol('sense_helm'), antiCol: df.getCol('antisense_helm')});
    expect(df.col('sense_helm+antisense_helm (oligo)') !== null, true, 'the combined column in the table');
    expect(col.semType, OLIGO_SEM_TYPE);
    expect(col.get(0), 'RNA1{m(G)[sp].m(A)[sp].m(C)p.m(U)p.m(G)p.m(A)p.m(A)p.m(U)p.m(A)p.m(U)p.m(A)p.m(A)p.m(A)p.m(C)p.m(U)p.m(U)p.m(G)[sp].m(U)[sp].m(G).[L3]}|RNA2{m(C)[sp].m(A)[sp].m(C)p.m(A)p.m(A)p.m(G)p.m(U)p.m(U)p.m(U)p.m(A)p.m(U)p.m(A)p.m(U)p.m(U)p.m(C)p.m(A)p.m(G)[sp].m(U)[sp].m(C)}$$$$');
    // aso-0034 is a single-strand ASO: no antisense, no second chain
    expect(col.get(33), 'RNA1{[lna](C)[sp].[lna](A)[sp].[lna](G)[sp].d(T)[sp].d(G)[sp].d(T)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(T)[sp].d(G)[sp].d(C)[sp].d(T)[sp].d(C)[sp].d(T)[sp].d(A)[sp].[lna](T)[sp].[lna](A)[sp].[lna](A)}$$$$');
  });
});
