import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {expect} from "@datagrok-libraries/utils/src/test";
import {_package} from "../package-test";

export async function requireText(name: string): Promise<string> {
  return await _package.files.readAsText(name);
}

export async function _testSearchSubstructure(csv: string, pattern: string, trueIndices: number[], params: any | null = null) {
  // await grok.functions.call('Chem:initChem');
  const df = DG.DataFrame.fromCsv(csv);
  const col = df.columns[0];
  const bitset: DG.BitSet = (params !== null) ?
    (await grok.chem.searchSubstructure(col, pattern, params)) :
    (await grok.chem.searchSubstructure(col, pattern));
  const bitsetString = bitset.toBinaryString();
  const bitsetArray = [...bitsetString];
  for (let k = 0; k < trueIndices.length; ++k) {
    expect(bitsetArray[trueIndices[k]] === '1', true);
    bitsetArray[trueIndices[k]] = '0';
  }
  for (let i = 0; i < bitsetArray.length; ++i)
    expect(bitsetArray[i] === '0', true);
}

export async function _testSearchSubstructureSARSmall(params: any | null = null) {
  const file = await requireText('sar_small.csv');

  // await grok.functions.call('Chem:initChem');
  const df = DG.DataFrame.fromCsv(file);
  const col = df.columns[0];
  const bitset: DG.BitSet = (params !== null) ?
    (await grok.chem.searchSubstructure(col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', params)) :
    (await grok.chem.searchSubstructure(col, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1'));
  const countDataframe = col.length;
  const countResult = bitset.trueCount;
  expect(countDataframe, 200);
  expect(countResult, 90);
}

export async function _testSearchSubstructureAllParameters(foo: any) {
  await foo({substructLibrary: false});
  await foo({substructLibrary: true});
  await foo({});
  await foo();
}
