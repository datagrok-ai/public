//name: AddMorganFP
//description: Adds a column with Morgan fingerprints
//language: javascript
//input: dataframe inDf
//output: dataframe outDf

let module = chem.getRdKitModule();
let molecules = inDf.columns[0];
const length = molecules.length;
let bitset = DG.BitSet.create(100);
bitset.setAll(true);
let morganFpCol = DG.Column.fromType(DG.TYPE.OBJECT, 'cache_morganFP', length);
for (let i = 0; i < length; ++i) {
  morganFpCol.set(i, bitset.toBinaryString());
}
outDf = DG.DataFrame.create(length);
outDf.columns.add(molecules);
outDf.columns.add(morganFpCol);