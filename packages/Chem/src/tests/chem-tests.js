async function testSubstructureSearch() {
  let t = grok.data.demo.getDemoTable('molecules');
  await grok.chem.searchSubstructure(t.col('smiles'), 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
}
