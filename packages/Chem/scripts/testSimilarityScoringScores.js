//name: testSimilarityScoringScores
//language: javascript
        
grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv').then(molecules =>
    grok.chem.similarityScoring(molecules.col('smiles'), 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1').then(similar =>
        grok.shell.addTableView(similar)
    )
);