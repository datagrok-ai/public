class ChemScriptsPackage extends DG.Package {

  //name: ActivityCliffs
  //tags: demo, chem, rdkit 
  //description: get activity cliffs
  //sample: chem/activity_cliffs.csv
  //top-menu: Chem | Activity Cliffs
  //input: dataframe df [Input data table]
  //input: column smiles {type:categorical; semType: Molecule} [Molecules, in SMILES format]
  //input: column activities
  //input: double similarity = 80 [Similarity cutoff]
  async activityCliffs(df, smiles, activities, similarity) {
    
    function renderLines (sp, pairs) {

      let np = pairs.rowCount;
      let ctx = sp.getInfo()['canvas'].getContext('2d');
    
      for (let i = 0; i < np; i++) {
        ctx.beginPath();
        ctx.strokeStyle = 'green';
        ctx.lineWidth = 1;

        let n1 = pairs.columns.byName('n1').get(i);
        let n2 = pairs.columns.byName('n2').get(i);

        let pointFrom = sp.worldToScreen(sp.dataFrame.get('~x_coord', n1), sp.dataFrame.get('~y_coord', n1));
        ctx.lineTo(pointFrom.x, pointFrom.y);
        let pointTo = sp.worldToScreen(sp.dataFrame.get('~x_coord', n2), sp.dataFrame.get('~y_coord', n2));
        ctx.lineTo(pointTo.x, pointTo.y);

        ctx.stroke();
      }
    }

    let f = await grok.functions.eval("ChemScripts:ActivityCliffsPy");
    let call = f.prepare({
      'data': df,
      'smiles': smiles,
      'activities': activities,
      'similarityValue': similarity
    });
    await call.call();

    var coords = call.getParamValue('output_coords');
    var pairs = call.getParamValue('output_pairs');

    df.columns.insert(coords.columns.byIndex(0));
    df.columns.insert(coords.columns.byIndex(1));
    df.columns.insert(coords.columns.byIndex(2));
    df.columns.insert(coords.columns.byIndex(3));
    
    df.columns.byName('x_coord').name = '~x_coord';
    df.columns.byName('y_coord').name = '~y_coord';
    df.columns.byName('size').name = '~size';
    df.columns.byName('n').name = '~n';

    let view = grok.shell.addTableView(df);
    let sp = view.addViewer(DG.Viewer.scatterPlot(df, {
      xColumnName: '~x_coord',
      yColumnName: '~y_coord',
      size: '~size'
    }));

    sp.props.showXSelector = false;
    sp.props.showYSelector = false;
    sp.props.markerMinSize = 5;
    sp.props.markerMaxSize = 30;
    sp.props.colorColumnName = 'activity';

    sp.onEvent('d4-before-draw-scene').subscribe(_ => renderLines(sp, pairs));
  }  
}