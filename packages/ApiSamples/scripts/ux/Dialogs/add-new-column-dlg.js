//Example "Add new column" dialog layout

let table = grok.data.testData('demog', 10000);
let grid = DG.Viewer.grid(table);
grid.root.style.height = '205px';
grid.columns.setVisible(['site']);

let addNewColumnDlg = ui.dialog('Add New Column');
let inputName = ui.stringInput('','');
inputName.input.placeholder = 'Name';
inputName.input.style.width = '100%';
inputName.input.style.maxWidth = '100%';

let inputType = ui.choiceInput('', 'Auto', ['Auto', 'Double', 'Int', 'String', 'Datetime', 'Bool']);
inputType.input.style.width = '100%';

let inputExpression = ui.textInput('','');
inputExpression.input.placeholder = 'Expression';
inputExpression.input.style.minHeight = '140px';
inputExpression.input.style.width = '100%';

let widgetColumnsFunc = null;
let widgetColumns = ui.div();
(async()=>{
widgetColumnsFunc = await DG.Func.byName('ColumnGridWidget').apply({df: table});
widgetColumnsFunc.root.lastChild.style.height = '350px';
widgetColumns.append(widgetColumnsFunc.root)
})();

let widgetFunctionsFunc = null;
let widgetFunctions = ui.div();
(async()=>{
widgetFunctionsFunc = await DG.Func.byName('FunctionsWidget').apply();
widgetFunctionsFunc.props.visibleTags = 'math,text,date,timespan,binning,logic,stats';
widgetFunctionsFunc.props.showSignature = true;
widgetFunctionsFunc.root.lastChild.style.height = '240px';
widgetFunctions.append(widgetFunctionsFunc.root)
})();

let addNewColumnLayout = 
    ui.div([
      ui.block50([
        ui.block([inputName], { style: { width: '65%' } }),
        ui.block([inputType], { style: { width: '35%' } }),
        ui.block([inputExpression]),
        ui.block([grid])
      ], { style: { paddingRight: '20px' } }),

      ui.block25([
        ui.block([widgetColumns])
      ], { style: { paddingRight: '20px' } }),

      ui.block25([
        ui.block([widgetFunctions])
      ])
    ], {style:{width:'700px', height:'400px'}});

addNewColumnDlg.add(addNewColumnLayout).show({x:550, y:100});