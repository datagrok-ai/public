//Example of various input-forms

let view = grok.shell.newView('Form types');

let basicLayout = ui.panel([
  ui.h1('Basic form'),
  ui.inputs([
    ui.stringInput('Name',''),
    ui.textInput('Description',''),
    ui.choiceInput('Type', 'Text', ['Text','Number','Bool']),
    ui.choiceInput('Suggested values', 'Query', ['Query']),
    ui.choiceInput('Query', 'List of countries', ['List of countries']),
    ui.intInput('Minimun',0),
    ui.intInput('Maximum',100),
    ui.boolInput('Add slicer', true),
    ui.buttonsInput([
      ui.bigButton('Apply'),
      ui.button('Reset')
    ])
  ]),
]);

let inlineLayout = ui.panel([
  ui.h1('Inline form'),
  ui.inputs([
    ui.divH([
      ui.choiceInput('Tables', 'Demog',['Demog']), 
      ui.choiceInput('', 'Crime',['Crime']),
      ui.button(ui.iconFA('folder-open'))
    ]),
    ui.divH([
      ui.stringInput('Key Columns', 'Age'), 
      ui.stringInput('', 'Age'), 
      ui.button(ui.iconFA('plus'))
    ]),
    ui.divH([ui.choiceInput('Link Type', 'row to filter', ['row to filter', 'row to row', 'filter to filter'])]),
    ui.divH([
      ui.buttonsInput([
        ui.bigButton('Apply'),
        ui.button('Reset')
      ])
    ])
  ])
]);
//Setting the same input width for all fields
$(inlineLayout)
  .find('input, select')
  .css('min-width','130px');


let columnLayout = ui.panel([
  ui.h1('Two column form'),
  ui.divH([
  ui.div([
    ui.inputs([
      ui.stringInput('Name',''),
      ui.textInput('Description',''),
      ui.choiceInput('Type', 'Text', ['Text','Number','Bool']),
      ui.choiceInput('Suggested values', 'Query', ['Query']),
      ui.choiceInput('Query', 'List of countries', ['List of countries']),
      ui.boolInput('Show assistance', false),
      ui.buttonsInput([
        ui.bigButton('Apply'),
        ui.button('Reset')
      ])
    ])
  ]),
  ui.div([
    ui.inputs([
      ui.choiceInput('Calculation type', 'Percent of total',['Percent of total']),
      ui.boolInput('Compute total', true),
      ui.choiceInput('Compute usign', 'Pane (down)',['Pane (down)']),
      ui.choiceInput('Compute type', 'Month of order',['Month of order']),
      ui.intInput('Minimun',0),
      ui.intInput('Maximum',100),
    ]),
  ]),
])
]);

let advancedFormAcc = ui.accordion();
let slider = ui.element('input');
slider.type = 'range';
slider.min = '0';
slider.max = '100';
let sliderLabel = ui.label('Value range');
sliderLabel.classList = 'ui-label ui-input-label';

let min = ui.intInput('Minimun',0);
let max = ui.intInput('Maximum',100);
min.type = 'number';
max.type = 'number';

advancedFormAcc.addPane('Advanced',()=>
  ui.inputs([
    ui.choiceInput('Suggested values', 'Query', ['Query']),
    ui.choiceInput('Query', 'List of countries', ['List of countries']),

    ui.boolInput('Add slicer', true),
    ui.div([
      sliderLabel,
      slider
    ], 'ui-input-root')
  ])
, true);

let advancedForm = ui.panel([
  ui.h1('Advanced form'),
  ui.inputs([
    ui.stringInput('Name',''),
    ui.textInput('Description',''),
    ui.choiceInput('Type', 'Text', ['Text','Number','Bool']),
  ]),
  advancedFormAcc.root,
  ui.inputs([
    ui.buttonsInput([
      ui.bigButton('Apply'),
      ui.button('Reset')
    ])
  ])
]);


view.append(ui.div([
  basicLayout,
  ui.div([],{style:{borderBottom:'1px solid var(--grey-2)'}}),
  inlineLayout,
  ui.div([],{style:{borderBottom:'1px solid var(--grey-2)'}}),
  columnLayout,
  ui.div([],{style:{borderBottom:'1px solid var(--grey-2)'}}),
  advancedForm
])
);