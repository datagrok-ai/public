//Example of Data Source view

let view = grok.shell.newView('Data Source View', [ui.div(['Select view from the toolbox: ', ui.iconFA('arrow-up')])]);
let viewList = [
'No view', 
'apps', 
'connections', 
'databases',
'emails', //empty
'files',
'forum',
'functions',
'groups',
'help',
'jobs',
'js',
'layouts',
'models',
'notebooks',
'packages',
'projects',
'queries',
'queryruns',
'repositories',
'script', //empty
'scripts',
'settings',
'sketch', //empty
'text',
'users',
'webservices',
'welcome' //show projects
];

let choiceView = ui.choiceInput('View', 'No view', viewList, ()=>{
    view.root.innerHTML = '';
    if (choiceView.value == 'No view')
        view.root.innerHTML = '';
    else
        view.append(DG.View.createByType(choiceView.value))
}); 

view.setRibbonPanels([
    [
        choiceView.root
    ],
]);