// Example of navigation between views and by URLs

//Creating new views
let mainView = grok.shell.newView('Main view', [grok.shell.info('Main view')]);
let firstView = grok.shell.newView('View 1', [ui.bigButton('Back to main view', ()=> {grok.shell.v = mainView; grok.shell.info('Main view');})]);
let secondView = grok.shell.newView('View 2', [ui.bigButton('Back to main view', ()=> {grok.shell.v = mainView; grok.shell.info('Main view');})]);

//Set the URL path for views
firstView.basePath = '/firstView';
secondView.basePath = '/secondView';

mainView.append(ui.div([
  //Navigate to another view
  ui.button('First view', ()=>{grok.shell.v = firstView; grok.shell.info('View 1');}),
  //Navigate to another view by URL
  ui.button('Second view', ()=>{grok.shell.route('/secondView'); grok.shell.info('View 2');}),
]));

// Set the current view
grok.shell.v = mainView; 