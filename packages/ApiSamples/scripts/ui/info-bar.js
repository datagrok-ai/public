// Creating yellow information bars in views

let appDescription = ui.info(
  'This app fetches weather forecasts and displays them on a map',
  'Weather Forecast');
grok.shell.newView('New App', [appDescription]);

/*
  These call styles are also available:
  ui.info('Text');
  ui.info(ui.divText('Text'));
  ui.info(ui.divText('Text'), 'Header');
  ui.info([ui.divText('Text1'), ui.divText('Text2')]);
  ui.info([ui.divText('Text1'), ui.divText('Text2')], 'Header');
*/