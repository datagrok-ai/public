// Creating yellow information bars in views

grok.shell.newView('Weather App', [
  ui.info(
    'This application fetches weather forecasts and displays them on a map',
    'Weather Forecast Application'
  ),
  ui.info(
    ui.div([ui.p('This application fetches weather forecasts'), ui.p('and displays them on a map')]),
    ui.h1('Weather Forecast Application'),
    false
  ),
  ui.divText('Main Application Area')
]);

/*
  These call styles are also available:
  ui.info('Text');
  ui.info(ui.divText('Text'));
  ui.info(ui.divText('Text'), 'Header');
  ui.info([ui.divText('Text1'), ui.divText('Text2')]);
  ui.info([ui.divText('Text1'), ui.divText('Text2')], 'Header');
*/