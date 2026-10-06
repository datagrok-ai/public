// Package autostart functions run 3 seconds after the app has started;
// wait for them before reading state they set up (menus, panels, settings)

await grok.shell.autostartsCompleted;
grok.shell.info('Every package autostart has run');
