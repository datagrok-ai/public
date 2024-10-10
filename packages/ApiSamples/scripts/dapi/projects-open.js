// find project, open it and close everything else
grok.dapi.projects
  .filter('demog')
  .first()
  .then(p => p.open({closeAll: true}));

// a shortcut method
grok.dapi.projects.open('SAR');