// find project, open it and close everything else
let project = await grok.dapi.projects
  .filter('demog')
  .first();
project.open({closeAll: true})

// a shortcut method
grok.dapi.projects.open('demog');