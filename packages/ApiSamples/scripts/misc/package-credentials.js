//tags: Package
//help-url: https://datagrok.ai/help/develop/develop#packages
// Read more about package credentials: https://datagrok.ai/help/develop/how-to/manage-credentials

let p = await grok.dapi.packages.filter('Api Samples').first();
p.getCredentials().then((c) => grok.shell.info(c ? c.openParameters : 'Credentials are not set.'));