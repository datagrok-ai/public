//Read more about package credentials: https://datagrok.ai/help/develop/js-api#credentials

let p = new DG.Package();
p.name = "ApiSamples";

p.getCredentials().then((c) => grok.shell.info(c.parameters['test']));