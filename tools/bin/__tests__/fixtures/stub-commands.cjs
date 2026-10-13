// Preloaded with `node -r` in front of bin/grok.js: every command module is replaced by a stub
// that prints the parsed options it was called with, so a test sees exactly what grok.js passes on.
const Module = require('module');

const load = Module._load;
Module._load = function(request, parent, isMain) {
  if (request.startsWith('./commands/') && request !== './commands/help') {
    return new Proxy({}, {get: () => (args) => {
      const {_, ...options} = args;
      process.stdout.write(JSON.stringify(options));
      return true;
    }});
  }
  return load.apply(this, arguments);
};
