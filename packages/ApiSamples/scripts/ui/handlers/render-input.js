// ObjectHandler.renderInput: a meta-provided input for an object,
// such as the user selector for a User (null when the handler does not define one)

let user = await grok.dapi.users.current();
let input = DG.ObjectHandler.forEntity(user).renderInput(user);
grok.shell.newView('renderInput', [input.root]);
