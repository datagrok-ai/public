let user = DG.User.create();
user.login = 'newlogin';
user.firstName = 'new';
user.lastName = 'user';

await grok.dapi.users.save(user);