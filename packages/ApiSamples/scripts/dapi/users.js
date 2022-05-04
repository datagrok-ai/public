let user = DG.User.create();
user.login = 'newlogin';
user.status = USER_STATUS.STATUS_NEW;
user.firstName = 'new';
user.lastName = 'user';

await grok.dapi.users.save(user);