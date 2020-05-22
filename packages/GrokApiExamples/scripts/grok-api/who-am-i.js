// prints current user name and picture

var view = grok.shell.newView('User', []);
grok.dapi.users.current().then(function(user) {
    var p = user.picture;
    var d = ui.div([ui.h1(user.firstName + ' ' + user.lastName), p]);
    d.style.width = '300px';
    d.style.textAlign = 'center';
    p.style.width = '300px';
    p.style.height = '300px';
    p.style.fontSize = '150px';
    view.root.appendChild(d);
});