let view = grok.newView('Usage');

var last = new Date();
last.setHours(last.getHours() - 5);

var timer = setInterval(function() {
    grok.query('EventsOnDate', {'date': ('after ' + last.toLocaleString().replace(',', ''))})
        .then((t) => {
            grok.addTableView(t);
        });
}, 5000);

