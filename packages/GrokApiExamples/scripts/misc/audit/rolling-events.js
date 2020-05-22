let view = grok.shell.newView('Usage');

var last = new Date();
last.setHours(last.getHours() - 5);

var timer = setInterval(function() {
    grok.data.query('UsageAnalysis:EventsOnDate', {'date': ('after ' + last.toLocaleString().replace(',', ''))})
        .then((t) => {
            grok.shell.addTableView(t);
        });
}, 5000);

