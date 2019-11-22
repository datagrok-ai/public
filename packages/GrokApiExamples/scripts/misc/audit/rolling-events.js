let view = gr.newView('Usage');

var last = new Date();
last.setHours(last.getHours() - 5);

var timer = setInterval(function() {
    gr.query('EventsOnDate', {'date': ('after ' + last.toLocaleString().replace(',', ''))})
        .then((t) => {
            gr.addTableView(t);
        });
}, 5000);

