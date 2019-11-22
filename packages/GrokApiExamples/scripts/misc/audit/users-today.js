let view = gr.newView('Usage');

gr.query('UniqueUsersByDate', {'date': 'today'})
    .then(t => {
        let ids = Array.from(t.getCol('id').values());
        gr.dapi.getEntities(ids).then((users) => view.append(ui.list(users)));
    });
