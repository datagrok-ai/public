let view = grok.newView('Usage');

grok.query('UniqueUsersByDate', {'date': 'today'})
    .then(t => {
        let ids = Array.from(t.getCol('id').values());
        grok.dapi.getEntities(ids).then((users) => view.append(ui.list(users)));
    });
