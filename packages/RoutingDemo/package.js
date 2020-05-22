class RoutingDemoPackage extends DG.Package {

    //name: Routing Demo
    //tags: app
    startApp(context) {
        let parser = document.createElement('a');
        parser.href = window.location;
        let pathSegments = parser.pathname.split('/');

        let t = grok.data.testData("cars");
        let view = grok.shell.addTableView(t);
        view.name = 'SDTM:LB:Preview';
        view.description = 'SDTM LB domain viewer';
        view.basePath = '/cars';
        let label = ui.h1();

        function setFilter(str) {
            label.innerText = str;
            view.path = '/' + str;
            for (let i = 0; i < t.rowCount; i++)
                t.filter.set(i, str === 'All' || t.get('Make', i) === str);
        }

        if (pathSegments.length > 4)
            setFilter(pathSegments[4]);
        else
            setFilter('All');

        view.toolbox = ui.div([
            label,
            ui.button('Honda', () => {
                setFilter('Honda');
            }),
            ui.button('Tesla', () => {
                setFilter('Tesla');
            })
        ], 'sdtm-controls,pure-form');
    }
}
