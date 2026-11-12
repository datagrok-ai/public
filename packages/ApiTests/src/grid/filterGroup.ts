import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { before, category, expect, settled, test, waitUntil } from '@datagrok-libraries/test/src/test';
import { filter } from 'wu';


category('FilterGroup', () => {
    DG.FilterGroup
    let tv: DG.TableView;
    let fg: DG.FilterGroup;
    before(async () => {
    });

    test('remove', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await settled(fg);

        //@ts-ignore
        let filters: Array<DG.Widget | DG.Filter> = fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM);
        let countOfFilters = filters.length;
        let countOfFiltersByName = fg.filters.length;

        expect(filters.length > 0, true);
        expect(fg.filters.length > 0, true);
        fg.remove(filters[0]);

        //@ts-ignore
        filters = fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM);
        expect(filters.length, countOfFilters - 1);
        expect(fg.filters.length, countOfFiltersByName - 1);
    });

    test('add', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await settled(fg);

        //@ts-ignore
        let filters: Array<DG.Widget | DG.Filter> = fg.getStates('age', DG.FILTER_TYPE.CATEGORICAL);
        let countOfFilters = filters.length;
        let countOfFiltersByName = fg.filters.length;
        fg.add({
            "type": DG.FILTER_TYPE.CATEGORICAL,
            "column": "age",
            "active": false,
            "filterOutMissingValues": false,
            "showMissingValuesOnly": false,
            "min": 10,
            "max": 20,
            "showHistogram": true,
            "showSlider": true,
            "showMinMax": false
        });

        //@ts-ignore
        filters = fg.getStates('age', DG.FILTER_TYPE.CATEGORICAL);
        expect(filters.length, countOfFilters + 1);
        expect(fg.filters.length, countOfFiltersByName + 1);
    });

    test('updateOrAdd', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await settled(fg);

        //@ts-ignore
        let filters: Array<DG.Widget | DG.Filter> = fg.getStates('height', DG.FILTER_TYPE.CATEGORICAL);
        let countOfFilters = filters.length;
        let countOfFiltersByName = fg.filters.length;

        fg.updateOrAdd({
            "type": DG.FILTER_TYPE.CATEGORICAL,
            "column": "height",
            "active": false,
            "filterOutMissingValues": false,
            "showMissingValuesOnly": false,
            "min": 10,
            "max": 20,
            "showHistogram": true,
            "showSlider": true,
            "showMinMax": false
        });

        //@ts-ignore
        filters = fg.getStates('height', DG.FILTER_TYPE.CATEGORICAL);
        expect(filters.length, countOfFilters + 1);
        expect(fg.filters.length, countOfFiltersByName + 1);

        //@ts-ignore
        filters = fg.getStates('height', DG.FILTER_TYPE.CATEGORICAL);
        countOfFilters = filters.length;
        countOfFiltersByName = fg.filters.length;

        fg.updateOrAdd({
            "type": DG.FILTER_TYPE.CATEGORICAL,
            "column": "height",
            "active": false,
            "filterOutMissingValues": false,
            "showMissingValuesOnly": false,
            "min": 10,
            "max": 20,
            "showHistogram": true,
            "showSlider": true,
            "showMinMax": false
        });

        //@ts-ignore
        filters = fg.getStates('height', DG.FILTER_TYPE.CATEGORICAL);
        expect(filters.length, countOfFilters);
        expect(fg.filters.length, fg.filters.length);
    });

    test('setEnabled', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await settled(fg);

        let filters: any[] = fg.getStates('height', DG.FILTER_TYPE.HISTOGRAM);
        fg.setEnabled(filters[0], false);
        let filtersData: any[] = fg.getStates('height', DG.FILTER_TYPE.HISTOGRAM);
        expect(filtersData[0]?.active, false);
    });

    test('setExpanded', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await settled(fg);

        fg.setExpanded(fg.filters[0], false);
        expect(fg.filters[0].root.style.display, 'none');
        fg.setExpanded(fg.filters[0], true);
        expect(fg.filters[0].root.style.display, '');
    });

    test('two filters on one column survive a layout round-trip', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await settled(fg);
        // the default panel already has a histogram on age; a categorical filter on the same column is a second one
        fg.add({type: DG.FILTER_TYPE.CATEGORICAL, column: 'age'});
        expect(fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM).length, 1);
        expect(fg.getStates('age', DG.FILTER_TYPE.CATEGORICAL).length, 1);

        const json = tv.saveLayout().toJson();
        const tv2 = grok.shell.addTableView(grok.data.demo.demog(100));
        tv2.loadLayout(DG.ViewLayout.fromJson(json));
        await waitUntil(() => tv2.viewers.some((v) => v.type === DG.VIEWER.FILTERS), 'No filter panel in the layout');
        const fg2 = tv2.getFiltersGroup();
        await settled(fg2);
        expect(fg2.getStates('age', DG.FILTER_TYPE.HISTOGRAM).length, 1);
        expect(fg2.getStates('age', DG.FILTER_TYPE.CATEGORICAL).length, 1);
    });
}, { clear: true });