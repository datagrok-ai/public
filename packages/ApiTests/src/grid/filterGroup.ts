import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { before, category, expect, test } from '@datagrok-libraries/utils/src/test';
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
        await DG.delay(50);

        //@ts-ignore
        let filters: Array<DG.Widget | DG.Filter> = fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM);
        let countOfFilters = filters.length;
        let countOfFiltersByName = fg.filters.length;

        await DG.delay(50);
        expect(filters.length > 0, true);
        expect(fg.filters.length > 0, true);
        fg.remove(filters[0]);

        //@ts-ignore
        filters = fg.getStates('age', DG.FILTER_TYPE.HISTOGRAM);
        await DG.delay(50);
        expect(filters.length, countOfFilters - 1);
        expect(fg.filters.length, countOfFiltersByName - 1);
    });

    test('add', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await DG.delay(50);

        //@ts-ignore
        let filters: Array<DG.Widget | DG.Filter> = fg.getStates('age', DG.FILTER_TYPE.CATEGORICAL);
        let countOfFilters = filters.length;
        let countOfFiltersByName = fg.filters.length;
        await DG.delay(50);
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
        await DG.delay(50);
        expect(filters.length, countOfFilters + 1);
        expect(fg.filters.length, countOfFiltersByName + 1);
    });

    test('updateOrAdd', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await DG.delay(50);

        //@ts-ignore
        let filters: Array<DG.Widget | DG.Filter> = fg.getStates('height', DG.FILTER_TYPE.CATEGORICAL);
        let countOfFilters = filters.length;
        let countOfFiltersByName = fg.filters.length;
        await DG.delay(50);

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
        await DG.delay(50);

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
        await DG.delay(50);

        let filters: any[] = fg.getStates('height', DG.FILTER_TYPE.HISTOGRAM);
        fg.setEnabled(filters[0], false);
        await DG.delay(50);
        let filtersData: any[] = fg.getStates('height', DG.FILTER_TYPE.HISTOGRAM);
        await DG.delay(50);
        expect(filtersData[0]?.active, false);
    });

    test('setExpanded', async () => {
        tv = grok.shell.addTableView(grok.data.demo.demog(100));
        fg = tv.getFiltersGroup();
        await DG.delay(50);

        fg.setExpanded(fg.filters[0], false);
        expect(fg.filters[0].root.style.display, 'none');
        await DG.delay(50);
        fg.setExpanded(fg.filters[0], true);
        expect(fg.filters[0].root.style.display, '');
    });
}, { clear: true });