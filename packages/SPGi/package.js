class SPGiPackage extends DG.Package {

    //tags: app
    startApp(context) {
        let view = DG.View.create();
        view.name = 'SPGi';
        view.description = 'SPGi projects view';

        let projectsView = DG.ProjectsView.create({});
        projectsView.permanentFilter = '#spgi';
        projectsView.searchValue = '';
        let header = ui.e('div', 'spgi-view-header');
        let content = ui.e('div', 'spgi-view-content');

        let title = ui.divText("SPGi Dashboard");
        title.classList.add('spgi-view-title');
        header.appendChild(title);
        header.style.backgroundImage = `url(${this.webRoot}images/background.png)`;
        content.appendChild(projectsView.root);
        view.root.appendChild(header);
        view.root.appendChild(content);
        view.root.classList.add('d4-flex-col');
        view.root.classList.add('grok-app-view');
        grok.shell.addView(view);

        grok.events.onViewAdded.subscribe((view) => {
            if (view.type === DG.VIEW_TYPE.TABLE_VIEW && view.name === 'Main')
                this.layoutMain(view);
        });

        grok.events.onProjectOpened.subscribe((p) => {
            if (p.name.startsWith('BFC'))
                grok.shell.v = grok.shell.getTableView('Main');
        });
    }


    //input: dataframe table
    prepare(table) {
        // Remove dummy columns
        if (table.columns.contains('RowID'))
            table.columns.remove('RowID');

        // Pack units into tags
        let names = table.columns.names();
        const factors = {
            'm': 1e-3,
            'u': 1e-6,
            'n': 1e-9,
            'p': 1e-12
        }
        for (let n = 0; n < table.columns.length; n++) {
            let name = names[n];
            for (let m = n + 1; m < table.columns.length; m++) {
                if (names[m].startsWith(name) && names[m] === `${name} Unit`) {
                    let col = table.columns.byName(name);
                    let units = table.columns.byName(names[m]);
                    if (units.categories.filter((c) => c !== '').length === 1) {
                        let unit = '';
                        for (let r = 0; r < col.length; r++) {
                            let currUnit = units.get(r);
                            if (currUnit !== '') {
                                // if (currUnit.length > 1) {
                                //     col.set(r, col.get(r) * factors[currUnit.substring(0, 1)]);
                                //     if (unit === '')
                                //         unit = currUnit.substring(1);
                                // } else {
                                     if (unit === '')
                                         unit = currUnit;
                                // }
                            }
                        }
                        col.setTag('units', unit);
                        table.columns.remove(names[m]);
                    }
                }
            }
        }

        // Date, multiline strings
        table.columns.names().forEach(function (name) {
            let col = table.columns.byName(name);

            if (name.endsWith('Date')) {
                if (col.type === DG.TYPE.DATE_TIME)
                    col.setTag('format', 'M/d/yyyy');
            }

            if (col.type === DG.TYPE.STRING)
                for (let r = 0; r < col.length; r++)
                    col.set(r, col.get(r).replace(/\\n/g, ' ').trim());
        });

        table.setTag('spgi', '');
    }


    //description: Prepares all opened tables
    prepareAll() {
        for (let table of grok.shell.tables)
            this.prepare(table);
    }


    //description: Generates layout for "Main" table
    layoutMain(view) {
        let table = view.table;
        const floatNone = 2.6789344063684636e-34;

        // DRCs ID to image URL
        if (table.columns.contains('Competition assay Curve ID')) {
            let col = table.col('Competition assay Curve ID')
            for (let n = 0; n < col.length; n++) {
                let value = col.get(n)
                if (value !== '') {
                    let idx = Math.abs(this.stringHashCode(value)) % 10;
                    col.set(n, `${this.webRoot}images/drcs/drc000${idx}.png`);
                }
            }
            col.semType = 'Image';
        }

        view.grid.columns.setOrder([
            'Unique_Identifier',
            'Structure',
            'Competition assay',
            'Competition assay Date',
            'Competition assay Curve ID',
            'IDEA',
            'Origin',
            'First Synthesis Date',
            'Lab Notebook Id',
            'CSF_Series',
            'CSF_Scaffold',
            'createdAt',
            'status',
            'tags',
            'Cellular assay 1',
            'Cellular assay 1 Date',
            'Cellular assay 1 Curve ID'
        ]);
        table.setTag(DG.TAGS.TOOLTIP, 'Structure\nUnique_Identifier\nauthor\nstatus');

        this.formatsMap = new Map();
        if (grok.shell.tableByName('SupData').d !== null) {
            let table = grok.shell.tableByName('SupData');
            for (let n = 0; n < table.rowCount; n++)
                if (table.get('Property Name', n) === 'STDF DartConditionalFormatting') {
                    let format = JSON.parse(table.get('Property Value', n))['conditionalFormatting'];
                    for (let option of format)
                        option.bgColor = parseInt(`ff${option.bgColor.substr(1)}`, 16);
                    this.formatsMap.set(table.get('Column Name', n), format);
                }
        }

        view.grid.onCellPrepare(function (gc) {
            if (!gc.isTableCell || !this.formatsMap.has(gc.gridColumn.name) || gc.cell.value === floatNone)
                return;
            gc.style.backColor = this.resolveBackColor(gc.gridColumn.name, gc.cell.value);
        }.bind(this));

        let col = table.columns.byName('CSF_Scaffold');
        if (col.semType !== 'Molecule') {
            const scaffoldsMap = new Map([
                ['Scaffold 1', 'C1CCNC1'],
                ['Scaffold 2', 'c1ccc(-c2nnc[nH]2)cc1'],
                ['Scaffold 3', 'c1ccc(-c2nncn2Cc2nccs2)cc1'],
                ['', '']
            ]);
            col.semType = 'Molecule';
            for (let n = 0; n < table.rowCount; n++)
                col.set(n, scaffoldsMap.get(col.get(n)));
            col.compact();
        }
    }


    //description: Resolves background color using format and value
    resolveBackColor(columnName, value) {
        let format = this.formatsMap.get(columnName);
        for (let option of format) {
            if (option.operator !== undefined && option.operand !== undefined) {
                if (option.operator === 'greaterThan')
                    if (value > option.operand)
                        return option.bgColor;
                else
                    console.log(option.operator);
            } else
                return option.bgColor;
        }
    }


    //input: string str
    //output: int hash
    stringHashCode(str) {
        let hash = 0;
        for (let i = 0; i < str.length; i++)
            hash = str.charCodeAt(i) + ((hash << 5) - hash);
        return hash;
    }


    //name: Samples Availability
    //description: Samples availability info panel
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    //condition: t.tags.contains("spgi")
    samplesAvailability(smiles) {
        let table = grok.shell.tableByName('Availability');
        let noInfo = new DG.Widget(ui.divText('Information not available'));
        if (table.d === null)
            return noInfo;
        let selection = table.selection.clone();
        selection.setAll(false);
        for (let n = 0; n < table.rowCount; n++)
            if (table.get('Structure', n) === smiles)
                selection.set(n, true);
        if (selection.trueCount === 0)
            return noInfo;
        let available = table.clone(
            selection, [
                'SMF Availability Amount',
                'SMF Availability Amount Unit',
                'SMF Availability Concentration',
                'SMF Availability Location'], false);
        available.col('SMF Availability Amount').name = 'Amount';
        available.col('SMF Availability Amount Unit').name = 'Unit';
        available.col('SMF Availability Concentration').name = 'Concentration';
        available.col('SMF Availability Location').name = 'Location';
        let unit = available.col('Unit');
        let amount = available.col('Amount');
        for (let n = 0; n < available.rowCount; n++)
            if (unit.get(n) !== '')
                unit.set(n, `${amount.get(n)} ${unit.get(n)}`);
        available.columns.remove('Amount');
        unit.name = 'Amount';
        available.selection.setAll(false);
        let grid = DG.Grid.create(available);
        grid.root.style.width = '400px';
        grid.root.style.height = '125px';
        let button = ui.bigButton('ORDER');
        button.style.marginBottom = '6px';
        return new DG.Widget(ui.div([button, ui.div([grid.root])]));
    }
}
