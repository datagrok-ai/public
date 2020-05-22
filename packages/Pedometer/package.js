class PedometerPackage extends DG.Package {

    //input: dataframe accel [Accelerometry data table]
    //input: column x {semType: Accelerometer-X} [X axis]
    //input: column y {semType: Accelerometer-Y} [Y axis]
    //input: column z {semType: Accelerometer-Z} [Z axis]
    //input: column timeOffset {semType: Time-Offset} [Time offset column]
    pedometer(accel, x, y, z, timeOffset) {
        let view = grok.shell.getTableView(accel.name);

        let viewer = view.markup();
        while (viewer.root.firstChild)
            viewer.root.removeChild(viewer.root.firstChild);
        let sampleRate = this._estimateSampleRate(timeOffset);
        viewer.root.appendChild(ui.divText('Sample rate: ' + sampleRate.toString() + ' Hz'));

        this._detectSteps(accel, x, y, z, sampleRate, true, result => {
            if (accel.columns.contains('steps'))
                accel.columns.remove('steps');

            let steps = accel.columns.add((new DG.DataFrame(result)).col('steps'), true);
            view.lineChart({
                xColumnName: timeOffset.name,
                yColumnNames: [x.name, y.name, z.name, steps.name]
            });

            let filtered = viewer.root.appendChild(ui.divText('Filtered steps: 0'));
            let selected = viewer.root.appendChild(ui.divText('Selected steps: 0'));

            function countSteps(index) {
                let count = 0;
                for (let n = -1; (n = index.findNext(n, true)) !== -1;)
                    if (steps.get(n) === 1)
                        count++;
                return count;
            }

            this.subs.push(accel.onFilterChanged.subscribe((_) => filtered.innerText = 'Filtered steps: ' + countSteps(accel.filter)));
            this.subs.push(accel.onSelectionChanged.subscribe((_) => selected.innerText = 'Selected steps: ' + countSteps(accel.selection)));
        });
    }

    //name: Step counter
    //description: Step counter info panel
    //tags: panel, widgets
    //input: dataframe accel {semType: Accelerometer} [Accelerometry data table]
    //output: widget result
    //condition: stepCounterCondition(t)
    stepCounter(accel) {
        let columns = accel.columns.toList();
        let x = columns.find((c) => c.semType === 'Accelerometer-X');
        let y = columns.find((c) => c.semType === 'Accelerometer-Y');
        let z = columns.find((c) => c.semType === 'Accelerometer-Z');
        let timeOffset = columns.find((c) => c.semType === 'Time-Offset');
        let sampleRate = this._estimateSampleRate(timeOffset);
        let root = ui.div();

        this._detectSteps(accel, x, y, z, sampleRate, false, result => {
            let steps = new DG.DataFrame(result).col('steps').stats.sum;
            while (root.firstChild)
                root.removeChild(root.firstChild);
            root.appendChild(ui.divText(steps.toString() + ' steps'));
        });

        root.appendChild(ui.loader());

        return new DG.Widget(root);
    }

    _estimateSampleRate(timeOffset) {
        return 1.0 / (timeOffset.get(1) - timeOffset.get(0));
    }

    _detectSteps(accel, x, y, z, sampleRate, showProgress, callback) {
        grok.functions.call('Pedometer:DetectSteps', {
            "accel": accel.d,
            "x": x.name,
            "y": y.name,
            "z": z.name,
            "sample_rate": sampleRate
        }, showProgress).then(callback);
    }
}
