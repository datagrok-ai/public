class PedometerPackage extends GrokPackage {

    init() {
        console.log('Pedometer package initialized.');
    }

    //tags: app
    startApp() {
    }

    //input: dataframe accel [Accelerometry data table]
    //input: column x {semType: Accelerometer-X} [X axis]
    //input: column y {semType: Accelerometer-Y} [Y axis]
    //input: column z {semType: Accelerometer-Z} [Z axis]
    pedometer(accel, x, y, z) {
        let view = gr.getTableView(accel.name);

        let viewer = view.markup();
        while (viewer.root.firstChild)
            viewer.root.removeChild(viewer.root.firstChild);
        let time = accel.col('time');
        let sampleRate = 1.0 / (time.get(1) - time.get(0));
        viewer.root.appendChild(ui.divText('Sample rate: ' + sampleRate.toString() + ' Hz'));

        gr.callFunc('Pedometer:DetectSteps', {
            "accel": accel.d,
            "x": x.name,
            "y": y.name,
            "z": z.name,
            "sample_rate": sampleRate
        }, true)
            .then(result => {
                let steps = accel.cols.add(new DataFrame(result).col('steps'), true);
                view.lineChart({
                    xColumnName: time.name,
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

                accel.onFilterChanged((_) => filtered.innerText = 'Filtered steps: ' + countSteps(accel.filter));
                accel.onSelectionChanged((_) => selected.innerText = 'Selected steps: ' + countSteps(accel.selection));
            });
    }
}
