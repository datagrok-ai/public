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

        let viewer = view.markup({content: ''});
        let sampleRate = ui.intInput('Sample rate, Hz', 32);
        viewer.root.appendChild(ui.inputs([sampleRate]));

        gr.callFunc('Demo:PythonScripts:DetectSteps', {"accel": accel, "x": x, "y": y, "z": z, "sample_rate": sampleRate.value})
            .then(result => {
                view.lineChart({
                    xColumnName: 'time',
                    yColumnNames: [x.name, y.name, z.name, 'steps']
                })
            });
    }
}
