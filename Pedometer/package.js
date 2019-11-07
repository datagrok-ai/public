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
        let sampleRate = 32
        // Sample rate input
        // Execute StepsDetector
        // Line Chart with columns
    }
}
