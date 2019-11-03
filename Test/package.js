class TestPackage extends GrokPackage {

    // Guaranteed to get called exactly once before the execution of any function below
    init() { console.log('Test package initialized.'); }

    //description: test code
    test() {
        let view = gr.addTableView(gr.testData('demog', 5000));
        view.grid.sort(['AGE', 'HEIGHT'], [true, false]);

        view.grid.onCellPrepare(function (gc) {
            gc.style.backColor = 0xFFFFFF00;
        });
    }
}
