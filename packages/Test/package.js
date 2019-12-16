class TestPackage extends GrokPackage {

    //description: test code
    test() {
        let view = gr.addTableView(gr.testData('demog', 5000));
        view.grid.sort(['AGE', 'HEIGHT'], [true, false]);

        view.grid.onCellPrepare(function (gc) {
            gc.style.backColor = 0xFFFFFF00;
        });
    }
}
