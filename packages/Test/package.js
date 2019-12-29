class TestPackage extends GrokPackage {

    //description: test code
    test() {
        let view = grok.addTableView(grok.testData('demog', 5000));
        view.grid.sort(['AGE', 'HEIGHT'], [true, false]);

        view.grid.onCellPrepare(function (gc) {
            gc.style.backColor = 0xFFFFFF00;
        });
    }

    //name: AsyncTest
    //output: string userInfo
    async asyncTest() {
        var userInfo = await fetch('https://api.github.com/users/skalkin');
        return userInfo;
    }
}
