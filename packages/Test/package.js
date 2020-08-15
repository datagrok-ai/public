class TestPackage extends DG.Package {

    //description: test code
    test() {
        let view = grok.shell.addTableView(grok.data.demo.demog());
        view.grid.sort(['AGE', 'HEIGHT'], [true, false]);

        view.grid.onCellPrepare(function (gc) {
            gc.style.backColor = 0xFFFFFF00;
        });
    }

    //name: AsyncTest
    //output: string userInfo
    async asyncTest() {
        let userInfo = await fetch('https://api.github.com/users/skalkin');
        return userInfo;
    }
}
