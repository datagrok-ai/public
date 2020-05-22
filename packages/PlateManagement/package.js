class PlateManagementApp extends DG.Package {
    start(context) {
         let plates = grok.data.testData('wells', 5000);
         let view = grok.shell.addTableView(plates);
         view.addViewer('shape map');
    }
}