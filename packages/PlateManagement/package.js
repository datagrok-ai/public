class PlateManagementApp extends grok.Package {
    start(context) {
         let plates = grok.testData('wells', 5000);
         let view = grok.addTableView(plates);
         view.addViewer('shape map');
    }
}