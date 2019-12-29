class PlateManagementApp extends GrokPackage {
    start(context) {
         let plates = grok.testData('wells', 5000);
         let view = grok.addTableView(plates);
         view.addViewer('shape map');
    }
}