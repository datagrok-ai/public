class PlateManagementApp extends GrokPackage {
    start(context) {
         let plates = gr.testData('wells', 5000);
         let view = gr.addTableView(plates);
         view.addViewer('shape map');
    }
}