class WeatherApp extends GrokPackage {
    start(context) {
         gr.addTableView(gr.testData('demog', 5000));
    }
}