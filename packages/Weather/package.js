class WeatherApp extends GrokPackage {
    start(context) {
         grok.addTableView(grok.testData('demog', 5000));
    }
}