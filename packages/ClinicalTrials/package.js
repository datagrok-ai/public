class ClinicalTrialsPackage extends GrokPackage {

    //tags: app
    start(context) {
         grok.balloon.info("Clinical Trials App launched.");

         grok.script('DbQuery(Demo:Aact:Aact, "countries", aggregations=["count(*)"], fields=["name"], groupByFields=["name"])')
           .then(t => grok.addTableView(t));
    }
}