class ClinicalTrialsApp extends GrokPackage {
    start(context) {
         gr.balloon.info("Clinical Trials App launched.");

         gr.script('DbQuery(Demo:Aact:Aact, "countries", aggregations=["count(*)"], fields=["name"], groupByFields=["name"])')
           .then(t => gr.addTableView(t));
    }
}