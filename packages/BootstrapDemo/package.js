class BootstrapDemoPackage extends GrokPackage {

    //name: Custom UI demo
    //tags: app
    //meta.icon-path: coreui.png
    //description: This is an app built to demonstrate how to make fully custom UI
    startApp() {
        var webRoot = this.webRoot;
        var x = new XMLHttpRequest();
        x.open("GET", this.webRoot + "index.html", true);
        x.onload = function (){
            document.body.className = 'app header-fixed sidebar-fixed aside-menu-fixed sidebar-lg-show  pace-done';
            document.body.innerHTML = x.responseText.replace(/#{WEB_ROOT}/g, webRoot);

            grok.openTable("e1792340-bd30-11e8-ed33-716dfb17fcbb").then(function (demog) {
                grok.balloon.info(demog.name);
                var view = grok.addTableView(demog);
                var scatterPlot = view.addViewer('scatter plot');
                document.body.querySelector("#scatterplot").style.height = '500px';
                document.body.querySelector("#scatterplot").appendChild(scatterPlot.root);
            } );
            //var demog  = grok.testData('demog', 5000);

            //here: rewrite links, rewrite paths, etc.
        }
        x.send(null);
    }


    //tags: test
    //meta.icon-path: package.png
    test() {
        alert(this.webRoot);
    }

}
