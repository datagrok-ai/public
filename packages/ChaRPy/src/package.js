/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

//name: toScript
//top-menu: convert | toScript
export async function toScript() {


    function strReplace(string, propsObj) {

        let toRemove = Object.keys(propsObj['look'])
        let toInsert = Object.values(propsObj['look'])

        var i;
        let newString;
        for (i = 0; i < toRemove.length; i++) {

            newString = string.replace(toRemove[i], toInsert[i]);
            if (toRemove[i] === 'showRegressionLine' || toInsert[i] === true) {
                newString = string.concat(` + geom_smooth(method=lm, se=TRUE)`);
            }

            string = newString;
        }
        return newString;
    }

    let view = grok.shell.addTableView(grok.data.demo.demog());
    let plot = view.scatterPlot({
        x: 'height',
        y: 'weight',
        size: 'age',
        color: 'race'
    });
    plot.setOptions({
        showRegressionLine: true,
        markerType: 'square'
    });

    let properties = JSON.parse(plot.getOptions());
    let rString = `require(ggplot2)
ggplot(data, aes(x=xColumnName, y=yColumnName, size=sizeColumnName, color=colorColumnName)) + geom_point() +
scale_size_continuous(range = c(markerMinSize,markerMaxSize))`;

    let rCode = strReplace(rString, properties)

    var input = document.createElement('TEXTAREA');
    input.value = rCode
    ui.dialog('Windows')
        .add(input)
        .onOK(() => { grok.shell.info('OK!'); })
        .show();
}


