import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoView} from './demo-app'
import {DEMO_APP_HIERARCHY} from './const';
import {demoApp} from '../package';


export class DemoAppWidget extends DG.Widget { 

    constructor() {
        super(ui.panel([], 'tutorial-widget'));
        
        const demoView = new DemoView();
        const funcs = demoView.funcs;
        
        grok.shell.dockManager.close(demoView.dockPanel);
        grok.shell.windows.showToolbox = true;
        grok.shell.windows.showRibbon = true;
        grok.shell.windows.showHelp = true;
        grok.shell.windows.showProperties = true;
        grok.shell.windows.help.syncCurrentObject = true;

        let acc = ui.accordion();
        DEMO_APP_HIERARCHY.children.forEach(item => {
            acc.addCountPane(item.name, () => {
                const root = ui.divV([], {style:{gap:'4px'}});
                funcs.forEach(element => {
                    if(element.category === item.name)
                        root.append(ui.link(element.name, () => {
                            demoApp();
                            demoView.tree.items.find(node => node.text === element.name)?.root.click();
                        }))
                })
                return root;
            }, () => {
                return funcs.filter(direction => direction.category === item.name).length;
            }, false)
        })
        this.root.append(acc.root);
    }

}
