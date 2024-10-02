import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoView} from './demo-app'
import {DEMO_APP_HIERARCHY} from './const';
import {demoApp} from '../package';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import { TreeViewGroup } from 'datagrok-api/dg';


export class DemoAppWidget extends DG.Widget { 

    constructor() {
        super(ui.panel([], 'tutorial-widget'));
        
        const demoView = new DemoView();
        const funcs = demoView.funcs;
        const searchInput = ui.input.search('',{value: ''});
        
        grok.shell.windows.showToolbox = true;
        grok.shell.windows.showRibbon = true;
        grok.shell.windows.showHelp = true;
        grok.shell.windows.showProperties = true;
        grok.shell.windows.help.syncCurrentObject = true;
       
        let tree = ui.tree();
        tree.root.classList.add('demo-app-widget');

        for (let i = 0; i < DEMO_APP_HIERARCHY.children.length; ++i) {
            const directionFuncs = funcs.filter((func) => {
              return (func.func.options[DG.FUNC_OPTIONS.DEMO_PATH] as string).includes(DEMO_APP_HIERARCHY.children[i].name);
            });
            
            for (let j = 0; j < directionFuncs.length; ++j) {
              const path = directionFuncs[j].path.split('|').map((s) => s.trim());
      
              if (path.length > 2) {
                let groupPath = path[0];
                let treePath = tree.getOrCreateGroup(path[0], {path: groupPath}, false);
                (treePath.root.firstElementChild as HTMLElement).dataset.name = path[0];
                for (let i = 1; i < path.length - 1; i++) {
                  groupPath += `/${path[i]}`;
                  treePath = treePath.getOrCreateGroup(path[i], {path: groupPath}, false);
                  (treePath.root.firstElementChild as HTMLElement).dataset.name = path[i];
                }
      
                const item = treePath.item(directionFuncs[j].name, {path: directionFuncs[j].path});
                item.root.onmouseover = (event) => {
                  const packageMessage = `Part of the ${directionFuncs[j].func.package.name === 'Tutorials' ?
                    'platform core' : `${directionFuncs[j].func.package.name} package`}`;
                  ui.tooltip.show(directionFuncs[j].func.description ?
                    ui.divV([directionFuncs[j].func.description, ui.element('br'), packageMessage]) :
                    ui.div(packageMessage), event.clientX, event.clientY);
                };
      
                item.root.onmouseout = (_) => {
                  ui.tooltip.hide();
                };
              } else {
                const folder = tree.getOrCreateGroup(directionFuncs[j].category, {path: path[0]}, false);
                (folder.root.firstElementChild as HTMLElement).dataset.name = directionFuncs[j].category;
                const item = folder.item(directionFuncs[j].name, {path: directionFuncs[j].path});
      
                item.root.onmouseover = (event) => {
                  const packageMessage = `Part of the ${directionFuncs[j].func.package.name === 'Tutorials' ?
                    'platform core' : `${directionFuncs[j].func.package.name} package`}`;
                  ui.tooltip.show(directionFuncs[j].func.description ?
                    ui.divV([directionFuncs[j].func.description, ui.element('br'), packageMessage]) :
                    ui.div(packageMessage), event.clientX, event.clientY);
                };
      
                item.root.onmouseout = (_) => {
                  ui.tooltip.hide();
                };
              }
            }
          }

          DG.debounce(tree.onSelectedNodeChanged, 300).subscribe(async (value) => {
            if (value.root.classList.contains('d4-tree-view-item')) {
                demoApp();
                demoView.tree.items.find(node => node.text === value.text)?.root.click();
            }
          });

          searchInput.onChanged.subscribe((value) => {
            
            const foundFuncs = funcs.filter((func) => {
                return func.name.toLowerCase().includes(value.toLowerCase()) ||
                  func.func.description.toLowerCase().includes(value.toLowerCase()) ||
                  func.keywords.toLowerCase().includes(value.toLowerCase())
              });
              
              const dom = tree.root.getElementsByClassName('d4-tree-view-node');
          
              for (let i = 0; i < dom.length; i++) {
                const item = dom[i] as HTMLElement;
                const foundFunc = foundFuncs.find((func) => func.name.toLowerCase() === item.innerText.toLowerCase());
                if (foundFunc) {
                  const foundFuncPath = foundFunc.path.split('|').map((s) => s.trim());
                  item.classList.remove('hidden');
                  if (item.classList.contains('d4-tree-view-item')) {
                    for (let i = foundFuncPath.length - 2; i >= 0; i--) {
                      const currentCategory = tree.root.querySelector(`[data-name="${foundFuncPath[i]}"]`);
                      currentCategory?.classList.remove('hidden');
                      let group = tree.items.find(item => item.text === foundFuncPath[i]) as TreeViewGroup;
                      group.expanded = true;
                    }
                  }
                }
                else if (item.innerText.toLowerCase().includes(value.toLowerCase())) {
                  item.classList.remove('hidden');
                }
                else
                  item.classList.add('hidden');

                if (value == '') {
                    for (let i = 0; i < tree.items.length; i++){
                        let group = tree.items[i] as TreeViewGroup;
                        group.expanded = false;
                    }
                }
              }
          });
      
          searchInput.input.onkeyup = (event) => {
            if (event.key === 'Escape')
              searchInput.fireChanged();
          };
      
          const closeIcon = searchInput.root.getElementsByClassName('ui-input-icon-right')[0] as HTMLElement;
          closeIcon.onclick = () => {
            searchInput.value = '';
            searchInput.fireChanged();
          };

          this.root.append(ui.divV([searchInput.root, tree.root]));
    }

}
