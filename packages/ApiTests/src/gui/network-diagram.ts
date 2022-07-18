import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Network Diagram', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('networkDiagram.visual', async () => {
    let networkDiagramIcon = document.getElementsByClassName('svg-network-diagram')[0] as HTMLElement;
    networkDiagramIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Network diagram');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Network diagram') {
            throw 'Network diagram viewer was not closed'
        }
    }
  }); 
  test('networkDiagram.api', async () => {
    let networkDiagram = v.addViewer(DG.VIEWER.NETWORK_DIAGRAM, {
        "node1ColumnName": "age",
        "node2ColumnName": "race",
        "edgeColorColumnName": "height",
        "edgeWidth": 0.15
      }); await delay(500);

    if (networkDiagram.props.node1ColumnName != 'age')
        throw 'Node1 column has not been set' 
    if (networkDiagram.props.node2ColumnName != 'race')
        throw 'Node2 column has not been set'     
    if (networkDiagram.props.edgeColorColumnName != 'height')
        throw 'Edge column has not been set'   
    if (networkDiagram.props.edgeWidth != 0.15)
        throw 'edgeWidth prooerty has not been set'    

    networkDiagram.setOptions({
        title: 'Test Network Diagram',
        node1Shape: 'box',
        suspendSimulation: true 
    }); await delay(500);

    if (networkDiagram.props.title != 'Test Network Diagram')
        throw 'Title property has not been set' 
    if (networkDiagram.props.node1Shape != 'box')
        throw 'node1Shape property has not been set to box' 
    if (!networkDiagram.props.suspendSimulation)
        throw 'suspendSimulation property has not been set to TRUE'
  });  
  test('networkDiagram.serialization', async () => {
    let project = DG.Project.create();
    project.name = 'Test project with Network Diagram'
    project.addChild(demog.getTableInfo());
    project.addChild(v.saveLayout());  
    await grok.dapi.layouts.save(v.saveLayout());
    await grok.dapi.tables.uploadDataFrame(demog);
    await grok.dapi.tables.save(demog.getTableInfo());
    await grok.dapi.projects.save(project);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with Network Diagram');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Network diagram');

    let networkDiagram:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Network diagram') {
            networkDiagram = Array.from(v.viewers)[i];
            break;
        }
    }

    if (networkDiagram!.props.node1ColumnName != 'age')
        throw 'Node1 column has not been deserialized' 
    if (networkDiagram!.props.node2ColumnName != 'race')
        throw 'Node2 column has not been deserialized'     
    if (networkDiagram!.props.edgeColorColumnName != 'height')
        throw 'Edge golumn has not been deserializedt'   
    if (networkDiagram!.props.edgeWidth != 0.15)
        throw 'edgeWidth prooerty has not been deserialized'   
    if (networkDiagram!.props.node1Shape != 'box')
        throw 'node1Shape property has not been deserialized' 
    if (!networkDiagram!.props.suspendSimulation)
        throw 'suspendSimulation property has not been deserialized'
    if (networkDiagram!.props.title != 'Test Network Diagram')
        throw 'tittle property has not been deserialized' 
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Network Diagram').first())
  }); 
});
