/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Network Diagram', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('networkDiagram.visual', async () => {
    const networkDiagramIcon = document.getElementsByClassName('svg-network-diagram')[0] as HTMLElement;
    networkDiagramIcon.click(); 
    await awaitCheck(() => document.querySelector('.d4-network-diagram') !== null, 'network diagram not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Network diagram');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null,
      'network diagram properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Network diagram viewer was not closed', 1000);
  });

  test('networkDiagram.api', async () => {
    const networkDiagram = v.addViewer(DG.VIEWER.NETWORK_DIAGRAM, {
      'node1ColumnName': 'age',
      'node2ColumnName': 'race',
      'edgeColorColumnName': 'height',
      'edgeWidth': 0.15,
    });
    await awaitCheck(() => document.querySelector('.d4-network-diagram') !== null, 'network diagram not found', 3000);

    if (networkDiagram.props.node1ColumnName != 'age')
      throw 'Node1 column has not been set'; 
    if (networkDiagram.props.node2ColumnName != 'race')
      throw 'Node2 column has not been set';     
    if (networkDiagram.props.edgeColorColumnName != 'height')
      throw 'Edge column has not been set';   
    if (networkDiagram.props.edgeWidth != 0.15)
      throw 'edgeWidth prooerty has not been set';    

    networkDiagram.setOptions({
      title: 'Test Network Diagram',
      node1Shape: 'box',
      suspendSimulation: true, 
    });

    await awaitCheck(() => networkDiagram.props.title === 'Test Network Diagram',
      'title property has not been set', 1000);
    if (networkDiagram.props.node1Shape != 'box')
      throw 'node1Shape property has not been set to box'; 
    if (!networkDiagram.props.suspendSimulation)
      throw 'suspendSimulation property has not been set to TRUE';
  });  

  // Does not work through Test Manager
  test('networkDiagram.serialization', async () => {
    await uploadProject('Test project with Network Diagram', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Network Diagram');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Network diagram');
    const networkDiagram = findViewer('Network diagram', v);

    if (networkDiagram!.props.node1ColumnName != 'age')
      throw 'Node1 column has not been deserialized'; 
    if (networkDiagram!.props.node2ColumnName != 'race')
      throw 'Node2 column has not been deserialized';     
    if (networkDiagram!.props.edgeColorColumnName != 'height')
      throw 'Edge golumn has not been deserializedt';   
    if (networkDiagram!.props.edgeWidth != 0.15)
      throw 'edgeWidth prooerty has not been deserialized';   
    if (networkDiagram!.props.node1Shape != 'box')
      throw 'node1Shape property has not been deserialized'; 
    if (!networkDiagram!.props.suspendSimulation)
      throw 'suspendSimulation property has not been deserialized';
    if (networkDiagram!.props.title != 'Test Network Diagram')
      throw 'tittle property has not been deserialized'; 
  }); 

  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Network Diagram').first());
  }); 
});
