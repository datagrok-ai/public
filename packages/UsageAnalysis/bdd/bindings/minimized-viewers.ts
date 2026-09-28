/* Minimized viewers: the title-bar icon that minimizes a viewer, the icon it becomes in the ribbon's
   Minimized panel, and the filter panel's reset icon. A minimized viewer leaves the view (it is no
   longer among the table view's viewers) and is shown, live, inside a super tooltip while its ribbon
   icon is hovered, so the viewer kind reaches the preview by itself. The minimize icon shows only
   while the viewer's title bar is hovered. */
import {element} from '@datagrok-libraries/bdd';

const PANEL = 'xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " panel-base ")][1]';
const minimizeIcon = (viewer: string): string =>
  `[name="viewer-${viewer}"] >> ${PANEL}//*[contains(@class, "panel-titlebar")]//*[@name="icon-window-minimize"]`;

element('scatter plot minimize icon', {selector: minimizeIcon('Scatter-plot')});
element('histogram minimize icon', {selector: minimizeIcon('Histogram')});
element('minimized scatter plot icon', {selector: '.grok-minimized-viewer:has([name="icon-scatter-plot"])'});
element('minimized histogram icon', {selector: '.grok-minimized-viewer:has([name="icon-histogram"])'});
element('filter panel reset icon', {selector: '[name="viewer-Filters"] [name="icon-arrow-rotate-left"]'});

