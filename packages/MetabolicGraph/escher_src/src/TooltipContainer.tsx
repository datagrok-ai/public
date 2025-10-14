/** @jsx h */
import {Ref} from 'preact';
import CallbackManager from './ts/CallbackManager';
import PlacedDiv from './ts/PlacedDiv';
import renderWrapper from './renderWrapper';
import _ from 'underscore';
import ZoomContainer from './ts/ZoomContainer';
import {EscherMap} from './ts/escherMap';
import Settings from './ts/Settings';
import {D3Selection, TooltipComponentProps} from './ts/types';

/**
 * Manage the tooltip that lives in a PlacedDiv.
 * @param selection
 * @param tooltipComponent
 * @param zoomContainer
 * @param map
 */
export default class TooltipContainer {
  div: D3Selection;
  tooltipRef: Ref<any>;
  zoomContainer: ZoomContainer;
  callbackManager: CallbackManager;
  map: EscherMap;
  settings: Settings;
  delayHideTimeout: any;
  currentTooltip: {type: string, id: string} | null;
  placedDiv: PlacedDiv;
  constructor(selection: D3Selection, TooltipComponent: any, zoomContainer: ZoomContainer, map: EscherMap, settings: Settings) {
    this.div = selection.append('div').attr('id', 'tooltip-container');
    this.tooltipRef = null;

    this.zoomContainer = zoomContainer;
    this.setUpZoomCallbacks(zoomContainer);

    // Create callback manager
    this.callbackManager = new CallbackManager();

    this.div.on('mouseover', this.cancelHideTooltip.bind(this));
    this.div.on('mouseleave', this.hide.bind(this));

    this.map = map;
    this.setUpMapCallbacks(map);

    this.settings = settings;

    this.delayHideTimeout = null;
    this.currentTooltip = null;

    renderWrapper(
      TooltipComponent,
      null,
      (passProps) => this.callbackManager.set('pass_props', passProps),
      this.div.node(),
      (instance) => { this.tooltipRef = instance; }
    );
    this.passProps({
      display: false,
      disableTooltips: () => this.disableTooltips()
    });
  }

  /**
   * Disable tooltips in the settings
   */
  disableTooltips() {
    this.settings.set('enable_tooltips', false);
    // draw to update tooltip settings
    this.map.draw_everything();
    this.hide();
    this.map.set_status(`Tooltips disabled. You can enable them again in the
                         settings menu.`, 3000);
  }

  /**
    * Function to pass props for the tooltips. Run without an argument to
    * rerender the component
    * @param {Object} props - Props that the tooltip will use
    */
  passProps(props: Partial<TooltipComponentProps> = {}) {
    this.callbackManager.run('pass_props', null, props);
  }

  /**
   * Sets up the appropriate callbacks to show the input
   * @param {object} map - map object
   */
  setUpMapCallbacks(map: EscherMap) {
    this.placedDiv = new PlacedDiv(this.div, map, undefined, false);

    // connect callbacks to show tooltip
    map.callback_manager.set('show_tooltip.tooltip_container', (type: string, d: {bigg_id: string, name: string, data_string: string, xPos: number, yPos: number, label_x: number, label_y: number}) => {
      this.show(type, d);
    });

    // callback to hide / delay hide tooltip
    map.callback_manager.set('hide_tooltip.tooltip_container', () => this.hide());
    map.callback_manager.set('delay_hide_tooltip.tooltip_container', () => this.delayHide());

    // Hides the tooltip when the canvas is touched
    map.sel.selectAll('.canvas-group').on('touchend', () => this.hide());
  }

  setUpZoomCallbacks(zoomContainer: ZoomContainer) {
    zoomContainer.callbackManager.set('zoom.tooltip_container', function() {
      if (this.is_visible())
        this.hide();
    }.bind(this));
    zoomContainer.callbackManager.set('go_to.tooltip_container', function() {
      if (this.is_visible())
        this.hide();
    }.bind(this));
  }

  /**
   * Return visibility of tooltip container.
   * @return {Boolean} Whether tooltip is visible.
   */
  is_visible() { // eslint-disable-line camelcase
    return this.placedDiv.is_visible();
  }

  /**
   * Show and place the input.
   * @param {string} type - 'reaction_label', 'node_label', or 'gene_label'
   * @param {Object} d - D3 data for DOM element
   */
  show(type: string, d: {bigg_id: string, name: string, data_string: string, xPos: number, yPos: number, label_x: number, label_y: number}) {
    // get rid of a lingering delayed hide
    this.cancelHideTooltip();

    if (_.contains(['reaction_label', 'node_label', 'gene_label', 'reaction_object', 'node_object'], type)) {
      // Use a default height if the ref hasn't been connected yet
      const tooltipSize = (this.tooltipRef !== null && 'get_size' in this.tooltipRef && typeof this.tooltipRef.get_size === 'function') ?
        this.tooltipRef.get_size() :
        {width: 270, height: 100};
      this.currentTooltip = {type, id: d[type.replace('_label', '_id').replace('_object', '_id')]};
      const windowTranslate = this.zoomContainer.windowTranslate;
      const windowScale = this.zoomContainer.windowScale;
      const mapSize = this.map !== null ? this.map.get_size() : {width: 1000, height: 1000};
      const offset = {x: 0, y: 0};
      const startPosX = type === 'reaction_object' ? d.xPos : d.label_x;
      const startPosY = type === 'reaction_object' ? d.yPos : d.label_y;
      const rightEdge = windowScale * startPosX + windowTranslate.x + tooltipSize.width;
      const bottomEdge = windowScale * startPosY + windowTranslate.y + tooltipSize.height;
      if (mapSize.width < 500) {
        if (rightEdge > mapSize.width)
          offset.x = -(rightEdge - mapSize.width) / windowScale;

        if (bottomEdge > mapSize.height - 74)
          offset.y = -(bottomEdge - mapSize.height + 77) / windowScale;
      } else {
        if (windowScale * startPosX + windowTranslate.x + 0.5 * tooltipSize.width > mapSize.width)
          offset.x = -tooltipSize.width / windowScale;
        else if (rightEdge > mapSize.width)
          offset.x = -(rightEdge - mapSize.width) / windowScale;

        if (windowScale * startPosY + windowTranslate.y + 0.5 * tooltipSize.height > mapSize.height - 45)
          offset.y = -(tooltipSize.height) / windowScale;
        else if (bottomEdge > mapSize.height - 45)
          offset.y = -(bottomEdge - mapSize.height + 47) / windowScale;
      }
      const coords = {x: startPosX + offset.x, y: startPosY + 10 + offset.y};
      this.placedDiv.place(coords);
      this.passProps({
        display: true,
        biggId: d.bigg_id,
        name: d.name,
        loc: coords,
        data: d.data_string,
        type: type.replace('_label', '').replace('node', 'metabolite').replace('_object', '')
      });
    } else
      throw new Error('Tooltip not supported for object type ' + type);
  }

  /**
   * Hide the input.
   */
  hide() {
    this.placedDiv.hide();
    this.currentTooltip = null;
  }

  /**
   * Hide the input after a short delay, so that mousing onto the tooltip does not
   * cause it to hide.
   */
  delayHide() {
    this.delayHideTimeout = setTimeout(() => this.hide(), 100);
  }

  cancelHideTooltip() {
    if (this.delayHideTimeout !== null)
      clearTimeout(this.delayHideTimeout);
  }
}
