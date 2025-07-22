/* eslint-disable no-invalid-this */
/* eslint-disable camelcase */
/**
 * PlacedDiv. A container to position an html div to match the coordinates of a
 * SVG element.
 */
import {EscherMap} from './escherMap';
import {Coord, D3Selection} from './types';
import * as utils from './utils';


export default class PlacedDiv {
  private div!: D3Selection;
  private map!: EscherMap;
  private displacement!: {x: number, y: number};
  private shouldReposition!: boolean;
  private visible!: boolean;

  constructor(div: D3Selection<any>, map: EscherMap, displacement = {x: 0, y: 0}, shouldReposition = true) {
    this.div = div;

    this.map = map;

    this.displacement = displacement;

    this.shouldReposition = shouldReposition;
    // begin hidden

    this.visible = true;

    this.hide();
  }

  is_visible(): boolean {
    return this.visible;
  }

  /**
   * Position the html div to match the given SVG coordinates.
   */
  place(coords: Coord) {
    // show the input

    this.div.style('display', null);

    // move the new input

    const window_translate = this.map.zoomContainer.windowTranslate;

    const window_scale = this.map.zoomContainer.windowScale;

    const map_size = this.map.get_size();

    /**
     * If shouldReposition is true, the div is placed so that it does not render
     * outside of the viewable area of the window. Math.max is used so that it
     * does not overflow to the left and top, Math.min is used so that it does not
     * overflow to the right or bottom. If the screen is tool small to show the
     * entire div, the div will overflow to the right and bottom.
     */

    if (this.shouldReposition) {
      const left = Math.max(20,
        Math.min(map_size.width - 270,
          (window_scale * coords.x + window_translate.x -

                                  this.displacement.x)));
      const top = Math.max(20,
        Math.min(map_size.height - 40,
          (window_scale * coords.y + window_translate.y -

                                  this.displacement.y)));

      this.div.style('position', 'absolute')

        .style('display', 'block')

        .style('left', `${left}px`)

        .style('top', `${top}px`);
    } else {
      this.div.style('position', 'absolute')
        .style('display', 'block')

        .style('left', `${window_scale * coords.x + window_translate.x - this.displacement.x}px`)

        .style('top', `${window_scale * coords.y + window_translate.y - this.displacement.y}px`);
    }

    this.visible = true;
  }

  /**
   * Hide the PlacedDiv.
   */
  hide() {
    if (this.visible) {
      this.div.style('display', 'none');

      this.visible = false;
    }
  }
}
