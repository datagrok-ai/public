import * as utils from './utils';
import PlacedDiv from './PlacedDiv';
import completely from './completely';
import DirectionArrow from './DirectionArrow';
import CobraModel from './CobraModel';
import _ from 'underscore';
// @ts-ignore
import {mouse as d3Mouse} from 'd3-selection';
import {EscherMap} from './escherMap';
import {CobraGene, CobraMetabolite, Coord, D3Selection, MapNode, MapReaction} from './types';
import ZoomContainer from './ZoomContainer';
import Settings from './Settings';

/**
 * BuildInput
 * @param selection - A d3 selection for the BuildInput.
 * @param map - A Map instance.
 * @param zoomContainer - A ZoomContainer instance.
 * @param settings - A Settings instance.
 */

let selectedOrderCounter = 0;

export default class BuildInput {
  placed_div: PlacedDiv;
  map: EscherMap;
  completely: any;
  direction_arrow: DirectionArrow;
  settings: Settings;
  target_coords: {x: number, y: number} | null;
  start_reaction_listener: boolean = false;
  zoomContainer: ZoomContainer;
  is_active: boolean = false;
  clear_escape: any;

  constructor(selection: D3Selection, map: EscherMap, zoomContainer: ZoomContainer, settings: Settings) {
    // set up container
    const newSel = selection.append('div').attr('id', 'rxn-input');
    this.placed_div = new PlacedDiv(newSel, map, {x: 240, y: 0});
    this.placed_div.hide();

    // set up complete.ly
    this.completely = completely(newSel.node()!, {backgroundColor: '#eee'});

    // close button
    newSel.append('button').attr('class', 'button input-close-button')
      .text('×')
      .on('mousedown', () => this.hideDropdown());

    // map
    this.map = map;
    // set up the reaction direction arrow
    const defaultAngle = 90; // degrees
    this.direction_arrow = new DirectionArrow(map.sel);
    this.direction_arrow.setRotation(defaultAngle);
    this.setUpMapCallbacks(map);

    // zoom container
    this.zoomContainer = zoomContainer;
    this.setUpZoomCallbacks(zoomContainer);

    // settings
    this.settings = settings;

    // toggle off
    this.toggle(false);
    this.target_coords = null;
  }

  setUpMapCallbacks(map: EscherMap) {
    // input
    map.callback_manager.set('select_metabolite_with_id.input', (selectedNode: MapNode, coords: Coord) => {
      if (this.is_active) {
        const hasModel = this.reload(selectedNode, coords, false);
        if (hasModel) this.showDropdown(coords);
      }
      this.map.unhighlightShortestPaths();
      this.hideTarget();
    });
    map.callback_manager.set('select_selectable.input', (count: number, selectedNode: MapNode, coords: Coord) => {
      this.hideTarget();
      if (selectedNode && typeof selectedNode === 'object')
        selectedNode.selectionOrder = selectedOrderCounter++;

      if (count === 1 && this.is_active && coords) {
        const hasModel = this.reload(selectedNode, coords, false);
        if (hasModel) this.showDropdown(coords);
      } else
        this.toggle(false);

      if (count === 1 && this.settings.get('highlight_same_nodes'))
        this.map.highlightSameNodes();
      this.map.selectionChanged();
      this.map.unhighlightShortestPaths();
    });
    map.callback_manager.set('deselect_nodes', () => {
      this.map.unhighlightShortestPaths();
      this.direction_arrow.hide();
      this.hideDropdown();
    });

    // svg export
    map.callback_manager.set('before_svg_export', () => {
      this.direction_arrow.hide();
      this.hideTarget();
    });
  }

  setUpZoomCallbacks(zoomContainer: ZoomContainer) {
    // TODO this is broken.
    // Should place either for selected or for location on zoom or pan.
    // zoomContainer.callbackManager.set('zoom_change.input', () => {
    //   if (this.is_active) {
    //     this.place_at_selected()
    //   }
    // })
  }

  is_visible() { // eslint-disable-line camelcase
    return this.placed_div.is_visible();
  }

  toggle(onOff: boolean) {
    if (onOff === undefined) this.is_active = !this.is_active;
    else this.is_active = onOff;
    if (this.is_active) {
      this.toggleStartReactionListener(true);
      let hasModelAndSelection = true;
      if (_.isNull(this.target_coords))
        hasModelAndSelection = this.reloadAtSelected() as boolean;
      else
        this.placed_div.place(this.target_coords);

      if (hasModelAndSelection) {
        this.showDropdown();
        this.map.set_status('Click on the canvas or an existing metabolite');
      }
      this.direction_arrow.show();
    } else {
      this.toggleStartReactionListener(false);
      this.hideDropdown();
      this.map.set_status(null);
      this.direction_arrow.hide();
    }
  }

  showDropdown(coords?: Coord) {
    // escape key
    this.clear_escape = this.map.key_manager
      .addEscapeListener(() => this.hideDropdown());
    // dropdown
    this.completely.input.blur();
    this.completely.repaint();
    this.completely.setText('');
    this.completely.input.focus();
  }

  hideDropdown() {
    // escape key
    if (this.clear_escape) this.clear_escape();
    this.clear_escape = null;
    // dropdown
    this.placed_div.hide();
    this.completely.input.blur();
    this.completely.hideDropDown();
  }

  place(coords: Coord) {
    this.placed_div.place(coords);
    this.direction_arrow.setLocation(coords);
    this.direction_arrow.show();
  }

  /**
   * Reload data for autocomplete box and redraw box at the first selected node.
   * @return {Boolean} Returns true if a model is present and a node is selected.
   */
  reloadAtSelected() {
    // get the selected node
    this.map.deselect_text_labels();
    const selectedNode = this.map.select_single_node();
    if (selectedNode == null) return false;
    const coords = {x: selectedNode.x, y: selectedNode.y};
    // reload the reaction input
    return this.reload(selectedNode, coords, false);
  }

  alreadyDrawn(biggId: string, reactions: { [key: string]: MapReaction }) {
    for (const drawnId in reactions) {
      if (reactions[drawnId].bigg_id === biggId)
        return true;
    }
    return false;
  }

  /**
   * Reload data for autocomplete box and redraw box at the new coordinates.
   * @param {} selectedNode -
   * @param {} coords -
   * @param {Boolean} startingFromScratch -
   * @return {Boolean} Returns true if a model is present.
   */
  reload(selectedNode: MapNode | null, coords: Coord, startingFromScratch: boolean) {
    // Try finding the selected node
    if (!startingFromScratch && !selectedNode) {
      console.error('No selected node and not starting from scratch');
      return;
    }

    this.place(coords);

    if (this.map.cobra_model == null) {
      this.completely.setText('Cannot add: No model.');
      // this.completely.repaint()
      return false;
    }

    // settings
    const showNames = this.settings.get('identifiers_on_map') === 'name';
    const allowDuplicates = this.settings.get('allow_building_duplicate_reactions');

    // Find selected
    const options = [];
    const cobraReactions = this.map.cobra_model.reactions;
    const cobraMetabolites = this.map.cobra_model.metabolites;
    const reactions = this.map.reactions;
    const hasDataOnReactions = this.map.has_data_on_reactions;
    const selectedMetName = (selectedNode ? (showNames ? selectedNode.name : selectedNode.bigg_id) : '');
    const boldMetsInStr = (str: string, mets: string[]) =>
      str.replace(new RegExp('(^| )(' + mets.join('|') + ')($| )', 'g'), '$1<b>$2</b>$3');

    // for reactions
    const reactionSuggestions: { [key: string]: boolean } = {};
    for (const biggId in cobraReactions) {
      const reaction = cobraReactions[biggId];
      const reactionName = reaction.name;
      const showReactionName = (showNames ? reactionName : biggId);

      // ignore drawn reactions
      if ((!allowDuplicates) && this.alreadyDrawn(biggId, reactions))
        continue;


      // check segments for match to selected metabolite
      for (const metBiggId in reaction.metabolites) {
        // if starting with a selected metabolite, check for that id
        if (startingFromScratch || metBiggId === selectedNode?.bigg_id) {
          // don't add suggestions twice
          if (biggId in reactionSuggestions) continue;

          // get the metabolite names or IDs
          let mets: { [key: string]: CobraMetabolite } = {};
          const showMetNames = [];
          let metId;
          if (showNames) {
            for (metId in reaction.metabolites) {
              const name = cobraMetabolites[metId].name;
              mets[name] = reaction.metabolites[metId];
              showMetNames.push(name);
            }
          } else {
            mets = utils.clone(reaction.metabolites);
            for (metId in reaction.metabolites)
              showMetNames.push(metId);
          }
          const showGeneNames = _.flatten(
            reaction.genes.map((g: CobraGene) => [g.name, g.biggId])
          );
          // get the reaction string
          const reactionString = CobraModel.build_reaction_string(Object.fromEntries(Object.entries(mets).map(([key, value]) => [key, value.coefficient])),
            reaction.reversibility);
          // make the matches list and filter out any missing entries (e.g.
          // missing gene names from model
          const matches = [showReactionName].concat(showMetNames).concat(showGeneNames).filter((x) => x);

          if (hasDataOnReactions) {
            options.push({
              reaction_data: reaction.data,
              html: '<b>' + showReactionName + '</b>' + ': ' + reaction.data_string,
              matches,
              id: biggId
            });
          } else {
            options.push({
              html: ('<b>' + showReactionName + '</b>' + '\t' +
                     boldMetsInStr(reactionString, [selectedMetName!])),
              matches,
              id: biggId
            });
          }
          reactionSuggestions[biggId!] = true;
        }
      }
    }

    // Generate the array of reactions to suggest and sort it
    const sortFn = hasDataOnReactions ?
      (x: any, y: any) => Math.abs(x.reaction_data) > Math.abs(y.reaction_data) ? -1 : 1 :
      (x: any, y: any) => x.html.toLowerCase() < y.html.toLowerCase() ? -1 : 1;

    // set up the box with data
    this.completely.options = options.sort(sortFn);

    // TODO test this behavior
    // if (strings_to_display.length==1) this.completely.setText(strings_to_display[0])
    // else this.completely.setText("")
    this.completely.setText('');

    const checkAndBuild = (id: string) => {
      if (id !== null) {
        // make sure the selected node exists, in case changes were made in the meantime
        if (startingFromScratch) {
          this.map.new_reaction_from_scratch(id,
            coords,
            this.direction_arrow.getRotation());
        } else {
          if (!(selectedNode!.node_id! in this.map.nodes)) {
            console.error('Selected node no longer exists');
            this.hideDropdown();
            return;
          }
          this.map.new_reaction_for_metabolite(id,
            selectedNode!.node_id!,
            this.direction_arrow.getRotation());
        }
      }
    };
    this.completely.onEnter = function(id: string) {
      this.setText('');
      this.onChange('');
      checkAndBuild(id);
    };

    return true;
  }

  /**
   * Toggle listening for a click to place a new reaction on the canvas.
   */
  toggleStartReactionListener(onOff: boolean) {
    if (onOff === undefined)
      this.start_reaction_listener = !this.start_reaction_listener;
    else if (this.start_reaction_listener === onOff)
      return;
    else
      this.start_reaction_listener = onOff;


    if (this.start_reaction_listener) {
      const node = this.map.sel.node();
      this.map.sel.on('click.start_reaction', () => {
        // TODO fix this hack
        if (this.direction_arrow.dragging) return;
        // reload the reaction input
        const coords = {
          x: d3Mouse(node)[0],
          y: d3Mouse(node)[1]
        };
        // unselect metabolites
        this.map.deselect_nodes();
        this.map.deselect_text_labels();
        // reload the reaction input
        const hasModel = this.reload(null, coords, true);
        if (hasModel) {
          // show the dropdown
          this.showDropdown(coords);
        }
        // generate the target symbol
        this.showTarget(this.map, coords);
      });
      this.map.sel.style('cursor', 'pointer');
    } else {
      this.map.sel.on('click.start_reaction', null);
      this.map.sel.style('cursor', null);
      this.hideTarget();
    }
  }

  hideTarget() {
    if (this.target_coords)
      this.map.sel.selectAll('.start-reaction-target').remove();

    this.target_coords = null;
  }

  showTarget(map: EscherMap, coords: Coord) {
    const s = map.sel.selectAll('.start-reaction-target').data([12, 5]);
    s.enter()
      .append('circle')
      .classed('start-reaction-target', true)
      .attr('r', function(d) { return d; })
      .style('stroke-width', 4)
      .merge(s as any)
      .style('visibility', 'visible')
      .attr('transform', 'translate(' + coords.x + ',' + coords.y + ')');
    this.target_coords = coords;
  }
}
