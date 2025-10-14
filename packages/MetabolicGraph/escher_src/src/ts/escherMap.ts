/* eslint-disable no-invalid-this */
/* eslint-disable camelcase */
import * as utils from './utils';
import Draw from './Draw';
import Behavior from './Behavior';
import Scale from './Scale';
import * as build from './build';
import UndoStack from './UndoStack';
import CallbackManager from './CallbackManager';
import KeyManager from './KeyManager';
import Canvas from './Canvas';
import * as dataStyles from './dataStyles';
import SearchIndex from './SearchIndex';

import bacon from 'baconjs';
import _ from 'underscore';
import {BaseType, select as d3Select} from 'd3-selection';
import {findOutGoing, findInGoing, getReactionPathSegments, findKthShortestPath, pathSeparatorString} from './algo';
import {Coord, D3Selection, D3SVGSelection, EscherMapDataType, IDType, MapData, MapNode, MapReaction, MapSegment, SettingsType, TextLabel} from './types';
import Settings from './Settings';
import ZoomContainer from './ZoomContainer';
import CobraModel from './CobraModel';

// interface ZoomContainer {
//   get_size(): { width: number; height: number };
//   goTo(zoom: number, pos: { x: number; y: number }): void;
//   windowScale: number;
//   windowTranslate: { x: number; y: number };
//   _goToSvg(scale: number, translate: { x: number; y: number }, callback: () => void): void;
// }


interface BezierDisplacement {
  reactionId: string;
  segmentId: string;
  bez: string;
  bezierId: string;
  displacement: { x: number; y: number };
}

function _on_array(fn: Function) {
  return function(array: Array<any>) { return fn(...array); };
}

/**
 * Map - Defines the metabolic map data, and manages drawing and building.
 * @param svg: The parent SVG container for the map.
 * @param css:
 * @param selection: A d3 selection for a node to place the map inside.
 * @param selection:
 * @param zoomContainer:
 * @param settings:
 * @param cobra_model:
 * @param canvas_size_and_loc:
 * @param enable_search:
 * @param map_name: (Optional, Default: 'new map')
 * @param map_id: (Optional, Default: A string of random characters.)
 * @param map_description: (Optional, Default: '')
 *
 * Callbacks
 * ---------
 *
 * map.callback_manager.run('set_status', null, status)
 * map.callback_manager.run('toggle_beziers', null, beziers_enabled)
 * map.callback_manager.run('select_metabolite_with_id', null, selected_node, coords)
 * map.callback_manager.run('select_selectable', null, node_count, selected_node, coords)
 * map.callback_manager.run('deselect_nodes')
 * map.callback_manager.run('select_text_label')
 * map.callback_manager.run('before_svg_export')
 * map.callback_manager.run('after_svg_export')
 * map.callback_manager.run('before_png_export')
 * map.callback_manager.run('after_png_export')
 * map.callback_manager.run('before_convert_map')
 * map.callback_manager.run('after_convert_map')
 * this.callback_manager.run('calc_data_stats__reaction', null, changed)
 * this.callback_manager.run('calc_data_stats__metabolite', null, changed)
 *
 */
export class EscherMap {
  callback_manager: CallbackManager;
  svg: D3SVGSelection;
  defs: D3Selection<SVGDefsElement>;
  canvas: Canvas;
  sel: D3Selection<Element>;
  zoomContainer: ZoomContainer;
  settings: Settings;
  cobra_model: CobraModel | null;
  reactions: { [key: string]: MapReaction };
  nodes: { [key: string]: MapNode };
  text_labels: { [key: string]: TextLabel };
  largest_ids: {
    reactions: number
    nodes: number
    segments: number
    text_labels: number
  };
  undo_stack: UndoStack;
  behavior: Behavior;
  draw: Draw;
  key_manager: KeyManager;
  enable_search: boolean;
  search_index: SearchIndex;
  map_name: string;
  map_id: string;
  map_description: string;
  beziers_enabled: boolean;
  has_data_on_reactions: boolean;
  has_data_on_nodes: boolean;
  imported_reaction_data: any;
  imported_metabolite_data: any;
  imported_gene_data: any;
  beziers: any;
  data_statistics: any;
  rotation_on: boolean;
  scale: any;
  private _status_timer?: NodeJS.Timeout | null | number;
  last_action: {function: keyof EscherMap, args: any[]} | null = null;
  shortestPathInfo: { hash: string; k: number; } | null = null;
  outGoingInfo: { hash: string; k: number; } | null = null;
  inGoingInfo: { hash: string; k: number; } | null = null;
  // save action that will be set from outside
  //todo: add types here

  constructor(svg: D3SVGSelection,
    css: string,
    selection: D3Selection<Element>,
    zoomContainer: ZoomContainer,
    settings: Settings,
    cobra_model: CobraModel | null,
    canvas_size_and_loc: { x: number, y: number, width: number, height: number } | null = null,
    enable_search: boolean = false,
    map_name: string = '',
    map_id: string = '',
    map_description: string = '') {
    if (canvas_size_and_loc == null) {
      const size = zoomContainer.get_size();
      canvas_size_and_loc = {
        x: -size.width,
        y: -size.height,
        width: size.width * 3,
        height: size.height * 3
      };
    }

    if (_.isUndefined(map_name) || map_name == null || map_name === '')
      map_name = 'new_map';
    else
      map_name = String(map_name);


    if (_.isUndefined(map_id) || map_id == null || map_id === '')
      map_id = utils.generate_map_id();
    else
      map_id = String(map_id);


    if (_.isUndefined(map_description) || map_description == null)
      map_description = '';
    else
      map_description = String(map_description);


    // set up the callbacks
    this.callback_manager = new CallbackManager();

    // set up the defs
    this.svg = svg;
    this.defs = utils.setup_defs(svg as D3Selection<any>, css);

    // make the canvas
    this.canvas = new Canvas(selection, canvas_size_and_loc);

    this.setup_containers(selection);
    this.sel = selection;
    this.zoomContainer = zoomContainer;

    this.settings = settings;

    // set the model AFTER loading the datasets
    this.cobra_model = cobra_model;

    this.largest_ids = {
      reactions: -1,
      nodes: -1,
      segments: -1,
      text_labels: -1
    };

    // make the undo/redo stack
    this.undo_stack = new UndoStack();

    // make a behavior object
    this.behavior = new Behavior(this, this.undo_stack);

    // draw manager
    this.draw = new Draw(this.behavior, this.settings, this);

    // make a key manager
    this.key_manager = new KeyManager();
    this.key_manager.settings = settings;
    this.key_manager.ctrlEqualsCmd = true;

    // make the search index
    this.enable_search = enable_search;
    this.search_index = new SearchIndex();

    // map properties
    this.map_name = map_name;
    this.map_id = map_id;
    this.map_description = map_description;

    // deal with the window
    const window_translate = {x: 0, y: 0};
    const window_scale = 1;

    // hide beziers
    this.beziers_enabled = false;

    // data
    this.has_data_on_reactions = false;
    this.has_data_on_nodes = false;
    this.imported_reaction_data = null;
    this.imported_metabolite_data = null;
    this.imported_gene_data = null;

    this.nodes = {};
    this.reactions = {};
    this.beziers = {};
    this.text_labels = {};

    // Update data with null to populate data-specific attributes. Also calculates
    // data stats for the first time.
    this.apply_reaction_data_to_map(null);
    this.apply_metabolite_data_to_map(null);
    this.apply_gene_data_to_map(null);

    // make the scales
    this.scale = new Scale();
    // initialize stats
    this.scale.connectToSettings(this.settings, this, this.get_data_statistics.bind(this));

    // rotation mode off
    this.rotation_on = false;
  }

  /**
   * Load a json map and add necessary fields for rendering.
   */
  static from_data(map_data: EscherMapDataType, svg: D3SVGSelection, css: string, selection: D3Selection, zoomContainer: ZoomContainer, settings: Settings,
    cobra_model: CobraModel | null, enable_search: boolean = false): EscherMap {
    const canvas = map_data[1].canvas;
    const map_name = map_data[0].map_name;
    const map_id = map_data[0].map_id;
    const map_description = (map_data[0].map_description.replace(/(\nLast Modified.*)+$/g, '') +
                           '\nLast Modified ' + new Date(Date.now()).toString());
    const map = new EscherMap(svg, css, selection, zoomContainer, settings, cobra_model,
      canvas, enable_search, map_name, map_id, map_description);

    map.reactions = map_data[1].reactions;
    map.nodes = map_data[1].nodes;
    map.text_labels = map_data[1].text_labels;

    for (const n_id in map.nodes) {
      const node = map.nodes[n_id];

      // clear all the connected segments
      node.connected_segments = [];

      //  populate the nodes search index.
      if (enable_search) {
        if (node.node_type !== 'metabolite') continue;
        map.search_index.insert('n' + n_id, {name: node.bigg_id,
          data: {type: 'metabolite',
            node_id: n_id}});
        map.search_index.insert('n_name' + n_id, {name: node.name,
          data: {type: 'metabolite',
            node_id: n_id}});
      }
    }

    // Propagate coefficients and reversibility, build the connected
    // segments, add bezier points, and populate the reaction search index.
    for (const r_id in map.reactions) {
      const reaction = map.reactions[r_id];

      // reaction search index
      if (enable_search) {
        map.search_index.insert('r' + r_id,
          {'name': reaction.bigg_id,
            'data': {type: 'reaction',
              reaction_id: r_id}});
        map.search_index.insert('r_name' + r_id,
          {'name': reaction.name,
            'data': {type: 'reaction',
              reaction_id: r_id}});
        for (const g_id in reaction.genes) {
          const gene = reaction.genes[g_id];
          map.search_index.insert('r' + r_id + '_g' + g_id,
            {'name': gene.bigg_id,
              'data': {type: 'reaction',
                reaction_id: r_id}});
          map.search_index.insert('r' + r_id + '_g_name' + g_id,
            {'name': gene.name,
              'data': {type: 'reaction',
                reaction_id: r_id}});
        }
      }

      // keep track of any bad segments
      const segments_to_delete = [];
      for (const s_id in reaction.segments) {
        const segment = reaction.segments[s_id];

        // propagate reversibility
        segment.reversibility = reaction.reversibility;

        // if there is an error with to_ or from_ nodes, remove this segment
        if (!(segment.from_node_id in map.nodes) || !(segment.to_node_id in map.nodes)) {
          console.warn('Bad node references in segment ' + s_id + '. Deleting segment.');
          segments_to_delete.push(s_id);
          continue;
        }

        const from_node = map.nodes[segment.from_node_id];
        const to_node = map.nodes[segment.to_node_id];

        // propagate coefficients
        reaction.metabolites.forEach(function(met) {
          if (met.bigg_id === from_node.bigg_id)
            segment.from_node_coefficient = met.coefficient;
          else if (met.bigg_id === to_node.bigg_id)
            segment.to_node_coefficient = met.coefficient;
        })

        // build connected segments
        ;[from_node, to_node].forEach(function(node) {
          node.connected_segments.push({segment_id: s_id,
            reaction_id: r_id});
        });

        // If the metabolite has no bezier points, then add them.
        const start = map.nodes[segment.from_node_id];
        const end = map.nodes[segment.to_node_id];
        if (start['node_type'] == 'metabolite' || end['node_type'] == 'metabolite') {
          const midpoint = utils.c_plus_c(start, utils.c_times_scalar(utils.c_minus_c(end, start)!, 0.5));
          if (segment.b1 == null) segment.b1 = midpoint;
          if (segment.b2 == null) segment.b2 = midpoint;
        }
      }
      // delete the bad segments
      segments_to_delete.forEach(function(s_id) {
        delete reaction.segments[s_id];
      });
    }

    // add text_labels to the search index
    if (enable_search) {
      for (const label_id in map.text_labels) {
        const label = map.text_labels[label_id];
        map.search_index.insert('l' + label_id, {'name': label.text,
          'data': {type: 'text_label',
            text_label_id: label_id}});
      }
    }

    // populate the beziers
    map.beziers = build.newBeziersForReactions(map.reactions);

    // get largest ids for adding new reactions, nodes, text labels, and
    // segments
    map.largest_ids.reactions = get_largest_id(map.reactions);
    map.largest_ids.nodes = get_largest_id(map.nodes);
    map.largest_ids.text_labels = get_largest_id(map.text_labels);

    let largest_segment_id = 0;
    for (const id in map.reactions) {
      largest_segment_id = get_largest_id(map.reactions[id].segments,
        largest_segment_id);
    }
    map.largest_ids.segments = largest_segment_id;

    // update data with null to populate data-specific attributes
    map.apply_reaction_data_to_map(null);
    map.apply_metabolite_data_to_map(null);
    map.apply_gene_data_to_map(null);

    return map;

    /**
     * Return the largest integer key in obj, or current_largest, whichever is
     * bigger.
     */
    function get_largest_id(obj: any, current_largest?: number) {
      if (_.isUndefined(current_largest)) current_largest = 0;
      if (_.isUndefined(obj)) return current_largest;
      return Math.max.apply(null, Object.keys(obj).map(function(x) {
        return parseInt(x);
      }).concat([current_largest]));
    }
  }

  // ---------------------------------------------------------------------
  // more setup
  // ---------------------------------------------------------------------

  setup_containers(sel: D3Selection) {
    sel.append('g')
      .attr('id', 'reactions');
    sel.append('g')
      .attr('id', 'nodes');
    sel.append('g')
      .attr('id', 'beziers');
    sel.append('g')
      .attr('id', 'text-labels');
  }

  reset_containers() {
    this.sel.select('#reactions')
      .selectAll('.reaction')
      .remove();
    this.sel.select('#nodes')
      .selectAll('.node')
      .remove();
    this.sel.select('#beziers')
      .selectAll('.bezier')
      .remove();
    this.sel.select('#text-labels')
      .selectAll('.text-label')
      .remove();
  }

  // -------------------------------------------------------------------------
  // Appearance
  // -------------------------------------------------------------------------

  /** Set the status of the map, with an optional expiration
      time. Rendering the status is taken care of by the Builder.

      Arguments
      ---------

      status: The status string.

      time: An optional time, in ms, after which the status is set to ''.

  */

  set_status(status: string | null, time?: number) {
    this.callback_manager.run('set_status', null, status);
    // clear any other timers on the status bar
    if (this._status_timer)
      clearTimeout(this._status_timer);
    this._status_timer = null;

    if (time !== undefined) {
      this._status_timer = setTimeout(() => {
        this.callback_manager.run('set_status', null, '');
      }, time);
    }
  }

  /**
   * Clear the map data
   */
  clearMapData() {
    this.reactions = {};
    this.beziers = {};
    this.nodes = {};
    this.text_labels = {};
    this.map_name = 'new_map';
    this.map_id = utils.generate_map_id();
    this.map_description = '';
  }

  has_cobra_model() {
    return (this.cobra_model !== null);
  }

  /**
   * Draw the all reactions, nodes, & text labels.
   */
  draw_everything() {
    this.draw_all_reactions(true, true); // also draw beziers
    this.draw_all_nodes(true);
    this.draw_all_text_labels();
    this.unhighlightShortestPaths();
  }

  /** Draw all reactions, and clear deleted reactions.

      Arguments
      ---------

      draw_beziers: (Boolean, default True) Whether to also draw the bezier
      control points.

      clear_deleted: (Optional, Default: true) Boolean, if true, then also
      clear deleted nodes.

  */
  draw_all_reactions(draw_beziers: boolean, clear_deleted: boolean) {
    if (_.isUndefined(draw_beziers)) draw_beziers = true;
    if (_.isUndefined(clear_deleted)) clear_deleted = true;

    // Draw all reactions.
    const reaction_ids = [];
    for (const reaction_id in this.reactions)
      reaction_ids.push(reaction_id);

    // If draw_beziers is true, just draw them all, rather than deciding
    // which ones to draw.
    this.draw_these_reactions(reaction_ids, false);
    if (draw_beziers && this.beziers_enabled)
      this.draw_all_beziers();

    // Clear all deleted reactions.
    if (clear_deleted)
      this.clear_deleted_reactions(draw_beziers);
  }

  /**
   * Draw specific reactions. Does nothing with exit selection. Use
   * clear_deleted_reactions to remove reactions from the DOM.
   * reactions_ids: An array of reaction_ids to update.
   * draw_beziers: (Boolean, default True) Whether to also draw the bezier control
   * points.
   */
  draw_these_reactions(reaction_ids: string[], draw_beziers?: boolean) {
    if (_.isUndefined(draw_beziers)) draw_beziers = true;

    // find reactions for reaction_ids
    const reaction_subset = utils.object_slice_for_ids_ref(this.reactions,
      reaction_ids);

    // function to update reactions
    const update_fn = (sel: D3Selection) => {
      return this.draw.update_reaction(sel, this.scale, this.cobra_model,
        this.nodes, this.defs,
        this.has_data_on_reactions);
    };

    // draw the reactions
    utils.draw_an_object(this.sel, '#reactions', '.reaction', reaction_subset,
      'reaction_id', this.draw.create_reaction.bind(this.draw),
      update_fn);

    if (draw_beziers) {
      // particular beziers to draw
      const bezier_ids = build.bezierIdsForReactionIds(reaction_subset);
      this.draw_these_beziers(bezier_ids);
    }
  }

  /**
   * Remove any reactions that are not in *this.reactions*.
   * draw_beziers: (Boolean, default True) Whether to also clear deleted bezier
   * control points.
   */
  clear_deleted_reactions(draw_beziers?: boolean) {
    if (_.isUndefined(draw_beziers)) draw_beziers = true;

    // Remove deleted reactions and segments
    utils.draw_an_object(
      this.sel, '#reactions', '.reaction', this.reactions, 'reaction_id', null,
      function(update_selection: D3Selection<any>) {
        // Draw segments
        utils.draw_a_nested_object(
          update_selection, '.segment-group', 'segments', 'segment_id', null,
          null, function(sel: any) { sel.remove(); }
        );
      },
      function(sel: D3Selection) {
        sel.remove();
      }
    );

    if (draw_beziers === true)
      this.clear_deleted_beziers();
  }

  /**
   * Draw all nodes, and clear deleted nodes.
   * @param clear_deleted: (Optional, Default: true) Boolean, if true, then also
   * @param clear deleted nodes.
   */
  draw_all_nodes(clear_deleted?: boolean) {
    if (clear_deleted == undefined) clear_deleted = true;

    const node_ids = [];
    for (const node_id in this.nodes)
      node_ids.push(node_id);

    this.draw_these_nodes(node_ids);

    // clear the deleted nodes
    if (clear_deleted)
      this.clear_deleted_nodes();
  }

  /** Draw specific nodes.

      Does nothing with exit selection. Use clear_deleted_nodes to remove
      nodes from the DOM.

      Arguments
      ---------

      nodes_ids: An array of node_ids to update.

  */
  draw_these_nodes(node_ids: IDType[]) {
    // find reactions for reaction_ids
    const node_subset = utils.object_slice_for_ids_ref(this.nodes, node_ids);

    // functions to create and update nodes
    const create_fn = (sel: D3Selection) => {
      return this.draw.create_node(sel,
        this.nodes,
        this.reactions);
    };
    const update_fn = (sel: D3Selection) => {
      return this.draw.update_node(sel,
        this.scale,
        this.has_data_on_nodes,
        this.behavior.selectableMousedown,
        this.behavior.selectableClick,
        this.behavior.nodeMouseover,
        this.behavior.nodeMouseout,
        this.behavior.selectableDrag,
        this.behavior.nodeLabelDrag);
    };

    // draw the nodes
    utils.draw_an_object(this.sel, '#nodes', '.node', node_subset, 'node_id',
      create_fn, update_fn);
  }

  /**
   * Remove any nodes that are not in *this.nodes*.
   */
  clear_deleted_nodes() {
    // Run remove for exit selection
    utils.draw_an_object(this.sel, '#nodes', '.node', this.nodes, 'node_id',
      null, null, function(sel: D3Selection) { sel.remove(); });
  }

  /**
   * Draw all text_labels.
   */
  draw_all_text_labels() {
    this.draw_these_text_labels(Object.keys(this.text_labels));

    // Clear all deleted text_labels
    this.clear_deleted_text_labels();
  }

  /**
   * Draw specific text_labels. Does nothing with exit selection. Use
   * clear_deleted_text_labels to remove text_labels from the DOM.
   * @param {Array} text_labels_ids - An array of text_label_ids to update.
   */
  draw_these_text_labels(text_label_ids: IDType[]) {
    // Find reactions for reaction_ids
    const text_label_subset = utils.object_slice_for_ids_ref(this.text_labels, text_label_ids);

    // Draw the text_labels
    utils.draw_an_object(this.sel, '#text-labels', '.text-label',
      text_label_subset, 'text_label_id',
      this.draw.create_text_label.bind(this.draw),
      this.draw.update_text_label.bind(this.draw));
  }

  /**
   * Remove any text_labels that are not in *this.text_labels*.
   */
  clear_deleted_text_labels() {
    utils.draw_an_object(this.sel, '#text-labels', '.text-label',
      this.text_labels, 'text_label_id', null, null,
      function(sel: D3Selection) { sel.remove(); });
  }

  /**
   * Draw all beziers, and clear deleted reactions.
   */
  draw_all_beziers() {
    const bezier_ids = [];
    for (const bezier_id in this.beziers)
      bezier_ids.push(bezier_id);

    this.draw_these_beziers(bezier_ids);

    // clear delete beziers
    this.clear_deleted_beziers();
  }

  draw_these_beziers(bezier_ids: IDType[]) {
    /** Draw specific beziers.

        Does nothing with exit selection. Use clear_deleted_beziers to remove
        beziers from the DOM.

        Arguments
        ---------

        beziers_ids: An array of bezier_ids to update.

    */
    // find reactions for reaction_ids
    const bezier_subset = utils.object_slice_for_ids_ref(this.beziers, bezier_ids);

    // function to update beziers
    const update_fn = (sel: D3Selection) => {
      return this.draw.update_bezier(sel,
        this.beziers_enabled,
        this.behavior.bezierDrag,
        this.behavior.bezierMouseover,
        this.behavior.bezierMouseout,
        this.nodes,
        this.reactions);
    };

    // draw the beziers
    utils.draw_an_object(this.sel, '#beziers', '.bezier', bezier_subset,
      'bezier_id', this.draw.create_bezier.bind(this.draw),
      update_fn);
  }

  clear_deleted_beziers() {
    /** Remove any beziers that are not in *this.beziers*.

     */
    // remove deleted
    utils.draw_an_object(this.sel, '#beziers', '.bezier', this.beziers,
      'bezier_id', null, null,
      function(sel: D3Selection) { sel.remove(); });
  }

  show_beziers() {
    this.toggle_beziers(true);
  }

  hide_beziers() {
    this.toggle_beziers(false);
  }

  toggle_beziers(on_off?: boolean) {
    if (_.isUndefined(on_off)) this.beziers_enabled = !this.beziers_enabled;
    else this.beziers_enabled = on_off;
    this.draw_all_beziers();
    this.callback_manager.run('toggle_beziers', null, this.beziers_enabled);
  }

  /**
   * Returns True if the scale has changed.
   * @param {Array} keys - (Optional) The keys in reactions to apply data to.
   */
  apply_reaction_data_to_map(data: { [key: string]: (number | string | null)[] }, keys?: any) {
    const styles = this.settings.get('reaction_styles');
    const compareStyle = this.settings.get('reaction_compare_style');
    const hasData = dataStyles.apply_reaction_data_to_reactions(
      this.reactions,
      data,
      styles,
      compareStyle,
      keys
    );
    this.has_data_on_reactions = hasData;
    this.imported_reaction_data = hasData ? data : null;

    return this.calc_data_stats('reaction');
  }

  /**
   * Returns True if the scale has changed.
   * @param {Array} keys - (Optional) The keys in nodes to apply data to.
   */
  apply_metabolite_data_to_map(data: { [key: string]: (number | string | null)[] }, keys?: any) {
    const styles = this.settings.get('metabolite_styles');
    const compare_style = this.settings.get('metabolite_compare_style');

    const has_data = dataStyles.apply_metabolite_data_to_nodes(this.nodes,
      data, styles,
      compare_style,
      keys);
    this.has_data_on_nodes = has_data;
    this.imported_metabolite_data = has_data ? data : null;

    return this.calc_data_stats('metabolite');
  }

  /**
   * Returns True if the scale has changed.
   * gene_data_obj: The gene data object, with the following style:
   * { reaction_id: { rule: 'rule_string', genes: { gene_id: value } } }
   * @param {Array} keys - (Optional) The keys in reactions to apply data to.
   */
  apply_gene_data_to_map(gene_data_obj: any, keys?: any) {
    const styles = this.settings.get('reaction_styles');
    const compare_style = this.settings.get('reaction_compare_style');
    const identifiers_on_map = this.settings.get('identifiers_on_map');
    const and_method_in_gene_reaction_rule = this.settings.get('and_method_in_gene_reaction_rule');

    const has_data = dataStyles.apply_gene_data_to_reactions(this.reactions, gene_data_obj,
      styles, identifiers_on_map,
      compare_style,
      and_method_in_gene_reaction_rule,
      keys);
    this.has_data_on_reactions = has_data;
    this.imported_gene_data = has_data ? gene_data_obj : null;

    return this.calc_data_stats('reaction');
  }

  // ------------------------------------------------
  // Data domains
  // ------------------------------------------------

  get_data_statistics() {
    return this.data_statistics;
  }

  /**
   * Returns True if the stats have changed.
   * @param {String} type - Either 'metabolite' or 'reaction'
   */
  calc_data_stats(type: string) {
    if (['reaction', 'metabolite'].indexOf(type) === -1)
      throw new Error('Bad type ' + type);


    // make the data structure
    if (!this.data_statistics) {
      this.data_statistics = {};
      this.data_statistics[type] = null;
    } else if (!(type in this.data_statistics))
      this.data_statistics[type] = null;


    // default min and max
    const vals: any[] = [];
    if (type === 'metabolite') {
      for (const node_id in this.nodes) {
        const node = this.nodes[node_id];
        // check number
        if (_.isUndefined(node.data))
          console.error('metabolite missing ');
        else if (node.data !== null)
          vals.push(node.data);
      }
    } else if (type == 'reaction') {
      for (const reaction_id in this.reactions) {
        const reaction = this.reactions[reaction_id];
        // check number
        if (_.isUndefined(reaction.data))
          console.error('reaction data missing ');
        else if (reaction.data !== null)
          vals.push(reaction.data);
      }
    }

    // If no vals, then set to null
    if (vals.length === 0) {
      const wasNull = this.data_statistics[type] == null;
      this.data_statistics[type] = null;
      if (type === 'reaction')
        this.callback_manager.run('calc_data_stats__reaction', null, !wasNull);
      else
        this.callback_manager.run('calc_data_stats__metabolite', null, !wasNull);

      return !wasNull; // true means did change which means was not null
    }

    // If there are vals and it was previously null, then create an object
    if (this.data_statistics[type] == null)
      this.data_statistics[type] = {};


    let same = true;

    // calculate these statistics
    const quartiles = utils.quartiles(vals);
    const funcs = [
      ['min', _on_array(Math.min)],
      ['max', _on_array(Math.max)],
      ['mean', utils.mean],
      ['Q1', function() { return quartiles[0]; }],
      ['median', function() { return quartiles[1]; }],
      ['Q3', function() { return quartiles[2]; }],
    ];
    funcs.forEach((ar) => {
      let new_val;
      const name = ar[0];
      if (vals.length === 0)
        new_val = null;
      else {
        const fn = ar[1] as (array: number[]) => number;
        new_val = fn(vals as number[]);
      }
      if (new_val != this.data_statistics[type][name as any])
        same = false;

      this.data_statistics[type][name as any] = new_val;
    });

    // Deal with max === min
    if (this.data_statistics[type]['min'] === this.data_statistics[type]['max'] &&
        this.data_statistics[type]['min'] !== null) {
      const min = this.data_statistics[type]['min'];
      const max = this.data_statistics[type]['max'];
      this.data_statistics[type]['min'] = min - 1 - (Math.abs(min) * 0.1);
      this.data_statistics[type]['max'] = max + 1 + (Math.abs(max) * 0.1);
    }

    if (type === 'reaction')
      this.callback_manager.run('calc_data_stats__reaction', null, !same);
    else
      this.callback_manager.run('calc_data_stats__metabolite', null, !same);

    return !same;
  }

  // ---------------------------------------------------------------------
  // Node interaction
  // ---------------------------------------------------------------------

  get_coords_for_node(node_id: IDType) {
    const node = this.nodes[node_id];
    const coords = node ? {x: node.x, y: node.y} : null;
    return coords;
  }

  get_selected_node_ids() {
    const selected_node_ids: IDType[] = [];
    this.sel.select('#nodes')
      .selectAll('.selected')
      .each((d: any) => { selected_node_ids.push(d.node_id); });
    return selected_node_ids;
  }

  getSelectedNodes() {
    const selected_nodes: {[key: IDType]: MapNode} = {};
    this.sel.select('#nodes')
      .selectAll('.selected')
      .each((d: any) => {
        selected_nodes[d.node_id] = this.nodes[d.node_id];
      });
    return selected_nodes;
  }

  get_selected_text_label_ids() {
    const selected_text_label_ids: IDType[] = [];
    this.sel.select('#text-labels')
      .selectAll('.selected')
      .each(function(d: any) { selected_text_label_ids.push(d.text_label_id); });
    return selected_text_label_ids;
  }

  get_selected_text_labels() {
    const selected_text_labels: Record<IDType, TextLabel> = {};
    this.sel.select('#text-labels')
      .selectAll('.selected')
      .each((d: any) => {
        selected_text_labels[d.text_label_id] = this.text_labels[d.text_label_id];
      });
    return selected_text_labels;
  }

  select_all() {
    /** Select all nodes and text labels.

     */
    this.sel.selectAll('#nodes,#text-labels')
      .selectAll('.node,.text-label')
      .classed('selected', true);
  }

  select_none() {
    /** Deselect all nodes and text labels.

     */
    this.sel.selectAll('.selected')
      .classed('selected', false);
  }

  invert_selection() {
    /** Invert selection of nodes and text labels.

     */
    const selection = this.sel.selectAll('#nodes,#text-labels')
      .selectAll('.node,.text-label');
    selection.classed('selected', function() {
      return !d3Select(this).classed('selected');
    });
  }

  select_metabolite_with_id(node_id: IDType) {
    /** Select a metabolite with the given id, and turn off the reaction
        target.

    */
    // deselect all text labels
    this.deselect_text_labels();

    const node_selection = this.sel.select('#nodes').selectAll('.node');
    let coords;
    let selected_node;
    node_selection.classed('selected', (d: any) => {
      const selected = String(d.node_id) == String(node_id);
      if (selected) {
        selected_node = d;
        coords = {x: d.x, y: d.y};
      }
      return selected;
    });
    this.sel.selectAll('.start-reaction-target').style('visibility', 'hidden');
    this.callback_manager.run('select_metabolite_with_id', null, selected_node, coords);
  }

  select_selectable(node: BaseType, _d: any, shift_key_on?: boolean) {
    /** Select a metabolite or text label, and manage the shift key. */
    shift_key_on = _.isUndefined(shift_key_on) ? false : shift_key_on;
    const classable_selection = this.sel.selectAll('#nodes,#text-labels')
      .selectAll('.node,.text-label');
    let classable_node;
    if (d3Select(node).attr('class').indexOf('text-label') == -1) {
      // node
      //@ts-ignore
      classable_node = node!.parentNode;
    } else {
      // text-label
      classable_node = node;
    }
    // toggle selection
    if (shift_key_on) {
      // toggle this node
      d3Select(classable_node)
        .classed('selected', !d3Select(classable_node).classed('selected'));
    } else {
      // unselect all other nodes, and select this one
      classable_selection.classed('selected', false);
      d3Select(classable_node).classed('selected', true);
    }
    // run the select_metabolite callback
    const selected_nodes = this.sel.select('#nodes').selectAll('.selected');
    let node_count = 0;
    let coords;
    let selected_node;
    selected_nodes.each(function(d: any) {
      selected_node = d;
      coords = {x: d.x, y: d.y};
      node_count++;
    });
    this.callback_manager.run('select_selectable', null, node_count, selected_node, coords);
  }

  /**
   * Unselect all but one selected node, and return the node. If no nodes are
   * selected, return null.
   */
  select_single_node(): MapNode | null {
    let out = null;
    const node_selection = this.sel.select('#nodes').selectAll('.selected');
    node_selection.classed('selected', function(d, i) {
      if (i === 0) {
        out = d;
        return true;
      } else
        return false;
    });
    return out;
  }

  deselect_nodes() {
    const node_selection = this.sel.select('#nodes').selectAll('.node');
    node_selection.classed('selected', false);
    this.callback_manager.run('deselect_nodes');
  }

  select_text_label(sel: D3Selection, d: any) {
    // deselect all nodes
    this.deselect_nodes();
    // Find the new selection. Ignore shift key and only allow single selection
    // for now.
    const text_label_selection = this.sel.select('#text-labels').selectAll('.text-label');
    text_label_selection.classed('selected', function(p) { return d === p; });
    const selected_text_labels = this.sel.select('#text-labels').selectAll('.selected');
    let coords;
    selected_text_labels.each(function(d: any) {
      coords = {x: d.x, y: d.y};
    });
    this.callback_manager.run('select_text_label');
  }

  deselect_text_labels() {
    const text_label_selection = this.sel.select('#text-labels').selectAll('.text-label');
    text_label_selection.classed('selected', false);
  }


  /**
   * Align selected nodes and/or reactions vertically. Undoable.
   */
  align_vertical() {
    return this._align(false);
  }

  /**
   * Align selected nodes and/or reactions horizontally. Undoable
   */
  align_horizontal() {
    return this._align(true);
  }

  /**
   * Generic function for aligning nodes.
   * @param {Boolean} isHorizontal - If true, align horizontal; else vertical.
   */
  _align(isHorizontal?: boolean) {
    const selected = this.getSelectedNodes();
    // Get markers and primary nodes
    const markersAndPrimary = _.pick(
      selected,
      (node) => node.node_type !== 'metabolite' || !!node.node_is_primary
    );

    // 1. generally operates on markers and primary nodes, bringing along
    // unconnected secondary nodes and beziers
    // 2. if only selected secondary nodes, then align those
    const alignByPrimary = Object.keys(markersAndPrimary).length > 0;
    const toAlign = alignByPrimary ? markersAndPrimary : selected;
    const keysToAlign = Object.keys(toAlign);

    // Get new x location
    const mean = keysToAlign.reduce((accum, val) => {
      return accum + (isHorizontal ? toAlign[val]!.y : toAlign[val]!.x);
    }, 0) / keysToAlign.length;

    // Align. Remember displacements for undo/redo.
    const displacements = _.pairs(toAlign).map(([nodeId, node]) => ({
      nodeId,
      displacement: isHorizontal ? {x: 0, y: mean - node!.y} : {x: mean - node!.x, y: 0}
    }));
    const bezierDisplacements: BezierDisplacement[] = [];
    const movedSecondaryNodes: { [key: string]: boolean } = {};

    // Align unconnected secondary metabolites
    if (alignByPrimary) {
      _.mapObject(toAlign, (node, nodeId) => {
        if (!node) return;
        node.connected_segments.map((segmentLink) => {
          // get each connected node
          const segmentId = segmentLink.segment_id;
          const reactionId = segmentLink.reaction_id;
          const segment = this.reactions[reactionId].segments[segmentId];
          const isToNode = segment.to_node_id === node.node_id;
          const otherNodeId = isToNode ? segment.from_node_id : segment.to_node_id;
          const otherNode = this.nodes[otherNodeId];
          if (!otherNode) return;
          const bez = isToNode ? 'b2' : 'b1';

          // align this side bezier if the other node is selected (and that node
          // will handle its bezier)
          if (otherNode.node_id! in selected && segment[bez]) {
            const bezierId = build.bezierIdForSegmentId(segmentId, bez);
            bezierDisplacements.push({
              reactionId,
              segmentId,
              bez,
              bezierId,
              displacement: (isHorizontal ?
                {x: 0, y: node.y - segment[bez].y} :
                {x: node.x - segment[bez].x, y: 0})
            });

            // If the node is secondary and unconnected, then move it along with
            // the current node.
            if (otherNode.node_type === 'metabolite' &&
                !otherNode.node_is_primary &&
                !(otherNodeId in movedSecondaryNodes)) {
              // If all the connected segments are connected to selected nodes, then move it
              const connected = otherNode.connected_segments.filter((segmentLink) => {
                const segment = this.reactions[reactionId].segments[segmentId];
                const isToNode = segment.to_node_id === otherNode.node_id;
                return isToNode ? segment.from_node_id in selected : segment.to_node_id in selected;
              });
              if (otherNode.connected_segments.length <= connected.length) {
                // then move it with the same displacement as the parent
                displacements.push({
                  nodeId: otherNodeId,
                  displacement: isHorizontal ? {x: 0, y: mean - node.y} : {x: mean - node.x, y: 0}
                });
                // remember not to move this again
                movedSecondaryNodes[otherNodeId] = true;
              }
            }
          }
        });
      });
    }

    // reusable movement function for aligning and undo/redo
    const _moveNodes = (disps: { nodeId: string; displacement: { x: number; y: number } }[], bezDisps: { reactionId: string; segmentId: string; bez: string; bezierId: string; displacement: { x: number; y: number } }[]) => {
      let reactionIds: string[] = [];
      disps.map((d) => {
        // TODO abstract this approach in a function because the alternative
        // (saving the node itself) causes bugs)
        const node = this.nodes[d.nodeId];
        const updated = build.moveNodeAndDependents(
          node,
          d.nodeId,
          this.reactions,
          this.beziers,
          d.displacement
        );
        reactionIds = utils.uniqueConcat([reactionIds, updated.reaction_ids]);
      });
      // move beziers
      bezDisps.map((d) => {
        const segment = this.reactions[d.reactionId].segments[d.segmentId];
        const bez = d.bez as keyof typeof segment;
        (segment[bez] as Coord) = utils.c_plus_c(segment[bez] as Coord, d.displacement as Coord)!;
        this.beziers[d.bezierId].x = segment[bez].x;
        this.beziers[d.bezierId].y = segment[bez].y;
      });

      this.draw_these_nodes(disps.map((d) => d.nodeId));
      this.draw_these_reactions(reactionIds, true); // and beziers
    };

    // undo /redo
    this.undo_stack.push(
      // undo
      () => {
        const reverse = (disps: ({ nodeId: string; displacement: { x: number; y: number } } & any)[]) => disps.map((d) => ({
          ...d,
          displacement: {x: -d.displacement.x, y: -d.displacement.y}
        }));
        _moveNodes(reverse(displacements), reverse(bezierDisplacements));
      },
      // redo
      () => {
        _moveNodes(displacements, bezierDisplacements);
      }
    ).do(); // do the first time

    // finish
    this.set_status(alignByPrimary ? 'Aligned reactions' : 'Aligned nodes', 3000);
  }


  // ---------------------------------------------------------------------
  // Delete
  // ---------------------------------------------------------------------

  /**
   * Delete the selected nodes and associated segments and reactions, and selected
   * labels. Undoable.
   */
  delete_selected() {
    const selected_nodes = this.getSelectedNodes();
    const selected_text_labels = this.get_selected_text_labels();
    if (Object.keys(selected_nodes).length >= 1 ||
        Object.keys(selected_text_labels).length >= 1)
      this.delete_selectable(selected_nodes, selected_text_labels, true);
  }

  /**
   * Delete the nodes and associated segments and reactions. Undoable.
   * selected_nodes: An object that is a subset of map.nodes.
   * selected_text_labels: An object that is a subset of map.text_labels.
   * should_draw: A boolean argument to determine whether to draw the changes to
   * the map.
   */
  delete_selectable(selected_nodes: { [key: string]: MapNode }, selected_text_labels: { [key: string]: TextLabel }, should_draw: boolean) {
    const out = this.segments_and_reactions_for_nodes(selected_nodes);
    let segment_objs_w_segments = out.segment_objs_w_segments; // TODO repeated values here
    let reactions = out.reactions;

    // copy nodes to undelete
    const saved_nodes = utils.clone(selected_nodes);
    const saved_segment_objs_w_segments = utils.clone(segment_objs_w_segments);
    const saved_reactions = utils.clone(reactions);
    const saved_text_labels = utils.clone(selected_text_labels);
    const delete_and_draw = (nodes: any, reactions: any, segment_objs: any,
      selected_text_labels: any) => {
      // delete nodes, segments, and reactions with no segments
      this.delete_node_data(Object.keys(selected_nodes));
      this.delete_segment_data(segment_objs); // also deletes beziers
      this.delete_reaction_data(Object.keys(reactions));
      this.delete_text_label_data(Object.keys(selected_text_labels));

      // apply the reaction and node data
      let changed_r_scale = false;
      let changed_m_scale = false;
      if (this.has_data_on_reactions)
        changed_r_scale = this.calc_data_stats('reaction');

      if (this.has_data_on_nodes)
        changed_m_scale = this.calc_data_stats('metabolite');


      // redraw
      if (should_draw) {
        if (changed_r_scale)
          this.draw_all_reactions(true, true);
        else
          this.clear_deleted_reactions(); // also clears segments and beziers
        if (changed_m_scale)
          this.draw_all_nodes(true);
        else
          this.clear_deleted_nodes();
        this.clear_deleted_text_labels();
      }
    };

    // delete
    delete_and_draw(selected_nodes, reactions, segment_objs_w_segments,
      selected_text_labels);

    // add to undo/redo stack
    this.undo_stack.push(() => {
      // undo
      // redraw the saved nodes, reactions, and segments

      this.extend_nodes(saved_nodes);
      this.extend_reactions(saved_reactions);
      const reaction_ids_to_draw = Object.keys(saved_reactions);
      for (const segment_id in saved_segment_objs_w_segments) {
        const segment_obj = saved_segment_objs_w_segments[segment_id];

        const segment = segment_obj.segment;
        this.reactions[segment_obj.reaction_id]
          .segments[segment_obj.segment_id] = segment;

        // updated connected nodes
        const to_from = [segment.from_node_id, segment.to_node_id];
        to_from.forEach((node_id) => {
          // not necessary for the deleted nodes
          if (node_id in saved_nodes) return;
          const node = this.nodes[node_id];
          node.connected_segments.push({reaction_id: segment_obj.reaction_id,
            segment_id: segment_obj.segment_id});
        });

        // extend the beziers
        const seg_id = segment_obj.segment_id;
        const r_id = segment_obj.reaction_id;
        const seg_o: any = {};
        seg_o[seg_id] = segment_obj.segment;
        utils.extend(this.beziers, build.newBeziersForSegments(seg_o, r_id));

        if (reaction_ids_to_draw.indexOf(segment_obj.reaction_id) === -1)
          reaction_ids_to_draw.push(segment_obj.reaction_id);
      }

      // Apply the reaction and node data. If the scale changes, redraw
      // everything.
      if (this.has_data_on_reactions) {
        const scale_changed = this.calc_data_stats('reaction');
        if (scale_changed) this.draw_all_reactions(true, false);
        else this.draw_these_reactions(reaction_ids_to_draw);
      } else
        if (should_draw) this.draw_these_reactions(reaction_ids_to_draw);

      if (this.has_data_on_nodes) {
        const scale_changed = this.calc_data_stats('metabolite');
        if (should_draw) {
          if (scale_changed) this.draw_all_nodes(false);
          else this.draw_these_nodes(Object.keys(saved_nodes));
        }
      } else
        if (should_draw) this.draw_these_nodes(Object.keys(saved_nodes));


      // redraw the saved text_labels
      utils.extend(this.text_labels, saved_text_labels);
      if (should_draw) this.draw_these_text_labels(Object.keys(saved_text_labels));
      // copy text_labels to re-delete
      selected_text_labels = utils.clone(saved_text_labels);

      // copy nodes to re-delete
      selected_nodes = utils.clone(saved_nodes);
      segment_objs_w_segments = utils.clone(saved_segment_objs_w_segments);
      reactions = utils.clone(saved_reactions);
    }, function() {
      // redo
      // clone the nodes and reactions, to redo this action later
      delete_and_draw(selected_nodes, reactions, segment_objs_w_segments,
        selected_text_labels);
    });
  }

  /**
   * Delete nodes, and remove from search index.
   */
  delete_node_data(nodeIds: IDType[]) {
    nodeIds.forEach((nodeId) => {
      if (this.enable_search && this.nodes[nodeId].node_type === 'metabolite') {
        const found = (this.search_index.remove('n' + nodeId) &&
                     this.search_index.remove('n_name' + nodeId));
        if (!found)
          console.warn('Could not find deleted metabolite in search index');
      }
      delete this.nodes[nodeId];
    });
  }

  /**
   * Delete segments, update connected_segments in nodes, and delete bezier
   * points.
   * @param {Object} segment_objs - Object with values like
   *                                { reaction_id: '123', segment_id: '456' }
   */
  delete_segment_data(segment_objs: { reaction_id: string; segment_id: string }[]) {
    for (const segment_id in segment_objs) {
      const segment_obj = segment_objs[segment_id];
      const reaction = this.reactions[segment_obj.reaction_id];

      // segment already deleted
      if (!(segment_obj.segment_id in reaction.segments)) return;

      const segment = reaction.segments[segment_obj.segment_id]
      // updated connected nodes
      ;[segment.from_node_id, segment.to_node_id].forEach((node_id) => {
        if (!(node_id in this.nodes)) return;
        const node = this.nodes[node_id];
        node.connected_segments = node.connected_segments.filter(function(so) {
          return so.segment_id != segment_obj.segment_id;
        });
      })

      // remove beziers
      ;['b1', 'b2'].forEach((bez) => {
        const bez_id = build.bezierIdForSegmentId(segment_obj.segment_id, bez);
        delete this.beziers[bez_id];
      });

      delete reaction.segments[segment_obj.segment_id];
    }
  }

  /**
   * Delete reactions, segments, and beziers, and remove reaction from search
   * index.
   */
  delete_reaction_data(reaction_ids: IDType[]) {
    reaction_ids.forEach((reaction_id) => {
      // remove beziers
      const reaction = this.reactions[reaction_id];
      for (const segment_id in reaction.segments) {
        ;['b1', 'b2'].forEach((bez) => {
          const bez_id = build.bezierIdForSegmentId(segment_id, bez);
          delete this.beziers[bez_id];
        });
      }
      // delete reaction
      delete this.reactions[reaction_id];
      // remove from search index
      const found = (this.search_index.remove('r' + reaction_id) &&
                   this.search_index.remove('r_name' + reaction_id));
      if (!found) {
        console.warn('Could not find deleted reaction ' +
                     reaction_id + ' in search index');
      }
      for (const g_id in reaction.genes) {
        const found = (this.search_index.remove('r' + reaction_id + '_g' + g_id) &&
                     this.search_index.remove('r' + reaction_id + '_g_name' + g_id));
        if (!found) {
          console.warn('Could not find deleted gene ' +
                       g_id + ' in search index');
        }
      }
    });
  }

  /**
   * Delete text labels for an array of IDs
   */
  delete_text_label_data(text_label_ids: IDType[]) {
    text_label_ids.forEach((text_label_id) => {
      // delete label
      delete this.text_labels[text_label_id];
      // remove from search index
      const found = this.search_index.remove('l' + text_label_id);
      if (!found)
        console.warn('Could not find deleted text label in search index');
    });
  }

  // ---------------------------------------------------------------------
  // Building
  // ---------------------------------------------------------------------

  _extend_and_draw_metabolite(new_nodes: { [key: string]: MapNode }, selected_node_id: string) {
    this.extend_nodes(new_nodes);
    const keys = [selected_node_id];
    if (this.has_data_on_nodes) {
      if (this.imported_metabolite_data == null)
        throw new Error('imported_metabolite_data should not be null');

      const scale_changed = this.apply_metabolite_data_to_map(this.imported_metabolite_data,
        keys);
      if (scale_changed)
        this.draw_all_nodes(false);
      else
        this.draw_these_nodes(keys);
    } else
      this.draw_these_nodes(keys);
  }

  /**
   * Draw a reaction on a blank canvas.
   * @param {String} starting_reaction - bigg_id for a reaction to draw.
   * @param {Coords} coords - coordinates to start drawing
   */
  new_reaction_from_scratch(starting_reaction: string, coords: { x: number; y: number }, direction: number) {
    // If there is no cobra model, error
    if (!this.cobra_model) {
      console.error('No CobraModel. Cannot build new reaction');
      return;
    }

    // Set reaction coordinates and angle. Be sure to clone the reaction.
    const cobra_reaction = utils.clone(this.cobra_model.reactions[starting_reaction]);

    // check for empty reactions
    if (_.size(cobra_reaction.metabolites) === 0)
      throw Error('No metabolites in reaction ' + cobra_reaction.bigg_id);


    // create the first node
    const reactant_ids = _.map(cobra_reaction.metabolites,
      (coeff: number, met_id: string) => [coeff, met_id])
      .filter((x) => Number(x[0]) < 0) // coeff < 0
      .map((x) => x[1]); // metabolite id
    // get the first reactant or else the first product
    const metaboliteId = reactant_ids.length > 0 ?
      reactant_ids[0] :
      Object.keys(cobra_reaction.metabolites)[0];
    const metabolite = this.cobra_model.metabolites[metaboliteId];
    const selected_node_id = String(++this.largest_ids.nodes);
    const label_d = build.getMetLabelLoc(utils.to_radians(direction), 0, 1,
      true, metaboliteId as IDType, 1);
    const selected_node: MapNode = {
      connected_segments: [],
      x: coords.x,
      y: coords.y,
      node_is_primary: true,
      label_x: coords.x + label_d.x,
      label_y: coords.y + label_d.y,
      name: metabolite.name,
      //@ts-ignore
      bigg_id: metaboliteId,
      node_type: 'metabolite'
    };
    let new_nodes: { [key: string]: MapNode } = {};
    new_nodes[selected_node_id] = selected_node;

    // draw
    this._extend_and_draw_metabolite(new_nodes, selected_node_id);

    // clone the nodes and reactions, to redo this action later
    const saved_nodes = utils.clone(new_nodes);

    // draw the reaction
    const out = this.new_reaction_for_metabolite(starting_reaction,
      selected_node_id,
      direction, false);
    const reaction_redo = out.redo;
    const reaction_undo = out.undo;

    // add to undo/redo stack
    this.undo_stack.push(() => {
      // Undo. First undo the reaction.
      reaction_undo();
      // Get the nodes to delete
      this.delete_node_data(Object.keys(new_nodes));
      // Save the nodes and reactions again, for redo
      new_nodes = utils.clone(saved_nodes);
      // Draw
      this.clear_deleted_nodes();
      // Deselect
      this.deselect_nodes();
    }, () => {
      // Redo. Clone the nodes and reactions, to redo this action later.
      this._extend_and_draw_metabolite(new_nodes, selected_node_id);
      // Now redo the reaction
      reaction_redo();
    });

    return;
  }

  /**
   * Add new nodes to data and search index.
   */
  extend_nodes(new_nodes: { [key: string]: MapNode }) {
    if (this.enable_search) {
      for (const node_id in new_nodes) {
        const node = new_nodes[node_id];
        if (node.node_type != 'metabolite')
          continue;
        this.search_index.insert('n' + node_id,
          {'name': node.bigg_id,
            'data': {type: 'metabolite',
              node_id: node_id}});
        this.search_index.insert('n_name' + node_id,
          {'name': node.name,
            'data': {type: 'metabolite',
              node_id: node_id}});
      }
    }
    utils.extend(this.nodes, new_nodes);
  }

  /**
   * Add new reactions to data and search index.
   */
  extend_reactions(new_reactions: { [key: string]: MapReaction }) {
    if (this.enable_search) {
      for (const r_id in new_reactions) {
        const reaction = new_reactions[r_id];
        this.search_index.insert('r' + r_id, {'name': reaction.bigg_id,
          'data': {type: 'reaction',
            reaction_id: r_id}});
        this.search_index.insert('r_name' + r_id, {'name': reaction.name,
          'data': {type: 'reaction',
            reaction_id: r_id}});
        for (const g_id in reaction.genes) {
          const gene = reaction.genes[g_id];
          this.search_index.insert('r' + r_id + '_g' + g_id,
            {'name': gene.bigg_id,
              'data': {type: 'reaction',
                reaction_id: r_id}});
          this.search_index.insert('r' + r_id + '_g_name' + g_id,
            {'name': gene.name,
              'data': {type: 'reaction',
                reaction_id: r_id}});
        }
      }
    }
    utils.extend(this.reactions, new_reactions);
  }


  _extend_and_draw_reaction(new_nodes: { [key: string]: MapNode }, new_reactions: { [key: string]: MapReaction }, new_beziers: { [key: string]: any },
    selected_node_id: string) {
    this.extend_reactions(new_reactions);
    utils.extend(this.beziers, new_beziers);
    // Remove the selected node so it can be updated
    this.delete_node_data([selected_node_id]); // TODO this is a hack. fix
    this.extend_nodes(new_nodes);

    // Apply the reaction and node data to the scales. If the scale changes,
    // redraw everything.
    const keys = Object.keys(new_reactions);
    if (this.has_data_on_reactions) {
      let scale_changed = false;
      if (this.imported_reaction_data) {
        scale_changed = this.apply_reaction_data_to_map(this.imported_reaction_data,
          keys);
      } else if (this.imported_gene_data)
        scale_changed = this.apply_gene_data_to_map(this.imported_gene_data, keys);
      else {
        throw new Error('imported_gene_data or imported_reaction_data should ' +
                        'not be null');
      }
      if (scale_changed)
        this.draw_all_reactions(true, false);
      else
        this.draw_these_reactions(keys);
    } else
      this.draw_these_reactions(keys);


    const node_keys = Object.keys(new_nodes);
    if (this.has_data_on_nodes) {
      if (this.imported_metabolite_data == null)
        throw new Error('imported_metabolite_data should not be null');

      const scale_changed = this.apply_metabolite_data_to_map(this.imported_metabolite_data,
        node_keys);
      if (scale_changed)
        this.draw_all_nodes(false);
      else
        this.draw_these_nodes(node_keys);
    } else
      this.draw_these_nodes(node_keys);


    // select new primary metabolite
    for (const node_id in new_nodes) {
      const node = new_nodes[node_id];
      if (node.node_is_primary && node_id != selected_node_id) {
        this.select_metabolite_with_id(node_id);
        const new_coords = {x: node.x, y: node.y};
        if (this.zoomContainer)
          this.zoomContainer.translateOffScreen(new_coords);
      }
    }
  }

  /**
   * Build a new reaction starting with selected_met. Undoable.
   * @param {String} reaction_bigg_id - The BiGG ID of the reaction to draw.
   * @param {String} selected_node_id - The ID of the node to begin drawing with.
   * @param {Number} direction - The direction to draw in.
   * @param {Boolean} [apply_undo_redo=true] - If true, then add to the undo
   * stack. Otherwise, just return the undo and redo functions.
   * @return An object of undo and redo functions:
   *   { undo: undo_function, redo: redo_function }
   */
  new_reaction_for_metabolite(reaction_bigg_id: IDType, selected_node_id: IDType,
    direction: number, apply_undo_redo?: boolean) {
    // default args
    if (apply_undo_redo === undefined) apply_undo_redo = true;

    // get the metabolite node
    const selected_node = this.nodes[selected_node_id];

    // Set reaction coordinates and angle. Be sure to copy the reaction
    // recursively.
    const cobra_reaction = this.cobra_model!.reactions[reaction_bigg_id];

    // build the new reaction
    const out = build.newReaction(
      reaction_bigg_id,
      cobra_reaction,
      this.cobra_model!.metabolites,
      selected_node_id,
      utils.clone(selected_node),
      this.largest_ids,
      this.settings.get('cofactors'),
      direction
    );
    let new_nodes = out.new_nodes;
    let new_reactions = out.new_reactions;
    let new_beziers = out.new_beziers;

    // Draw
    this._extend_and_draw_reaction(new_nodes, new_reactions,
      new_beziers, selected_node_id);

    // Clone the nodes and reactions, to redo this action later
    const saved_nodes = utils.clone(new_nodes);
    const saved_reactions = utils.clone(new_reactions);
    const saved_beziers = utils.clone(new_beziers);

    // Add to undo/redo stack
    const undo_fn = () => {
      // Undo. Get the nodes to delete.
      delete new_nodes[selected_node_id];
      this.delete_node_data(Object.keys(new_nodes));
      this.delete_reaction_data(Object.keys(new_reactions)); // also deletes beziers
      this.select_metabolite_with_id(selected_node_id);
      // Save the nodes and reactions again, for redo
      new_nodes = utils.clone(saved_nodes);
      new_reactions = utils.clone(saved_reactions);
      new_beziers = utils.clone(saved_beziers);
      // Draw
      if (this.has_data_on_reactions) {
        const scale_changed = this.calc_data_stats('reaction');
        if (scale_changed)
          this.draw_all_reactions(true, true);
        else {
          // Also clears segments and beziers
          this.clear_deleted_reactions(true);
        }
      } else {
        // Also clears segments and beziers
        this.clear_deleted_reactions(true);
      }
      if (this.has_data_on_nodes) {
        const scaleChanged = this.calc_data_stats('metabolite');
        if (scaleChanged)
          this.draw_all_nodes(true);
        else
          this.clear_deleted_nodes();
      } else
        this.clear_deleted_nodes();
    };
    const redo_fn = () => {
      // redo
      // clone the nodes and reactions, to redo this action later
      this._extend_and_draw_reaction(new_nodes, new_reactions,
        new_beziers, selected_node_id);
    };

    if (apply_undo_redo)
      this.undo_stack.push(undo_fn, redo_fn);


    return {undo: undo_fn, redo: redo_fn};
  }

  cycle_primary_node() {
    const selected_nodes = this.getSelectedNodes();
    if (_.isEmpty(selected_nodes)) return;
    // Get the first node
    const node_id = Object.keys(selected_nodes)[0];
    const node = selected_nodes[node_id];
    const reactions = this.reactions;
    const nodes = this.nodes;
    // make the other reactants or products secondary
    // 1. Get the connected anchor nodes for the node
    const connected_anchor_ids: IDType[] = [];
    let reactions_to_draw: IDType[] = [];
    nodes[node_id].connected_segments.forEach(function(segment_info) {
      reactions_to_draw = [segment_info.reaction_id];
      let segment: MapSegment;
      try {
        segment = reactions[segment_info.reaction_id].segments[segment_info.segment_id];
        if (segment === undefined) throw new Error('undefined segment');
      } catch (e: any) {
        console.warn('Could not find connected segment ' + segment_info.segment_id);
        return;
      }
      connected_anchor_ids.push(segment.from_node_id == node_id ?
        segment.to_node_id : segment.from_node_id);
    });
    // can only be connected to one anchor
    if (connected_anchor_ids.length != 1) {
      console.error('Only connected nodes with a single reaction can be selected');
      return;
    }
    const connected_anchor_id = connected_anchor_ids[0];
    // 2. find nodes connected to the anchor that are metabolites
    const related_node_ids = [node_id];
    nodes[connected_anchor_id].connected_segments.forEach(function(segment_info) { // deterministic order
      let segment: MapSegment;
      try {
        segment = reactions[segment_info.reaction_id].segments[segment_info.segment_id];
        if (segment === undefined) throw new Error('undefined segment');
      } catch (e) {
        console.warn('Could not find connected segment ' + segment_info.segment_id);
        return;
      }
      const conn_met_id = segment.from_node_id == connected_anchor_id ? segment.to_node_id : segment.from_node_id;
      const conn_node = nodes[conn_met_id];
      if (conn_node.node_type == 'metabolite' && conn_met_id != node_id)
        related_node_ids.push(String(conn_met_id));
    });
    // 3. make sure they only have 1 reaction connection, and check if
    // they match the other selected nodes
    for (let i = 0; i < related_node_ids.length; i++) {
      if (nodes[related_node_ids[i]].connected_segments.length > 1) {
        console.error('Only connected nodes with a single reaction can be selected');
        return;
      }
    }
    for (const a_selected_node_id in selected_nodes) {
      if (a_selected_node_id != node_id && related_node_ids.indexOf(a_selected_node_id) == -1) {
        console.warn('Selected nodes are not on the same reaction');
        return;
      }
    }
    // 4. change the primary node, and change coords, label coords, and beziers
    const nodes_to_draw: IDType[] = [];
    let last_i = related_node_ids.length - 1;
    const last_node = nodes[related_node_ids[last_i]];
    let last_is_primary = last_node.node_is_primary;
    let last_coords = {x: last_node.x, y: last_node.y,
      label_x: last_node.label_x, label_y: last_node.label_y};
    if (last_node.connected_segments.length > 1)
      console.warn('Too many connected segments for node ' + last_node.node_id);
    const last_segment_info = last_node.connected_segments[0]; // guaranteed above to have only one
    let last_segment: MapSegment;
    try {
      last_segment = reactions[last_segment_info.reaction_id].segments[last_segment_info.segment_id];
      if (last_segment === undefined) throw new Error('undefined segment');
    } catch (e) {
      console.error('Could not find connected segment ' + last_segment_info.segment_id);
      return;
    }
    let last_bezier = {b1: last_segment.b1, b2: last_segment.b2};
    let primary_node_id: IDType | undefined;
    related_node_ids.forEach(function(related_node_id) {
      const node = nodes[related_node_id];
      const this_is_primary = node.node_is_primary;
      const these_coords = {x: node.x, y: node.y,
        label_x: node.label_x, label_y: node.label_y};
      const this_segment_info = node.connected_segments[0];
      const this_segment = reactions[this_segment_info.reaction_id].segments[this_segment_info.segment_id];
      const this_bezier = {b1: this_segment.b1, b2: this_segment.b2};
      node.node_is_primary = last_is_primary;
      node.x = last_coords.x; node.y = last_coords.y;
      node.label_x = last_coords.label_x; node.label_y = last_coords.label_y;
      this_segment.b1 = last_bezier.b1; this_segment.b2 = last_bezier.b2;
      last_is_primary = this_is_primary;
      last_coords = these_coords;
      last_bezier = this_bezier;
      if (node.node_is_primary) primary_node_id = related_node_id;
      nodes_to_draw.push(related_node_id);
    });
    // 5. cycle the connected_segments array so the next time, it cycles differently
    const old_connected_segments = nodes[connected_anchor_id].connected_segments;
    last_i = old_connected_segments.length - 1;
    const new_connected_segments = [old_connected_segments[last_i]];
    old_connected_segments.forEach(function(segment, i) {
      if (last_i == i) return;
      new_connected_segments.push(segment);
    });
    nodes[connected_anchor_id].connected_segments = new_connected_segments;
    // 6. draw the nodes
    this.draw_these_nodes(nodes_to_draw);
    this.draw_these_reactions(reactions_to_draw);
    // 7. select the primary node
    if (primary_node_id)
      this.select_metabolite_with_id(primary_node_id);

    return;
  }

  /**
   * Toggle the primary/secondary status of each selected node. Undoable.
   */
  toggle_selected_node_primary() {
    const selected_node_ids = this.get_selected_node_ids();
    const go = (ids: IDType[]) => {
      const nodes_to_draw: { [key: string]: MapNode } = {};
      const hide_secondary_metabolites = this.settings.get('hide_secondary_metabolites');
      ids.forEach((id) => {
        if (!(id in this.nodes)) {
          console.warn('Could not find node: ' + id);
          return;
        }
        const node = this.nodes[id];
        if (node.node_type == 'metabolite') {
          node.node_is_primary = !node.node_is_primary;
          nodes_to_draw[id] = node;
        }
      });
      // draw the nodes
      this.draw_these_nodes(Object.keys(nodes_to_draw));
      // draw associated reactions
      if (hide_secondary_metabolites) {
        const out = this.segments_and_reactions_for_nodes(nodes_to_draw);
        const reaction_ids_to_draw_o: { [key: string]: boolean } = {};
        for (const id in out.segment_objs_w_segments) {
          const r_id = out.segment_objs_w_segments[id].reaction_id;
          reaction_ids_to_draw_o[r_id] = true;
        }
        this.draw_these_reactions(Object.keys(reaction_ids_to_draw_o));
      }
    };

    // go
    go(selected_node_ids);

    // add to the undo stack
    this.undo_stack.push(function() {
      go(selected_node_ids);
    }, function() {
      go(selected_node_ids);
    });
  }

  segments_and_reactions_for_nodes(nodes: { [key: string]: MapNode }) {
    /** Get segments and reactions that should be deleted with node deletions

     */
    const segment_objs_w_segments: { [key: string]: { reaction_id: string; segment_id: string; segment: MapSegment } } = {};
    const these_reactions: { [key: string]: MapReaction } = {};
    const segment_ids_for_reactions: { [key: string]: string[] } = {};
    const reactions = this.reactions;
    // for each node
    for (const node_id in nodes) {
      const node = nodes[node_id];
      // find associated segments and reactions
      node.connected_segments.forEach(function(segment_obj) {
        let segment;
        try {
          segment = reactions[segment_obj.reaction_id].segments[segment_obj.segment_id];
          if (segment === undefined) throw new Error('undefined segment');
        } catch (e) {
          console.warn('Could not find connected segments for node');
          return;
        }
        const segment_obj_w_segment = utils.clone(segment_obj);
        //@ts-ignore
        segment_obj_w_segment['segment'] = utils.clone(segment);
        //@ts-ignore
        segment_objs_w_segments[segment_obj.segment_id] = segment_obj_w_segment;
        if (!(segment_obj.reaction_id in segment_ids_for_reactions))
          segment_ids_for_reactions[segment_obj.reaction_id] = [];
        segment_ids_for_reactions[segment_obj.reaction_id].push(segment_obj.segment_id);
      });
    }
    // find the reactions that should be deleted because they have no segments left
    for (const reaction_id in segment_ids_for_reactions) {
      const reaction = reactions[reaction_id];
      const these_ids = segment_ids_for_reactions[reaction_id];
      let has = true;
      for (const segment_id in reaction.segments)
        if (these_ids.indexOf(segment_id) == -1) has = false;

      if (has) these_reactions[reaction_id] = reaction;
    }
    return {segment_objs_w_segments: segment_objs_w_segments, reactions: these_reactions};
  }

  add_label_to_search_index(id: IDType, text: string) {
    this.search_index.insert('l' + id, {
      name: text,
      data: {type: 'text_label', text_label_id: id}
    });
  }

  new_text_label(coords: { x: number; y: number }, text: string) {
    // Make an label
    const out = build.newTextLabel(this.largest_ids, text, coords);
    this.text_labels[out.id] = out.label;
    this.draw_these_text_labels([out.id]);
    // Add to the search index
    if (text !== '')
      this.add_label_to_search_index(out.id, text);

    return out.id;
  }

  /**
   * Edit a text label. Undoable.
   * @param {} text_label_id -
   * @param {String} new_value -
   * @param {Boolean} should_draw -
   * @param {Boolean} [is_new=false] - If true, then the text label is all new, so
   * delete it on undo and create it again on redo.
   */
  edit_text_label(text_label_id: IDType, new_value: string, should_draw: boolean, is_new?: boolean) {
    if (_.isUndefined(is_new)) is_new = false;

    if (new_value === '')
      throw new Error('Should not be called for empty string');


    const edit_and_draw = (new_val: string, should_draw: boolean) => {
      // set the new value
      const label = this.text_labels[text_label_id];
      label.text = new_val;
      if (should_draw)
        this.draw_these_text_labels([text_label_id]);

      // Update in the search index
      const record_id = 'l' + text_label_id;
      const found = this.search_index.remove(record_id);
      if (!is_new && !found)
        console.warn('Could not find modified text label in search index');

      this.search_index.insert(record_id, {
        name: new_val,
        data: {type: 'text_label', text_label_id: text_label_id}
      });
    };

    // Save old value
    const saved_label = utils.clone(this.text_labels[text_label_id]);

    // Edit the label
    edit_and_draw(new_value, should_draw);

    // Add to undo stack
    this.undo_stack.push(() => {
      if (is_new) {
        this.delete_text_label_data([text_label_id]);
        this.clear_deleted_text_labels();
      } else
        edit_and_draw(saved_label.text, true);
    }, () => {
      if (is_new) {
        this.text_labels[text_label_id] = utils.clone(saved_label);
        this.text_labels[text_label_id].text = new_value;
        this.draw_these_text_labels([text_label_id]);
        this.add_label_to_search_index(text_label_id, new_value);
      } else
        edit_and_draw(new_value, true);
    });
  }

  // -------------------------------------------------------------------------
  // Zoom
  // -------------------------------------------------------------------------

  /**
   * Zoom to fit all the nodes. Returns error if one is raised.
   * @param {} margin - optional argument to set the margins as a fraction of
   * height.
   */
  zoom_extent_nodes(margin?: number) {
    this._zoom_extent(margin, 'nodes');
  }

  /**
   * Zoom to fit the canvas. Returns error if one is raised.
   * @param {} margin - optional argument to set the margins as a fraction of
   * height.
   */
  zoom_extent_canvas(margin?: number) {
    this._zoom_extent(margin, 'canvas');
  }

  /**
   * Zoom to fit the canvas or all the nodes. Returns error if one is raised.
   * @param {} margin - optional argument to set the margins.
   * @param {} mode - Values are 'nodes', 'canvas'.
   */
  _zoom_extent(margin?: number, mode?: 'nodes' | 'canvas') {
    // Optional args
    if (_.isUndefined(margin)) margin = (mode === 'nodes' ? 0.2 : 0);
    if (_.isUndefined(mode)) mode = 'canvas';

    let new_zoom;
    let new_pos;
    const size = this.get_size();
    // scale margin to window size
    margin = margin * size.height;

    if (mode === 'nodes') {
      // get the extent of the nodes
      const min: { x: number | null; y: number | null } = {x: null, y: null}; // TODO make infinity?
      const max: { x: number | null; y: number | null } = {x: null, y: null};
      for (const node_id in this.nodes) {
        const node = this.nodes[node_id];
        if (min.x == null) min.x = node.x;
        if (min.y == null) min.y = node.y;
        if (max.x == null) max.x = node.x;
        if (max.y == null) max.y = node.y;

        min.x = Math.min(min.x, node.x);
        min.y = Math.min(min.y, node.y);
        max.x = Math.max(max.x, node.x);
        max.y = Math.max(max.y, node.y);
      }
      // set the zoom
      new_zoom = Math.min((size.width - margin * 2) / (max.x! - min.x!),
        (size.height - margin * 2) / (max.y! - min.y!));
      new_pos = {x: - (min.x! * new_zoom) + margin + ((size.width - margin * 2 - (max.x! - min.x!) * new_zoom) / 2),
        y: - (min.y! * new_zoom) + margin + ((size.height - margin * 2 - (max.y! - min.y!) * new_zoom) / 2)};
    } else if (mode == 'canvas') {
      // center the canvas
      new_zoom = Math.min((size.width - margin * 2) / (this.canvas.width),
        (size.height - margin * 2) / (this.canvas.height));
      new_pos = {x: - (this.canvas.x * new_zoom) + margin + ((size.width - margin * 2 - this.canvas.width * new_zoom) / 2),
        y: - (this.canvas.y * new_zoom) + margin + ((size.height - margin * 2 - this.canvas.height * new_zoom) / 2)};
    } else
      return console.error('Did not recognize mode');

    this.zoomContainer.goTo(new_zoom, new_pos);
    return null;
  }

  get_size() {
    return this.zoomContainer.get_size();
  }

  zoom_to_reaction(reaction_id: IDType) {
    const reaction = this.reactions[reaction_id];
    const new_zoom = 0.5;
    const size = this.get_size();
    const new_pos = {x: - reaction.label_x * new_zoom + size.width / 2,
      y: - reaction.label_y * new_zoom + size.height / 2};
    this.zoomContainer.goTo(new_zoom, new_pos);
  }

  zoom_to_node(node_id: IDType) {
    const node = this.nodes[node_id];
    const new_zoom = 0.5;
    const size = this.get_size();
    const new_pos = {
      x: - node.label_x * new_zoom + size.width / 2,
      y: - node.label_y * new_zoom + size.height / 2,
    };
    this.zoomContainer.goTo(new_zoom, new_pos);
  }

  zoom_to_text_label(text_label_id: IDType) {
    const text_label = this.text_labels[text_label_id];
    const new_zoom = 0.5;
    const size = this.get_size();
    const new_pos = {
      x: - text_label.x * new_zoom + size.width / 2,
      y: - text_label.y * new_zoom + size.height / 2,
    };
    this.zoomContainer.goTo(new_zoom, new_pos);
  }

  highlight_reaction(reaction_id: IDType) {
    this.highlight(this.sel.selectAll('#r' + reaction_id).selectAll('text'));
  }

  highlight_node(node_id: IDType) {
    this.highlight(this.sel.selectAll('#n' + node_id).selectAll('text'));
  }

  highlight_text_label(text_label_id: IDType) {
    this.highlight(this.sel.selectAll('#l' + text_label_id).selectAll('text'));
  }

  highlight(sel?: D3Selection | null) {
    this.sel.selectAll('.highlight')
      .classed('highlight', false);
    if (sel !== null)
      sel.classed('highlight', true);
  }

  // -------------------------------------------------------------------------
  // IO
  // -------------------------------------------------------------------------

  save() {
    utils.download_json(this.map_for_export(), this.map_name);
  }

  map_for_export() {
    const out = [{map_name: this.map_name,
      map_id: this.map_id,
      map_description: this.map_description,
      homepage: 'https://escher.github.io',
      schema: 'https://escher.github.io/escher/jsonschema/1-0-0#'
    },
    {reactions: utils.clone(this.reactions),
      nodes: utils.clone(this.nodes),
      text_labels: utils.clone(this.text_labels),
      canvas: this.canvas.sizeAndLocation()}
    ];

    // remove extra data
    for (const r_id in out[1].reactions) {
      const reaction = out[1].reactions[r_id];
      const new_reaction: { [key: string]: any } = {};
      const attrs = ['name', 'bigg_id', 'reversibility', 'label_x', 'label_y',
        'gene_reaction_rule', 'genes', 'metabolites'];
      attrs.forEach(function(attr) {
        new_reaction[attr] = reaction[attr];
      });
      new_reaction['segments'] = {};
      for (const s_id in reaction.segments) {
        const segment = reaction.segments[s_id];
        const new_segment: { [key: string]: any } = {};
        const attrs = ['from_node_id', 'to_node_id', 'b1', 'b2'];
        attrs.forEach(function(attr) {
          new_segment[attr] = segment[attr as keyof MapSegment];
        });
        new_reaction['segments'][s_id] = new_segment;
      }
      out[1].reactions[r_id] = new_reaction as MapReaction;
    }
    for (const n_id in out[1].nodes) {
      const node = out[1].nodes[n_id];
      const new_node: { [key: string]: any } = {};
      let attrs: string[] = [];
      if (node.node_type === 'metabolite') {
        attrs = ['node_type', 'x', 'y', 'bigg_id', 'name', 'label_x', 'label_y',
          'node_is_primary'];
      } else
        attrs = ['node_type', 'x', 'y'];

      attrs.forEach(function(attr) {
        new_node[attr] = node[attr];
      });
      out[1].nodes[n_id] = new_node as MapNode;
    }
    for (const t_id in out[1].text_labels) {
      const text_label = out[1].text_labels[t_id];
      const new_text_label: { [key: string]: any } = {};
      const attrs = ['x', 'y', 'text'];
      attrs.forEach(function(attr) {
        new_text_label[attr] = text_label[attr as keyof TextLabel];
      });
      out[1].text_labels[t_id] = new_text_label as TextLabel;
    }
    // canvas
    const canvas_el = out[1].canvas;
    const new_canvas_el: { [key: string]: any } = {};
    const attrs = ['x', 'y', 'width', 'height'];
    attrs.forEach(function(attr) {
      new_canvas_el[attr] = canvas_el![attr as keyof typeof canvas_el];
    });
    out[1].canvas = new_canvas_el as typeof canvas_el;

    return out;
  }

  /**
   * Rescale the canvas and save as svg/png.
   */
  saveMap(callbackBefore: string, callbackAfter: string, mapType: 'svg' | 'png') {
    // Run the before callback
    this.callback_manager.run(callbackBefore);

    // Turn off zoom and translate so that illustrator likes the map
    const windowScale = this.zoomContainer.windowScale;
    const windowTranslate = this.zoomContainer.windowTranslate;
    const canvasSizeAndLoc = this.canvas.sizeAndLocation();
    const mouseNodeSizeAndTrans = {
      w: this.canvas.mouseNode.attr('width'),
      h: this.canvas.mouseNode.attr('height'),
      transform: this.canvas.mouseNode.attr('transform')
    };

    this.zoomContainer._goToSvg(
      1.0,
      {x: -canvasSizeAndLoc.x, y: -canvasSizeAndLoc.y},
      () => {
        this.svg.attr('width', canvasSizeAndLoc.width);
        this.svg.attr('height', canvasSizeAndLoc.height);
        this.canvas.mouseNode.attr('width', '0px');
        this.canvas.mouseNode.attr('height', '0px');
        this.canvas.mouseNode.attr('transform', null);

        // hide the segment control points
        const hidden_sel = this.sel.selectAll('.multimarker-circle,.midmarker-circle,#canvas,.bezier,#rotation-center,.direction-arrow,.start-reaction-target')
          .style('visibility', 'hidden');

        // do the export
        if ( mapType === 'svg')
          utils.downloadSvg('saved_map', this.svg as D3Selection<any>, true);
        else if (mapType === 'png')
          utils.downloadPng('saved_map', this.svg);


        // revert everything
        this.zoomContainer._goToSvg(windowScale, windowTranslate, () => {
          this.svg.attr('width', null);
          this.svg.attr('height', null);
          this.canvas.mouseNode.attr('width', mouseNodeSizeAndTrans.w);
          this.canvas.mouseNode.attr('height', mouseNodeSizeAndTrans.h);
          this.canvas.mouseNode.attr('transform', mouseNodeSizeAndTrans.transform);
          // unhide the segment control points
          hidden_sel.style('visibility', null);

          // run the after callback
          this.callback_manager.run(callbackAfter);
        });
      }
    );
  }

  save_svg() {
    this.saveMap('before_svg_export', 'after_svg_export', 'svg');
  }

  save_png() {
    this.saveMap('before_png_export', 'after_png_export', 'png');
  }

  /**
   * Assign the descriptive names and gene_reaction_rules from the model to the
   * map. If no map is loaded, then throw an Error. If some reactions are not in
   * the model, then warn in the status.
   */
  convert_map() {
    // Run the before callback
    this.callback_manager.run('before_convert_map');

    // Check the model
    if (!this.has_cobra_model())
      throw Error('No COBRA model loaded.');

    const model = this.cobra_model;

    // IDs for reactions and metabolites not found in the model
    const reactions_not_found: { [key: string]: boolean } = {};
    const reaction_attrs = ['name', 'gene_reaction_rule', 'genes'];
    const met_nodes_not_found: { [key: string]: boolean } = {};
    const metabolite_attrs = ['name'];
    let found = false;
    // convert reactions
    for (const reaction_id in this.reactions) {
      const reaction = this.reactions[reaction_id];
      found = false;
      // find in cobra model
      for (const model_reaction_id in model!.reactions) {
        const modelReaction = model!.reactions[model_reaction_id];
        if (modelReaction.bigg_id == reaction.bigg_id) {
          reaction_attrs.forEach((attr) => {
            reaction[attr] = modelReaction[attr];
          });
          // check matching metabolites & warn if not matching. if matching, check
          // & reverse
          let matches = true;
          let looksReversed = null;
          for (const metId in modelReaction.metabolites) {
            const modelCoeff = modelReaction.metabolites[metId];
            const mapMet = _.find(reaction.metabolites, (x) => x.bigg_id === metId);
            if (mapMet === undefined) {
              matches = false;
              break;
            }
            const mapCoeff = mapMet.coefficient;
            // keep track of whether any of these are reversed
            if (looksReversed == null)
              looksReversed = (modelCoeff > 0) !== (mapCoeff > 0);

            // make sure there are no mismatches in reversibility direction
            if ((looksReversed === true && ((modelCoeff > 0) === (mapCoeff > 0))) ||
                (looksReversed === false && ((modelCoeff > 0) !== (mapCoeff > 0)))) {
              matches = false;
              break;
            }
          }
          if (looksReversed && matches) {
            // looks reversed with not mismatches, then reverse the direction
            reaction.metabolites.forEach((met) => {
              met.coefficient = -met.coefficient;
            });
            // propagate changes into segments
            for (const segmentId in reaction.segments) {
              const segment = reaction.segments[segmentId];

              // propagate reversibility
              segment.reversibility = reaction.reversibility;

              const from_node = this.nodes[segment.from_node_id];
              const to_node = this.nodes[segment.to_node_id];

              // propagate coefficients
              reaction.metabolites.forEach((met) => {
                if (met.bigg_id === from_node.bigg_id)
                  segment.from_node_coefficient = met.coefficient;
                else if (met.bigg_id === to_node.bigg_id)
                  segment.to_node_coefficient = met.coefficient;
              });
            }
          }
          if (!matches) {
            console.warn(`Metabolites for ${modelReaction.bigg_id} are different in model and map. Could
 not check and fix direction.`);
            break;
          }
          found = true;
        }
      }
      if (!found) reactions_not_found[reaction_id] = true;
    }
    // convert metabolites
    for (const node_id in this.nodes) {
      const node = this.nodes[node_id];
      // only look at metabolites
      if (node.node_type != 'metabolite') continue;
      found = false;
      // find in cobra model
      for (const model_metabolite_id in model!.metabolites) {
        const model_metabolite = model!.metabolites[model_metabolite_id];
        if (model_metabolite.bigg_id == node.bigg_id) {
          metabolite_attrs.forEach(function(attr) {
            node[attr] = model_metabolite[attr];
          });
          found = true;
        }
      }
      if (!found)
        met_nodes_not_found[node_id] = true;
    }

    // status
    const n_reactions_not_found = Object.keys(reactions_not_found).length;
    const n_met_nodes_not_found = Object.keys(met_nodes_not_found).length;
    const status_delay = 10000;
    if (n_reactions_not_found === 0 && n_met_nodes_not_found === 0)
      this.set_status('Successfully converted attributes.', status_delay);
    else if (n_met_nodes_not_found === 0) {
      this.set_status('Converted attributes, but count not find ' + n_reactions_not_found +
                      ' reactions in the model.', status_delay);
      this.settings.set('highlight_missing', true);
    } else if (n_reactions_not_found === 0) {
      this.set_status('Converted attributes, but count not find ' + n_met_nodes_not_found +
                      ' metabolites in the model.', status_delay);
      this.settings.set('highlight_missing', true);
    } else {
      this.set_status('Converted attributes, but count not find ' + n_reactions_not_found +
                      ' reactions and ' + n_met_nodes_not_found + ' metabolites in the model.',
      status_delay);
      this.settings.set('highlight_missing', true);
    }

    // redraw
    this.draw_everything();

    // run the after callback
    this.callback_manager.run('after_convert_map');
  }

  mergeSelectedNodes() {
    try {
      const selNodes = this.getSelectedNodes();
      const nodeObjs = Object.values(selNodes);
      const nodeKeys = Object.keys(selNodes);
      if (Object.values(selNodes).length != 2)
        return;
      const firstNodeId = nodeKeys[0];
      const secondNodeId = nodeKeys[1];
      if (nodeObjs.every((a) => a.bigg_id) && nodeObjs[0].bigg_id == nodeObjs[1].bigg_id) {
        this.nodes[firstNodeId].connected_segments.forEach((seg) => {
          this.nodes[secondNodeId].connected_segments.push({...seg});
        });
        this.nodes[firstNodeId] && (delete this.nodes[firstNodeId]);
        Object.values(this.reactions).forEach((r) => {
          Object.values(r.segments).forEach((seg) => {
            if (seg.from_node_id == firstNodeId)
              seg.from_node_id = secondNodeId;
            if (seg.to_node_id == firstNodeId)
              seg.to_node_id = secondNodeId;
          });
        });

        this.sel.selectAll('.selected').classed('selected', false);
        this.draw_everything();
      }
    } catch (e) {
      console.error(e);
    }
  }

  highlightSameNodes() {
    try {
      const selNodes = Object.values(this.getSelectedNodes());
      if (selNodes.length === 1 && selNodes[0].node_type == 'metabolite' && selNodes[0].bigg_id) {
        const nodeBiggId = selNodes[0].bigg_id;
        const applicableNodeIds = Object.entries(this.nodes).filter(([nodeId, node]) => node.node_type == 'metabolite' && node.bigg_id == nodeBiggId).map(([nodeId, _]) => nodeId);
        this.sel.selectAll(applicableNodeIds.map((id) => `#n${id}`).join(',')).classed('selected', true);

        // Object.values(this.nodes).forEach((n) => {
        //   if (n.bigg_id === this.nodes[selNodeId].bigg_id) {
        //     this.sel.selectAll('#n' + n.node_id).selectAll('text').classed('highlight', true);
        //   }
        // })
      }
    } catch (e) {
      console.error(e);
    }
  }

  findAndHighlightShortest(inc = 0) {
    this.last_action = null;

    const nodes = Array.from(new Set(Object.values(this.getSelectedNodes())
      .sort((node1, node2) => ((!!node1.selectionOrder ? node1.selectionOrder : 0) - (!!node2.selectionOrder ? node2.selectionOrder : 0)))
      .map((a) => a.bigg_id)));
    if (nodes.length !== 2)
      return;

    let path;
    let k;
    try {
      const reactionHash = `r${nodes[0]}_r${nodes[1]}`;
      k = typeof inc === 'string' ? Number.parseInt(inc) : inc == 0 ? 1 : Math.max(this.shortestPathInfo ? (this.shortestPathInfo.hash === reactionHash ? this.shortestPathInfo.k + inc : 1) : 1, 1);
      if (k === 0)
        k = 1;
      path = findKthShortestPath(this, nodes[0], nodes[1], k);
      this.shortestPathInfo = {hash: reactionHash, k};
      this.outGoingInfo = null;
      this.inGoingInfo = null;
    } catch (e) {
      console.error(`Error finding ${k}th shortest path`, e);
      return;
    }
    this.unhighlightShortestPaths();
    //const kthPath = findKthShortestPath(this, nodes[0], nodes[1], 20);
    if (!path)
      return;
    for (let i = 0; i < path.length - 1; i++) {
      const reactionId = path[i].reaction_id;
      const firtNodeBiggId = path[i].met;
      const secondNodeBiggId = path[i + 1].met;
      if (!reactionId || !firtNodeBiggId || !secondNodeBiggId)
        continue;
      const segments = getReactionPathSegments(this, firtNodeBiggId, secondNodeBiggId, reactionId);
      segments?.forEach((seg) => {
        this.sel.select(`#s${seg}.segment-group`).classed('shortest-path', true);
      });
    }

    // path.filter(p => p.reaction_id != null).forEach(p => {
    //   this.sel.select(`#r${p.reaction_id}.reaction`).classed('shortest-path', true)
    // })
    this.last_action = {
      function: 'findAndHighlightShortest',
      args: [k.toString()]
    };
  }

  findKthOutOrIngoingReactions(inc: string | number, outGoingFlag = true) {
    this.last_action = null;
    const nodes = Array.from(new Set(Object.values(this.getSelectedNodes()).map((a) => a.bigg_id)));
    if (nodes.length !== 1)
      return;
    let k: number;

    const biggId = nodes[0];
    this.unhighlightShortestPaths();

    const getSavedInfo = () => outGoingFlag ? this.outGoingInfo : this.inGoingInfo;
    const setSavedInfo = (info: { hash: string; k: number; }) => outGoingFlag ? (this.outGoingInfo = info) : (this.inGoingInfo = info);
    const calcFunc = outGoingFlag ? findOutGoing : findInGoing;
    // If you want to set k explicitely, pass inc as string
    k = typeof inc === 'string' ? Number.parseInt(inc) : inc == 0 ? 1 : Math.max(getSavedInfo() ? (getSavedInfo()!.hash === biggId ? getSavedInfo()!.k + inc : 1) : 1, 1);
    if (k === 0)
      k = 1;
    this.outGoingInfo = null;
    this.inGoingInfo = null;
    this.shortestPathInfo = null;
    setSavedInfo({hash: biggId, k});
    const visitedMetabolites = new Set<IDType>();
    const visitedReactions = new Set<string>();
    const newMetabolites = new Set<IDType>();
    // outGoing.forEach((o) => {
    //   visitedMetabolites.add(o[1]);
    //   visitedReactions.add(o[0]);
    // });
    const sign = (num: number) => num >= 0 ? 1 : -1;

    const queue: IDType[] = [biggId];
    for (let iter = 0; iter < k; iter++) {
      newMetabolites.clear();
      while (queue.length > 0) {
        const currentId = queue.shift() as IDType;
        if (visitedMetabolites.has(currentId))
          continue;
        visitedMetabolites.add(currentId);
        const outGoing = calcFunc(this, currentId).map((a) => a.split(pathSeparatorString));
        outGoing.forEach((o) => {
          const reactionId = o[0];
          const toId = o[1];
          if (!visitedReactions.has([currentId, o[0], o[1]].join(pathSeparatorString)) && !visitedMetabolites.has(toId)) {
            newMetabolites.add(toId);
            // what we need to do here is to go through every reaction and find ones which have from and to metabolites in correct order
            const reactionIds: IDType[] = [];
            Object.values(this.reactions).forEach((reaction) => {
              const reactionMets = reaction.metabolites;
              const fromMet = reactionMets.find((m) => m.bigg_id === currentId);
              const toMet = reactionMets.find((m) => m.bigg_id === toId);
              if (fromMet && toMet && (reaction.reversibility ||
                (outGoingFlag && sign(fromMet.coefficient) === -1 && sign(toMet.coefficient) === 1) ||
                (!outGoingFlag && sign(fromMet.coefficient) === 1 && sign(toMet.coefficient) === -1)
              ))
                reactionIds.push(reaction.reaction_id);
            });
            reactionIds.forEach((rId) => {
              const segments = getReactionPathSegments(this, currentId, toId, rId);
              segments && segments.forEach((seg) => {
                const nthShortest = Math.max(Math.min(k, 4), 1);
                this.sel.select(`#s${seg}.segment-group`).classed('shortest-path', true);
              //.classed(`shortest-path-${nthShortest}`, true);
              });
            });
          }
        });

        outGoing.forEach((o) => {
          visitedReactions.add([currentId, o[0], o[1]].join(pathSeparatorString));
        });
      }
      Array.from(newMetabolites).forEach((m) => queue.push(m));
    }

    this.last_action = {
      function: 'findKthOutOrIngoingReactions',
      args: [k.toString(), outGoingFlag]
    };
  }

  select_nodes(nodeIds: IDType[]) {
    nodeIds.filter((id) => id && id.toString()).forEach((id) => this.sel.selectAll('#n' + id.toString()).classed('selected', true));
  }

  selectionChanged() {
    // handle different things when selection changes
    const element = document.getElementById('selection-info-container');
    const prefix = 'Selected: ';
    if (!element)
      return;
    const nodes = Array.from(new Set(Object.values(this.getSelectedNodes())
      .sort((node1, node2) => ((!!node1.selectionOrder ? node1.selectionOrder : 0) - (!!node2.selectionOrder ? node2.selectionOrder : 0)))
      .map((a) => a.bigg_id)));

    if (nodes.length === 0) {
      element.innerHTML = '';
      return;
    }
    if (nodes.length === 1) {
      element.innerHTML = `<h3>${prefix}${nodes[0]}</h3>`;
      return;
    }
    if (nodes.length === 2) {
      element.innerHTML = `<h3>${prefix}${nodes[0]} - ${nodes[1]}</h3>`;
      return;
    }

    element.innerHTML = ``;
  }

  unhighlightShortestPaths() {
    this.sel.selectAll('.shortest-path').classed('shortest-path', false);
    this.sel.selectAll('.shortest-path-1').classed('shortest-path-1', false);
    this.sel.selectAll('.shortest-path-2').classed('shortest-path-2', false);
    this.sel.selectAll('.shortest-path-3').classed('shortest-path-3', false);
    this.sel.selectAll('.shortest-path-4').classed('shortest-path-4', false);
    this.last_action = null;
  }
}
