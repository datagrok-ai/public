/* eslint-disable camelcase */
/**
 * Draw. Manages creating, updating, and removing objects during d3 data
 * binding.
 *
 * Arguments
 * ---------
 *
 * behavior: An escher.Behavior object.
 * settings: An escher.Settings object.
 *
 * Callbacks
 * ---------
 *
 * draw.callback_manager.run('create_membrane', draw, enter_selection)
 * draw.callback_manager.run('update_membrane', draw, update_selection)
 * draw.callback_manager.run('create_reaction', draw, enter_selection)
 * draw.callback_manager.run('update_reaction', draw, update_selection)
 * draw.callback_manager.run('create_reaction_label', draw, enter_selection)
 * draw.callback_manager.run('update_reaction_label', draw, update_selection)
 * draw.callback_manager.run('create_segment', draw, enter_selection)
 * draw.callback_manager.run('update_segment', draw, update_selection)
 * draw.callback_manager.run('create_bezier', draw, enter_selection)
 * draw.callback_manager.run('update_bezier', draw, update_selection)
 * draw.callback_manager.run('create_node', draw, enter_selection)
 * draw.callback_manager.run('update_node', draw, update_selection)
 * draw.callback_manager.run('create_text_label', draw, enter_selection)
 * draw.callback_manager.run('update_text_label', draw, update_selection)
 *
 */

import CallbackManager from './CallbackManager';
import * as utils from './utils';
import * as dataStyles from './dataStyles';
import {format as d3_format} from 'd3-format';
import Behavior from './Behavior';
import {EscherMap} from './escherMap';
import Settings from './Settings';
import {Coord, D3Selection} from './types';
import Scale from './Scale';

export default class Draw {
  behavior: Behavior;
  settings: Settings;
  map: EscherMap;
  callback_manager: CallbackManager;

  constructor(behavior: Behavior, settings: Settings, map: EscherMap) {
    this.behavior = behavior;
    this.settings = settings;
    this.map = map;
    this.callback_manager = new CallbackManager();
  }

  /**
 * Create membranes in the enter_selection.
 * @param {} enter_selection - The D3 enter selection.
 * @returns {} The selection of the new nodes.
 */
  create_membrane(enter_selection: D3Selection) {
    const rect = enter_selection
      .append('rect')
      .attr('class', 'membrane');

    this.callback_manager.run('create_membrane', this, enter_selection);

    return rect;
  }

  /**
 * Update the membrane
 */
  update_membrane(update_selection: D3Selection) {
    update_selection
      .attr('width', function(d: any) { return d.width; })
      .attr('height', function(d: any) { return d.height; })
      .attr('transform', function(d: any) { return 'translate(' + d.x + ',' + d.y + ')'; })
      .style('stroke-width', function(d: any) { return 10; })
      .attr('rx', function(d: any) { return 20; })
      .attr('ry', function(d: any) { return 20; });

    this.callback_manager.run('update_membrane', this, update_selection);
  }

  /**
 * Create reactions in the enter_selection.
 * @param {} enter_selection - The D3 enter selection.
 * @returns {} The selection of the new nodes.
 */
  create_reaction(enter_selection: D3Selection) {
  // attributes for new reaction group
    const group = enter_selection.append('g')
      .attr('id', function(d: any) { return 'r' + d.reaction_id; })
      .attr('class', 'reaction');
    this.create_reaction_label(group);

    this.callback_manager.run('create_reaction', this, enter_selection);

    return group;
  }

  /**
 * Run on the update selection for reactions.
 * update_selection: The D3.js update selection.
 * scale: A Scale object.
 * cobra_model: A CobraModel object.
 * drawn_nodes: The nodes object (e.g. Map.nodes).
 * defs: The defs object generated by utils.setup_defs() (e.g. Map.defs).
 * has_data_on_reactions: Boolean to determine whether data needs to be drawn.
 */
  update_reaction(update_selection: D3Selection, scale: any, cobra_model: any, drawn_nodes: any,
    defs: any, has_data_on_reactions: boolean) {
  // Update reaction label
    update_selection.select('.reaction-label-group')
      .call((sel: D3Selection<any>) => {
        return this.update_reaction_label(sel, has_data_on_reactions);
      });

    // draw segments
    utils.draw_a_nested_object(update_selection, '.segment-group', 'segments', 'segment_id',
      this.create_segment.bind(this),
      (sel: D3Selection<any>) => {
        return this.update_segment(sel, scale, cobra_model,
          drawn_nodes, defs,
          has_data_on_reactions);
      },
      (sel: D3Selection<any>) => {
        sel.remove();
      });

    // run the callback
    this.callback_manager.run('update_reaction', this, update_selection);
  }

  /**
 * Draw reaction labels in the enter selection.
 * @param {} enter_selection - The D3 enter selection.
 * @returns {} The selection of the new nodes.
 */
  create_reaction_label(enter_selection: D3Selection<any>, tool?: any) {
    const group = enter_selection
      .append('g')
      .attr('class', 'reaction-label-group');
    group.append('text').attr('class', 'reaction-label label');
    group.append('g').attr('class', 'all-genes-label-group');

    this.callback_manager.run('create_reaction_label', this, enter_selection);

    return group;
  }

  /**
 * Run on the update selection for reaction labels.
 * @param {D3 Selection} update_selection - The D3.js update selection.
 * @param {Boolean} has_data_on_reactions - Whether data needs to be drawn.
 */
  update_reaction_label(update_selection: D3Selection<any>, has_data_on_reactions: boolean) {
    const decimal_format = d3_format('.4g');
    const identifiers_on_map = this.settings.get('identifiers_on_map');
    const reaction_data_styles = this.settings.get('reaction_styles');
    const show_gene_reaction_rules = this.settings.get('show_gene_reaction_rules');
    const hide_all_labels = this.settings.get('hide_all_labels');
    const gene_font_size = this.settings.get('gene_font_size');
    const reactionLabelMouseover = this.behavior.reactionLabelMouseover;
    const reactionLabelMouseout = this.behavior.reactionLabelMouseout;
    const reactionLabelTouch = this.behavior.reactionLabelTouch;
    const geneLabelMouseover = this.behavior.geneLabelMouseover;
    const geneLabelMouseout = this.behavior.geneLabelMouseout;
    const geneLabelTouch = this.behavior.geneLabelTouch;

    // label location
    update_selection
      .attr('transform', function(d: any) {
        return 'translate(' + d.label_x + ',' + d.label_y + ')';
      })
      .call(this.behavior.turnOffDrag)
      .call(this.behavior.reactionLabelDrag);

    // update label visibility
    const label = update_selection.select('.reaction-label')
      .attr('visibility', hide_all_labels ? 'hidden' : 'visible');

    if (!hide_all_labels) {
      label
        .text(function(d: any) {
          let t = d[identifiers_on_map];
          if (has_data_on_reactions &&
            reaction_data_styles.indexOf('text') !== -1)
            t += ' ' + d.data_string;

          return t;
        })
        .on('mouseover', reactionLabelMouseover)
        .on('mouseout', reactionLabelMouseout)
        .on('touchend', reactionLabelTouch);
    }

    const add_gene_height = function(y: number, i: number) {
      return y + (gene_font_size * 1.5 * (i + 1));
    };

    // gene label
    const all_genes_g = update_selection.select('.all-genes-label-group')
      .selectAll('.gene-label-group')
      .data(function(d: any) {
        const show_gene_string = ('gene_string' in d &&
                              d.gene_string !== null &&
                              show_gene_reaction_rules &&
                              (!hide_all_labels) &&
                              reaction_data_styles.indexOf('text') !== -1);
        const show_gene_reaction_rule = ('gene_reaction_rule' in d &&
                                     d.gene_reaction_rule !== null &&
                                     show_gene_reaction_rules &&
                                     (!hide_all_labels));
        if (show_gene_string) {
        // TODO do we ever use gene_string?
          console.warn('Showing gene_string. See TODO in source.');
          return d.gene_string;
        } else if (show_gene_reaction_rule) {
        // make the gene string with no data
          const sd = dataStyles.gene_string_for_data(d.gene_reaction_rule, null,
            d.genes, null,
            identifiers_on_map, null);
          // add coords for tooltip
          sd.forEach(function(td: any, i: number) {
            td.label_x = d.label_x;
            td.label_y = add_gene_height(d.label_y, i);
          });
          return sd;
        } else
          return [];
      });

    // enter
    const gene_g = all_genes_g.enter()
      .append('g')
      .attr('class', 'gene-label-group');
    gene_g.append('text')
      .attr('class', 'gene-label')
      .style('font-size', gene_font_size + 'px');

    // update
    const gene_update = gene_g.merge(all_genes_g as D3Selection<any>);
    gene_update.attr('transform', function(d: any, i: number) {
      return 'translate(0, ' + add_gene_height(0, i) + ')';
    });
    // update text
    gene_update
      .select('text')
      .text((d: any) => d.text)
      .on('mouseover', geneLabelMouseover)
      .on('mouseout', geneLabelMouseout)
      .on('touchend', geneLabelTouch);

    // exit
    all_genes_g.exit().remove();

    this.callback_manager.run('update_reaction_label', this, update_selection);
  }

  /**
 * Create segments in the enter_selection.
 * @param {} enter_selection - The D3 enter selection.
 * @returns {} The selection of the new nodes.
 */
  create_segment(enter_selection: D3Selection) {
  // create segments
    const g = enter_selection
      .append('g')
      .attr('class', 'segment-group')
      .attr('id', function(d: any) { return 's' + d.segment_id; });

    // create reaction arrow
    g.append('path')
      .attr('class', 'segment');

    g.append('g')
      .attr('class', 'arrowheads');

    g.append('g')
      .attr('class', 'stoichiometry-labels');

    this.callback_manager.run('create_segment', this, enter_selection);

    return g;
  }

  /**
 * Update segments in update selection.
 * @param {} -
 * @param {} -
 * @param {} -
 * @param {} -
 * @param {} -
 * @param {} -
 * @return {}
 */
  update_segment(update_selection: D3Selection, scale: any, cobra_model: any,
    drawn_nodes: any, defs: any, has_data_on_reactions: boolean) {
    const reaction_data_styles = this.settings.get('reaction_styles');
    const should_size = (has_data_on_reactions && reaction_data_styles.indexOf('size') !== -1);
    const should_color = (has_data_on_reactions && reaction_data_styles.indexOf('color') !== -1);
    const no_data_size = this.settings.get('reaction_no_data_size');
    const no_data_color = this.settings.get('reaction_no_data_color');

    // update segment attributes
    const highlight_missing = this.settings.get('highlight_missing');
    const hide_secondary_metabolites = this.settings.get('hide_secondary_metabolites');
    const primary_r = this.settings.get('primary_metabolite_radius');
    const secondary_r = this.settings.get('secondary_metabolite_radius');

    const objectMouseover = this.behavior.reactionObjectMouseover;
    const objectMouseout = this.behavior.reactionObjectMouseout;

    const get_arrow_size = function(data: any, should_size: boolean) {
      let width = 20;
      let height = 13;
      if (should_size) {
        height = (data == null ? no_data_size : scale.reaction_size(data));
        // check for nan
        if (isNaN(height))
          height = no_data_size;

        width = height * 2;
      }
      return {width: width, height: height};
    };
    const get_disp = function(arrow_size: any, reversibility: boolean, coefficient: number, node_is_primary: boolean) {
      const arrow_height = ((reversibility || coefficient > 0) ?
        arrow_size.height : 0);
      const r = node_is_primary ? primary_r : secondary_r;
      return r + arrow_height + 10;
    };

    // update arrows
    update_selection
      .selectAll('.segment')
      .datum(function() {
      // Concatenate the segment data with the reaction data from its parent node
      // @ts-ignore
      // eslint-disable-next-line no-invalid-this
        return Object.assign({}, this.parentNode.__data__, this.parentNode.parentNode.__data__);
      })
      .style('visibility', function(d: any) {
        const start = drawn_nodes[d.from_node_id];
        const end = drawn_nodes[d.to_node_id];
        if (hide_secondary_metabolites &&
          ((end['node_type'] === 'metabolite' && !end.node_is_primary) ||
           (start['node_type'] === 'metabolite' && !start.node_is_primary)))
          return 'hidden';

        return null;
      })
      .attr('d', function(d: any) {
        if (d.from_node_id == null || d.to_node_id == null)
          return null;

        let start = drawn_nodes[d.from_node_id];
        let end = drawn_nodes[d.to_node_id];
        const b1 = d.b1;
        const b2 = d.b2;
        // if metabolite, then displace the arrow
        if (start['node_type'] === 'metabolite') {
          const arrow_size = get_arrow_size(d.data, should_size);
          const disp = get_disp(arrow_size, d.reversibility,
            d.from_node_coefficient,
            start.node_is_primary);
          const direction = (b1 == null) ? end : b1;
          start = displacedCoords(disp, start, direction, 'start');
        }
        if (end['node_type'] == 'metabolite') {
          const arrow_size = get_arrow_size(d.data, should_size);
          const disp = get_disp(arrow_size, d.reversibility,
            d.to_node_coefficient,
            end.node_is_primary);
          const direction = (b2 == null) ? start : b2;
          end = displacedCoords(disp, direction, end, 'end');
        }
        let curve = ('M' + start.x + ',' + start.y + ' ');
        if (b1 !== null && b2 !== null) {
          curve += ('C' + b1.x + ',' + b1.y + ' ' +
                  b2.x + ',' + b2.y + ' ');
        }
        curve += (end.x + ',' + end.y);
        return curve;
      })
      .style('stroke', function(d: any) {
        // @ts-ignore
        // eslint-disable-next-line no-invalid-this
        const reaction_id = this.parentNode.parentNode.__data__.bigg_id;
        const show_missing = (highlight_missing &&
                          cobra_model !== null &&
                          !(reaction_id in cobra_model.reactions));
        if (show_missing)
          return 'red';

        if (should_color) {
          const f = d.data;
          return f == null ? no_data_color : scale.reaction_color(f);
        }
        return null;
      })
      .style('stroke-width', function(d) {
        if (should_size) {
          const f = d.data;
          return f == null ? no_data_size : scale.reaction_size(f);
        } else
          return null;
      })
      .attr('pointer-events', 'visibleStroke')
      .on('mouseover', objectMouseover)
      .on('mouseout', objectMouseout);

    // new arrowheads
    const arrowheads = update_selection.select('.arrowheads')
      .selectAll('.arrowhead')
      .data(function(d: any) {
        const arrowheads: any[] = [];
        const start = drawn_nodes[d.from_node_id];
        const b1 = d.b1;
        const end = drawn_nodes[d.to_node_id];
        const b2 = d.b2;
        // hide_secondary_metabolites option
        if (hide_secondary_metabolites &&
          ((end['node_type'] === 'metabolite' && !end.node_is_primary) ||
           (start['node_type'] === 'metabolite' && !start.node_is_primary)))
          return arrowheads;


        if (start.node_type === 'metabolite' &&
          (d.reversibility || d.from_node_coefficient > 0)) {
          const arrow_size = get_arrow_size(d.data, should_size);
          const disp = get_disp(arrow_size, d.reversibility,
            d.from_node_coefficient,
            start.node_is_primary);
          const direction = (b1 == null) ? end : b1;
          const rotation = utils.to_degrees(utils.get_angle([start, direction])) + 90;
          const loc = displacedCoords(disp, start, direction, 'start');
          arrowheads.push({
            data: d.data,
            x: loc!.x,
            y: loc!.y,
            size: arrow_size,
            rotation: rotation,
            show_arrowhead_flux: (((d.from_node_coefficient < 0) === d.reverse_flux) || d.data === 0)
          });
        }

        if (end.node_type === 'metabolite' &&
          (d.reversibility || d.to_node_coefficient > 0)) {
          const arrow_size = get_arrow_size(d.data, should_size);
          const disp = get_disp(arrow_size, d.reversibility,
            d.to_node_coefficient,
            end.node_is_primary);
          const direction = (b2 == null) ? start : b2;
          const rotation = utils.to_degrees(utils.get_angle([end, direction])) + 90;
          const loc = displacedCoords(disp, direction, end, 'end');
          arrowheads.push({
            data: d.data,
            x: loc!.x,
            y: loc!.y,
            size: arrow_size,
            rotation: rotation,
            show_arrowhead_flux: (((d.to_node_coefficient < 0) === d.reverse_flux) || d.data === 0)
          });
        }

        if (d.unconnected_segment_with_arrow) {
          const arrow_size = get_arrow_size(d.data, should_size);
          const direction = end;
          const rotation = utils.to_degrees(utils.get_angle([start, direction])) + 90;
          arrowheads.push({
            data: d.data,
            x: start.x,
            y: start.y,
            size: arrow_size,
            rotation: rotation,
            show_arrowhead_flux: (((d.to_node_coefficient < 0) === d.reverse_flux) || d.data === 0)
          });
        }

        return arrowheads;
      });
    arrowheads.enter().append('path')
      .classed('arrowhead', true)
    // update arrowheads
    // @ts-ignore
      .merge(arrowheads)
      .attr('d', function(d) {
        return ('M' + [-d.size.width / 2, 0] +
              ' L' + [0, d.size.height] +
              ' L' + [d.size.width / 2, 0] + ' Z');
      }).attr('transform', function(d) {
        return 'translate(' + d.x + ',' + d.y + ')rotate(' + d.rotation + ')';
      }).style('fill', function(d) {
        if (should_color) {
          if (d.show_arrowhead_flux) {
          // show the flux
            const f = d.data;
            return f == null ? no_data_color : scale.reaction_color(f);
          } else {
          // if the arrowhead is not filled because it is reversed
            return '#FFFFFF';
          }
        }
        // default fill color
        return null;
      }).style('stroke', function(d) {
        if (should_color) {
        // show the flux color in the stroke whether or not the fill is present
          const f = d.data;
          return f == null ? no_data_color : scale.reaction_color(f);
        }
        // default stroke color
        return null;
      });
    // remove
    arrowheads.exit().remove();

    // new stoichiometry labels
    const stoichiometry_labels = update_selection.select('.stoichiometry-labels')
      .selectAll('.stoichiometry-label')
      .data(function(d: any) {
        const labels: any[] = [];
        const start = drawn_nodes[d.from_node_id];
        const b1 = d.b1;
        const end = drawn_nodes[d.to_node_id];
        const b2 = d.b2;
        const disp_factor = 1.5;

        // hide_secondary_metabolites option
        if (hide_secondary_metabolites &&
          ((end['node_type'] == 'metabolite' && !end.node_is_primary) ||
           (start['node_type'] == 'metabolite' && !start.node_is_primary)))
          return labels;


        if (start.node_type === 'metabolite' && (Math.abs(d.from_node_coefficient) != 1)) {
          const arrow_size = get_arrow_size(d.data, should_size);
          const disp = disp_factor * get_disp(arrow_size, false, 0, end.node_is_primary);
          let direction = (b1 == null) ? end : b1;
          direction = utils.c_plus_c(direction, utils.rotate_coords(direction, 0.5, start));
          let loc = displacedCoords(disp, start, direction, 'start');
          loc = utils.c_plus_c(loc, {x: 0, y: 7})!;
          labels.push({
            coefficient: Math.abs(d.from_node_coefficient),
            x: loc.x,
            y: loc.y,
            data: d.data,
          });
        }

        if (end.node_type === 'metabolite' && (Math.abs(d.to_node_coefficient) !== 1)) {
          const arrow_size = get_arrow_size(d.data, should_size);
          const disp = disp_factor * get_disp(arrow_size, false, 0, end.node_is_primary);
          let direction = (b2 == null) ? start : b2;
          direction = utils.c_plus_c(direction,
            utils.rotate_coords(direction, 0.5, end));
          let loc = displacedCoords(disp, direction, end, 'end');
          loc = utils.c_plus_c(loc, {x: 0, y: 7})!;
          labels.push({
            coefficient: Math.abs(d.to_node_coefficient),
            x: loc.x,
            y: loc.y,
            data: d.data,
          });
        }
        return labels;
      });

    // add labels
    stoichiometry_labels.enter()
      .append('text')
      .attr('class', 'stoichiometry-label')
      .attr('text-anchor', 'middle')
    // update stoichiometry_labels
    // @ts-ignore
      .merge(stoichiometry_labels)
      .attr('transform', function(d) {
        return 'translate(' + d.x + ',' + d.y + ')';
      })
      .text(function(d) {
        return d.coefficient;
      })
      .style('fill', function(d) {
        if (should_color) {
        // show the flux color
          const f = d.data;
          return f == null ? no_data_color : scale.reaction_color(f);
        }
        // default segment color
        return null;
      });

    // remove
    stoichiometry_labels.exit().remove();

    this.callback_manager.run('update_segment', this, update_selection);
  }

  /**
 * Create beziers in the enter_selection.
 * @param {} enter_selection - The D3 enter selection.
 * @returns {} The selection of the new nodes.
 */
  create_bezier(enter_selection: D3Selection) {
    const g = enter_selection.append('g')
      .attr('id', function(d: any) { return d.bezier_id; })
      .attr('class', function(d: any) { return 'bezier'; });
    g.append('path')
      .attr('class', 'connect-line');
    g.append('circle')
      .attr('class', function(d: any) { return 'bezier-circle ' + d.bezier; })
      .style('stroke-width', String(1) + 'px')
      .attr('r', String(7) + 'px');

    this.callback_manager.run('create_bezier', this, enter_selection);

    return g;
  }

  /**
 * Update beziers in update_selection.
 */
  update_bezier(update_selection: D3Selection, show_beziers: boolean, drag_behavior: any,
    mouseover: any, mouseout: any, drawn_nodes: any, drawn_reactions: any) {
    const hide_secondary_metabolites = this.settings.get('hide_secondary_metabolites');

    if (!show_beziers) {
      update_selection.attr('visibility', 'hidden');
      return;
    } else
      update_selection.attr('visibility', 'visible');


    // hide secondary
    update_selection
      .style('visibility', function(d: any) {
        const seg_data = drawn_reactions[d.reaction_id].segments[d.segment_id];
        const start = drawn_nodes[seg_data.from_node_id];
        const end = drawn_nodes[seg_data.to_node_id];
        if (hide_secondary_metabolites &&
          ((end['node_type'] === 'metabolite' && !end.node_is_primary) ||
           (start['node_type'] === 'metabolite' && !start.node_is_primary)))
          return 'hidden';

        return null;
      });

    // Draw bezier points
    update_selection.select('.bezier-circle')
      // @ts-ignore
      .call(this.behavior.turnOffDrag)
      // @ts-ignore
      .call(drag_behavior)
      .on('mouseover', mouseover)
      .on('mouseout', mouseout)
      .attr('transform', function(d: any) {
        if (d.x == null || d.y == null) return '';
        return 'translate(' + d.x + ',' + d.y + ')';
      });

    // Update bezier line
    update_selection
      .select('.connect-line')
      .attr('d', function(d: any) {
        const segment_d = drawn_reactions[d.reaction_id].segments[d.segment_id];
        const node = d.bezier === 'b1' ?
          drawn_nodes[segment_d.from_node_id] :
          drawn_nodes[segment_d.to_node_id];
        if (d.x == null || d.y == null || node.x == null || node.y == null)
          return '';

        return 'M' + d.x + ', ' + d.y + ' ' + node.x + ',' + node.y;
      });

    this.callback_manager.run('update_bezier', this, update_selection);
  }

  /**
 * Create nodes in the enter_selection.
 * @param {} enter_selection - The D3 enter selection.
 * @param {} drawn_nodes - The nodes object (e.g. Map.nodes).
 * @param {} drawn_reactions - The reactions object (e.g. Map.reactions).
 * @returns {} The selection of the new nodes.
 */
  create_node(enter_selection: D3Selection, drawn_nodes: any, drawn_reactions: any) {
  // create nodes
    const g = enter_selection
      .append('g')
      .attr('class', 'node')
      .attr('id', function(d: any) { return 'n' + d.node_id; });

    // create metabolite circle and label
    g.append('circle')
      .attr('class', function(d: any) {
        let c = 'node-circle';
        if (d.node_type !== null)
          c += (' ' + d.node_type + '-circle');
        return c;
      });

    // labels
    const metabolite_groups = g.filter(function(d: any) {
      return d.node_type === 'metabolite';
    });

    metabolite_groups.append('text')
      .attr('class', 'node-label label');

    this.callback_manager.run('create_node', this, enter_selection);

    return g;
  }

  /**
 * Run on the update selection for nodes.
 * @param {D3 Selection} update_selection - The D3.js update selection.
 * @param {Scale} scale - A Scale object.
 * @param {Boolean} has_data_on_nodes - Boolean to determine whether data needs to be drawn.
 * @param {Function} mousedown_fn - A function to call on mousedown for a node.
 * @param {Function} click_fn - A function to call on click for a node.
 * @param {Function} mouseover_fn - A function to call on mouseover for a node.
 * @param {Function} mouseout_fn - A function to call on mouseout for a node.
 * @param {D3 Behavior} drag_behavior - The D3.js drag behavior object for the nodes.
 * @param {D3 Behavior} label_drag_behavior - The D3.js drag behavior object for the node labels.
 */
  update_node(update_selection: D3Selection<any>, scale: Scale, has_data_on_nodes: boolean,
    mousedown_fn: any, click_fn: any, mouseover_fn: any, mouseout_fn: any,
    drag_behavior: any, label_drag_behavior: any) {
  // update circle and label location
    const hide_secondary_metabolites = this.settings.get('hide_secondary_metabolites');
    const primary_r = this.settings.get('primary_metabolite_radius');
    const secondary_r = this.settings.get('secondary_metabolite_radius');
    const marker_r = this.settings.get('marker_radius');
    const hide_all_labels = this.settings.get('hide_all_labels');
    const identifiers_on_map = this.settings.get('identifiers_on_map');
    const metabolite_data_styles = this.settings.get('metabolite_styles');
    const no_data_style = {color: this.settings.get('metabolite_no_data_color'),
      size: this.settings.get('metabolite_no_data_size')};
    const labelMouseover = this.behavior.nodeLabelMouseover;
    const labelMouseout = this.behavior.nodeLabelMouseout;
    const labelTouch = this.behavior.nodeLabelTouch;
    const objectMouseover = this.behavior.nodeObjectMouseover;
    const objectMouseout = this.behavior.nodeObjectMouseout;

    const mg = update_selection
      .select('.node-circle')
      .attr('transform', function(d: any) {
        return 'translate(' + d.x + ',' + d.y + ')';
      })
      .style('visibility', function(d: any) {
        return hideNode(d, hide_secondary_metabolites) ? 'hidden' : null;
      })
      .attr('r', function(d: any) {
        if (d.node_type === 'metabolite') {
          const should_scale = (has_data_on_nodes &&
                            metabolite_data_styles.indexOf('size') !== -1);
          if (should_scale) {
            const f = d.data;
            return f == null ? no_data_style['size'] : scale.metabolite_size(f);
          } else
            return d.node_is_primary ? primary_r : secondary_r;
        }
        // midmarkers and multimarkers
        return marker_r;
      })
      .style('fill', function(d: any) {
        if (d.node_type === 'metabolite') {
          const should_color_data = (has_data_on_nodes &&
                                 metabolite_data_styles.indexOf('color') !== -1);
          if (should_color_data) {
            const f = d.data;
            return f == null ? no_data_style['color'] : scale.metabolite_color(f);
          } else
            return null;
        }
        // midmarkers and multimarkers
        return null;
      })
      // @ts-ignore
      .call(this.behavior.turnOffDrag)
      .call(drag_behavior)
      .on('mousedown', mousedown_fn)
      .on('click', click_fn)
      .on('mouseover', objectMouseover)
      .on('mouseout', objectMouseout);

    // update node label visibility
    const node_label = update_selection
      .select('.node-label')
      .attr('visibility', hide_all_labels ? 'hidden' : 'visible');
    if (!hide_all_labels) {
      node_label
        .style('visibility', function(d: any) {
          return hideNode(d, hide_secondary_metabolites) ? 'hidden' : null;
        })
        .attr('transform', function(d: any) {
          return 'translate(' + d.label_x + ',' + d.label_y + ')';
        })
        .text(function(d: any) {
          let t = d[identifiers_on_map];
          if (has_data_on_nodes && metabolite_data_styles.indexOf('text') !== -1)
            t += ' ' + d.data_string;
          return t;
        })
        // @ts-ignore
        .call(this.behavior.turnOffDrag)
        .call(label_drag_behavior)
        .on('mouseover', labelMouseover)
        .on('mouseout', labelMouseout)
        .on('touchend', labelTouch);
    }

    this.callback_manager.run('update_node', this, update_selection);

    function hideNode(d: any, hide_secondary_metabolites: boolean) {
      return (d.node_type === 'metabolite' &&
            hide_secondary_metabolites &&
            !d.node_is_primary);
    }
  }

  /**
 * Create text labels in the enter_selection.
 * @param {} enter_selection - The D3 enter selection.
 * @returns {} The selection of the new nodes.
 */
  create_text_label(enter_selection: D3Selection) {
    const g = enter_selection.append('g')
      .attr('id', function(d: any) { return 'l' + d.text_label_id; })
      .attr('class', 'text-label');
    g.append('text')
      .attr('class', 'label');

    this.callback_manager.run('create_text_label', this, enter_selection);

    return g;
  }

  update_text_label(update_selection: D3Selection<any>) {
    const mousedown = this.behavior.textLabelMousedown;
    const click = this.behavior.textLabelClick;
    const turnOffDrag = this.behavior.turnOffDrag;
    const drag = this.behavior.selectableDrag;

    update_selection
      .select('.label')
      .text(function(d: any) { return d.text; })
      .attr('transform', function(d: any) {
        return 'translate(' + d.x + ',' + d.y + ')';
      })
      .on('mousedown', mousedown)
      .on('click', click)
      .call(turnOffDrag as any)
      .call(drag as any);

    this.callback_manager.run('update_text_label', this, update_selection);
  }
}

export function displacedCoords(reactionArrowDisplacement: number, start: Coord, end: Coord, displace: string) {
  const length = reactionArrowDisplacement;
  const hyp = utils.distance(start, end);
  if (!length || !hyp) {
    console.warn('No space for displacement');
    return {x: start.x, y: start.y};
  }
  if (displace === 'start') {
    return {
      x: start.x + length * (end.x - start.x) / hyp,
      y: start.y + length * (end.y - start.y) / hyp
    };
  } else if (displace === 'end') {
    return {
      x: end.x - length * (end.x - start.x) / hyp,
      y: end.y - length * (end.y - start.y) / hyp
    };
  } else
    console.error('bad displace value: ' + displace);
}
