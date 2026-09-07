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

/** Parse '#rgb', '#rrggbb', 'rgb(...)' or 'rgba(...)' into [r, g, b]. */
function parseColor(color: string): [number, number, number] | null {
  const rgbMatch = /^rgba?\(([^)]+)\)/.exec(color);
  if (rgbMatch) {
    const parts = rgbMatch[1].split(',').map((x) => parseFloat(x));
    if (parts.length >= 3 && !parts.slice(0, 3).some(isNaN))
      return [parts[0], parts[1], parts[2]];
    return null;
  }
  const hexMatch = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(color.trim());
  if (hexMatch) {
    const h = hexMatch[1];
    if (h.length === 3) {
      return [parseInt(h[0] + h[0], 16), parseInt(h[1] + h[1], 16),
        parseInt(h[2] + h[2], 16)];
    }
    return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16),
      parseInt(h.slice(4, 6), 16)];
  }
  return null;
}

/**
 * Pick a flow-dot color that stays visible on the segment it rides on:
 * a lightened version of the edge color on dark edges, a darkened one on
 * pale edges.
 */
function flowDotColor(edgeColor: string): string {
  const rgb = parseColor(edgeColor);
  if (!rgb) return '#ffffff';
  const [r, g, b] = rgb;
  const luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b;
  const mix = (x: number, target: number, t: number) =>
    Math.round(x + (target - x) * t);
  return luminance < 140 ?
    `rgb(${mix(r, 255, 0.7)}, ${mix(g, 255, 0.7)}, ${mix(b, 255, 0.7)})` :
    `rgb(${mix(r, 0, 0.45)}, ${mix(g, 0, 0.45)}, ${mix(b, 0, 0.45)})`;
}

/** On-screen edge width above which the flux dots get their drop-shadow glow.
 * Compositing the shadow across a whole zoomed-out map lags the browser, and
 * it is only appreciable up close anyway. */
const FLOW_SHADOW_MIN_EDGE_PX = 13;

type FlowAnimated = SVGPathElement & {
  __flowAnim?: {anim: Animation, key: string, segWidth: number, halo: string}
};

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
    // toggle the flux-dot glow as the zoom level changes
    this.map.zoomContainer?.callbackManager.set('zoom_change.flux-flow',
      () => this.updateFlowShadows());
  }

  /**
   * Show the flux-dot drop-shadow glow only when zoomed in close enough for
   * it to matter; compositing it across a zoomed-out map lags the browser.
   */
  updateFlowShadows() {
    const scale = this.map.zoomContainer?.windowScale ?? 1;
    const svgNode = this.map.svg?.node() as Element | null;
    if (!svgNode) return;
    svgNode.querySelectorAll('.segment-flow').forEach((node) => {
      const el = node as FlowAnimated;
      if (!el.__flowAnim) return;
      el.style.filter = el.__flowAnim.segWidth * scale >= FLOW_SHADOW_MIN_EDGE_PX ?
        el.__flowAnim.halo : '';
    });
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

    // flux-direction animation overlay (dots drifting along the segment)
    g.append('path')
      .attr('class', 'segment-flow');

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

    // Flux-direction animation overlay: translucent dots drifting along each
    // segment in the direction the flux actually flows (see .segment-flow CSS).
    // Shown only when reaction data is loaded and the segment carries a
    // nonzero flux; speed tracks the flux magnitude relative to the scale.
    // Driven by the Web Animations API rather than CSS keyframes: SVG stroke
    // properties reject unitless var()/calc() keyframe values in some browsers,
    // which silently kills the animation.
    const animate_flux = this.settings.get('animate_flux') &&
      typeof Element !== 'undefined' && 'animate' in Element.prototype &&
      !(typeof window !== 'undefined' &&
        window.matchMedia?.('(prefers-reduced-motion: reduce)').matches);
    let max_abs_domain = 0;
    if (has_data_on_reactions && animate_flux) {
      const domain = scale.reaction_color.domain() as number[];
      max_abs_domain = domain.reduce((m, x) => Math.max(m, Math.abs(x)), 0);
    }
    const window_scale = this.map?.zoomContainer?.windowScale ?? 1;
    update_selection
      .select('.segment-flow')
      .each(function(d: any) {
        // @ts-ignore
        // eslint-disable-next-line no-invalid-this
        const el = this as FlowAnimated;
        const flux = d.data;
        const segEl = el.parentElement?.querySelector('.segment') as SVGPathElement | null;
        const show = Boolean(animate_flux && has_data_on_reactions && segEl &&
          flux != null && !isNaN(flux) && flux !== 0 &&
          segEl.style.visibility !== 'hidden');
        if (!show) {
          if (el.__flowAnim) {
            el.__flowAnim.anim.cancel();
            delete el.__flowAnim;
          }
          el.style.display = 'none';
          el.removeAttribute('d');
          return;
        }

        // Does forward (positive) flux run from -> to along this path?
        let forwardAlong: boolean;
        if (d.to_node_coefficient != null) {
          forwardAlong = d.to_node_coefficient > 0; // products receive flux, reactants emit it
        } else if (d.from_node_coefficient != null) {
          forwardAlong = d.from_node_coefficient < 0;
        } else {
          // marker-to-marker segment: which side of the reaction the non-center
          // marker sits on is only recorded on its sibling metabolite segments
          const reaction = (el.parentNode!.parentNode as any).__data__;
          const fromNode = drawn_nodes[d.from_node_id];
          const markerId = (fromNode && fromNode.node_type === 'midmarker') ?
            d.to_node_id : d.from_node_id;
          const markerIsFrom = markerId === d.from_node_id;
          let side = 0; // <0: reactant side, >0: product side
          const segments = reaction ? reaction.segments : {};
          for (const sid in segments) {
            const s = segments[sid];
            if (s.from_node_id === markerId && s.to_node_coefficient != null) {
              side = s.to_node_coefficient;
              break;
            }
            if (s.to_node_id === markerId && s.from_node_coefficient != null) {
              side = s.from_node_coefficient;
              break;
            }
          }
          if (side === 0) {
            // exchange/demand reactions: this marker has no metabolite segments
            // at all, so it sits opposite every metabolite-bearing segment
            for (const sid in segments) {
              const s = segments[sid];
              const coef = s.to_node_coefficient != null ?
                s.to_node_coefficient : s.from_node_coefficient;
              if (coef != null) {
                side = -coef;
                break;
              }
            }
          }
          // a reactant-side marker feeds the center; a product-side marker drains it
          forwardAlong = side > 0 ? !markerIsFrom : markerIsFrom;
        }
        // d.data holds |flux| when the 'abs' style is on, so the sign lives in reverse_flux
        const along = d.reverse_flux ? !forwardAlong : forwardAlong;

        // dots sized to sit inside the segment stroke
        let segWidth = 10; // default .segment stroke-width from the embedded css
        if (should_size) {
          segWidth = flux == null ? no_data_size : scale.reaction_size(flux);
          if (isNaN(segWidth)) segWidth = no_data_size;
        }
        const dotWidth = Math.max(2, segWidth * 0.8);
        const period = Math.max(20, dotWidth * 3);

        // one period per 0.3s at the top of the scale, easing down to 1.5s
        let duration = 0.6;
        if (max_abs_domain > 0) {
          const rel = Math.min(1, Math.abs(flux) / max_abs_domain);
          duration = rel <= 0 ? 1.5 : Math.min(1.5, Math.max(0.3, 0.3 / Math.sqrt(rel)));
        }

        el.style.display = 'inline';
        el.setAttribute('d', segEl!.getAttribute('d') ?? '');
        el.style.strokeWidth = dotWidth + 'px';
        // dots take a lightened/darkened shade of this segment's own color,
        // with a soft same-color halo around each bead
        const dotColor = flowDotColor(getComputedStyle(segEl!).stroke);
        el.style.stroke = dotColor;
        // stacking shadows intensifies the glow; a single larger radius only dilutes it.
        // Applied only at high zoom (see updateFlowShadows) — it lags the browser otherwise.
        const haloOne = `drop-shadow(0 0 ${(dotWidth * 0.6).toFixed(1)}px ${dotColor})`;
        const halo = `${haloOne} ${haloOne}`;
        el.style.filter = segWidth * window_scale >= FLOW_SHADOW_MIN_EDGE_PX ? halo : '';
        // a dot every `period` units: a tiny dash with round caps renders as a circle
        const dash = 0.1;
        el.style.strokeDasharray = dash + 'px ' + period + 'px';

        // (re)start the animation only when its parameters actually changed,
        // so redraws (e.g. time-course steps) don't visibly reset the phase
        const key = [period, duration.toFixed(2), along].join('|');
        if (!el.__flowAnim || el.__flowAnim.key !== key) {
          if (el.__flowAnim)
            el.__flowAnim.anim.cancel();
          // decreasing dashoffset slides the dots from the path start toward its end
          const anim = el.animate(
            [
              {strokeDashoffset: (period + dash) + 'px'},
              {strokeDashoffset: '0px'}
            ],
            {
              duration: duration * 1000,
              iterations: Infinity,
              direction: along ? 'normal' : 'reverse',
              easing: 'linear'
            }
          );
          el.__flowAnim = {anim, key, segWidth, halo};
        } else {
          el.__flowAnim.segWidth = segWidth;
          el.__flowAnim.halo = halo;
        }
      });

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
