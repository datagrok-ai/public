/* eslint-disable camelcase */
/**
 * For documentation of this class, see docs/javascript_api.rst
 */

/** @jsx h */
import * as utils from './ts/utils';
import BuildInput from './ts/BuildInput';
import ZoomContainer from './ts/ZoomContainer';
import {EscherMap} from './ts/escherMap';
import CobraModel from './ts/CobraModel';
import Brush from './ts/Brush';
import CallbackManager from './ts/CallbackManager';
import Settings from './ts/Settings';
import TextEditInput from './ts/TextEditInput';
import * as dataStyles from './ts/dataStyles';
import renderWrapper from './renderWrapper';
import SettingsMenu from './SettingsMenu';
import MenuBar from './MenuBar';
import SearchBar from './SearchBar';
import ButtonPanel from './ButtonPanel';
import TooltipContainer from './TooltipContainer';
import DistributionTooltipComponent from './DistributionTooltip';
import _ from 'underscore';
import {
  select as d3Select,
  selection as d3Selection
} from 'd3-selection';

// Include custom font set for icons
import '../icons/css/fontello.css';

// Include GUI CSS normally with webpack
import './css/Builder.css';

// Import CSS as a string to embed. This also works from lib because css/src get
// uploaded to NPM.
// @ts-ignore
import builderEmbed from '!!raw-loader!./css/Builder-embed.css';
import {CobraModelData, D3Selection, MapData, ReactionSamplingDistribution, SamplingFunctionResult, SettingsType} from './ts/types';
import {Ref} from 'preact';

export type BuilderType = Builder;
export type BuilderConstructor = typeof Builder;
export type BuilderSaveState = ReturnType<typeof Builder.prototype.getSavingState>

class Builder {
  reaction_sampling_distribution?: ReactionSamplingDistribution | null = null;

  map_data: MapData | null = null;
  model_data: CobraModelData | null = null;
  embeddedCss: string | null = null;
  selection: D3Selection | null = null;
  menu_div: D3Selection | null = null;
  button_div: D3Selection | null = null;
  search_bar_div: D3Selection | null = null;
  searchBarRef: D3Selection | null = null;
  semanticOptions: any | null = null;
  mode: string | null = null;
  settings: Settings | null = null;
  isFullScreen: boolean | null = null;
  savedFullScreenSettings: any | null = null;
  savedFullScreenParent: any | null = null;
  clearFullScreenEscape: any | null = null;
  callback_manager: CallbackManager | null = null;
  has_custom_reaction_styles: boolean;
  map: EscherMap | null = null;
  mapToolsContainer: D3Selection | null = null;
  zoom_container: ZoomContainer | null = null;
  tooltip_container: TooltipContainer | null = null;
  cobra_model: CobraModel | null = null;
  build_input: BuildInput;
  text_edit_input: TextEditInput;
  brush: Brush;
  settingsMenuRef: Ref<SettingsMenu> | null = null;
  menuBarRef: Ref<MenuBar> | null = null;
  buttonPanelRef: Ref<ButtonPanel> | null = null;
  update_model_timer: any;
  status_bar: any;

  constructor(mapData: MapData, modelData: CobraModelData | null, embeddedCss: string, selection: D3Selection, options: Partial<SettingsType>) {
    // Defaults
    if (!selection)
      selection = d3Select('body').append('div');
    else if (selection instanceof d3Selection) {
      // D3 V4 selection
    } else if ('node' in selection) {
      // If user passes in a selection from an different d3 version/instance,
      // then reselect.
      selection = d3Select(selection.node());
    } else {
      // HTML Element
      selection = d3Select(selection);
    }
    if (!options)
      options = {};

    if (!embeddedCss)
      embeddedCss = builderEmbed;


    this.reaction_sampling_distribution = {
      lower_bound: -25,
      upper_bound: 25,
      data: new Map()
    };

    this.map_data = mapData;
    this.model_data = modelData;
    this.embeddedCss = embeddedCss;
    this.selection = selection;
    this.menu_div = null;
    this.button_div = null;
    this.search_bar_div = null;
    this.searchBarRef = null;
    this.semanticOptions = null;
    this.mode = 'zoom';

    // apply this object as data for the selection
    this.selection.datum(this);
    // @ts-ignore
    this.selection.__builder__ = this;

    // Remember if the user provided a custom value for reaction_styles
    this.has_custom_reaction_styles = Boolean(options.reaction_styles);

    // set defaults
    const optionsWithDefaults: SettingsType = utils.set_options(options, {
      // view options
      menu: 'all',
      pathFindingDisabled: false,
      saveAction: null,
      loadAction: null,
      scroll_behavior: 'pan',
      use_3d_transform: false,
      highlight_same_nodes: false,
      enable_editing: true,
      enable_keys: true,
      enable_search: true,
      fill_screen: false,
      zoom_to_element: null,
      full_screen_button: false,
      ignore_bootstrap: false,
      disabled_buttons: null,
      semantic_zoom: null,
      // map, model, and styles
      starting_reaction: null,
      never_ask_before_quit: false,
      unique_map_id: null, // deprecated
      primary_metabolite_radius: 20,
      secondary_metabolite_radius: 10,
      marker_radius: 5,
      gene_font_size: 18,
      hide_secondary_metabolites: false,
      show_gene_reaction_rules: false,
      hide_all_labels: false,
      canvas_size_and_loc: null,
      // applied data
      // reaction
      reaction_data: null,
      reaction_styles: ['color', 'size', 'text', 'abs'],
      reaction_compare_style: 'log2_fold',
      reaction_scale: null,
      reaction_scale_preset: 'GaBuGeRd',
      reaction_no_data_color: '#dcdcdc',
      reaction_no_data_size: 8,
      // samplingFunction
      samplingFunction: null,
      // gene
      gene_data: null,
      and_method_in_gene_reaction_rule: 'mean',
      // metabolite
      metabolite_data: null,
      metabolite_styles: ['color', 'size', 'text'],
      metabolite_compare_style: 'log2_fold',
      metabolite_scale: null,
      metabolite_scale_preset: 'WhYlRd',
      metabolite_no_data_color: '#ffffff',
      metabolite_no_data_size: 10,
      // View and build options
      identifiers_on_map: 'bigg_id',
      highlight_missing: false,
      allow_building_duplicate_reactions: false,
      cofactors: [
        'atp', 'adp', 'nad', 'nadh', 'nadp', 'nadph', 'gtp', 'gdp', 'h', 'coa',
        'ump', 'h2o', 'ppi'
      ],
      // Extensions
      tooltip_component: DistributionTooltipComponent,
      enable_tooltips: ['label'],
      enable_keys_with_tooltip: true,
      // Callbacks
      first_load_callback: null
    } as SettingsType, {
      primary_metabolite_radius: true,
      secondary_metabolite_radius: true,
      marker_radius: true,
      gene_font_size: true,
      reaction_no_data_size: true,
      metabolite_no_data_size: true
    });

    // Check the location
    if (utils.check_for_parent_tag(this.selection, 'svg')) {
      throw new Error('Builder cannot be placed within an svg node ' +
                      'because UI elements are html-based.');
    }

    // The options that are erased when the settings menu is canceled
    const conditional = [
      'identifiers_on_map',
      'scroll_behavior',
      'hide_secondary_metabolites',
      'show_gene_reaction_rules',
      'hide_all_labels',
      'allow_building_duplicate_reactions',
      'highlight_missing',
      'enable_tooltips',
      'reaction_scale_preset',
      'reaction_no_data_color',
      'reaction_no_data_size',
      'reaction_scale',
      'reaction_styles',
      'reaction_compare_style',
      'and_method_in_gene_reaction_rule',
      'metabolite_scale_preset',
      'metabolite_scale',
      'metabolite_styles',
      'metabolite_compare_style',
      'metabolite_no_data_color',
      'metabolite_no_data_size',
      'samplingFunction',
    ];
    // this.options and this.settings used to have different functions, but now
    // they are aliases
    this.settings = new Settings(optionsWithDefaults, conditional as (keyof SettingsType)[]);

    // Warn if full/fill screen options conflict
    if (this.settings.get('fill_screen') && this.settings.get('full_screen_button')) {
      this.settings.set('full_screen_button', false);
      console.warn('The option full_screen_button has no effect when fill_screen is true');
    }

    // force full screen for fill_screen option
    this.isFullScreen = false;
    if (this.settings.get('fill_screen')) {
      d3Select('html').classed('fill-screen', true);
      d3Select('body').classed('fill-screen', true);
      this.selection.classed('fill-screen-div', true);
      this.isFullScreen = true;
    }
    this.savedFullScreenSettings = null;
    this.savedFullScreenParent = null;
    this.clearFullScreenEscape = null;

    // Set up this callback manager
    this.callback_manager = new CallbackManager();
    const firstLoadCallback = this.settings.get('first_load_callback');
    this.callback_manager.set('first_load', () => {
      if (firstLoadCallback !== null)
        firstLoadCallback(this);

      let oldModel = JSON.parse(JSON.stringify(this.model_data)); // need to keep this for resetting
      const getTooltipProps = () => {
        const sliderChangeFunc = (a: number[], b: string) => {
          if (!a || !b || !Array.isArray(a) || typeof b !== 'string' || !this.model_data || !this.model_data.reactions)
            return;
          const lower_bound = a[0];
          const upper_bound = a[1];
          const id = b;
          for (let i = 0; i < this.model_data.reactions.length; i++) {
            if (this.model_data.reactions[i].id === id) {
              this.model_data.reactions[i].lower_bound = lower_bound;
              this.model_data.reactions[i].upper_bound = upper_bound;
              break;
            }
          }
        };

        const debouncedSliderAction = _.debounce(sliderChangeFunc, 300);


        return {map: this.map, model: this.model_data, sliderChange: debouncedSliderAction, lowerRange: -25,
          upperRange: 25, step: 0.01, oldModel: oldModel, objectives: {}, samplingFunction: this.settings.get('samplingFunction') ? this.sampleReactions.bind(this) : null,
          fluxDistributions: this.reaction_sampling_distribution, reactionData: this.settings.get('reaction_data')
        };
      };
      // Set up the tooltip container


      this.tooltip_container.passProps(getTooltipProps()); // reactionData misssing
      this.callback_manager.set('load_map', (mpData) => {
        this.set_reaction_data(null);
        oldModel = JSON.parse(JSON.stringify(this.model_data));
        this.tooltip_container.passProps(getTooltipProps());
      });
      this.callback_manager.set('load_model', (modelData) => {
        this.set_reaction_data(null);
        oldModel = JSON.parse(JSON.stringify(this.model_data));
        this.tooltip_container.passProps(getTooltipProps());
      });
    });


    // Set up the zoom container
    this.zoom_container = new ZoomContainer(this.selection,
      this.settings.get('scroll_behavior'),
      this.settings.get('use_3d_transform'));
    // Zoom container status changes
    // this.zoom_container.callbackManager.set('svg_start', () => {
    //   if (this.map) this.map.set_status('Drawing ...')
    // })
    // this.zoom_container.callbackManager.set('svg_finish', () => {
    //   if (this.map) this.map.set_status('')
    // })
    this.zoom_container.callbackManager.set('zoom_change', () => {
      if (this.settings.get('semantic_zoom')) {
        const scale = this.zoom_container.windowScale;
        const optionObject = this.settings.get('semantic_zoom')
          .sort((a, b) => a.zoomLevel - b.zoomLevel)
          .find((a) => a.zoomLevel > scale);
        if (optionObject) {
          let didChange = false;
          _.mapObject(optionObject.options, (value, key) => {
            if (this.settings.get(key) !== value) {
              this.settings.set(key, value);
              didChange = true;
            }
          });
          if (didChange) this._updateData(false, true);
        }
      }
    });
    this.settings.streams.use_3d_transform.onValue((val) => {
      this.zoom_container.setUse3dTransform(val);
    });
    this.settings.streams.scroll_behavior.onValue((val) => {
      this.zoom_container.setScrollBehavior(val);
    });

    // Reactive tooltip settings
    this.settings.streams.enable_tooltips.onValue((val) => {
      this._updateTooltipSetting(val);
    });

    // Make a container for other map-related tools that will be reset on map load
    // TODO only create these once in the Builder constructor
    this.mapToolsContainer = this.selection.append('div')
      .attr('class', 'map-tools-container');

    // Status in both modes
    this._createStatus(this.selection);

    // Load the model, map, and update data in both
    this.load_model(this.model_data, false);

    // Append the bars and menu divs to the document
    const s = this.selection
      .append('div').attr('class', 'search-menu-container')
      .append('div').attr('class', 'search-menu-container-inline');
    this.menu_div = s.append('div');
    this.search_bar_div = s.append('div');
    this.button_div = this.selection.append('div');

    // Need to defer map loading to let webpack CSS load properly. Hack:
    // Delaying 50ms to make sure the css calculations on map size take
    // place.
    _.delay(() => {
      this.load_map(this.map_data, false);

      const messageFn = this._reactionCheckAddAbs();
      this._updateData(true, true);

      // Setting callbacks. TODO enable atomic updates. Right now, every time the
      // menu closes, everything is drawn.
      this.settings.statusBus.onValue((x) => {
        if (x === 'accept') {
          this._updateData(true, true, ['reaction', 'metabolite'], false);
          if (this.zoom_container !== null) {
            // TODO make this automatic
            const newBehavior = this.settings.get('scroll_behavior');
            this.zoom_container.setScrollBehavior(newBehavior);
          }
          if (this.map !== null) {
            this.map.draw_all_nodes(false);
            this.map.draw_all_reactions(true, false);
            this.map.select_none();
          }
        }
      });

      if (messageFn !== null) setTimeout(messageFn, 500);

      // Finally run callback
      _.defer(() => this.callback_manager.run('first_load', this));
    }, 50);
  }

  sampleReactions() {
    const samplingFunction = this.settings.get('samplingFunction');
    if (samplingFunction) {
      //  function should return same signature as this.reaction_sampling_distribution
      const data = samplingFunction(this.model_data);

      const dataReceived = (d: SamplingFunctionResult) => {
        if (!d)
          throw new Error('Sampling function did not return any data');
        this.reaction_sampling_distribution.lower_bound = d.lower_bound;
        this.reaction_sampling_distribution.upper_bound = d.upper_bound;
        this.reaction_sampling_distribution.data.clear();
        for (const [key, value] of d.data)
          this.reaction_sampling_distribution.data.set(key, value);
      };
      const isPromise = data instanceof Promise;
      if (isPromise)
        data.then(dataReceived);
      else
        dataReceived(data);
    }
  }

  // builder.options is deprecated
  get options() {
    throw new Error('builder.options is deprecated. Use builder.settings.get() ' +
                    'and builder.settings.set() instead.');
  }
  set options(_) {
    throw new Error('builder.options is deprecated. Use builder.settings.get() ' +
                    'and builder.settings.set() instead.');
  }

  /**
   * For documentation of this function, see docs/javascript_api.rst.
   */
  load_model(modelData, shouldUpdateData = true) { // eslint-disable-line camelcase
    this.model_data = modelData;
    this.reaction_sampling_distribution.data.clear();

    // Check the cobra model
    if (_.isNull(modelData))
      this.cobra_model = null;
    else
      this.cobra_model = CobraModel.from_cobra_json(modelData);


    if (this.map) {
      this.map.cobra_model = this.cobra_model;
      if (shouldUpdateData)
        this._updateData(true, false);

      if (this.settings.get('highlight_missing'))
        this.map.draw_all_reactions(false, false);
    }

    this.callback_manager.run('load_model', null, modelData, shouldUpdateData);
  }

  /**
   * For documentation of this function, see docs/javascript_api.rst
   */
  load_map(mapData, shouldUpdateData = true) { // eslint-disable-line camelcase
    // Store map options that might be changed by semantic_zoom function
    const tempSemanticOptions = {};
    if (this.settings.get('semantic_zoom')) {
      for (const level of this.settings.get('semantic_zoom')) {
        Object.keys(level.options).map((option) => {
          if (tempSemanticOptions[option] === undefined)
            tempSemanticOptions[option] = this.settings.get(option);
        });
      }
      this.semanticOptions = Object.assign({}, tempSemanticOptions);
    }

    // remove the old map and related divs
    utils.remove_child_nodes(this.zoom_container.zoomedSel);
    utils.remove_child_nodes(this.mapToolsContainer);

    const zoomedSel = this.zoom_container.zoomedSel;
    const svg = this.zoom_container.svg;

    // remove the old map side effects
    if (this.map)
      this.map.key_manager.toggle(false);


    if (mapData !== null) {
      // import map
      this.map = EscherMap.from_data(mapData,
        svg,
        this.embeddedCss,
        zoomedSel,
        this.zoom_container,
        this.settings,
        this.cobra_model,
        this.settings.get('enable_search'));
    } else {
      // new map
      this.map = new EscherMap(svg,
        this.embeddedCss,
        zoomedSel,
        this.zoom_container,
        this.settings,
        this.cobra_model,
        this.settings.get('canvas_size_and_loc'),
        this.settings.get('enable_search'));
    }

    // Connect status bar
    this._setupStatus(this.map);
    this.map.set_status('Loading map ...');

    // Connect tooltips
    this._updateTooltipSetting(this.settings.get('enable_tooltips'));

    // Set the data for the map
    if (shouldUpdateData)
      this._updateData(false, true);


    // Set up the reaction input with complete.ly
    this.build_input = new BuildInput(this.mapToolsContainer, this.map,
      this.zoom_container, this.settings);

    // Set up the text edit input
    this.text_edit_input = new TextEditInput(this.mapToolsContainer, this.map,
      this.zoom_container);

    // Set up the Brush
    this.brush = new Brush(zoomedSel, false, this.map, '.canvas-group');
    // reset brush when canvas resizes in brush mode
    this.map.canvas.callbackManager.set('resize', () => {
      if (this.mode === 'brush') this.brush.toggle(true);
    });

    // Set up menus
    this.setUpSettingsMenu(this.mapToolsContainer);
    this.setUpButtonPanel(this.mapToolsContainer);

    // share a parent container for menu bar and search bar
    const sel = this.mapToolsContainer
      .append('div').attr('class', 'search-menu-container')
      .append('div').attr('class', 'search-menu-container-inline');
    this.setUpMenuBar(sel);
    this.setUpSearchBar(sel);

    // Set up the tooltip container
    this.tooltip_container = new TooltipContainer(
      this.mapToolsContainer,
      this.settings.get('tooltip_component'),
      this.zoom_container,
      this.map,
      this.settings
    );

    // Set up key manager
    this.map.key_manager.assignedKeys = this.getKeys();
    // Tell the key manager about the reaction input and search bar
    this.map.key_manager.inputList = [
      this.build_input,
      this.searchBarRef,
      () => this.settingsMenuRef,
      this.text_edit_input
    ];
    if (!this.settings.get('enable_keys_with_tooltip'))
      this.map.key_manager.inputList.push(this.tooltip_container);

    // Make sure the key manager remembers all those changes
    this.map.key_manager.update();
    // Turn it on/off
    this.map.key_manager.toggle(this.settings.get('enable_keys'));
    this.settings.streams.enable_keys.onValue((val) => {
      // get keys given latest settings
      this.map.key_manager.toggle(val);
    });

    // Disable clears
    const newDisabledButtons = this.settings.get('disabled_buttons') || [];
    if (!this.settings.get('reaction_data'))
      newDisabledButtons.push('Clear reaction data');

    if (!this.settings.get('gene_data'))
      newDisabledButtons.push('Clear gene data');

    if (!this.settings.get('metabolite_data'))
      newDisabledButtons.push('Clear metabolite data');

    if (!this.settings.get('enable_search'))
      newDisabledButtons.push('Find');

    if (!this.settings.get('enable_editing'))
      newDisabledButtons.push('Show control points');

    this.settings.set('disabled_buttons', newDisabledButtons);

    // Set up selection box
    if (this.settings.get('zoom_to_element')) {
      const type = this.settings.get('zoom_to_element').type;
      const elementId = this.settings.get('zoom_to_element').id;
      if (_.isUndefined(type) || ['reaction', 'node'].indexOf(type) === -1)
        throw new Error('zoom_to_element type must be "reaction" or "node"');

      if (_.isUndefined(elementId))
        throw new Error('zoom_to_element must include id');

      if (type === 'reaction')
        this.map.zoom_to_reaction(elementId);
      else if (type === 'node')
        this.map.zoom_to_node(elementId);
    } else if (mapData)
      this.map.zoom_extent_canvas();
    else {
      if (this.settings.get('starting_reaction') && this.cobra_model !== null) {
        // Draw default reaction if no map is provided
        const size = this.zoom_container.get_size();
        const startCoords = {x: size.width / 2, y: size.height / 4};
        this.map.new_reaction_from_scratch(this.settings.get('starting_reaction'),
          startCoords, 90);
        this.map.zoom_extent_nodes();
      } else
        this.map.zoom_extent_canvas();
    }

    // Start in zoom mode for builder, view mode for viewer
    if (this.settings.get('enable_editing'))
      this.zoom_mode();
    else
      this.view_mode();

    // when enabled_editing changes, go to view mode
    this.settings.streams.enable_editing.onValue((val) => {
      if (val) this.zoom_mode();
      else this.view_mode();
    });

    // confirm before leaving the page
    if (this.settings.get('enable_editing'))
      this._setupConfirmBeforeExit();


    // draw
    this.map.draw_everything();

    this.map.set_status('');

    this.callback_manager.run('load_map', null, mapData, shouldUpdateData);
  }

  /**
   * Function to pass props for the settings menu. Run without an argument to
   * rerender the component
   * @param {Object} props - Props that the settings menu will use
   */
  passPropsSettingsMenu(props = {}) {
    this.map.callback_manager.run('pass_props_settings_menu', null, props);
  }

  /**
   * Initialize the settings menu
   */
  setUpSettingsMenu(sel) {
    this.settingsMenuRef = null;
    renderWrapper(
      SettingsMenu,
      (instance: Ref<SettingsMenu>) => { this.settingsMenuRef = instance; },
      (passProps) => this.map.callback_manager.set('pass_props_settings_menu', passProps),
      sel.append('div').node()
    );
    this.passPropsSettingsMenu({
      display: false,
      settings: this.settings,
      map: this.map
    });

    // redraw menu when settings change
    _.mapObject(this.settings.streams, (stream, key) => {
      stream.onValue((value) => {
        this.passPropsSettingsMenu();
      });
    });

    // recalculate data when switching to/from absolute value
    this.settings.streams.reaction_styles
      .map((x) => _.contains(x, 'abs'))
      .skipDuplicates()
      .onValue(() => this._updateData(false, true));
    this.settings.streams.metabolite_styles
      .map((x) => _.contains(x, 'abs'))
      .skipDuplicates()
      .onValue(() => this._updateData(false, true));
  }

  /**
   * Function to pass props for the menu bar
   * @param {Object} props - Props that the menu bar will use
   */
  passPropsMenuBar(props = {}) {
    this.map.callback_manager.run('pass_props_menu_bar', null, props);
  }

  /**
   * Initialize the menu bar
   * @param {D3 Selection} sel - The d3 selection to render in.
   */
  setUpMenuBar(sel) {
    this.menuBarRef = null;
    renderWrapper(
      MenuBar,
      (instance) => { this.menuBarRef = instance; },
      (passProps) => this.map.callback_manager.set('pass_props_menu_bar', passProps),
      sel.append('div').node()
    );
    this.passPropsMenuBar({
      display: this.settings.get('menu') === 'all',
      settings: this.settings,
      sel: this.selection,
      mode: this.mode,
      map: this.map,
      saveMap: () => {
        // Revert options changed by semanticZoom to their original values if option is active
        if (this.semanticOptions) {
          Object.entries(this.semanticOptions).map(([key, value]) => {
            this.settings.set(key, value);
          });
          this._updateData();
        }
        this.map.save();
      },
      loadMap: (file) => this.load_map(file),
      assignKeyLoadMap: (fn) => {
        // connect the key for this input
        this.map.key_manager.assignedKeys.load_map.fn = fn;
      },
      save_svg: () => this.map.save_svg(),
      save_png: () => this.map.save_png(),
      clear_map: () => { this.clear_map(); },
      loadModel: (file) => this.load_model(file, true),
      assignKeyLoadModel: (fn) => {
        // connect the key for this input
        this.map.key_manager.assignedKeys.load_model.fn = fn;
      },
      clearModel: () => {
        this.load_model(null);
        this.callback_manager.run('clear_model');
      },
      updateRules: () => this.map.convert_map(),
      setReactionData: (d: { [key: string]: (number | string | null) } | null) => this.set_reaction_data(d),
      clearReactionData: () => this.set_reaction_data(null),
      setGeneData: (d) => this.set_gene_data(d),
      clearGeneData: () => this.set_gene_data(null, true),
      setMetaboliteData: (d) => this.set_metabolite_data(d),
      clearMetaboliteData: (d) => this.set_metabolite_data(null),
      setMode: (mode: string) => this._setMode(mode),
      deleteSelected: () => this.map.delete_selected(),
      undo: () => this.map.undo_stack.undo(),
      redo: () => this.map.undo_stack.redo(),
      align_vertical: () => this.map.align_vertical(),
      align_horizontal: () => this.map.align_horizontal(),
      togglePrimary: () => this.map.toggle_selected_node_primary(),
      cyclePrimary: () => this.map.cycle_primary_node(),
      selectAll: () => this.map.select_all(),
      selectNone: () => this.map.select_none(),
      invertSelection: () => this.map.invert_selection(),
      zoom_in: () => this.zoom_container.zoom_in(),
      zoom_out: () => this.zoom_container.zoom_out(),
      zoomExtentNodes: () => this.map.zoom_extent_nodes(),
      zoomExtentCanvas: () => this.map.zoom_extent_canvas(),
      full_screen: () => this.full_screen(),
      search: () => this.passPropsSearchBar({display: true}),
      toggleBeziers: () => this.map.toggle_beziers(),
      renderSettingsMenu: () => this.passPropsSettingsMenu({display: true})
    });

    // redraw when beziers change
    this.map.callback_manager.set('toggle_beziers', () => {
      this.passPropsMenuBar();
    });

    // redraw when disabledButtons change
    this.settings.streams.disabled_buttons.onValue((value) => {
      this.passPropsMenuBar();
    });

    // redraw when mode changes
    this.callback_manager.set('set_mode', (mode) => {
      this.passPropsMenuBar({mode});
    });

    // redraw when menu option changes
    this.settings.streams.menu.onValue((menu) => {
      this.passPropsMenuBar({display: menu === 'all'});
    });

    this.settings.streams.saveAction?.onValue((value) => {
      this.passPropsMenuBar();
    });
    this.settings.streams.loadAction?.onValue((value) => {
      this.passPropsMenuBar();
    });
    this.settings.streams.pathFindingDisabled?.onValue((value) => {
      this.passPropsButtonPanel();
    });

    // redraw when full screen button changes
    this.settings.streams.full_screen_button.onValue((value) => {
      this.passPropsMenuBar();
    });
  }

  /**
   * Function to pass props for the search bar
   * @param {Object} props - Props that the search bar will use
   */
  passPropsSearchBar(props = {}) {
    this.map.callback_manager.run('pass_props_search_bar', null, props);
  }

  /**
   * Initialize the search bar
   * @param {D3 Selection} sel - The d3 selection to render in.
   */
  setUpSearchBar(sel: D3Selection) {
    this.searchBarRef = null;
    renderWrapper(
      SearchBar,
      (instance) => { this.searchBarRef = instance; },
      (passProps) => this.map.callback_manager.set('pass_props_search_bar', passProps),
      sel.append('div').node()
    );
    this.passPropsSearchBar({
      display: false,
      searchIndex: this.map.search_index,
      map: this.map
    });
  }

  /**
   * Function to pass props for the button panel
   * @param {Object} props - Props that the tooltip will use
   */
  passPropsButtonPanel(props = {}) {
    this.map.callback_manager.run('pass_props_button_panel', null, props);
  }

  /**
   * Initialize the button panel
   */
  setUpButtonPanel(sel: D3Selection) {
    renderWrapper(
      ButtonPanel,
      null,
      (passProps) => this.map.callback_manager.set('pass_props_button_panel', passProps),
      sel.append('div').node()
    );
    this.passPropsButtonPanel({
      display: _.contains(['all', 'zoom'], this.settings.get('menu')),
      mode: this.mode,
      settings: this.settings,
      setMode: (mode: string) => this._setMode(mode),
      zoomContainer: this.zoom_container,
      map: this.map,
      buildInput: this.build_input,
      full_screen: () => this.full_screen(),
    });

    // redraw when mode changes
    this.callback_manager.set('set_mode', (mode) => {
      this.passPropsButtonPanel({mode});
    });

    // redraw when full screen button changes
    this.settings.streams.full_screen_button.onValue((value) => {
      this.passPropsButtonPanel();
    });
  }

  /**
   * Set the mode
   */
  _setMode(mode: string) {
    this.mode = mode;

    // input
    this.build_input.toggle(mode === 'build');
    this.build_input.direction_arrow.toggle(mode === 'build');
    // brush
    this.brush.toggle(mode === 'brush');
    // zoom
    this.zoom_container.togglePanDrag(mode === 'zoom' || mode === 'view');
    // resize canvas
    this.map.canvas.toggleResize(mode !== 'view');

    // Behavior. Be careful of the order becuase rotation and
    // toggle_selectable_drag both use Behavior.selectableDrag.
    if (mode === 'rotate') {
      this.map.behavior.toggleSelectableDrag(false); // before toggle_rotation_mode
      this.map.behavior.toggleRotationMode(true); // XX
    } else {
      this.map.behavior.toggleRotationMode(mode === 'rotate'); // before toggleSelectableDrag
      this.map.behavior.toggleSelectableDrag(mode === 'brush'); // XX
    }
    this.map.behavior.toggleSelectableClick(mode === 'build' || mode === 'brush'); // XX
    this.map.behavior.toggleLabelDrag(mode === 'brush'); // XX
    this.map.behavior.toggleTextLabelEdit(mode === 'text'); // XX
    this.map.behavior.toggleBezierDrag(mode === 'brush'); // XX

    // edit selections
    if (mode === 'view' || mode === 'text')
      this.map.select_none();

    if (mode === 'rotate')
      this.map.deselect_text_labels();


    this.map.draw_everything();
    // what's not allowing me to delete this? XX above

    // callback
    this.callback_manager.run('set_mode', null, mode);
  }

  /** For documentation of this function, see docs/javascript_api.rst. */
  view_mode() { // eslint-disable-line camelcase
    this.callback_manager.run('view_mode');
    this._setMode('view');
  }

  /** For documentation of this function, see docs/javascript_api.rst. */
  build_mode() { // eslint-disable-line camelcase
    this.callback_manager.run('build_mode');
    this._setMode('build');
  }

  /** For documentation of this function, see docs/javascript_api.rst. */
  brush_mode() { // eslint-disable-line camelcase
    this.callback_manager.run('brush_mode');
    this._setMode('brush');
  }

  /** For documentation of this function, see docs/javascript_api.rst. */
  zoom_mode() { // eslint-disable-line camelcase
    this.callback_manager.run('zoom_mode');
    this._setMode('zoom');
  }

  /** For documentation of this function, see docs/javascript_api.rst. */
  rotate_mode() { // eslint-disable-line camelcase
    this.callback_manager.run('rotate_mode');
    this._setMode('rotate');
  }

  /** For documentation of this function, see docs/javascript_api.rst. */
  text_mode() { // eslint-disable-line camelcase
    this.callback_manager.run('text_mode');
    this._setMode('text');
  }

  _reactionCheckAddAbs() {
    const currStyle = this.settings.get('reaction_styles');
    if (
      this.settings.get('reaction_data') &&
      !this.has_custom_reaction_styles &&
      !_.contains(currStyle, 'abs')
    ) {
      this.settings.set('reaction_styles', currStyle.concat('abs'));
      return () => {
        this.map.set_status('Visualizing absolute value of reaction data. ' +
                            'Change this option in Settings.', 5000);
      };
    }
    return null;
  }

  /**
   * For documentation of this function, see docs/javascript_api.rst.
   */
  set_reaction_data(data: { [key: string]: (number | string | null) } | null) { // eslint-disable-line camelcase
    this.settings.set('reaction_data', data);

    // clear gene data
    if (data)
      this.settings.options.gene_data = null;


    const messageFn = this._reactionCheckAddAbs();

    this._updateData(true, true, ['reaction']);

    if (messageFn) messageFn();
    else this.map.set_status('');

    const disabledButtons = this.settings.get('disabled_buttons') || [];
    const buttonName = 'Clear reaction data';
    const geneButtonName = 'Clear gene data';
    const index = disabledButtons.indexOf(buttonName);
    if (data && index !== -1) {
      disabledButtons.splice(index, 1);
      const gInd = disabledButtons.indexOf(geneButtonName);
      if (gInd === -1) disabledButtons.push(geneButtonName);
      this.settings.set('disabled_buttons', disabledButtons);
    } else if (!data && index === -1) {
      disabledButtons.push(buttonName);
      this.settings.set('disabled_buttons', disabledButtons);
    }

    this.tooltip_container.passProps({reactionData: data}); // pass new reaction props.
  }

  /**
   * For documentation of this function, see docs/javascript_api.rst.
   */
  set_gene_data(data: any, clearGeneReactionRules: boolean = false) { // eslint-disable-line camelcase
    this.settings.set('gene_data', data);

    if (clearGeneReactionRules)
      this.settings.set('show_gene_reaction_rules', false);


    // clear reaction data; show gene reaction rules
    if (data) {
      this.settings.options.reaction_data = null;
      this.settings.set('show_gene_reaction_rules', true);
    }

    this._updateData(true, true, ['reaction']);
    this.map.set_status('');

    const disabledButtons = this.settings.get('disabled_buttons') || [];
    const index = disabledButtons.indexOf('Clear gene data');
    const buttonName = 'Clear gene data';
    const reactionButtonName = 'Clear reaction data';
    if (index > -1 && data) {
      disabledButtons.splice(index, 1);
      const rInd = disabledButtons.indexOf('Clear reaction data');
      if (rInd === -1) disabledButtons.push(reactionButtonName);
      this.settings.set('disabled_buttons', disabledButtons);
    } else if (index === -1 && !data) {
      disabledButtons.push(buttonName);
      this.settings.set('disabled_buttons', disabledButtons);
    }
  }

  /**
   * For documentation of this function, see docs/javascript_api.rst.
   */
  set_metabolite_data(data: { [key: string]: (number | string | null)[] } | null) { // eslint-disable-line camelcase
    this.settings.set('metabolite_data', data);

    this._updateData(true, true, ['metabolite']);
    this.map.set_status('');

    const disabledButtons = this.settings.get('disabled_buttons') || [];
    const buttonName = 'Clear metabolite data';
    const index = disabledButtons.indexOf(buttonName);
    if (index > -1 && data) {
      disabledButtons.splice(index, 1);
      this.settings.set('disabled_buttons', disabledButtons);
    } else if (index === -1 && !data) {
      disabledButtons.push(buttonName);
      this.settings.set('disabled_buttons', disabledButtons);
    }
  }

  _makeGeneDataObject(geneData: any, cobraModel: CobraModel, map: EscherMap) {
    const allReactions = {};
    if (cobraModel !== null)
      utils.extend(allReactions, cobraModel.reactions);

    // extend, overwrite
    if (map !== null)
      utils.extend(allReactions, map.reactions, true);


    // this object has reaction keys and values containing associated genes
    return dataStyles.importAndCheck(geneData, 'gene_data', allReactions);
  }

  /**
   * Clear the map
   */
  clear_map() { // eslint-disable-line camelcase
    this.callback_manager.run('clear_map');
    this.map.clearMapData();
    this._updateData(true, true, ['reaction', 'metabolite'], false);
    this.map.draw_everything();
  }

  /**
   * Set data and settings for the model.
   * update_model: (Boolean) Update data for the model.
   * update_map: (Boolean) Update data for the map.
   * kind: (Optional, Default: all) An array defining which data is being updated
   * that can include any of: ['reaction', 'metabolite'].
   * should_draw: (Optional, Default: true) Whether to redraw the update sections
   * of the map.
   */
  _updateData(
    updateModel = false,
    updateMap = false,
    kind: string[] = ['reaction', 'metabolite'],
    shouldDraw = true
  ) {
    const updateReactionData = _.contains(kind, 'reaction');
    const updateMetaboliteData = _.contains(kind, 'metabolite');
    let metaboliteDataObject: { [key: string]: (number | string | null)[] } | null = null;
    let reactionDataObject: { [key: string]: (number | string | null)[] } | null = null;
    let geneDataObject: { [key: string]: { [key: string]: (number | string | null)[] } } | null = null;

    // -------------------
    // First map, and draw
    // -------------------

    // metabolite data
    if (updateMetaboliteData && updateMap && this.map !== null) {
      metaboliteDataObject = dataStyles.importAndCheck(this.settings.get('metabolite_data'),
        'metabolite_data');
      this.map.apply_metabolite_data_to_map(metaboliteDataObject);
      if (shouldDraw)
        this.map.draw_all_nodes(false);
    }

    // reaction data
    if (updateReactionData) {
      if (this.settings.get('reaction_data') && updateMap && this.map !== null) {
        reactionDataObject = dataStyles.importAndCheck(this.settings.get('reaction_data'),
          'reaction_data');
        this.map.apply_reaction_data_to_map(reactionDataObject);
        if (shouldDraw)
          this.map.draw_all_reactions(false, false);
      } else if (this.settings.get('gene_data') && updateMap && this.map !== null) {
        geneDataObject = this._makeGeneDataObject(this.settings.get('gene_data'),
          this.cobra_model, this.map);
        this.map.apply_gene_data_to_map(geneDataObject);
        if (shouldDraw)
          this.map.draw_all_reactions(false, false);
      } else if (updateMap && this.map !== null) {
        // clear the data
        this.map.apply_reaction_data_to_map(null);
        if (shouldDraw)
          this.map.draw_all_reactions(false, false);
      }
    }

    // ----------------------------------------------------------------
    // Then the model, after drawing. Delay by 5ms so the the map draws
    // first.
    // ----------------------------------------------------------------

    // If this function runs again, cancel the previous model update
    if (this.update_model_timer)
      clearTimeout(this.update_model_timer);


    const delay = 5;
    this.update_model_timer = setTimeout(() => {
      // metabolite_data
      if (updateMetaboliteData && updateModel && this.cobra_model !== null) {
        // if we haven't already made this
        if (!metaboliteDataObject) {
          metaboliteDataObject = dataStyles.importAndCheck(this.settings.get('metabolite_data'),
            'metabolite_data');
        }
        this.cobra_model.apply_metabolite_data(metaboliteDataObject,
          this.settings.get('metabolite_styles'),
          this.settings.get('metabolite_compare_style'));
      }

      // reaction data
      if (updateReactionData) {
        if (this.settings.get('reaction_data') && updateModel && this.cobra_model !== null) {
          // if we haven't already made this
          if (!reactionDataObject) {
            reactionDataObject = dataStyles.importAndCheck(this.settings.get('reaction_data'),
              'reaction_data');
          }
          this.cobra_model.apply_reaction_data(reactionDataObject,
            this.settings.get('reaction_styles'),
            this.settings.get('reaction_compare_style'));
        } else if (this.settings.get('gene_data') && updateModel && this.cobra_model !== null) {
          if (!geneDataObject) {
            geneDataObject = this._makeGeneDataObject(this.settings.get('gene_data'),
              this.cobra_model, this.map);
          }
          this.cobra_model.apply_gene_data(geneDataObject,
            this.settings.get('reaction_styles'),
            this.settings.get('identifiers_on_map'),
            this.settings.get('reaction_compare_style'),
            this.settings.get('and_method_in_gene_reaction_rule'));
        } else if (updateModel && this.cobra_model !== null) {
          // clear the data
          this.cobra_model.apply_reaction_data(null,
            this.settings.get('reaction_styles'),
            this.settings.get('reaction_compare_style'));
        }
      }

      // callback
      this.callback_manager.run('update_data', null, updateModel, updateMap,
        kind, shouldDraw);
    }, delay);
  }

  _createStatus(selection: D3Selection) {
    this.status_bar = selection.append('div').attr('id', 'status');
  }

  _setupStatus(map: EscherMap) {
    map.callback_manager.set('set_status', (status) => this.status_bar.html(status));
  }

  _updateTooltipSetting(setting: string | string[]) {
    this.map.behavior.toggleLabelMouseover(setting && setting.includes('label'));
    this.map.behavior.toggleObjectMouseover(setting && setting.includes('object'));
  }

  /**
   * Define keyboard shortcuts
   */
  getKeys() {
    const map = this.map;
    const zoom_container = this.zoom_container; // eslint-disable-line camelcase
    return {
      save: {
        key: 'ctrl+s',
        target: map,
        fn: map.save
      },
      save_svg: {
        key: 'ctrl+shift+s',
        target: map,
        fn: map.save_svg
      },
      save_png: {
        key: 'ctrl+shift+p',
        target: map,
        fn: map.save_png
      },
      load_map: {
        key: 'ctrl+o',
        fn: null // defined by button
      },
      convert_map: {
        target: map,
        fn: map.convert_map
      },
      load_model: {
        key: 'ctrl+m',
        fn: null // defined by button
      },
      zoom_in_ctrl: {
        key: 'ctrl+=',
        target: zoom_container,
        fn: zoom_container.zoom_in
      },
      zoom_in: {
        key: '=',
        target: zoom_container,
        fn: zoom_container.zoom_in,
        ignoreWithInput: true
      },
      zoom_out_ctrl: {
        key: 'ctrl+-',
        target: zoom_container,
        fn: zoom_container.zoom_out
      },
      zoom_out: {
        key: '-',
        target: zoom_container,
        fn: zoom_container.zoom_out,
        ignoreWithInput: true
      },
      extent_nodes_ctrl: {
        key: 'ctrl+0',
        target: map,
        fn: map.zoom_extent_nodes
      },
      extent_nodes: {
        key: '0',
        target: map,
        fn: map.zoom_extent_nodes,
        ignoreWithInput: true
      },
      extent_canvas_ctrl: {
        key: 'ctrl+1',
        target: map,
        fn: map.zoom_extent_canvas
      },
      extent_canvas: {
        key: '1',
        target: map,
        fn: map.zoom_extent_canvas,
        ignoreWithInput: true
      },
      view_mode: {
        target: this,
        fn: this.view_mode,
        ignoreWithInput: true
      },
      show_settings_ctrl: {
        key: 'ctrl+,',
        fn: () => this.passPropsSettingsMenu({display: true})
      },
      show_settings: {
        key: ',',
        fn: () => this.passPropsSettingsMenu({display: true}),
        ignoreWithInput: true
      },
      build_mode: {
        key: 'n',
        target: this,
        fn: this.build_mode,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      zoom_mode: {
        key: 'z',
        target: this,
        fn: this.zoom_mode,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      brush_mode: {
        key: 'v',
        target: this,
        fn: this.brush_mode,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      rotate_mode: {
        key: 'r',
        target: this,
        fn: this.rotate_mode,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      text_mode: {
        key: 't',
        target: this,
        fn: this.text_mode,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      toggle_beziers: {
        key: 'b',
        target: map,
        fn: map.toggle_beziers,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      delete_ctrl: {
        key: 'ctrl+backspace',
        target: map,
        fn: map.delete_selected,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      delete: {
        key: 'backspace',
        target: map,
        fn: map.delete_selected,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      delete_del: {
        key: 'del',
        target: map,
        fn: map.delete_selected,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      align_vertical: {
        key: 'alt+l',
        target: map,
        fn: map.align_vertical
      },
      align_horizontal: {
        key: 'shift+alt+l',
        target: map,
        fn: map.align_horizontal
      },
      toggle_primary: {
        key: 'p',
        target: map,
        fn: map.toggle_selected_node_primary,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      cycle_primary: {
        key: 'c',
        target: map,
        fn: map.cycle_primary_node,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      direction_arrow_right: {
        key: 'right',
        target: this.build_input.direction_arrow,
        fn: this.build_input.direction_arrow.right,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      direction_arrow_down: {
        key: 'down',
        target: this.build_input.direction_arrow,
        fn: this.build_input.direction_arrow.down,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      direction_arrow_left: {
        key: 'left',
        target: this.build_input.direction_arrow,
        fn: this.build_input.direction_arrow.left,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      direction_arrow_up: {
        key: 'up',
        target: this.build_input.direction_arrow,
        fn: this.build_input.direction_arrow.up,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      undo: {
        key: 'ctrl+z',
        target: map.undo_stack,
        fn: map.undo_stack.undo,
        requires: 'enable_editing'
      },
      redo: {
        key: 'ctrl+shift+z',
        target: map.undo_stack,
        fn: map.undo_stack.redo,
        requires: 'enable_editing'
      },
      select_all: {
        key: 'ctrl+a',
        target: map,
        fn: map.select_all,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      select_none: {
        key: 'ctrl+shift+a',
        target: map,
        fn: map.select_none,
        ignoreWithInput: true,
        requires: 'enable_editing'
      },
      invert_selection: {
        target: map,
        fn: map.invert_selection,
        requires: 'enable_editing'
      },
      search_ctrl: {
        key: 'ctrl+f',
        fn: () => this.passPropsSearchBar({display: true}),
        requires: 'enable_search'
      },
      search: {
        key: 'f',
        fn: () => this.passPropsSearchBar({display: true}),
        ignoreWithInput: true,
        requires: 'enable_search'
      }
    };
  }

  /**
   * Ask if the user wants to exit the page (to avoid unplanned refresh).
   */
  _setupConfirmBeforeExit() {
    window.onbeforeunload = (_) => this.settings.get('never_ask_before_quit') ?
      null :
      'You will lose any unsaved changes.';
  }

  /**
   * Toggle full screen mode.
   */
  full_screen() { // eslint-disable-line camelcase
    // these settings can update in full screen if provided
    const fullScreenSettings = [
      'menu',
      'scroll_behavior',
      'enable_editing',
      'enable_keys',
      'enable_tooltips'
    ];

    if (this.isFullScreen) {
      d3Select('html').classed('fill-screen', false);
      d3Select('body').classed('fill-screen', false);
      this.selection.classed('fill-screen-div', false);
      this.isFullScreen = false;

      // clear escape listener
      if (this.clearFullScreenEscape) {
        this.clearFullScreenEscape();
        this.clearFullScreenEscape = null;
      }

      // hack for full screen in jupyterlab / notebook
      if (this.savedFullScreenParent) {
        const parentNode = this.savedFullScreenParent.node();
        parentNode.insertBefore(this.selection.remove().node(), parentNode.firstChild);
        this.savedFullScreenParent = null;
      }

      // apply the saved settings
      if (this.savedFullScreenSettings !== null) {
        _.mapObject(this.savedFullScreenSettings, (v, k) => {
          this.settings.set(k, v);
        });
      }
      this.savedFullScreenSettings = null;
    } else {
      // save current settings
      const fullScreenButton = this.settings.get('full_screen_button');
      if (_.isObject(fullScreenButton)) {
        this.savedFullScreenSettings = (
          _.chain(fullScreenButton)
            .pairs()
            .map(([k, v]) => {
              if (_.contains(fullScreenSettings, k)) {
                const currentSetting = this.settings.get(k);
                this.settings.set(k, v);
                return [k, currentSetting];
              } else {
                console.warn(`${k} not recognized as an option for full_screen_button`);
                return [null, null];
              }
            })
            .filter(([k, v]) => k)
            .object()
            .value()
        );
      }

      d3Select('html').classed('fill-screen', true);
      d3Select('body').classed('fill-screen', true);
      this.selection.classed('fill-screen-div', true);
      this.isFullScreen = true;

      // hack for full screen in jupyterlab
      this.savedFullScreenParent = d3Select(this.selection.node().parentNode as Element);
      const bodyNode = d3Select('body').node() as Element;
      bodyNode.insertBefore(this.selection.remove().node(), bodyNode.firstChild);

      // set escape listener
      this.clearFullScreenEscape = this.map.key_manager.addEscapeListener(
        () => this.full_screen()
      );
    }
    this.map.zoom_extent_canvas();
    this.passPropsButtonPanel({isFullScreen: this.isFullScreen});
    this.passPropsMenuBar({isFullScreen: this.isFullScreen});
  }

  getSavingState() {
    return {
      type: 'Datagrok metabolic network',
      selectedNodes: Object.keys(this.map.getSelectedNodes()),
      lastAction: this.map.last_action,
      mapJson: this.map.map_for_export(),
      cobraModel: this.model_data,
      description: '',
    };
  }

  // here, state is stringified JSON
  loadSavingState(stateJSON: string | null) {
    if (!stateJSON)
      return;
    const state: ReturnType<typeof this.getSavingState> = JSON.parse(stateJSON);
    if (!state.mapJson)
      return;
    if (state.cobraModel)
      this.load_model(state.cobraModel, true);

    this.load_map(state.mapJson, true);

    if (state.selectedNodes)
      this.map.select_nodes(state.selectedNodes);


    if (state.lastAction && state.lastAction.function && state.lastAction.args && state.lastAction.function in this.map) {
      if (state.lastAction.args.length == 1)
        this.map[state.lastAction.function](state.lastAction.args[0]);
      else if (state.lastAction.args.length == 2)
        this.map[state.lastAction.function](state.lastAction.args[0], state.lastAction.args[1]);
    }

    this.map.selectionChanged && this.map.selectionChanged();
  }
}

export default utils.class_with_optional_new(Builder);
