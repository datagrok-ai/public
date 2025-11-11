import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {USE_AUTO_SELECTION, PARETO_COLOR_CODING, USE_PARETO_AXES} from './defs';

export class ParetoFrontViewer extends DG.JsViewer {
  private title: string;
  private showTitle: boolean;
  private description: string;
  private descriptionPosition: string;
  private descriptionVisibilityMode: string;
  private minimizeColumnNames: string[];
  private maximizeColumnNames: string[];
  private labelColumnNames: string[];
  private useAutoSelection: boolean = true;
  private xColumnName: string;
  private yColumnName: string;
  private useParetoAxes: boolean;
  private legendVisibility: string;
  private legendPosition: string;
  private colorColumnName: string;
  private paretoColorCoding: boolean = true;
  private initialized = false;

  get type(): string {
    return 'ParetoFrontViewer';
  }

  constructor() {
    super();

    this.title = this.string('title', 'Pareto front');
    this.showTitle = this.bool('showTitle', false, {category: 'Description'});

    this.description = this.string('description');

    this.descriptionPosition = this.string('descriptionPosition', 'Top', {choices: ['Left', 'Right', 'Top', 'Bottom']});

    this.descriptionVisibilityMode = this.string('descriptionVisibilityMode', 'Auto', {
      choices: ['Auto', 'Always', 'Never']},
    );

    this.minimizeColumnNames = this.addProperty('minimizeColumnNames', DG.TYPE.COLUMN_LIST, null, {
      columnTypeFilter: DG.TYPE.NUMERICAL,
      category: 'Data',
      description: 'Columns with features to be minimized during Pareto optimization.',
    });

    this.maximizeColumnNames = this.addProperty('maximizeColumnNames', DG.TYPE.COLUMN_LIST, null, {
      columnTypeFilter: DG.TYPE.NUMERICAL,
      category: 'Data',
      description: 'Columns with features to be maximized during Pareto optimization.',
    });

    this.xColumnName = this.string('xAxisColumnName', null, {
      category: 'Axes',
      description: 'A column to be used on the X axis of the scatter plot.',
      nullable: false,
    });

    this.yColumnName = this.string('yAxisColumnName', null, {
      category: 'Axes',
      description: 'A column to be used on the Y axis of the scatter plot.',
      nullable: false,
    });

    this.useParetoAxes = this.bool('useParetoAxes', USE_PARETO_AXES, {
      defaultValue: USE_PARETO_AXES,
      category: 'Axes',
      description: 'Use optimized variables as scatter plot axes.',
    });

    this.labelColumnNames = this.addProperty('labelColumnsColumnNames', DG.TYPE.COLUMN_LIST, null, {
      category: 'Labels',
      description: 'Label columns to show next to the markers.',
    });

    this.useAutoSelection = this.bool('useAutoSelection', USE_AUTO_SELECTION, {
      category: 'Labels',
      description: 'Automatically select the most relevant columns for the legend.',
      defaultValue: USE_AUTO_SELECTION,
    });

    this.legendVisibility = this.string('legendVisibility', 'Auto', {choices: ['Auto', 'Always', 'Never']});
    this.legendPosition = this.string('legendPosition', 'Top', {
      choices: ['Auto', 'Left', 'Right', 'Top', 'Bottom', 'RightTop', 'RightBottom', 'LeftTop', 'LeftBottom'],
    });

    this.colorColumnName = this.string('colorColumnName', null, {
      category: 'Color',
      description: 'A column to be used for color-coding.',
      nullable: false,
    });

    this.paretoColorCoding = this.bool('paretoColorCoding', PARETO_COLOR_CODING, {
      category: 'Color',
      description: 'Color markers based on whether a point is Pareto optimal or not.',
      defaultValue: PARETO_COLOR_CODING,
    });
  } // constructor

  private init() {
    this.setOptions({
      useParetoAxes: USE_PARETO_AXES,
      useAutoSelection: USE_AUTO_SELECTION,
      paretoColorCoding: PARETO_COLOR_CODING,
    });

    console.log(this);
  }

  onTableAttached() {
    this.init();

    // Stream subscriptions
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));

    this.render();
  }

  // Cancel subscriptions when the viewer is detached
  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
    super.detach();
  }

  // Override to handle property changes
  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);
    if (this.initialized)
      this.render();
  }

  render(computeData = true) {
    if (computeData) {

    }
  } // render
}
