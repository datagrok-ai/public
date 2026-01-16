import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {DesirabilityProfile, PropertyDesirability} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {
  MPO_SCORE_CHANGED_EVENT,
  MpoProfileEditor,
} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';
import {MpoContextPanel} from './mpo-context-panel';

export class MpoProfileCreateView {
  readonly view: DG.View;
  private df: DG.DataFrame | null = null;
  private editor: MpoProfileEditor;
  private profile: DesirabilityProfile;

  private editorDiv: HTMLDivElement;
  private main: HTMLDivElement;

  private methodInput: DG.ChoiceInput<string | null> | undefined;

  private mpoContextPanel: MpoContextPanel | null = null;

  constructor(existingProfile?: DesirabilityProfile, showMethod: boolean = true) {
    this.view = DG.View.create();
    this.view.name = existingProfile ? 'Edit MPO Profile' : 'Create MPO Profile';

    // default profile
    this.profile = existingProfile ??{
      name: '',
      description: '',
      properties: {
        'Property 1': {weight: 1, min: 0, max: 1, line: []},
        'Property 2': {weight: 1, min: 0, max: 1, line: []},
        'Property 3': {weight: 1, min: 0, max: 1, line: []},
      },
    };

    this.editor = new MpoProfileEditor(undefined, true);
    this.editor.setProfile(this.profile);

    this.editorDiv = ui.divV([this.editor.root], {
      style: {
        flex: '1 1 0',
        width: '100%',
        height: '100%',
        minWidth: '0',
        minHeight: '0',
        overflow: 'auto',
      },
    });

    const controls = this.buildControls(showMethod);

    this.main = ui.divV([controls], {
      style: {
        width: '100%',
        height: '100%',
        display: 'flex',
        flexDirection: 'column',
        minWidth: '0',
        minHeight: '0',
      },
    });

    this.view.root.append(this.main);
    this.attachLayout();
    this.listenForProfileChanges();
    // grok.shell.addView(this.view);
  }

  private buildControls(showMethod: boolean = true): HTMLElement {
    const controls: DG.InputBase[] = [];
    this.methodInput = ui.input.choice('Method', {
      items: ['Manual', 'Probabilistic'],
      value: 'Manual',
      nullable: false,
      onValueChanged: () => this.attachLayout(), // update layout when method changes
    });

    const datasetInput = ui.input.table('Dataset', {
      nullable: true,
      onValueChanged: (df) => {
        this.df = df;
        this.attachLayout();
      },
    });

    if (showMethod)
      controls.push(this.methodInput);
    controls.push(datasetInput);

    return ui.divV([
      ui.h1('New MPO Profile'),
      ui.form(controls),
    ], {
      style: {
        gap: '10px',
        padding: '5px',
        flex: '0 0 auto',
      },
    });
  }

  // private async attachLayout() {
  //   while (this.main.children.length > 1)
  //     this.main.removeChild(this.main.lastChild!);

  //   // Determine profile
  //   if (this.methodInput?.value === 'Manual' && this.df) {
  //     // create profile for DataFrame
  //     this.editor.dataFrame = this.df;
  //     this.profile = this.createProfileForDf(this.df);
  //   } else if (!this.df) {
  //     // default profile
  //     this.profile = {
  //       name: '',
  //       description: '',
  //       properties: {
  //         'Property 1': {weight: 1, min: 0, max: 1, line: []},
  //         'Property 2': {weight: 1, min: 0, max: 1, line: []},
  //         'Property 3': {weight: 1, min: 0, max: 1, line: []},
  //       },
  //     };
  //   }

  //   this.editor.setProfile(this.profile);

  //   if (!this.df) {
  //     this.main.append(this.editorDiv);
  //     return;
  //   }

  //   if (!this.mpoContextPanel)
  //     this.mpoContextPanel = new MpoContextPanel(this.df);

  //   await grok.data.detectSemanticTypes(this.df);

  //   const gridHost = ui.div([], {
  //     style: {
  //       flex: '1 1 0',
  //       width: '100%',
  //       height: '100%',
  //       minWidth: '0',
  //       minHeight: '0',
  //       overflow: 'hidden',
  //       display: 'flex',
  //       flexDirection: 'column',
  //     },
  //   });

  //   const gridPlot = this.df.plot.grid();
  //   gridPlot.root.style.flex = '1 1 0';
  //   gridPlot.root.style.width = '100%';
  //   gridPlot.root.style.height = '100%';
  //   gridHost.append(gridPlot.root);

  //   // editor box flex 65%, grid box flex 35%
  //   const editorBox = ui.box(this.editorDiv, {style: {flex: '1 1 65%', minWidth: '0', minHeight: '0', height: '100%'}});
  //   const gridBox = ui.box(gridHost, {style: {flex: '1 1 35%', minWidth: '0', minHeight: '0', height: '100%'}});

  //   const split = ui.splitH([editorBox, gridBox], {}, true);

  //   Object.assign(split.style, {
  //     flex: '1 1 0',
  //     width: '100%',
  //     height: '100%',
  //     minWidth: '0',
  //     minHeight: '0',
  //   });

  //   if (this.df) {
  //     await this.mpoContextPanel?.render(
  //       this.profile!,
  //       this.editor.columnMapping,
  //       'Average',
  //     );
  //   }

  //   this.main.append(split);
  // }

  private async attachLayout() {
  // Show spinner
    ui.setUpdateIndicator(this.view.root, true, 'Updating layout...');

    try {
    // Clear previous layout (except controls)
      while (this.main.children.length > 1)
        this.main.removeChild(this.main.lastChild!);

      // Determine profile
      if (this.methodInput?.value === 'Manual' && this.df) {
        this.editor.dataFrame = this.df;
        this.editor.design = true;
        if (!this.profile)
          this.profile = this.createProfileForDf(this.df);
      } else if (!this.df && !this.profile) {
        if (!this.profile) {
          this.profile = {
            name: '',
            description: '',
            properties: {
              'Property 1': {weight: 1, min: 0, max: 1, line: []},
              'Property 2': {weight: 1, min: 0, max: 1, line: []},
              'Property 3': {weight: 1, min: 0, max: 1, line: []},
            },
          };
        }
      }

      this.editor.setProfile(this.profile);

      if (!this.df) {
        this.main.append(this.editorDiv);
        return;
      }

      if (!this.mpoContextPanel)
        this.mpoContextPanel = new MpoContextPanel(this.df);

      await grok.data.detectSemanticTypes(this.df);

      const gridHost = ui.div([], {
        style: {
          flex: '1 1 0',
          width: '100%',
          height: '100%',
          minWidth: '0',
          minHeight: '0',
          overflow: 'hidden',
          display: 'flex',
          flexDirection: 'column',
        },
      });

      const gridPlot = this.df.plot.grid();
      gridPlot.root.style.flex = '1 1 0';
      gridPlot.root.style.width = '100%';
      gridPlot.root.style.height = '100%';
      gridHost.append(gridPlot.root);

      const editorBox = ui.box(this.editorDiv, {style: {flex: '1 1 65%', minWidth: '0', minHeight: '0', height: '100%'}});
      const gridBox = ui.box(gridHost, {style: {flex: '1 1 35%', minWidth: '0', minHeight: '0', height: '100%'}});

      const split = ui.splitH([editorBox, gridBox], {}, true);

      Object.assign(split.style, {
        flex: '1 1 0',
        width: '100%',
        height: '100%',
        minWidth: '0',
        minHeight: '0',
      });

      if (this.df) {
        await this.mpoContextPanel?.render(
        this.profile!,
        this.editor.columnMapping,
        'Average',
        );
      }

      this.main.append(split);
    } finally {
    // Hide spinner
      ui.setUpdateIndicator(this.view.root, false);
    }
  }

  /** Creates a basic profile for the DataFrame columns */
  private createProfileForDf(df: DG.DataFrame): DesirabilityProfile {
    const props: {[key: string]: PropertyDesirability} = {};
    for (const col of df.columns.numerical)
      props[col.name] = {weight: 1, min: col.min, max: col.max, line: []}; // simple default, can enhance later

    return {
      name: '',
      description: '',
      properties: props,
    };
  }

  private async listenForProfileChanges(): Promise<void> {
    grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      if (this.df && this.profile && this.mpoContextPanel) {
        await this.mpoContextPanel.render(
          this.profile,
          this.editor.columnMapping,
          'Average',
        );
      }
    });
  }
}
