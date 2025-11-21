import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

//@ts-ignore
import '../../css/ai.css';
import {Plan, Step} from './interfaces';

function renderHelm(value: string): HTMLElement {
  return ui.wait(async () => {
    //@ts-ignore
    const helmInput = await ui.input.helmAsync('helm', {
      editable: false,
    });
    helmInput.setStringValue(value);
    await DG.delay(200); // wait for proper sizing
    helmInput.getInput().addEventListener('click', () => {
      grok.shell.o = helmInput.getValue();
    });
    helmInput.getInput().addEventListener('dblclick', () => {
      helmInput.showEditorDialog();
    });

    helmInput.getInput().style.width = '100%';
    helmInput.getInput().style.setProperty('height', '300px', 'important');
    return helmInput.getInput();
  });
}


export class AssistantRenderer {
  static renderPlan(plan: Plan): HTMLElement {
    const container = ui.divV([], 'ai-execute-plan');
    if (!plan?.steps?.length) {
      container.append(ui.divText('No steps generated.'));
      return container;
    }

    const analysis = ui.divV(
      (plan.analysis ?? []).map((line: string) => ui.divText(line))
    );

    const stepsAccordion = ui.accordion();
    plan.steps.map((s: Step, i: number) =>
      stepsAccordion.addPane(
        `Step ${i + 1}: ${s.function}`,
        () => {
          return ui.divV([
            ui.divText(`Action: ${s.action}`),
            ui.divText(`Inputs: ${JSON.stringify(s.inputs, null, 2)}`),
            ui.divText(`Outputs: ${JSON.stringify(s.outputs)}`),
          ]);
        },
        false
      )
    );

    container.append(ui.h2('Plan'), analysis, stepsAccordion.root);
    return container;
  }

  static renderResult(result: any): HTMLElement {
    const container = ui.div([ui.h2('Result')], 'ai-execute-result');
    const val = result?.finalResult?.value ?? result?.finalResult ?? result;
    const meta = result?.finalResult?.meta;

    if (meta?.propertyType === 'graphics')
      container.append(renderGraphics(val));
    else if (meta?.propertyType === 'widget')
      container.append(val.root);
    else if (typeof val === 'string' && meta?.semType === DG.SEMTYPE.MOLECULE)
      container.append(grok.chem.drawMolecule(val, 200, 300));
    else if (typeof val === 'string' && meta?.semType === DG.SEMTYPE.MACROMOLECULE &&
      (meta?.units?.toLowerCase() === 'helm' || meta?.options?.units?.toLowerCase() === 'helm')
    )
      container.append(renderHelm(val));
    else if (meta?.propertyType === 'dataframe') {
      const grid = val.plot.grid().root;
      grid.style.width = '100%';
      container.append(grid);
    } else
      container.append(ui.divText(String(val)));

    return container;
  }
}

function renderGraphics(intermediateResult: HTMLElement): HTMLElement {
  if (!intermediateResult?.style?.backgroundImage)
    return ui.divText('No graphics available.');

  try {
    const backgroundImageUrl = intermediateResult.style.backgroundImage;
    let base64Data = backgroundImageUrl.split('data:image/png;base64,')[1];
    base64Data = base64Data.slice(0, -2);
    return ui.image(`data:image/png;base64,${base64Data}`, 200, 300);
  } catch (err) {
    console.error('Failed to render graphics output:', err);
    return ui.divText('Failed to render graphics.');
  }
}
