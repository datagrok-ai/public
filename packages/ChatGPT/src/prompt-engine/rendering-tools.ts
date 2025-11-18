import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

//@ts-ignore
import '../css/chatgpt.css';

export class AssistantRenderer {
  static renderPlan(plan: any): HTMLElement {
    const container = ui.divV([], 'chatgpt-ask-ai-plan');
    if (!plan?.steps?.length) {
      container.append(ui.divText('No steps generated.'));
      return container;
    }

    const analysis = ui.divV(
      (plan.analysis ?? []).map((line: string) =>
        ui.divText(line)
      )
    );

    const steps = ui.list(
      plan.steps.map((s: any, i: number) =>
        ui.divV([
          ui.h3(`Step ${i + 1}: ${s.function}`),
          ui.divText(`Action: ${s.action}`),
          ui.divText(`Inputs: ${JSON.stringify(s.inputs, null, 2)}`),
          ui.divText(`Outputs: ${JSON.stringify(s.outputs)}`),
        ], 'chatgpt-ask-ai-plan-step')
      )
    );

    container.append(ui.h2('Plan'), analysis, steps);
    return container;
  }

  static renderResult(result: any): HTMLElement {
    const container = ui.div([ui.h2('Result')], 'chatgpt-ask-ai-result');
    const val = result?.finalResult?.value ?? result?.finalResult ?? result;
    const meta = result?.finalResult?.meta;

    if (meta?.propertyType === 'graphics')
      container.append(renderGraphics(val));
    else if (meta?.propertyType === 'widget')
      container.append(val.root);
    else if (typeof val === 'string' && meta?.semType === DG.SEMTYPE.MOLECULE)
      container.append(grok.chem.drawMolecule(val, 200, 300));
    else if (meta?.propertyType === 'dataframe')
      container.append(val.plot.grid().root);
    else
      container.append(ui.divText(String(val)));

    return container;
  }
}

function renderGraphics(intermediateResult: any): HTMLElement {
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
