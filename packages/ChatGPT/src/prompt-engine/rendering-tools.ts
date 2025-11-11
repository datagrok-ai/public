import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../css/ai.css';

export class AssistantRenderer {
  static renderPlan(plan: any): HTMLElement {
    const container = ui.divV([], 'assistant-plan');
    if (!plan?.steps?.length) {
      container.append(ui.divText('No steps generated.'));
      return container;
    }

    const analysis = ui.divV(
      (plan.analysis ?? []).map((line: string) =>
        ui.divText(line, { classes: 'assistant-plan-analysis' })
      )
    );

    const steps = ui.list(
      plan.steps.map((s: any, i: number) =>
        ui.divV([
          ui.h3(`Step ${i + 1}: ${s.function}`),
          ui.divText(`Action: ${s.action}`),
          ui.divText(`Inputs: ${JSON.stringify(s.inputs, null, 2)}`),
          ui.divText(`Outputs: ${JSON.stringify(s.outputs)}`),
        ], 'assistant-plan-step')
      )
    );

    container.append(ui.h2('Plan'), analysis, steps);
    return container;
  }
  
  static renderResult(result: any): HTMLElement {
    const container = ui.div([], 'assistant-result');
    const val = result?.finalResult?.value ?? result?.finalResult ?? result;
    const meta = result?.finalResult?.meta;
    
    if (meta?.propertyType === 'graphics') {
      container.append(renderGraphics(val));
    } else if (meta?.propertyType === 'widget') {
      container.append(val.root);
    } else if (typeof val === 'object') {
      container.append(ui.divText(JSON.stringify(val, null, 2)));
    } else {
      container.append(ui.divText(String(val)));
    }
    
    return container;
  }

  static renderFullOutput(plan: any, result: any): DG.Widget {
    const root = ui.divV([], { classes: 'assistant-widget' });
    root.append(
      this.renderPlan(plan),
      ui.h2('Result'),
      this.renderResult(result)
    );
    return new DG.Widget(root);
  }
}

function renderGraphics(intermediateResult: any): HTMLElement {
    if (!intermediateResult?.style?.backgroundImage) {
        return ui.divText('No graphics available.');
    }
    
    try {
        // let base64Data = intermediateResult.style.backgroundImage;
        // base64Data = base64Data.replace(/^url\("data:image\/png;base64,/, '').replace(/"\)$/, '');
        // return ui.image(`data:image/png;base64,${base64Data}`, 200, 300);
      const backgroundImageUrl = intermediateResult.style.backgroundImage;
      let base64Data = backgroundImageUrl.split('data:image/png;base64,')[1];
      base64Data = base64Data.slice(0, -2);
      return ui.image(`data:image/png;base64,${base64Data}`, 200, 300);
    } catch (err) {
      console.error('Failed to render graphics output:', err);
      return ui.divText('Failed to render graphics.');
    }
}