import {ChoiceInput, Control, Form, TextInput, Wizard, button, component, notify, span, tableFromMap}
  from '@datagrok-libraries/u2';

export function newCampaignWizard(create: (name: string, template: string, source: string) => Promise<void>): Wizard {
  const name = new TextInput({label: 'Name', required: true, placeholder: 'e.g. BTK follow-up'});
  const template = new ChoiceInput({label: 'Template', items: ['Kinase hit-to-lead', 'GPCR screen'], value: 'Kinase hit-to-lead'});
  const source = new ChoiceInput({label: 'Molecules', items: ['From a file', 'From a query', 'Draw them'], value: 'From a file'});
  const general = new Form().addAll([name, template]);
  return new Wizard({
    steps: [
      {id: 'general', title: 'General', content: general.root, canProceed: general.problem},
      {id: 'molecules', title: 'Molecules', content: () => new Form().add(source).root},
      {id: 'review', title: 'Review', finishText: 'Create', content: () => tableFromMap({
        'Name': span(name.signal), 'Template': span(template.signal), 'Molecules': span(source.signal)})},
    ],
    onFinish: () => create(name.value, template.value!, source.value!),
  });
}

// Live demo
async function createCampaign(name: string): Promise<void> {
  await new Promise((resolve) => setTimeout(resolve, 800));
  notify.info(`Campaign ${name} created`);
}

export function render(host: HTMLElement): Control {
  const demo = component(() => {
    const wizard = newCampaignWizard(createCampaign);
    return button('New campaign...', () => wizard.openInDialog('New campaign', {width: 520}), {primary: true});
  }, {name: 'new-campaign'});
  host.append(demo.root);
  return demo;
}
