import {Control, Form, TextInput, computed, span} from '@datagrok-libraries/u2';

type Campaign = {code: string, name: string};

class CampaignsPanel extends Control {
  constructor(load: () => Promise<Campaign[]>) {
    super();
    this.root.append(span('Loading…'));
    void this.show(load);
  }

  private async show(load: () => Promise<Campaign[]>): Promise<void> {
    const campaigns = await load();
    if (this.scope.isDisposed)
      return;
    this.build(() => {
      const names = campaigns.map((c) => new TextInput({label: c.code, value: c.name}));
      const renamed = computed(() => names.filter((input, i) => input.value !== campaigns[i].name).length);
      this.root.replaceChildren(new Form().addAll(names).root, span(computed(() => `${renamed.value} renamed`)));
    });
  }
}

export function campaignsPanel(load: () => Promise<Campaign[]>): Control {
  return new CampaignsPanel(load);
}
