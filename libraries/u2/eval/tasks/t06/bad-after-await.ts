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
    const names = campaigns.map((c) => new TextInput({label: c.code, value: c.name}));
    const renamed = computed(() => names.filter((input, i) => input.value !== campaigns[i].name).length);
    const status = span('');
    this.effect(() => status.textContent = `${renamed.value} renamed`);
    this.root.replaceChildren(new Form().addAll(names).root, status);
  }
}

export function campaignsPanel(load: () => Promise<Campaign[]>): Control {
  return new CampaignsPanel(load);
}
