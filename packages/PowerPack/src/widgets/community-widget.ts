/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class CommunityWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.panel([], 'welcome-community-widget'));

    grok.dapi.fetchProxy('https://community.datagrok.ai/latest.json', undefined, 48 * 60 * 60)
      .then((r) => r.json())
      .then((json) => {
        let links = json.topic_list.topics
          .slice(0, 10)
          .map((t: any) => ui.link(t.title, `https://community.datagrok.ai/t/${t.slug}`), 'Open in new tab');
        this.root.appendChild(ui.divV(links));
      });

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Community');
  }
}

