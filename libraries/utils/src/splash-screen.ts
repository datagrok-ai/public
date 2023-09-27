import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subscription} from 'rxjs';
import Timeout = NodeJS.Timeout;

export class SplashScreenView extends DG.ViewBase {
  private loader: HTMLElement = ui.div();
  private state: 'load' | 'none' | 'fail' = 'none';
  private sub?: Subscription;
  private timeout?: Timeout;

  constructor(name?: string) {
    super();
    this.name = name ?? 'Loading app';
    this.box = true;
  }

  private init() {
    this.root.innerHTML = '';
    this.loader = ui.div([ui.h2(this.name, {style: appNameCss}), ui.loader()],
      {classes: 'grok-wait', style: loaderCss});
    this.root.appendChild(this.loader);
  }

  add() {
    switch (this.state) {
    case 'load': return;
    case 'fail': this.close();
    }
    //@ts-ignore
    clearTimeout(this.timeout);
    this.init();
    this.sub = grok.events.onViewRemoving.subscribe((e) => {
      //@ts-ignore
      if (e.args.view === this)
        e.preventDefault();
    });
    grok.shell.addView(this);
    this.state = 'load';
    //@ts-ignore
    this.timeout = setTimeout(() => {
      if (this.state !== 'load') return;
      this.loader.querySelector('.grok-loader')?.remove();
      this.loader.appendChild(ui.divText('Cannot load app', {style: {
        ...appNameCss,
        top: '49%',
        color: 'var(--failure)',
        fontSize: '1.5em',
      }}));
      this.sub?.unsubscribe();
      this.state = 'fail';
    }, 30000);
  }

  remove() {
    if (this.state === 'none') return;
    //@ts-ignore
    clearTimeout(this.timeout);
    this.sub?.unsubscribe();
    this.close();
    this.state = 'none';
  }
}

// Move to css file

const appNameCss = {
  top: '44%',
  left: '0',
  right: '0',
  position: 'absolute',
  fontSize: '1.4em',
  width: '300px',
  marginLeft: 'auto',
  marginRight: 'auto',
  textAlign: 'center',
  color: 'var(--blue-1)',
  wordBreak: 'break-all'
};

const loaderCss = {
  inset: 0,
  position: 'absolute',
  zIndex: 1,
  backdropFilter: 'blur(5px)',
  backgroundColor: 'rgba(255, 255, 255, 0.93)',
  transition: '2s'
};
