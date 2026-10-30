import {AppHeader, Control, component, notify} from '@datagrok-libraries/u2';

export function hitDesignHeader(): AppHeader {
  return new AppHeader({
    icon: 'flask',
    iconSize: 96,
    title: 'Hit Design',
    subtitle: 'Design, score and submit molecules with your team',
    items: [
      'Sketch molecules or bring them from a file or a query',
      'Score them with your template\'s functions as you go',
      {text: 'Browse the demo campaign', onClick: () => notify.info('Opening the demo campaign')},
    ],
    learnMore: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README.md',
    bottomLine: true,
  });
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = component(() => hitDesignHeader());
  demo.root.classList.add('u2-adaptive');
  host.append(demo.root);
  return demo;
}
