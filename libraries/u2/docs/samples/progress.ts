import {Card, Control, ProgressBar, button, component, computed, notify, signal} from '@datagrok-libraries/u2';

export function scoringJob(total: number): Control {
  return component((self) => {
    const done = signal(0);
    const running = signal(false);
    const progress = new ProgressBar({showPercent: true, description: computed(() => `${done.value} of ${total} compounds`)});
    let stop = () => {};
    const step = (to: number) => {
      done.value = to;
      progress.value = to / total;
    };
    const finish = (message: string) => {
      stop();
      running.value = false;
      notify.info(message);
    };
    const run = () => {
      step(0);
      running.value = true;
      stop = self.setInterval(() => {
        step(Math.min(total, done.value + 50));
        if (done.value === total)
          finish(`Scored ${total} compounds`);
      }, 100);
    };
    return new Card({icon: 'cogs', title: 'Score compounds', subtitle: computed(() => running.value ? 'Running' : 'Idle'),
      children: [progress],
      footer: [
        button('Run', run, {primary: true, disabled: computed(() => running.value ? 'Already running' : null)}),
        button('Cancel', () => finish(`Cancelled after ${done.value} compounds`),
          {disabled: computed(() => running.value ? null : 'Nothing is running')}),
      ]});
  }, {name: 'scoring-job'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = scoringJob(2000);
  host.append(demo.root);
  return demo;
}
