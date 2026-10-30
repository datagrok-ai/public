import {AsyncPager, Control, VirtualList, button, component, computed, divH, on, show, span, statusLine}
  from '@datagrok-libraries/u2';
import './pager.css';

interface Sample {
  id: string;
  plate: string;
}

const PAGE = 50;

export function sampleList(fetchPage: (page: number) => Promise<Sample[]>, count: () => Promise<number>): Control {
  return component(() => {
    const pager = new AsyncPager<Sample>({fetchPage, count, pageSize: PAGE});
    const list = new VirtualList<Sample>({itemHeight: 28, keyOf: (s) => s.id, render: (s) => span(`${s.id} · ${s.plate}`)});
    list.root.classList.add('sample-list');
    list.setItems(pager.items);
    on(list.root, 'scroll', () => {
      if (list.root.scrollTop + list.root.clientHeight > list.root.scrollHeight - 10 * list.itemHeight)
        pager.loadMore();
    });
    const failed = computed(() => pager.state.value === 'error');
    const status = statusLine(computed(() => failed.value ? 'A page did not load' : null), computed(() =>
      pager.state.value === 'loading' ? 'Loading...' : `${pager.items.value.length} of ${pager.total.value ?? '…'} samples`));
    const retry = divH([]);
    show(retry, failed, () => button('Retry', () => pager.loadMore()));
    pager.reset();
    return [list, divH([status, retry])];
  }, {name: 'sample-list'});
}

// Live demo
let calls = 0;

async function fetchSamples(page: number): Promise<Sample[]> {
  await new Promise((resolve) => setTimeout(resolve, 400));
  if (++calls === 4)
    throw new Error('timeout');
  const from = page * PAGE;
  return Array.from({length: Math.max(0, Math.min(PAGE, 420 - from))},
    (_, i) => ({id: `S-${from + i + 1}`, plate: `Plate ${1 + Math.floor((from + i) / 96)}`}));
}

export function render(host: HTMLElement): Control {
  const demo = sampleList(fetchSamples, async () => 420);
  host.append(demo.root);
  return demo;
}
