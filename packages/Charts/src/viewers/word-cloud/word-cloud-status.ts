import * as DG from 'datagrok-api/dg';
import type {WordCloudViewer} from './word-cloud-viewer';

type Box = {x: number, y: number, width: number, height: number};

type LaidWord = {name: string, box: Box, font: string};

type Readings = {[name: string]: number | string | boolean};

/** Where the layout put every word it managed to draw, in CSS px of the chart canvas. The series'
 * graphic elements are the only record the layout leaves, and a word it dropped has none — so an
 * absent `word` area is the honest answer for a word that was counted but not drawn. */
function laidOutWords(chart: any): LaidWord[] {
  const data = chart?.getModel?.()?.getSeriesByIndex?.(0)?.getData?.();
  if (!data)
    return [];
  const words: LaidWord[] = [];
  for (let i = 0; i < data.count(); i++) {
    const el = data.getItemGraphicEl(i);
    if (!el)
      continue;
    const rect = el.getBoundingRect().clone();
    const transform = el.getComputedTransform();
    if (transform)
      rect.applyTransform(transform);
    const style = el.style ?? {};
    words.push({
      name: String(data.getName(i)),
      box: {x: rect.x, y: rect.y, width: rect.width, height: rect.height},
      font: `${style.fontWeight ?? ''} ${style.fontFamily ?? ''}`.trim(),
    });
  }
  return words;
}

/** How many words the layout has placed — the count the status reports, and what tells a render
 * pass that the picture is there: `setOption` returns before the layout runs, so a canvas is not
 * evidence of a cloud. */
export function laidOutWordCount(chart: any): number {
  return laidOutWords(chart).length;
}

/** What the word cloud shows, for automation: the laid-out words as hit areas in CSS px of the
 * chart canvas, and the readings a test compares. While the viewer shows a message instead of a
 * cloud it reports that message and nothing the previous frame drew — `render` returns before
 * re-creating the chart, so `chart` still holds the geometry of a cloud that is no longer on
 * screen. */
export function wordCloudStatus(v: WordCloudViewer): DG.IWidgetStatus & {values: Readings} {
  const error = v.renderError;
  const canvas = error === null ? v.root.querySelector('canvas') as HTMLCanvasElement | null : null;
  const hitAreas: {[name: string]: Box} = {};
  const values: Readings = {};

  values['column'] = v.wordColumnName;
  values['rows shown'] = v.filter.trueCount;

  if (canvas) {
    const origin = canvas.getBoundingClientRect();
    const words = laidOutWords(v.chart);
    for (const word of words)
      hitAreas[`word "${word.name}"`] = word.box;
    hitAreas['view'] = {x: 0, y: 0, width: origin.width, height: origin.height};
    values['words'] = words.length;
    values['word names'] = words.map((word) => word.name).join(', ');
    if (words.length > 0)
      values['font'] = words[0].font;
    for (const [name, count] of v.wordCounts)
      values[`rows of word "${name}"`] = count;
  }

  return {
    parts: canvas ? {root: v.root, canvas} : {root: v.root},
    hitAreas, values, shortcuts: {}, events: [], description: null, error,
  };
}
