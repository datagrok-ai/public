import * as ui from 'datagrok-api/ui';

/** What {@link nglStageOnPage} uses of an NGL stage (`NGL.Stage`, the platform's /js/common/ngl_viewer/ngl.js). */
export interface NglStageLike {
  tasks: {onZeroOnce(callback: () => void): void};
  mouseObserver: {dispose(): void};
  viewer: {renderer: {forceContextLoss(): void}};
  dispose(): void;
}

/** How often the hosts of the stages on the page are checked. */
const RELEASE_CHECK_MS = 1000;

/** A stage made in `host`, and what releases it. */
interface StageOnPage {
  host: HTMLElement;
  release: () => void;
}

const stagesOnPage = new Set<StageOnPage>();
let releaseTimer: ReturnType<typeof setInterval> | null = null;

/**
 * An NGL stage in `host` while `host` is on the page: `open(host)` makes it once `host` is attached, and the stage is
 * released ({@link releaseNglStage}) and `onRelease` runs once `host` has left the page. A page has few WebGL contexts
 * (Chrome keeps 16 and drops the one used least recently): the panels' stages were never released, so after a few
 * molecules the panel on screen lost its context and showed an empty white box. A panel built but never shown makes
 * no stage, and a host does not come back (the context panel makes a new panel for each molecule).
 */
export async function nglStageOnPage<S extends NglStageLike>(host: HTMLElement, open: (host: HTMLElement) => S,
  onRelease?: () => void): Promise<S> {
  await ui.tools.waitForElementInDom(host);
  const stage = open(host);
  stagesOnPage.add({host, release: () => {
    releaseNglStage(stage);
    onRelease?.();
  }});
  releaseTimer ??= setInterval(releaseStagesOffPage, RELEASE_CHECK_MS);
  return stage;
}

function releaseStagesOffPage(): void {
  for (const entry of stagesOnPage) {
    if (!entry.host.isConnected) {
      stagesOnPage.delete(entry);
      entry.release();
    }
  }
  if (stagesOnPage.size === 0 && releaseTimer !== null) {
    clearInterval(releaseTimer);
    releaseTimer = null;
  }
}

/**
 * Releases an NGL stage and what it holds outside itself, which `dispose()` keeps: its WebGL context (until garbage
 * collection, and three.js asks the browser to restore it when it is dropped) and its mouse observer's listeners on
 * the document. Done once the stage is idle: NGL builds a representation in a worker after `loadFile` has resolved,
 * and a stage disposed before then throws (BiostructureViewer's ngl-ui.ts).
 */
export function releaseNglStage(stage: NglStageLike): void {
  stage.tasks.onZeroOnce(() => {
    stage.mouseObserver?.dispose();
    stage.dispose(); // removes three.js's context listeners first: the lost context is not restored
    stage.viewer.renderer?.forceContextLoss();
  });
}
