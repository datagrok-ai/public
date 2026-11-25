import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Track} from '@datagrok-libraries/tutorials/src/track';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {TutorialRunner} from './tutorial-runner';


const trackPath = (track: Track) => `${Tutorial.APP_PATH}/${track.name.replaceAll(' ', '')}`;

/** The app's view: the tracks (or one of them) with their progress and cards; a card starts its tutorial. */
export function tutorialsView(tracks: Track[], track?: Track): DG.View {
  const view = DG.View.create(Tutorial.VIEW_CLASS);
  view.name = Tutorial.APP_TITLE;
  view.path = track ? trackPath(track) : Tutorial.APP_PATH;
  view.append(ui.div((track ? [track] : tracks).map((t) => new TutorialRunner(t).root), 'tutorials-view-tracks'));
  return view;
}

/** The Browse preview of a tutorial: what it teaches and whether it is done; only Start runs it. */
export function tutorialView(tutorial: Tutorial, track: Track): DG.View {
  const view = DG.View.create('tutorials-tutorial-view');
  view.name = tutorial.name;
  view.path = trackPath(track);
  const status = ui.divText('', 'tutorials-tutorial-view-status');
  tutorial.updateStatus().then(() => status.textContent = tutorial.status ? 'Completed' : '');
  view.append(ui.divV([
    ui.image(tutorial.imageUrl, 360, 240),
    ui.h1(tutorial.name),
    ui.divText(track.name, 'tutorials-tutorial-view-track'),
    ui.divText(tutorial.description),
    status,
    ui.div([ui.bigButton('Start', () => TutorialRunner.run(tutorial, track))]),
  ]));
  return view;
}
