import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TutorialRunner, TutorialSubstitute} from './tutorial-runner';
import {chem} from './tracks/chem';
import {eda} from './tracks/eda';
import {da} from './tracks/data-access';
import {ml} from './tracks/ml';
import {TutorialWidget} from './widget';
import '../css/tutorial.css';
import {Track} from '@datagrok-libraries/tutorials/src/track';
import {DemoView} from './demo-app/demo-app';


export const _package = new DG.Package();

const tracks: Track[] = [];

//name: Tutorials
//tags: app
//top-menu: Help | Tutorials @Toolbox Help | Tutorials
export function trackOverview() {
  const tutorialRunners = tracks.map((track) => new TutorialRunner(track));
  const root = ui.div([
    ...tutorialRunners.map((runner) => runner.root),
    ui.panel([], {id: 'tutorial-child-node', style: {paddingTop: '10px'}}),
  ], 'tutorials-root');

  grok.shell.dockManager.dock(root, DG.DOCK_TYPE.RIGHT, null, 'Tutorials', 0.3);
  setPath(window.location.pathname, tutorialRunners);
}

//output: widget tutorial
export function tutorialWidget(): DG.Widget {
  return new TutorialWidget(...tracks.map((track) => new TutorialRunner(track)));
}

//tags: init
export async function tutorialsInit() {
  const properties = await _package.getProperties();
  setProperties(properties as unknown as { [propertyName: string]: boolean });

  const tutorialFuncs = DG.Func.find({tags: ['tutorial']});
  const trackFuncs = DG.Func.find({tags: ['track']});
  const defaultIcon = `${_package.webRoot}package.png`;
  const tutorialIconPaths: { [packageName: string]: { [tutorialName: string]: {
    tutorial: TutorialSubstitute,
    imageUrl: string,
  } } } = {};

  for (const track of tracks) {
    for (const tutorial of track.tutorials) {
      if (tutorial.imageUrl === '')
        tutorial.imageUrl = `${_package.webRoot}images/${tutorial.name.toLowerCase().replace(/ /g, '-')}.png`;
    }
  }

  for (const func of trackFuncs) {
    const t = tracks.find((t) => t.name === func.options['name']);
    if (!t)
      tracks.push(new Track(func.options['name'], [], func.helpUrl ?? ''));
    else
      console.error(`Tutorials: Couldn't add a track. Track ${func.options['name']} already exists.`);
  }

  for (const func of tutorialFuncs) {
    const trackName = func.options['track'] ?? func.package.friendlyName;
    const tutorial = new TutorialSubstitute(func.options['name'], func.description, defaultIcon, func);

    if (func.options['icon']) {
      tutorialIconPaths[func.package.name] = {[func.options['name']]: {
        tutorial,
        imageUrl: func.options['icon'],
      }};
    }

    const track = tracks.find((t) => t.name === trackName);
    if (track) {
      const t = track.tutorials.find((t) => t.name === tutorial.name);
      if (!t)
        track.tutorials.push(tutorial);
      else {
        console.error(`Tutorials: Couldn't add a tutorial to track ${track.name}. ` +
          `Tutorial ${tutorial.name} already exists.`);
      }
    } else
      tracks.push(new Track(trackName, [tutorial], ''));
  }

  grok.events.onPackageLoaded.subscribe((p) => {
    if (p.name in tutorialIconPaths) {
      for (const tutorial in tutorialIconPaths[p.name]) {
        if (tutorialIconPaths[p.name].hasOwnProperty(tutorial)) {
          const data = tutorialIconPaths[p.name][tutorial];
          data.tutorial.imageUrl = `${data.tutorial.func.package.webRoot}${data.imageUrl}`;
        }
      }
    }
  });
}

//name: Demo
//tags: app
//description: Interactive demo of major Datagrok capabilities
export function demoApp() {
  const pathSegments = window.location.pathname.split('/');
  const demoView = new DemoView();
  grok.shell.addView(demoView);

  if (pathSegments.length > 4) {
    const category = pathSegments[4];
    if (pathSegments[pathSegments.length - 1] === category) {
      demoView.nodeView(category);
      return;
    }

    const viewerName = pathSegments[5].replaceAll('%20', ' ');
    const f = DemoView.findDemoFunc(`${category} | ${viewerName}`);
    if (f) {
      const viewPath = `${category}/${viewerName}`;
      demoView.startDemoFunc(f, viewPath);
    }
  }
}

function setProperties(properties: { [propertyName: string]: boolean }): void {
  if (tracks.length > 0)
    return;

  const registry: { [propertyName: string]: Track } = {
    'dataAnalysisTrack': eda,
    'machineLearningTrack': ml,
    'cheminformaticsTrack': chem,
    'dataAccessTrack': da,
  };

  for (const property in properties) {
    if (property in registry && properties[property] === true)
      tracks.push(registry[property]);
  }
}

function setPath(path: string, tutorialRunners: TutorialRunner[]): void {
  const pathParts = path.split('/');
  const removeSpaces = (s: string) => s.replaceAll(' ', '');
  const trackShortNames: {[key: string]: Track} = {
    'eda': eda,
    'ml': ml,
    'chem': chem,
    'access': da,
  };

  if (pathParts.length !== 5)
    return;

  const [trackName, tutorialName] = pathParts.slice(3);
  let track: Track | null = null;
  let trackIdx: number | null = null;
  if (trackName in trackShortNames) {
    track = trackShortNames[trackName];
    const _idx = tracks.findIndex((t) => t === track);
    trackIdx = _idx === -1 ? null : _idx;
  } else {
    track = tracks.find((t, idx) => {
      if (removeSpaces(t.name) === trackName) {
        trackIdx = idx;
        return true;
      }
      return false;
    }) ?? null;
  }
  const tutorial = track?.tutorials?.find((t) => removeSpaces(t.name) === tutorialName);
  if (tutorial && trackIdx != null)
    tutorialRunners[trackIdx].run(tutorial);
}
