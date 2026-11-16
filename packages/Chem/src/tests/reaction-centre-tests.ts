/**
 * Category `Reactions: centre`: the reaction centre Chem's pictures colour, and the legend that names its colours. The
 * Reaction Sketcher's default example (five mapped steps) was once drawn with whole rings red and blue: a route's later
 * step was coloured with the step before's map numbers (crux-core decision `route-step-maps`). Every value here
 * colours what a table written by hand from the mapped text says (crux-core scheme-draw tests/data/centre-steps.jsonl),
 * as Crux counts the bonds its picture colours (`ChemReactionScheme.get_centre_counts`: made, broken, changed on each
 * side); and the Reaction panel (so the app, and full screen) names the colours under a picture that has them, and only
 * then.
 */
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {ChemicalEngine} from '../engine/chem-engine';
import {ChemWasm, cruxModuleSync} from '../engine/crux-module';
import {cruxCentreCounts} from '../rendering/crux-reaction-drawing';
import {openReactionFullScreen, reactionPanel} from '../widgets/reaction-panel';
import {CENTRE_COLOURS, LEGEND_CLASS} from '../widgets/reaction-centre-legend';
import {defaultExample, examples} from '../utils/reactions/reaction-document';
import {inEngine} from './engine-utils';
import {loadFileAsText} from './utils';

/** What each picture colours, `made broken changed` (the table; a changed bond counts on each side). */
const DEFAULT_EXAMPLE = '4 6 0';
/** The default example's five steps drawn alone: amide coupling, Boc removal, SNAr, Suzuki, reductive amination. */
const DEFAULT_STEPS = ['1 1 0', '0 1 0', '1 1 0', '1 2 0', '1 1 0'];
/** The Examples menu in its order: paracetamol, Diels-Alder, amide template, Suzuki, Suzuki template, Mitsunobu, Ugi,
 * isotopes, imatinib. */
const EXAMPLES = ['0 0 0', '2 0 8', '1 1 0', '1 2 0', '1 2 0', '0 0 0', '0 0 0', '0 0 0', '0 0 0'];
/** files/tests/reactions-crux.csv's mapped rows; every other row (the routes too) is unmapped: `0 0 0`. */
const CRUX_CSV: {[id: string]: string} = {
  'dgdef-suzuki-coupling': '1 2 0', 'dgdef-heck-reaction': '1 1 0', 'dgdef-sonogashira-coupling': '0 0 0',
  'dgdef-grignard-ketone': '1 1 2', 'dgdef-wittig-reaction': '1 2 0', 'dgdef-amide-coupling': '1 1 0',
  'spec-01m': '1 1 0', 'spec-03m': '1 1 0', 'deli-DEL006-step1': '0 0 0', 'deli-DEL006-step2': '0 1 0',
  'deli-DEL006-step3': '1 0 0', 'deli-DEL006-step4': '1 0 0', 'own-smarts-09': '0 0 0', 'own-smarts-10': '1 2 0',
  'own-smarts-11': '0 0 0', 'own-smarts-05': '1 2 0', 'own-smarts-01': '1 1 0', 'own-smarts-03': '0 0 0',
};
const GREEN_RED = `${CENTRE_COLOURS[0]} ${CENTRE_COLOURS[1]}`;

/** What Crux's picture of `text` colours, `made broken changed`. */
function counts(text: string): string {
  return (cruxCentreCounts(text) ?? ['no counts']).join(' ');
}

/** The legends in `root`, in order: the colours each names (`data-colours`), '' for a hidden one. */
function legends(root: HTMLElement): string[] {
  return Array.from(root.querySelectorAll(`.${LEGEND_CLASS}`))
    .map((l) => (l as HTMLElement).style.display === 'none' ? '' : (l as HTMLElement).dataset.colours ?? '');
}

/** A CSV's rows by its header, the values quoted as Datagrok's CSV reader reads them. */
function rows(csv: string): {[column: string]: string}[] {
  const df = DG.DataFrame.fromCsv(csv);
  return Array.from({length: df.rowCount}, (_, i) =>
    Object.fromEntries(df.columns.toList().map((c) => [c.name, c.get(i) as string])));
}

category('Reactions: centre', () => {
  test('the app\'s examples colour what the table says, each step of the default example too', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const example = defaultExample().text;
      expect(counts(example), DEFAULT_EXAMPLE, 'the default example: a red bond a step, a 2nd in the Suzuki step, ' +
        'a green one in each step but the Boc removal, no blue (whole rings were red and blue before)');
      const scheme = new (cruxModuleSync() as unknown as ChemWasm).ChemReactionScheme(example, '');
      try {
        expect(scheme.get_steps().map(counts).join(', '), DEFAULT_STEPS.join(', '), 'its steps alone');
      } finally {
        scheme.free();
      }
      expect(examples().map((e) => counts(e.text)).join(', '), EXAMPLES.join(', '), 'the Examples menu');
    });
  });

  test('Chem\'s test reactions colour what the table says; unmapped ones nothing', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const crux = rows(await loadFileAsText('tests/reactions-crux.csv'));
      expect(crux.length, 41, 'the rows of reactions-crux.csv (110 lines: its RXN blocks span lines)');
      for (const r of crux)
        expect(counts(r['reaction']), CRUX_CSV[r['id']] ?? '0 0 0', r['id']);
      // the file's last line is blank, which Datagrok reads as an empty row: no reaction, so no counts (the engine
      // refuses an empty text, "not a reaction")
      const tests = rows(await loadFileAsText('test-reactions.csv'));
      expect(counts(''), 'no counts', 'an empty text');
      for (const r of tests.filter((t) => (t['reaction'] ?? '').trim() !== ''))
        expect(counts(r['reaction']), '0 0 0', r['description']);
      expect(tests.filter((t) => (t['reaction'] ?? '').trim() !== '').length, 15, 'test-reactions.csv\'s reactions');
    });
  });

  // a route's legend sat at the end, under its steps (and a second under the route), while a single reaction's sat
  // right under its picture: one legend now, right under the route's picture, its steps sharing it, in the panel and
  // in full screen
  test('the Reaction panel names the colours once, right under its picture, a route\'s steps sharing it; full ' +
    'screen too', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      const example = defaultExample().text;
      const underPicture = (root: HTMLElement) => {
        const legend = root.querySelector(`.${LEGEND_CLASS}`);
        return !!legend && legend.previousElementSibling === root.querySelector('canvas');
      };
      const panel = reactionPanel(example);
      try {
        expect(legends(panel.root).join(' | '), GREEN_RED, 'one legend: the route\'s, its steps sharing it');
        expect(underPicture(panel.root), true, 'right under the route\'s picture');
      } finally {
        panel.detach();
      }
      const scheme = new (cruxModuleSync() as unknown as ChemWasm).ChemReactionScheme(example, '');
      const steps = scheme.get_steps();
      scheme.free();
      const dialog = openReactionFullScreen(example, steps);
      try {
        expect(legends(dialog.root).join(' | '), GREEN_RED, 'full screen: one legend');
        expect(underPicture(dialog.root), true, 'full screen: right under the route\'s picture');
      } finally {
        dialog.close();
      }
      // the Diels-Alder example: green and blue, no red
      const dielsAlder = reactionPanel(examples()[1].text);
      try {
        expect(legends(dielsAlder.root).join(' | '), `${CENTRE_COLOURS[0]} ${CENTRE_COLOURS[2]}`, 'Diels-Alder');
      } finally {
        dielsAlder.detach();
      }
    });
  });

  test('a picture that colours nothing has no legend', async () => {
    await inEngine(ChemicalEngine.Crux, async () => {
      // unmapped: a route (its picture and its steps), a template mapped only on its R groups
      for (const text of [examples()[0].text, '[*:1]C(=O)O.[*:2]N>>[*:1]C(=O)N[*:2]']) {
        const panel = reactionPanel(text);
        try {
          expect(legends(panel.root).every((l) => l === ''), true, `${text}: ${legends(panel.root).join(' | ')}`);
        } finally {
          panel.detach();
        }
      }
    });
  });
});
