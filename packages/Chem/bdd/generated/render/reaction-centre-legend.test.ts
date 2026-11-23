/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/reaction-centre-legend.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.reactions.centre-legend]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-editor.js';
import '../../bindings/reactions-app-text.js';
import '../../bindings/sketcher-app.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {engineGate} from '../../bindings/engine.js';
import {appFullScreen, appRouteLegend, fullScreenRouteLegend, pngColours, previewLegend} from '../../bindings/reaction-centre-legend.js';
import {choose, exportFromRibbon, openAppOnDefault, save} from '../../bindings/reaction-export.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {fileDownloaded, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Every drawing that colours a reaction's centre names its colours in a legend under it", () => {
  const session = feature(test, "features/render/reaction-centre-legend.feature", import.meta.url);
  test("Every drawing that colours a reaction's centre names its colours in a legend under it", {tag: ["@journey", "@realizes:chem.reactions.centre-legend", "@eval-b7"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await run.scenario("The Reactions app's preview and its full-screen view name the default route's colours", async () => {
      await session.step(17, "Given user opens Chem's Reactions app on its default example", () => openAppOnDefault(page));
      await session.step(18, "Then the Reactions app should name the centre's colours once, right under its route: \"bond made, bond broken\"", () => appRouteLegend(page, "bond made, bond broken"));
      await session.step(19, "When user clicks the Full screen icon of the Reactions app's preview", () => appFullScreen(page));
      await session.step(20, "Then the full-screen Reaction dialog should name the centre's colours once, right under its route: \"bond made, bond broken\"", () => fullScreenRouteLegend(page, "bond made, bond broken"));
      await session.step(21, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A PNG exported with the centre's colours carries their legend", async () => {
      await session.step(24, "Given user opens Chem's Reactions app on its default example", () => openAppOnDefault(page));
      await session.step(25, "And user watches downloads", () => watchDownloads(page));
      await session.step(26, "When user clicks Export... in the Reactions app's ribbon", () => exportFromRibbon(page));
      await session.step(27, "And user chooses \"Screen\", \"double column (17.8 cm)\" and \"PNG 300 dpi\" in the Export panel", () => choose(page, "Screen", "double column (17.8 cm)", "PNG 300 dpi"));
      await session.step(28, "Then the Export panel's preview should carry the centre's legend", () => previewLegend(page));
      await session.step(29, "When user clicks Save in the Export panel", () => save(page));
      await session.step(30, "Then a file \"reaction.png\" should have been downloaded", () => fileDownloaded(page, "reaction.png"));
      await session.step(31, "And the saved PNG should have the centre's green and red in it", () => pngColours(page));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
