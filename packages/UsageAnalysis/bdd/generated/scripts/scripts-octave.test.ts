/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/scripts/scripts-octave.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.scripts]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {scriptResult, scriptResultListed} from '../../bindings/scripts.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableOpen} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, scriptsView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Octave template counts the cells of its sample table", () => {
  const session = feature(test, "features/scripts/scripts-octave.feature", import.meta.url);
  test("The Octave template counts the cells of its sample table", {tag: ["@journey", "@full-stand", "@serial", "@realizes:views.scripts"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens the Scripts view", () => scriptsView(page));
    await run.scenario("The Octave template runs with cars", async () => {
      await session.step(25, "When user clicks on New button", () => clickOn(page, el("New button")));
      await session.step(26, "And user picks \"Octave Script...\" from the open menu", () => pickFromOpenMenu(page, "Octave Script..."));
      await session.step(27, "Then the \"Template\" view should be current", () => viewIsCurrent(page, "Template"));
      await session.step(28, "And code editor should contain the text \"#language: octave\"", () => shouldContainText(page, el("code editor"), "#language: octave"));
      await session.step(29, "When user clicks on \"Open script sample table\" icon", () => clickOn(page, el("\"Open script sample table\" icon")));
      await session.step(30, "Then table \"cars\" should be open", () => tableOpen(page, "cars"));
      await session.step(31, "When user clicks on \"Run script (F5)\" icon", () => clickOn(page, el("\"Run script (F5)\" icon")));
      await session.step(32, "Then \"Template\" dialog should be visible", () => shouldBe(page, el("\"Template\" dialog"), "visible"));
      await session.step(33, "When user selects \"cars\" in Table input in \"Template\" dialog", () => selectIn(page, "cars", el("Table input in \"Template\" dialog")));
      await session.step(34, "And user clicks on OK button in \"Template\" dialog", () => clickOn(page, el("OK button in \"Template\" dialog")));
      await session.step(35, "Then the \"Template\" dialog should close", () => dialogCloses(page, "Template"));
      await session.step(36, "And the script results should list \"count\"", () => scriptResultListed(page, "count"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
      await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("Octave counts the cells of cars its own way", async () => {
      await session.step(41, "Then the script results should show \"count\" as \"527\"", () => scriptResult(page, "count", "527"));
    });
    run.finish();
  });
});
