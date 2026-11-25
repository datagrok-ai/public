/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/reactions-app-text.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.reactions.app]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '../../bindings/reaction-editor.js';
import '../../bindings/reaction-export.js';
import '../../bindings/sketcher-app.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {cruxSpot} from '../../bindings/crux.js';
import {engineGate} from '../../bindings/engine.js';
import {appDrawingEmpty, asksForArrow, clickAppClear, dragOnCanvas, openAppEmpty, openAppExample, previewAsks, textIs, viewWords} from '../../bindings/reactions-app-text.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Reaction Sketcher's Text follows its drawing: a lone molecule is no error", () => {
  const session = feature(test, "features/sketcher/reactions-app-text.feature", import.meta.url);
  test("The Reaction Sketcher's Text follows its drawing: a lone molecule is no error", {tag: ["@journey", "@realizes:chem.reactions.app", "@eval-b7"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await run.scenario("Clear canvas as the app opens empties the drawing and the Text; the preview asks for a reaction", async () => {
      await session.step(17, "Given user opens Chem's Reactions app on the default example", () => openAppExample(page));
      await session.step(18, "When user clicks on the Reactions app's Clear canvas", () => clickAppClear(page));
      await session.step(19, "Then the Reactions app's drawing should be empty", () => appDrawingEmpty(page));
      await session.step(20, "And the Reactions app's Text should be \"\"", () => textIs(page, ""));
      await session.step(21, "And the Reactions app's preview should ask for a reaction", () => previewAsks(page));
      await session.step(22, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Ethane drawn alone is no error; an arrow after it makes it a reaction", async () => {
      await session.step(25, "Given user opens Chem's Reactions app on an empty reaction", () => openAppEmpty(page));
      await session.step(26, "When user clicks on crux single bond tool", () => clickOn(page, el("crux single bond tool")));
      await session.step(27, "And user clicks on Crux canvas 30% across and 50% down", () => cruxSpot(page, 30, 50));
      await session.step(28, "Then the Reactions app's Text should be \"CC\"", () => textIs(page, "CC"));
      await session.step(29, "And the Reactions app's preview should say that a reaction needs an arrow", () => asksForArrow(page));
      await session.step(30, "And no text of the Reactions app should contain \"RDKit\", \"Crux\", \"M  END\" or \"confidence\"", () => viewWords(page, "RDKit", "Crux", "M  END", "confidence"));
      await session.step(31, "When user clicks on crux reaction arrow tool", () => clickOn(page, el("crux reaction arrow tool")));
      await session.step(32, "And user drags on Crux canvas from 50% to 75% across, 50% down", () => dragOnCanvas(page, 50, 75, 50));
      await session.step(33, "Then the Reactions app's Text should be \"CC>>\"", () => textIs(page, "CC>>"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
