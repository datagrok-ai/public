/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/reactions-app-help.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.reactions.app]
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
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
import {engineGate} from '../../bindings/engine.js';
import {clickFullScreen, dialogHelp, hoverFullScreen, openAppExample, tooltipSays} from '../../bindings/reactions-app-text.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {notShowsText, showsText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Reaction Sketcher's help: its preview's buttons say what they are, and its dialogs open Chem's Reactions page", () => {
  const session = feature(test, "features/sketcher/reactions-app-help.feature", import.meta.url);
  test("The Reaction Sketcher's help: its preview's buttons say what they are, and its dialogs open Chem's Reactions page", {tag: ["@journey", "@realizes:chem.reactions.app", "@eval-b7"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await run.scenario("The full-screen button names itself, and the full-screen dialog's help loads", async () => {
      await session.step(13, "Given user opens Chem's Reactions app on the default example", () => openAppExample(page));
      await session.step(14, "When user hovers over the Reactions app's full-screen button", () => hoverFullScreen(page));
      await session.step(15, "Then a tooltip should say \"Full screen\"", () => tooltipSays(page, "Full screen"));
      await session.step(16, "When user clicks the Reactions app's full-screen button", () => clickFullScreen(page));
      await session.step(17, "And user opens the help of the \"Reaction\" dialog", () => dialogHelp(page, "Reaction"));
      await session.step(18, "Then help panel should show the text \"The Reaction panel\"", () => showsText(page, el("help panel"), "The Reaction panel"));
      await session.step(19, "And help panel should not show the text \"Help unavailable\"", () => notShowsText(page, el("help panel"), "Help unavailable"));
      await session.step(20, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
