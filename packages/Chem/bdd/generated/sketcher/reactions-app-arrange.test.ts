/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/reactions-app-arrange.feature
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
import {openAppExample} from '../../bindings/reactions-app-text.js';
import {appArrange, appTextAsRemembered, rememberAppText} from '../../bindings/reactions-sketch-stack.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Arrange in the Reaction Sketcher lays out the drawing and leaves the Text exactly as written", () => {
  const session = feature(test, "features/sketcher/reactions-app-arrange.feature", import.meta.url);
  test("Arrange in the Reaction Sketcher lays out the drawing and leaves the Text exactly as written", {tag: ["@journey", "@realizes:chem.reactions.app", "@eval-b7"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And the Chemical Engine is \"Crux\"", () => engineGate(page, "Crux"));
    await run.scenario("The default route arranged: the Text as remembered", async () => {
      await session.step(13, "Given user opens Chem's Reactions app on the default example", () => openAppExample(page));
      await session.step(14, "When user remembers the Reactions app's Text", () => rememberAppText(page));
      await session.step(15, "And user chooses \"Arrange\" in the Clean Up and Layout menu of the Reactions app's sketcher", () => appArrange(page, "Arrange"));
      await session.step(16, "Then the Reactions app's Text should be as remembered", () => appTextAsRemembered(page));
      await session.step(17, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
