/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/platform/viewers-from-a-package.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../bindings/demo.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openToolbox, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Platform behaviour from a package", () => {
  const session = feature(test, "features/platform/viewers-from-a-package.feature", import.meta.url);
  test("A viewer from the toolbox", {tag: ["@platform", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(9, "When user opens toolbox", () => openToolbox(page));
    await session.step(10, "And user clicks on scatter plot icon on toolbox", () => clickOn(page, el("scatter plot icon on toolbox")));
    await session.step(11, "Then scatter plot viewer should be added to the open tableview", () => viewerAdded(page, "scatter plot"));
  });
});
