/* ---
generated: features/platform/toolbox-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openToolbox, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Adding viewers from the toolbox", () => {
  const session = feature(test);
  test("Scatter plot from the toolbox icon", {tag: ["@platform", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user is logged in", () => loggedIn(page));
    await test.step("And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await test.step("When user opens toolbox", () => openToolbox(page));
    await test.step("And user clicks on scatter plot icon on toolbox", () => clickOn(page, el("scatter plot icon on toolbox")));
    await test.step("Then scatter plot viewer should be added to the open tableview", () => viewerAdded(page, "scatter plot"));
  });
  test("Other viewers the same way [viewer=histogram]", {tag: ["@platform", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user is logged in", () => loggedIn(page));
    await test.step("And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await test.step("When user opens toolbox", () => openToolbox(page));
    await test.step("And user clicks on histogram icon on toolbox", () => clickOn(page, el("histogram icon on toolbox")));
    await test.step("Then histogram viewer should be added to the open tableview", () => viewerAdded(page, "histogram"));
  });
  test("Other viewers the same way [viewer=bar chart]", {tag: ["@platform", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user is logged in", () => loggedIn(page));
    await test.step("And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await test.step("When user opens toolbox", () => openToolbox(page));
    await test.step("And user clicks on bar chart icon on toolbox", () => clickOn(page, el("bar chart icon on toolbox")));
    await test.step("Then bar chart viewer should be added to the open tableview", () => viewerAdded(page, "bar chart"));
  });
});
