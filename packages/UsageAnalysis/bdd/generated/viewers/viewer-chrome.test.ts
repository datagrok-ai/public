/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/viewer-chrome.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.chrome]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Viewer title and description", () => {
  const session = feature(test, "features/viewers/viewer-chrome.feature", import.meta.url);
  test("bar chart shows and clears its title and description [viewer=bar chart]", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Given user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(14, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Show Title","true"],["Title","Demographics"]]));
    await session.step(17, "Then title of bar chart viewer should have text \"Demographics\"", () => shouldHaveText(page, el("title of bar chart viewer"), "Demographics"));
    await session.step(18, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Description","By race"],["Description Visibility Mode","Always"]]));
    await session.step(21, "Then description of bar chart viewer should have text \"By race\"", () => shouldHaveText(page, el("description of bar chart viewer"), "By race"));
    await session.step(22, "When user sets \"Description Position\" property of bar chart viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("bar chart viewer"), "Bottom"));
    await session.step(23, "Then description of bar chart viewer should be visible", () => shouldBe(page, el("description of bar chart viewer"), "visible"));
    await session.step(24, "When user sets \"Description Position\" property of bar chart viewer to \"Left\"", () => setProperty(page, "Description Position", el("bar chart viewer"), "Left"));
    await session.step(25, "Then description of bar chart viewer should be visible", () => shouldBe(page, el("description of bar chart viewer"), "visible"));
    await session.step(26, "When user sets \"Description Position\" property of bar chart viewer to \"Right\"", () => setProperty(page, "Description Position", el("bar chart viewer"), "Right"));
    await session.step(27, "Then description of bar chart viewer should be visible", () => shouldBe(page, el("description of bar chart viewer"), "visible"));
    await session.step(28, "When user sets \"Description Visibility Mode\" property of bar chart viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("bar chart viewer"), "Never"));
    await session.step(29, "Then description of bar chart viewer should be absent", () => shouldBe(page, el("description of bar chart viewer"), "absent"));
    await session.step(30, "When user sets properties of bar chart viewer:", () => setProperties(page, el("bar chart viewer"), [["Show Title","false"],["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]));
    await session.step(36, "Then description of bar chart viewer should be absent", () => shouldBe(page, el("description of bar chart viewer"), "absent"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("box plot shows and clears its title and description [viewer=box plot]", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Given user adds a box plot viewer", () => addViewer(page, "box plot"));
    await session.step(14, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Title","true"],["Title","Demographics"]]));
    await session.step(17, "Then title of box plot viewer should have text \"Demographics\"", () => shouldHaveText(page, el("title of box plot viewer"), "Demographics"));
    await session.step(18, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Description","By race"],["Description Visibility Mode","Always"]]));
    await session.step(21, "Then description of box plot viewer should have text \"By race\"", () => shouldHaveText(page, el("description of box plot viewer"), "By race"));
    await session.step(22, "When user sets \"Description Position\" property of box plot viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("box plot viewer"), "Bottom"));
    await session.step(23, "Then description of box plot viewer should be visible", () => shouldBe(page, el("description of box plot viewer"), "visible"));
    await session.step(24, "When user sets \"Description Position\" property of box plot viewer to \"Left\"", () => setProperty(page, "Description Position", el("box plot viewer"), "Left"));
    await session.step(25, "Then description of box plot viewer should be visible", () => shouldBe(page, el("description of box plot viewer"), "visible"));
    await session.step(26, "When user sets \"Description Position\" property of box plot viewer to \"Right\"", () => setProperty(page, "Description Position", el("box plot viewer"), "Right"));
    await session.step(27, "Then description of box plot viewer should be visible", () => shouldBe(page, el("description of box plot viewer"), "visible"));
    await session.step(28, "When user sets \"Description Visibility Mode\" property of box plot viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("box plot viewer"), "Never"));
    await session.step(29, "Then description of box plot viewer should be absent", () => shouldBe(page, el("description of box plot viewer"), "absent"));
    await session.step(30, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Show Title","false"],["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]));
    await session.step(36, "Then description of box plot viewer should be absent", () => shouldBe(page, el("description of box plot viewer"), "absent"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("histogram shows and clears its title and description [viewer=histogram]", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Given user adds a histogram viewer", () => addViewer(page, "histogram"));
    await session.step(14, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Show Title","true"],["Title","Demographics"]]));
    await session.step(17, "Then title of histogram viewer should have text \"Demographics\"", () => shouldHaveText(page, el("title of histogram viewer"), "Demographics"));
    await session.step(18, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Description","By race"],["Description Visibility Mode","Always"]]));
    await session.step(21, "Then description of histogram viewer should have text \"By race\"", () => shouldHaveText(page, el("description of histogram viewer"), "By race"));
    await session.step(22, "When user sets \"Description Position\" property of histogram viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("histogram viewer"), "Bottom"));
    await session.step(23, "Then description of histogram viewer should be visible", () => shouldBe(page, el("description of histogram viewer"), "visible"));
    await session.step(24, "When user sets \"Description Position\" property of histogram viewer to \"Left\"", () => setProperty(page, "Description Position", el("histogram viewer"), "Left"));
    await session.step(25, "Then description of histogram viewer should be visible", () => shouldBe(page, el("description of histogram viewer"), "visible"));
    await session.step(26, "When user sets \"Description Position\" property of histogram viewer to \"Right\"", () => setProperty(page, "Description Position", el("histogram viewer"), "Right"));
    await session.step(27, "Then description of histogram viewer should be visible", () => shouldBe(page, el("description of histogram viewer"), "visible"));
    await session.step(28, "When user sets \"Description Visibility Mode\" property of histogram viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("histogram viewer"), "Never"));
    await session.step(29, "Then description of histogram viewer should be absent", () => shouldBe(page, el("description of histogram viewer"), "absent"));
    await session.step(30, "When user sets properties of histogram viewer:", () => setProperties(page, el("histogram viewer"), [["Show Title","false"],["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]));
    await session.step(36, "Then description of histogram viewer should be absent", () => shouldBe(page, el("description of histogram viewer"), "absent"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("pc plot shows and clears its title and description [viewer=pc plot]", {tag: ["@viewers", "@realizes:viewers.chrome"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "Given user adds a pc plot viewer", () => addViewer(page, "pc plot"));
    await session.step(14, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Title","true"],["Title","Demographics"]]));
    await session.step(17, "Then title of pc plot viewer should have text \"Demographics\"", () => shouldHaveText(page, el("title of pc plot viewer"), "Demographics"));
    await session.step(18, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Description","By race"],["Description Visibility Mode","Always"]]));
    await session.step(21, "Then description of pc plot viewer should have text \"By race\"", () => shouldHaveText(page, el("description of pc plot viewer"), "By race"));
    await session.step(22, "When user sets \"Description Position\" property of pc plot viewer to \"Bottom\"", () => setProperty(page, "Description Position", el("pc plot viewer"), "Bottom"));
    await session.step(23, "Then description of pc plot viewer should be visible", () => shouldBe(page, el("description of pc plot viewer"), "visible"));
    await session.step(24, "When user sets \"Description Position\" property of pc plot viewer to \"Left\"", () => setProperty(page, "Description Position", el("pc plot viewer"), "Left"));
    await session.step(25, "Then description of pc plot viewer should be visible", () => shouldBe(page, el("description of pc plot viewer"), "visible"));
    await session.step(26, "When user sets \"Description Position\" property of pc plot viewer to \"Right\"", () => setProperty(page, "Description Position", el("pc plot viewer"), "Right"));
    await session.step(27, "Then description of pc plot viewer should be visible", () => shouldBe(page, el("description of pc plot viewer"), "visible"));
    await session.step(28, "When user sets \"Description Visibility Mode\" property of pc plot viewer to \"Never\"", () => setProperty(page, "Description Visibility Mode", el("pc plot viewer"), "Never"));
    await session.step(29, "Then description of pc plot viewer should be absent", () => shouldBe(page, el("description of pc plot viewer"), "absent"));
    await session.step(30, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Title","false"],["Title",""],["Description",""],["Description Visibility Mode","Auto"],["Description Position","Top"]]));
    await session.step(36, "Then description of pc plot viewer should be absent", () => shouldBe(page, el("description of pc plot viewer"), "absent"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
});
