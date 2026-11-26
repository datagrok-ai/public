/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/01-getting-started/02-open-your-data.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, expand, pasteInto, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Basics", () => {
  const session = feature(test, "features/howto/01-getting-started/02-open-your-data.feature", import.meta.url);
  test("Open your data", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Browse: the files and databases on the server");
    await session.step(19, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")), undefined, "Files > Demo: sample data that comes with Datagrok");
    await session.step(20, "And user expands Files---Demo tree node inside browse tree", () => expand(page, el("Files---Demo tree node inside browse tree")));
    await session.step(22, "And user clicks on Files---Demo---cars.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")), undefined, "One click previews a file");
    await session.step(24, "And user double-clicks on Files---Demo---cars.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")), undefined, "A double-click opens it as a table");
    await session.step(25, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(26, "And the table should have 30 rows", () => rowCount(page, 30));
    await session.step(28, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "A file on your computer: the open-file icon of Browse");
    await session.step(29, "And user uploads \"features/howto/01-getting-started/files/assay-results.csv\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "features/howto/01-getting-started/files/assay-results.csv", el("\"Open local file\" icon inside browse toolbar")));
    await session.step(30, "Then the \"assay-results\" view should be current", () => viewIsCurrent(page, "assay-results"));
    await session.step(31, "And the table should have 24 rows", () => rowCount(page, 24));
    await session.step(33, "When user drops the \"features/howto/01-getting-started/files/compound-yields.csv\" file of the project onto grid", () => dropFile(page, "features/howto/01-getting-started/files/compound-yields.csv", el("grid")), undefined, "Or drag a file from your computer onto the window");
    await session.step(34, "Then the \"compound-yields\" view should be current", () => viewIsCurrent(page, "compound-yields"));
    await session.step(35, "And the table should have 12 rows", () => rowCount(page, 12));
    await session.step(37, "When user clicks on browse tab", () => clickOn(page, el("browse tab")), undefined, "Copied cells or text: the Open text icon");
    await session.step(38, "And user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(39, "Then the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(41, "When user clicks on code editor", () => clickOn(page, el("code editor")), undefined, "Paste them in; the preview below follows");
    await session.step(42, "And user pastes \"compound,yield\\nCMP-201,72.5\\nCMP-202,65.0\\nCMP-203,81.2\\nCMP-204,58.9\" into code editor", () => pasteInto(page, "compound,yield\\nCMP-201,72.5\\nCMP-202,65.0\\nCMP-203,81.2\\nCMP-204,58.9", el("code editor")));
    await session.step(44, "And user clicks on DONE button", () => clickOn(page, el("DONE button")), undefined, "DONE opens it as a table");
    await session.step(45, "Then the table should have 4 rows", () => rowCount(page, 4));
    await session.step(47, "When user clicks on \"cars\" view tab", () => clickOn(page, el("\"cars\" view tab")), undefined, "Each table has its own tab: click one to go back");
    await session.step(48, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
  });
});
