jest.setMock('selenium-webdriver', webdriver);
// This file was generated using Selenium IDE
const tests = require("./commons.js");
global.Key = require('selenium-webdriver').Key;
global.URL = require('url').URL;
global.BASE_URL = configuration.baseUrl || 'https://dev.datagrok.ai/';
let vars = {};
jest.setTimeout(4350000);
describe("build-query-join-test", () => {
  it("build-query-join-test", async () => {
    await tests["build-query-join-test"](driver, vars);
    expect(true).toBeTruthy();
  });
});
beforeEach(() => {
  vars = {};
});
afterEach(async () => (cleanup()));