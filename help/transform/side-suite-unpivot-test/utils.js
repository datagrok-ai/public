async function waitForWindow(driver, handles, timeout) {
  await driver.sleep(timeout);
  const hndls = await driver.getAllWindowHandles();
  if (hndls.length > handles.length) {
    return hndls.find(h => (!handles.includes(h)));
  }
  throw new Error("New window did not appear before timeout");
}

module.exports = {
  waitForWindow
};