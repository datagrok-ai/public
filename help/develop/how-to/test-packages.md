<!-- TITLE: Test packages -->
# Package testing
To make sure that a package works correctly it needs to be properly tested. Each package should include a bunch of unit tests responsible for different aspects of the package such as UI or logic underneath. And developers should be able to easily run tests at any time (while developing a package or to perform regression testing).

## Test manager

'Test manager' is a tool within Datagrok platform which provides a convenient interface to select and run package unit tests with further results exploration. 'Test manager' itself is a part of DevTools package.

To start 'Test manager' go to top menu Tools -> Dev -> Test manager

![Test manager start](test-mngr-start.png)

After starting the tool you will see a list of all package tests divided by package name. Inside each package the tests are divided by category.
Using check boxes you can choose which tests you want to run at a time. You can choose either whole package or required category or exact tests inside a category. After all required tests are selected click on `RUN` button at the top.

![Tests list](test-mngr-tests-list.png)

Failed tests are marked with red cancel sign while passed test are marked with green tick mark.

![Tests list](test-mngr-results.png)

You can get more detailed information on tests results by clicking on a test/category/package name. Information will be shown on a property panel.

![Tests property panel](test-mngr-property-panel.png)


