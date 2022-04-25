<!-- TITLE: Test packages -->
# Package testing
To make sure that a package works correctly it need to properly tested. Each package should include a bunch of unit tests responsible for different aspects of te package such as UI or logic underneath. And developer should be able to easily run tests at any time (while developing a package or to perform regression testing).

## Test manager

'Test manager' is a tool within Datagrok platform which provides a convenient interface to select and run package unit tests with further esults exploration.

To start 'Test manager' go to top menu Tools -> Dev -> Test manager
![Test manager start](test-manager-start.png)

After starting the tool you will see a list of all package tests divided by package name. Inside each package tests are divided by category.
Using checkboxes you can choose which tests you want to run at a time. You can choose either whole package or required category or exact tests inside a category. After all required tests are selected click on `RUN` button at the top.
![Tests list](test-manager-tests-list.png)

Failed tests will be marked with red cancel sign while passed test will be marked with green tick mark.
![Tests list](test-manager-results.png)

You can get more detailed information on tests results by clicking on a test/category/package name. Information will be shown on a property panel.
![Tests property panel](test-manager-property-panel.png)

