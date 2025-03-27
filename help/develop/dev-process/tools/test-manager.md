---
title: "Test Manager"
---


# Test Manager Tour

Test manager is a tool that can be useful in test reproducing and debugging. It can be used to run any package test. To open `Test Manager`, you should run it from `Browser` using the next path: `Apps` > `Admin` > `Test Manager`. 

![TM open](tm-open.gif)

### Test Manager Options

Test Manager options:
| Option          | Description                                                                                   |
|-----------------|-----------------------------------------------------------------------------------------------|
| **Debug**       | Enables debug point before Test execution (Useless without **Browser Inspector (F12)**)       |
| **Benchmark**   | Runs test in Benchmark mode                                                                   |
| **Run Skipped** | Enables skipped test execution                                                                |

Context Menu options:
| Option                | Description                                                                             |
|-----------------------|-----------------------------------------------------------------------------------------|
| **Run**               | Execute test                                                                            |
| **Profile**           | Profile test                                                                            |
| **Copy**              | Copy name to the clipboard                                                              |
| **Copy URL**          | Copy the test's link to the clipboard                                                   |
| **Copy URL(running)** | Copy the test's link to the clipboard. Opening this link would invoke test              |
| **Run Force**         | Execute test even if it is marked as skipped                                            |

## Tests Execution

You can easily navigate in the Test Manager hierarchy to specify test or category execution. To invoke a test or category you can use the `Run` button or **Right Click** > `Run`. After the test run results appear in the **Context Panel**, you can easily see their tests' statuses and more data(`memoryDelta`, `ms`, `widgetsDelta`, `owner`) that can be useful. 

![TM open](tm-test-execution.gif)

### Tests Debugging

The test can be easily debugged by the Test Manager. To do this, you have to select the option `Debug` under the Search input and open **Browser Inspector(F12)**. The next test execution would be stopped on test invocation so to debug it you would need just Step Into it. 

![TM open](tm-debug-usage.gif)

### Tests Profiling 

Test Manager also offers the ability to profile tests and get more information about the test's execution. You should open **Browser Inspector(F12)** and **Right Click** > `Profile` to profile the test. After that, you would be able to see stages of the test's execution in the `Performance` tab of Inspector.

![TM open](tm-profiling.gif)

