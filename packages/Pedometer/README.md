# Pedometer

Pedometer is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform.
It showcases the ability to develop plugins that include scientific computations and a rich UI
that leverages the platform's capabilities.

Here are the files of particular interest:

* [scripts/detect_steps.py](https://github.com/datagrok-ai/public/blob/master/packages/Pedometer/scripts/detect_steps.py)
  : a Python function that detects steps based on raw accelerometer data.
* [package.js](https://github.com/datagrok-ai/public/blob/master/packages/Pedometer/package.js)
  : a JavaScript function that accepts columns containing `x`, `y`, and `z` acceleration components, calculates
  position of steps by invoking `detect_steps`, visualizes the raw data on a chart, and customizes
  chart's interactivity.

## Detecting steps

Let's take a detailed look at the header of the
[scripts/detect_steps.py](https://github.com/datagrok-ai/public/blob/master/packages/Pedometer/scripts/detect_steps.py) file.
It addition to containing information such as name and
description, it instructs the platform what the input and output parameters are, including their data types
and semantic types:

```
# name: detect steps
# description: detects positions of steps based on accelerometer data
# language: python
# tags: template, demo, accelerometer
# sample: sensors/accelerometer.csv
# input: dataframe accel [accelerometry data table]
# input: column x {semtype: accelerometer-x} [x axis]
# input: column y {semtype: accelerometer-y} [y axis]
# input: column z {semtype: accelerometer-z} [z axis]
# input: double sample_rate = 32 [sample rate, in hz]
# output: dataframe steps {action:join(accel)} [steps positions]
```

## Javascript

```js
//input: dataframe accel [Accelerometry data table]
//input: column x {semType: Accelerometer-X} [X axis]
//input: column y {semType: Accelerometer-Y} [Y axis]
//input: column z {semType: Accelerometer-Z} [Z axis]
//input: column timeOffset {semType: Time-Offset} [Time offset column]
pedometer(accel, x, y, z, timeOffset)
```

See also:

* [Grok API](https://datagrok.ai/help/develop/js-api)
* [Packages](https://datagrok.ai/help/develop/develop#packages)
* [Scripting](https://datagrok.ai/help/compute/scripting)
* [Info panels](https://datagrok.ai/help/discover/info-panels)
* [Semantic types](https://datagrok.ai/help/discover/semantic-types)
