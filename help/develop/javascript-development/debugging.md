# Debugging

This document discusses debugging of your packages for Datagrok in two popular IDEs &mdash; Visual Studio Code and 
JetBrains WebStorm.

## One-click debugging with Visual Studio Code

To configure VS Code for Datagrok development, run the `grok create` command with the `--ide=vscode` flag:

```shell
grok create my-datagrok-package --ide=vscode
```

**IMPORTANT**: The `--ide` flag only works on Windows. If you're using Linux or macOS, refer to the
[Creating a VS Code Configuration for Linux or macOS section].

The Grok CLI will create an additional `.vscode` folder with some useful configurations. Here are two more steps to
configure VS Code for debugging:

* In VS Code, open **Run* > **Install additional debuggers...**, find and install **Debugger for Chrome**, and then 
  restart VS Code.
* In VS Code, select **Activate View** > **Appearance** > **Show activity bar** to bring the Run button into the view. 
  You can also use Ctrl+Shift+D to run the application.

After that, in VS Code, click the **Run and Debug** icon > **Run and Debug**. The Grok CLI will prepare a webpack
package locally, showing errors if any, deploy it to the [configured Datagrok instance], and then run the browser in
debug mode.

The first time you run debugging, you need to enter your Datagrok credentials into Chrome. After you entered credentials,
close Chrome and restart debugging to see locals and stack traces as usual.

**See more**: [Debugging with VS Code video]

## Debugging with JetBrains IDEs

### Debugging with a Shell Script

To configure debugging with a shell script on WebStorm:

1. Select **Run** > **Edit Configurations...**.

![WebStorm: Edit Configurations](webstorm-debugging-01.png)

2. In the **Run/Debug Configurations** dialog, select `+`, and then select **Shell Script**.

![WebStorm: Adding a Shell Script configuration](webstorm-debugging-02.png)

3. Add the name for your configuration, and then add the `/c call webpack && grok publish dev && call echo` script.

![WebStorm: Shell Script configuration content](webstorm-debugging-03.png)

This script publishes your package to the Dev server. You can configure this or other scripts to publish the package to
the public server or run the package locally in a Docker container.

**See more**: [Datagrok configuration](./datagrok-configuration.md)

### JavaScript debug configuration

To configure JS Debug on Webstorm:

1. Select **Run** > **Edit Configurations...**.

![WebStorm: Edit Configurations](webstorm-debugging-01.png)

2. In the **Run/Debug Configurations** dialog, select `+`, and then select **JS Debug**.

3. Add a configuration name, and then add the previously created shell script in the **Before launch** section.

![WebStorm: JavaScript Debug configuration content](webstorm-debugging-05.png)

## Troubleshooting

Sometimes JavaScript debugging from the WebStorm IDE becomes impossible after it initially worked. This is a [known 
issue in WebStorm]. To fix this issue, remove these two files and restart the IDE:

* `%USERPROFILE%\AppData\Roaming\JetBrains\WebStorm2020.1\options\web-browsers.xml`
* `%USERPROFILE%\AppData\Roaming\JetBrains\WebStorm2020.1\options\other.xml`

Additional information:

* Debugging JavaScript in WebStorm: [video 1], [video 2]
* IntelliJ IDEA JavaScript debugging: [link]

## Source-Based Packages

Deploying such package locates it to the Datagrok host URI (such as `https://dev.datagrok.ai`) under
`api > packages/published/flies > <PACKAGE_NAME>/<VERSION>/_/<DIGIT>`, where you'd set breakpoints.

## Troubleshooting Debugging

1. Make sure that the webpack `devtool` configuration is set to `inline-source-map` for development. If this option 
   isn't set, the source maps aren't generated, and the IDE will not get source code locations.
2. Make sure the required plugins and debuggers for Chrome are installed in your IDE.

[Creating a VS Code Configuration Manually section]: #creating-a-vs-code-configuration-for-debugging-manually
[Debugging with VS Code video]: https://youtu.be/zVVmlRorpjg?list=PLIRnAn2pMh3kvsE5apYXqX0I9bk257_eY&t=871
[known issue in WebStorm]: https://intellij-support.jetbrains.com/hc/en-us/community/posts/360009567459-Webstorm-2020-2-1-Remote-Debugging-do-not-work
[known issue]: https://youtrack.jetbrains.com/issue/IDEA-229467
[JetBrains IDE Support plugin is no longer required]: https://intellij-support.jetbrains.com/hc/en-us/community/posts/360010507240-where-is-JETBRAINS-IDE-SUPPORT-chrome-extension-it-cant-be-found-anywhere-now-on-the-internet
[video 1]: https://www.youtube.com/watch?v=Qcqnmle6Wu8
[video 2]: https://www.youtube.com/watch?v=YNNDMpoGV0w
[link]: https://www.jetbrains.com/help/idea/debugging-javascript-in-chrome.html
[configured Datagrok instance]: datagrok-configuration.md