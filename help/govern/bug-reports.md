---
title: "Feedback"
sidebar_position: 5
---

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
```

Datagrok provides several ways to report issues, send feedback, and share diagnostic
information:
- **User-driven reports**:
  - **Report an error**: Use this when unexpected behavior occurs and handled or unhandled exceptions appear.
  
  <details>
  <summary>How to report</summary>
  <Tabs groupId="exceptions">
  
  <TabItem value="handled" label="Handled exception">
  
  Handled exception appear as a red message box in the upper right corner. Hover over it and click the **Report** icon to open the report dialog.
  
  ![Reporting of handled exception](img/handled_exception.gif)
  
  </TabItem>
  
  <TabItem value="unhandled" label="Unhandled exception">

  Unhandled exception indicates internal platform error and is indicated by a **red exclamation icon** on the **Sidebar**. Click the icon to open the report dialog.

  ![Reporting of unhandled exception](img/unhandled_exception.png)

  </TabItem>

  </Tabs>
  </details>

  - **Submit feedback**: Use this to raise non-exception issues, such as UX concerns, questions, or feature requests, 
    or as an alternative way to report an error.

  <details>
  <summary>How to submit</summary>

  Click the **Datagrok Logo Menu** on the **Sidebar**, select **Help > Feedback**, write your feedback, check **Include logs** if needed, and submit your message.

  ![Submit feedback](img/feedback.gif)

  </details>

- **Auto-reporting system**: Automatically log internal exceptions when enabled (see [Configuring error reporting system](#configuring-error-reporting-system)).
  Logs are stored in [Usage Analysis](audit/usage-analysis.md) and are never sent externally.
  Auto-reports help track recurring errors, monitor instance stability, and support in-depth debugging.

## Report an error

The **Report an error** dialog contains multiple tabs with diagnostic information from recent user actions. 
Before submitting, you can review, edit, or remove data — you remain in full control of what is included and sent.
No diagnostic data leaves your environment unless you explicitly allow it.

- **SUMMARY**: 
  - **Description**: Add context to help reproduce the issue or clarify expected behavior.  
  - **Screenshot**: Hover over it and click **REMOVE** to exclude sensitive image, if needed.  
  - **Email** checkbox: Select to send the report to the specified email addresses 
    (see [Configuring error reporting system](#configuring-error-reporting-system)).
    If unchecked, the report not sent externally, is stored in the system, and can be viewed via the **Usage Analysis** application.

- **DETAILS**, **ERRORS**, **TIMELINES**, **CONSOLE**, **CLIENT**, **SERVER**, **SERVICES**: Include technical traces, execution logs, and runtime diagnostics to assist in root cause analysis.

- **TABLES**: Select which open tables’ data and metadata to include. The table’s title, row count, and column count are always included by default.

:::note Note
If the reporting system is not configured, users can save the report as a JSON file via the **Save as json** icon and share it with the support team.

:::

![Report an error](img/report-an-error.gif)

## Configuring error reporting system

Administrators can configure error reporting under **Settings > Admin > Error reporting**:

- **Report email:** Specify one or more email addresses (comma-separated) to receive user-submitted reports. `feedback@datagrok.ai` is recommended.
- **Auto report errors:** Enable to automatic logging of internal exceptions.
  Automatic reports are **never sent externally**, are stored in the system, and can be viewed via the **Usage Analysis** application.

![Configuring error reporting system](img/сonfiguring-reporting-system.gif)
