---
title: FAQ
sidebar_position: 2
---

## Is my data private?

**Q:** What happens to my data when I open a local file in Datagrok?

**A:** When you open a local file in Datagrok (like dragging and dropping a file to your browser), you can analyze it without saving. This data stays in your browser's memory and isn't sent to the server unless you run resource-intensive server-side computations. Your data is gone when you close the browser tab. To save your work, you need to upload it to the server. Note that uploading data does not make it accessible to others. Your data stays private and visible to you only until you explicitly share it. Learn how to [save](../concepts/project/project.md#saving-entities-to-projects) and [share](../navigation/basic-tasks/basic-tasks.md#share) data.

**Q:** What data or telemetry is sent back to Datagrok? What egress ports/protocols are used?

**A:** Datagrok does not send anything without user permissions. If user reports error and explicitly checks "Send report back to Datagrok", email is sent to feedback@datagrok.ai. Datagrok can download images from Docker Hub or packages from NPM.