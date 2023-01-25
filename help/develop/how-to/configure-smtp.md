---
title: "SMTP configuration"
---

Datagrok supports [Mailgun email delivery platform](https://www.mailgun.com/) and different SMTP servers,
including [Amazon SES](https://aws.amazon.com/ses/). Configure a local SMTP server or use a cloud solution based on your
needs.

To configure email delivery for Datagrok:

1) Go to the Datagrok Settings section 'Admin'
2) Configure Sender Email from which users will get emails
3) Set Web Root, which is the root URL of the Datagrok instance
4) Set Api Root, which is the same as the root URL with the `/api` suffix
5) Email Service:
    1) Set Mailgun if you use their integration
        1) Configure Mailgun Domain and Mailgun key got from the Mailgun interface
    2) In all other cases, set SMTP
        1) Configure SMTP server address/DNS name. If you want to use the host SMTP server with dockerized Datagrok
           set `host.docker.internal`
        2) Set SMTP server port
        3) Set the SMTP user and password
        4) Use SMTP anonymous mode if you want to use an anonymous SMTP server
        5) Check SMTP secured to use [SMTPS](https://en.wikipedia.org/wiki/SMTPS)

![SMTP configuration](../img/smtp.png)
