---
title: "Configure SMTP"
sidebar_position: 2
description: Configure Mailgun, Amazon SES, or an SMTP server so Datagrok can send signup, confirmation, and password-reset emails.
keywords:
  - smtp server setup
  - mailgun integration
  - amazon ses
  - email delivery
  - sender email address
---

Datagrok supports [Mailgun email delivery platform](https://www.mailgun.com/), [Amazon SES](https://aws.amazon.com/ses/),
different SMTP servers, and an unauthenticated SMTP relay (Datagrok mailer). Configure a local SMTP server or use a cloud solution based on your needs.

To configure email delivery for Datagrok:

1. Go to the Datagrok Settings section 'Admin'  
2. Configure Sender Email from which users will get emails  
3. Set Web Root, which is the root URL of the Datagrok instance  
4. Set Api Root, which is the same as the root URL with the `/api` suffix  
5. Email Service:  
    - Set Mailgun if you use their integration
        - Configure Mailgun Domain and Mailgun key got from the Mailgun interface
    - Set SES to send through Amazon SES without SMTP credentials. Datagrok calls the SES v2 API with the
      IAM role of its pod or instance. That role needs the `ses:SendEmail` permission, and Sender Email must be
      a verified SES identity
        - Set AWS Region to the region of the SES identity. If empty, Datagrok uses the `AWS_REGION` environment variable
        - Optionally, set Configuration set to an SES configuration set for event tracking
    - Set Datagrok mailer to relay email through an SMTP relay without authentication, such as a mail
      container in your deployment. Set Mailer URL to the relay's host name. Datagrok connects to it on port 25
      without TLS
    - In all other cases, set SMTP
        - Configure SMTP server address/DNS name. If you want to use the host SMTP server with dockerized Datagrok
           set `host.docker.internal`
        - Set the SMTP server port. In most cases is 25
        - Set the SMTP user and password
        - Use SMTP anonymous mode if you want to use an anonymous SMTP server
        - Check SMTP is secured to use [SMTPS](https://en.wikipedia.org/wiki/SMTPS)

![SMTP configuration](smtp.png)
