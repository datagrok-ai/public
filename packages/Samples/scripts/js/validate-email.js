//name: validateEmail
//help-url: https://datagrok.ai/help/datagrok/concepts/functions/func-params-annotation
//language: javascript
//input: string email
//output: string invalidReason
invalidReason = /^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(email) ? null : 'Invalid email';
