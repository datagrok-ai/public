// First uncomment this errorneous line in the script:
// 2 ++ 2;
// Hit save the script like that, you should now see an error.
// Then comment this bad line above and run the script again.
// The alert will contain the error's text.

alert(await grok.shell.lastError);