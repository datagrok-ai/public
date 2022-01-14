cd ../libraries
forfiles /S /M package.json /C "cmd /c npm link"
cd ../js-api
npm link
cd ../packages
forfiles /S /M package.json /C "cmd /c npm run link-all"