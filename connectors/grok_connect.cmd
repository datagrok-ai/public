:: Base project directory and executable
SET GROK_CONNECT_DIR=grok_connect
SET GROK_CONNECT=grok_connect-1.0.3.jar

:: Remove target
del /s /q %GROK_CONNECT_DIR%\target

:: Get dart test dependencies
cd ..\ddt
call pub get

cd ..\grok_shared
call pub get

cd ..\ddt\bin\serialization_test
call pub get

cd ..\..\..\public\connectors

:: Build project
call mvn package

:: Run connector server with shared libraries
call java -Xmx4g -classpath %GROK_CONNECT_DIR%\lib\*;%GROK_CONNECT_DIR%\target\%GROK_CONNECT% grok_connect.GrokConnect
pause
