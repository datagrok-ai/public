set DATAGROK_VERSION=1.14.2
set GROK_CONNECT_VERSION=2.0.13

rem set DATAGROK_VERSION=bleeding-edge

docker-compose -f docker/localhost.docker-compose.mlb.yaml ^
  --project-name datagrok ^
  --profile datagrok ^
  --profile db ^
  --profile grok_connect ^
  --profile scripting ^
  --profile grok_spawner ^
  up -d