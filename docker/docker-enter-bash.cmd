set SVC=%1
if "%SVC%"=="" set SVC=jkg_python
docker exec -u root -it datagrok-%SVC%-1 bash
