FROM tomcat:8.0.48-jre8

RUN sed -i 's/port="8080"/port="8081"/' $CATALINA_HOME/conf/server.xml

WORKDIR $CATALINA_HOME/tmp

ENV HELM2_WEBSERVICE_VERSION=2.2.0

RUN curl -o $CATALINA_HOME/webapps/WebService.war "https://oss.sonatype.org/service/local/repositories/releases/content/org/pistoiaalliance/helm/helm2-webservice/${HELM2_WEBSERVICE_VERSION}/helm2-webservice-${HELM2_WEBSERVICE_VERSION}.war"

RUN cp -r $CATALINA_HOME/tmp $CATALINA_HOME/webapps

CMD ["catalina.sh", "run"]
