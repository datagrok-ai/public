package grok_connect.utils;

import com.amazonaws.auth.AWSCredentials;
import com.amazonaws.auth.AWSCredentialsProvider;
import com.simba.athena.amazonaws.auth.BasicAWSCredentials;
import com.simba.athena.amazonaws.auth.BasicSessionCredentials;

public class CustomAthenaSessionCredentialsProvider implements AWSCredentialsProvider {
    private final AWSCredentials credentials;

    public CustomAthenaSessionCredentialsProvider(String awsAccessKey, String awsSecretKey) {
        this.credentials = new BasicAWSCredentials(awsAccessKey, awsSecretKey);
    }

    public CustomAthenaSessionCredentialsProvider(String awsAccessKey, String awsSecretKey, String sessionToken) {
        credentials = new BasicSessionCredentials(awsAccessKey, awsSecretKey, sessionToken);
    }


    @Override
    public AWSCredentials getCredentials() {
        return credentials;
    }

    @Override
    public void refresh() {

    }
}
