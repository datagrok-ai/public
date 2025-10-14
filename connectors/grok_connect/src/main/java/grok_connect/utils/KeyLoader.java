package grok_connect.utils;

import org.bouncycastle.asn1.pkcs.PrivateKeyInfo;
import org.bouncycastle.jce.provider.BouncyCastleProvider;
import org.bouncycastle.openssl.PEMKeyPair;
import org.bouncycastle.openssl.PEMParser;
import org.bouncycastle.openssl.jcajce.JcaPEMKeyConverter;
import org.bouncycastle.openssl.jcajce.JceOpenSSLPKCS8DecryptorProviderBuilder;
import org.bouncycastle.operator.InputDecryptorProvider;
import org.bouncycastle.pkcs.PKCS8EncryptedPrivateKeyInfo;

import javax.crypto.Cipher;
import javax.crypto.EncryptedPrivateKeyInfo;
import javax.crypto.SecretKeyFactory;
import javax.crypto.spec.PBEKeySpec;
import java.io.StringReader;
import java.security.Key;
import java.security.KeyFactory;
import java.security.PrivateKey;
import java.security.Security;
import java.security.spec.PKCS8EncodedKeySpec;
import java.util.Base64;

public class KeyLoader {

    /**
     *
     * @param fileContent string content of .pem or .der file. If it's binary .der then we assume that it's BASE64 encoded
     * @param passphrase used to decrypt key, can be null
     * @return PrivateKey
     */
    public static PrivateKey get(String fileContent, String passphrase) throws Exception {
        if (Security.getProvider(BouncyCastleProvider.PROVIDER_NAME) == null)
            Security.addProvider(new BouncyCastleProvider());

        if (fileContent.contains("-----BEGIN")) {
            // Treat as PEM
            PEMParser pemParser = new PEMParser(new StringReader(fileContent));
            Object pemObject = pemParser.readObject();
            pemParser.close();

            PrivateKeyInfo privateKeyInfo;
            if (pemObject instanceof PKCS8EncryptedPrivateKeyInfo) {
                PKCS8EncryptedPrivateKeyInfo encryptedPrivateKeyInfo = (PKCS8EncryptedPrivateKeyInfo) pemObject;
                InputDecryptorProvider pkcs8Prov = new JceOpenSSLPKCS8DecryptorProviderBuilder().build(passphrase.toCharArray());
                privateKeyInfo = encryptedPrivateKeyInfo.decryptPrivateKeyInfo(pkcs8Prov);
            }
            else if (pemObject instanceof PrivateKeyInfo)
                privateKeyInfo = (PrivateKeyInfo) pemObject;
            else if (pemObject instanceof PEMKeyPair) {
                PEMKeyPair pemKeyPair = (PEMKeyPair) pemObject;
                return new JcaPEMKeyConverter()
                        .setProvider(BouncyCastleProvider.PROVIDER_NAME)
                        .getPrivateKey(pemKeyPair.getPrivateKeyInfo());
            }
            else
                throw new IllegalArgumentException("Unsupported PEM object");

            return new JcaPEMKeyConverter().setProvider(BouncyCastleProvider.PROVIDER_NAME).getPrivateKey(privateKeyInfo);
        }
        else {
            // Assume base64-encoded DER
            byte[] decoded = Base64.getDecoder().decode(fileContent);
            EncryptedPrivateKeyInfo encryptedInfo;
            try {
                encryptedInfo = new EncryptedPrivateKeyInfo(decoded);
            } catch (Exception e) {
                // If not encrypted, parse directly
                return KeyFactory.getInstance("RSA").generatePrivate(new PKCS8EncodedKeySpec(decoded));
            }

            Cipher cipher = Cipher.getInstance(encryptedInfo.getAlgName());
            PBEKeySpec pbeKeySpec = new PBEKeySpec(passphrase.toCharArray());
            SecretKeyFactory keyFactory = SecretKeyFactory.getInstance(encryptedInfo.getAlgName());
            Key pbeKey = keyFactory.generateSecret(pbeKeySpec);
            cipher.init(Cipher.DECRYPT_MODE, pbeKey, encryptedInfo.getAlgParameters());

            byte[] decryptedKey = encryptedInfo.getKeySpec(cipher).getEncoded();
            return KeyFactory.getInstance("RSA").generatePrivate(new PKCS8EncodedKeySpec(decryptedKey));
        }
    }
}
