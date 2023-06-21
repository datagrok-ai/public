package grok_connect.providers.utils;

import java.lang.annotation.Documented;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Documented
@Target(ElementType.METHOD)
@Retention(RetentionPolicy.RUNTIME)
public @interface Sql {
    /**
     *Path to init sql script
     */
    String path();

    /**
     *Path to sql script to restore data after method is proceeded
     */
    String restorePath() default "";
}
