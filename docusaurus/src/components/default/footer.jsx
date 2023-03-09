import styles from './styles.module.css';
import React from 'react';

export default function footer (){
    return (
        <section>
            <div className="contact-us-wrapper">
                <div className="contact-us-button contact-us">request a demo</div>
            </div>
            <footer className="py-3 px-4 text-light footer-landing">
                <div className="container-fluid">
                    <div className="row">
                        <div className="col-md-6"></div>
                        <div className="col-md-6 text-right">
                            Datagrok, Inc<br/>
                                <a className="text-light" href="mailto:info@datagrok.ai">info@datagrok.ai</a>
                        </div>
                    </div>
                </div>
            </footer>
        </section>
    );
}