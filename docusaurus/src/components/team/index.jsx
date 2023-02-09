import styles from './styles.module.css';
import React from 'react';

export function LayoutHeader ({children, meta}){
    return(
        <section className="container-fluid bg-light py-5">
                <div className="container">
                    <div className="row">
                        <div className="row py-5 my-3 align-items-center">
                            <div className="col-lg-7 mb-4 px-5">
                                <p className="display-5 text-secondary">{meta.description}</p>
                                <a href="https://datagrok.ai/company/careers" className="mt-3 btn btn-primary px-">JOIN OUR TEAM</a>
                            </div>
                            <div className="col-lg-5 my-3">
                                <iframe width="100%" height="280" className="rounded pt-2" src={meta.video}></iframe>
                            </div>
                        </div>
                    </div>
                </div>
        </section>
    );
}
export function LayoutBody ({children, meta}){
    return(

        <div className="row justify-content-md-center text-center">
            <div className="row py-3 align-items-center">
                <div className="col-lg-12">
                    {children}
                </div>
            </div>
        </div>
    )
}