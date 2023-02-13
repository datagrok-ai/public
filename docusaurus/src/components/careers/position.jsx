import styles from './styles.module.css';
import React from 'react';

export default function profileCard ({children, meta}){
    return(
        <div className={"col-lg-6 col-md-6 col-sm-6 col-6 my-3"}>
            <div className="card hover-shadow border-0 border-0 border-light shadow-sm h-100">
                <a className="btn link" data-toggle="modal" data-target={"#"+meta.id}>
                    <div className="card-body text-left pt-0">
                        <div className="mt-3">
                            <div className="h5 mt-3 pb-2">{meta.name}</div>
                            <div className="d-flex flex-row">
                                <span className={'badge badge-pill mr-2 '+ styles.badgedg }>{meta.type}</span>
                                <span className={'badge badge-pill badge-info '+ styles.badgedg }>{meta.location}</span>
                            </div>
                        </div>
                    </div>
                </a>
            </div>
            <div className="modal fade" tabIndex="-1" id={meta.id} aria-hidden="true">
                <div className={"modal-dialog modal-dialog-scrollable "+ styles.dialogSize}>
                    <div className={"modal-content "}>
                        <div className="modal-body text-center ">
                            <button type="button" className="close" data-dismiss="modal" aria-label="Close">
                                <i className="fal fa-times fa-lg text-secondary"></i>
                            </button>
                            <section className="container-fluid">
                                <div className="container">
                                    <div className="row text-left">
                                        <div className='col-lg-12 col-md-12'>
                                                <h1 className="">{meta.name}</h1>
                                                <div className="d-flex flex-row">
                                                    <span className={'badge badge-pill mr-2 '+ styles.badgedg }>{meta.type}</span>
                                                    <span className={'badge badge-pill badge-info '+ styles.badgedg }>{meta.location}</span>
                                                </div>
                                                <br></br>                                        </div>
                                        <div className="col-lg-12 col-md-12 text-left py-3">
                                        {children}
                                        <br></br>
                                        <a className='btn btn-primary px-4' href={"mailto:hr@datagrok.ai?subject="+meta.name}>Apply</a>
                                        </div>
                                    </div>
                                </div>
                            </section>    
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}