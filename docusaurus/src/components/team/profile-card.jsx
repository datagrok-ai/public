import styles from './styles.module.css';
import React from 'react';

export default function profileCard ({children, meta}){
    return(
        <div className={"col-lg-3 col-md-3 col-sm-6 col-6 my-3 " + styles.card}>
            <div className="card hover-shadow border-0 border-0 border-light shadow-sm h-100">
                <a className="btn link" data-toggle="modal" data-target={"#"+meta.id}>
                    <img src={meta.photo} className="rounded-circle mt-4" height="96" width="96"></img>
                    <div className="card-body pt-0">
                        <div className="my-3">
                            <div className="h6 text-dark mt-3 mb-1">{meta.name}</div>
                            <div className={"font-size-sm mb-2 " + styles.greyColor}>{meta.position}</div>
                        </div>
                    </div>
                </a>
            </div>
            <div className="modal fade" tabIndex="-1" id={meta.id} aria-hidden="true">
                <div className={"modal-dialog modal-dialog-scrollable " + styles.dialogSize}>
                    <div className="modal-content">
                        <div className="modal-body text-center ">
                            <button type="button" className="close" data-dismiss="modal" aria-label="Close">
                                <i className="far fa-times fa-xs text-secondary"></i>
                            </button>
                            <img src={meta.photo} className="rounded-circle mt-4 " height="124" width="124"></img>
                                <h4 className="modal-title mt-2" id="exampleModalLabel">{meta.name}</h4>
                                <p className="font-size-sm" style={{color:'var(--blue-2)'}}>{meta.position}</p>
                                {meta.linkedin != undefined ? <a href={meta.linkedin} target="_blank"><i className="fab fa-linkedin-in fa-lg mx-1"></i></a>:null}
                                {meta.github != undefined ? <a href={meta.github} target="_blank"><i className="fab fa-github fa-lg mx-1"></i></a>:null}
                                <div className="text-left p-3">
                                    <div className="text-secondary mb-3">
                                     {children}
                                    </div>
                                </div>

                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}