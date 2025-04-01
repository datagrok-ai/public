import Markdown from 'markdown-to-jsx';
import React, { useEffect, useState, Component } from 'react';

export default function Hero({children, color}){
    return(
        <section className='container-fluid p-5'>
            <div className='section' style={{backgroundColor: color != undefined ? color : 'white'}} >
            <div className="container hero-container pt-5 animate__animated animate__fadeIn" >
                <div className="row">
                    <div className="col-lg-12 col-sm-10 mx-auto text-center animate__animated animate__fadeIn section-content">
                        {children}
                    </div>
                </div>
            </div>
            </div>
        </section>
    );
}