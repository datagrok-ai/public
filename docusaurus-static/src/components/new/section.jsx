import Markdown from 'markdown-to-jsx';
import React from 'react';

export default function Section({children, color}){
    return(
        <section className='container-fluid py-5 px-0' style={{backgroundColor: color != undefined ? color : 'white'}}>
                <div className="container p-5 animate__animated animate__fadeIn" >
                    <div className="row">
                        <div className="col-lg-12 col-sm-10 mx-auto text-center section-content">
                            {children}
                        </div>
                    </div>
                </div>
        </section>
    )
}