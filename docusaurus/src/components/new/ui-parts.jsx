import Markdown from 'markdown-to-jsx';
import React, { useEffect, useState, Component } from 'react';

export function Button({label, type, size}){
    let btnClass = "btn brand-btn m-2 px-4 ";
    if (type != undefined){
        if (type == "outline") btnClass+="brand-btn-outline "
        if (type == "primary") btnClass+="brand-btn-primary "
        if (type == "accent") btnClass+="brand-btn-accent "
        if (type == "") btnClass+="btn-outline "
    } else {
        btnClass+="btn-outline"
    }
    size != undefined ? btnClass+= " btn-"+size : ""
    return(<a className={btnClass} type="submit" href="">{label}</a>);
  }
  

export function Image({src, frame}){
    if (frame)
        return (
            <div className="browser shadow-sm shadow my-5 rounded animate__animated animate__fadeInUp">
                <div className="browser-header">
                    <span className="dot"></span>
                    <span className="dot"></span>
                    <span className="dot"></span>
                </div>
                <img src={src} style={{borderTop: '1px solid var(--ifm-color-gray-200)'}}></img>
            </div>
        )
    else 
        return (
            <img src={src}></img>
        )

}

export function Tab({name, title, description, link, slides, indicators}){
    const list = slides.length ? Array.from(slides) : [slides];   
    const [activeIndex, setActiveIndex] = useState(0);
    return (
        <div className='row'>
            <div className="col-lg-6 col-xs-10 text-left">
                <div style={{ background: '#F9FAFB', borderTopLeftRadius: '0.5rem', borderBottomLeftRadius: '0.5rem', borderBottomRightRadius: '0.5rem', padding: '1.5rem', marginBottom: '1rem', boxShadow: '0 1px 2px rgba(0,0,0,0.05)' }}>
                    <div className='h4' style={{fontWeight: 600}}>{title}</div>
                    <p style={{color: "var(--ifm-color-gray-700)"}}>{description}</p>
                    <ul className='nav flex-column accesscarousel py-3'>
                        {list.map((item, index) => (
                            <li className='nav-item' key={index}>
                                <a
                                    className={index === activeIndex ? 'nav-link active' : 'nav-link'}
                                    data-target={`#accesscarousel-${getId(String(name)).toLocaleLowerCase()}`}
                                    data-slide-to={index}
                                    onClick={() => setActiveIndex(index)}
                                >
                                    {Array.isArray(item.text) ? (
                                        <ul style={{
                                            paddingLeft: "1.25rem",
                                            margin: 0,
                                            listStyleType: "disc",
                                            color: "#374151",
                                            fontSize: "1rem",
                                            lineHeight: "1.6"
                                        }}>
                                            {item.text.map((line, i) => (
                                                <li key={i} style={{ marginBottom: "0.5rem" }}>{line}</li>
                                            ))}
                                        </ul>
                                    ) : (
                                        <span style={{ color: "#374151", fontSize: "1rem", lineHeight: "1.6" }}>{item.text}</span>
                                    )}
                                </a>
                            </li>
                        ))}
                    </ul>
                    {link !== undefined && (
                        <a href={link} className='tab-body-link'>Learn more </a>
                    )}
                </div>
            </div>
            <div className="col-lg-6 col-xs-10">
                <div
                    id={`accesscarousel-${getId(String(name)).toLocaleLowerCase()}`}
                    className="accesscarousel carousel slide carousel-fade"
                    data-ride="carousel"
                    data-interval="false"
                >
                    <ol className={indicators !== undefined ? "carousel-indicators" : "d-none"}>
                        {list.map((item, index) => (
                            <li
                                data-target={`#accesscarousel-${getId(String(name)).toLocaleLowerCase()}`}
                                data-slide-to={index}
                                className={index === 0 ? "active" : ""}
                            ></li>
                        ))}
                    </ol>
                    <div className="carousel-inner">
                        {list.map((item, index) => (
                            <div className={index === 0 ? "carousel-item active" : "carousel-item"} key={index}>
                                {item.details ?? <Image src={item.image} />}
                            </div>
                        ))}
                    </div>
                </div>
            </div>
        </div>
    );
}

export function getId(val) {
    return val.replaceAll(' ', '-')
}

export function TabGroup({children}){
    const tabs = children.length ? Array.from(children) : [children];
    const id = Math.random();
    return (
        <>
            <ul className="nav nav-tabs nav-pills my-3 mt-5 features-nav-tab" id={`features-nav-tab-${id}`} role="tablist">
                {tabs.map((item, index) => (
                    <li className='nav-item' key={index}>
                        <a
                            className={index === 0 ? 'brand-btn btn-outline btn nav-link active' : 'brand-btn btn-outline btn nav-link'}
                            id={getId(item.props.name)}
                            data-toggle="tab"
                            role="tab"
                            aria-controls={`${getId(item.props.name)}-content`}
                            href={`#${getId(item.props.name)}-content`}
                        >
                            {item.props.name}
                        </a>
                    </li>
                ))}
            </ul>
            <div className="tab-content text-left container py-4" id={`features-nav-tab-content-${id}`}>
                {tabs.map((item, index) => (
                    <div
                        className={index === 0 ? "tab-pane fade show active" : "tab-pane fade"}
                        id={`${getId(item.props.name)}-content`}
                        role="tabpanel"
                        key={index}
                        aria-labelledby={item.props.name}
                    >
                        <Tab
                            name={item.props.name}
                            title={item.props.title}
                            description={item.props.description}
                            slides={item.props.slides}
                            link={item.props.link}
                        />
                    </div>
                ))}
            </div>
        </>
    );
}

export function Link({label, href, arrow}){
    return(<a href={href} className={arrow != undefined ? 'slider-body-link' : 'body-link'}>{label} </a>)
}

export function Slider({name, children, link, slides, indicators}){
    const list = slides.length ? Array.from(slides) : [slides];   
    return (
            <div className='row'>
                <div className="col-lg-6 col-xs-10 text-left d-flex flex-column align-self-center">
                    <div className='d-block'>
                    {children}
                    {
                        link != undefined 
                            ? <a href={link.href} className='slider-body-link'>{link.label} </a> 
                            : null
                    }
                    </div>
                </div>
                <div className="col-lg-6 col-xs-10 d-flex flex-column align-self-center">
                    <div id={"accesscarousel-"+String(name).toLocaleLowerCase()} className="accesscarousel carousel p-3 m-3 slide carousel-fade" data-ride="carousel" data-interval="false">
                    <ol className={indicators != undefined ? "carousel-indicators": "d-none"}>{list.map((item, index)=>{return(<li data-target={"#accesscarousel-"+String(name).toLocaleLowerCase()} data-slide-to={index} className={index == 0 ? "active": ''}></li>)})}</ol>
                        <div className="carousel-inner">
                            { list.map((item, index) => {
                                return (
                                    <div className={index == 0 ? "carousel-item active" : "carousel-item"} key={index}>
                                        <Image src={item.image}/>
                                    </div>
                                )
                            })}
                        </div>
                    </div>
                
                </div> 
            </div> 
    )
}

export function Col({children, color, col}){
    return(
        <div className={col+ " px-2"} >
            <div className='d-flex text-left p-3 pb-4 h-100 flex-column justify-content-center'style={{backgroundColor: color != undefined ? color : 'white', borderRadius:"8px"}}>
            {children}
        </div>
        </div>
    )
}

export function Row({children, color}){
    return(
        <div className="container p-4 section" style={{backgroundColor: color != undefined ? color : null}}>
            <div className="row">
                {children}
            </div>
        </div>
    )
}


