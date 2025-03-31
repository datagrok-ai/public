import Markdown from 'markdown-to-jsx';
import React, { useEffect, useState, Component } from 'react';


const navItems = [
    {
        name: "Solutions",
        href: "#",
        items: [
            {
                name: "Cheminformatics",
                href: "#",
                iconFA: 'hexagon',
            }, 
            {
                name: "Bioinformatics",
                href: "#",
                iconFA: 'dna',
            },
            {
                name: "Data Science",
                href: "#",
                iconFA: 'chart-bar'
            }, 
            {
                name: "Enterprise IT",
                href: "#",
                iconFA: "building"
            }, 
            ]
    },
    {
        name: "Company",
        href: "#",
        items: [
            {
                name: "Careers",
                href: "#",
                iconFA: "briefcase"
            }, 
            {
                name: "Team",
                href: "#",
                iconFA: "users"
            },
            {
                name: "Customer Stories",
                href: "#",
                iconFA: "comment-alt-smile"
            }]
    },
    {
        name: "Resources",
        href: "#",
        items: [
            {
                name: "Wiki",
                href: "#",
                iconFA: "book"
            }, 
            {
                name: "JS-API",
                href: "#",
                iconFA: "cog"
            }, 
            {
                name: "Community",
                href: "#",
                iconFA: "comments-alt"
            }, 
            {
                name: "Youtube",
                href: "#",
                iconFA: "play"
            }, 
            {
                name: "User Meetings",
                href: "#",
                iconFA: "users-class"
            }]
    }
];

let result = [];

function RenderDropDownItem(data){
    return (
        <>
        <li key={data.keyindex} className={data.page == data.name ? "nav-item dropdown active": "nav-item dropdown"}>
            <a className="nav-link dropdown-toggle mx-2" role="button"  id={"nav-"+data.name.trim()} data-toggle="dropdown" href={data.href}> {data.name} </a>
            <div className="dropdown-menu shadow p-3 ml-2 animate__animated animate__fadeIn animate__faster" aria-labelledby={"nav-"+data.name.trim()}>
                {data.items.map((item, index) => { 
                    return (<a className={data.page == item.name ? "active dropdown-item p-2 pl-2 pr-3" : "dropdown-item p-2 pl-2 pr-3" } href={item.href} key={index}><i className={item.iconFA != undefined ? "far fa-fw mr-2 fa-"+item.iconFA: ''}></i> {item.name} </a>)
                })}
            </div> 
        </li>
        </>
    );
}

function RenderItem(data){
    return (
        <>
        <li key={data.keyindex} className={data.page == data.name ? "nav-item active" : "nav-item"}>
            <a className="nav-link mx-2" href={data.href}><i className={data.iconFA != undefined ? "far fa-fw mr-2 fa-"+data.iconFA: ''}></i> {data.name} </a>
        </li>
        </>
    );
}

function NavItems(page){
    result = [];
    navItems.forEach((navItem, index)=>{
        if (navItem.items != undefined && navItem.items.length != 0)
             result.push(<RenderDropDownItem key={index} page={page.current} href={navItem.href} name={navItem.name} items={navItem.items}/>)
        else
             result.push(<RenderItem key={index} iconFA={navItem.iconFA} page={page.current} href={navItem.href} name={navItem.name} />)
     });
     return (<>{result.map((navitems) => { return (navitems) })}</>)
}

export default function navbar(current){
    return (
        <>
        <nav className="navbar fixed-top navbar-expand-lg px-5 animate__animated animate__fadeInDown">
            <a className="navbar-brand mb-0 h1 animate__animated" href="https://datagrok.ai/" 
            onMouseOver={(e)=>{e.target.classList.add('animate__pulse')}}
            onMouseLeave={(e)=>{e.target.classList.remove('animate__pulse')}} >datagrok</a>
            <div className='ml-auto mb-xs-4 d-lg-none'>
                    <button className="btn brand-btn btn-outline mr-2" type="submit">Get A Demo</button>
                    <button className="btn brand-btn brand-btn-primary" type="submit">Try Now</button>
            </div>
            <button className="btn navbar-toggler collapsed ml-4" type="button" data-toggle="collapse" data-target="#navbar" aria-controls="navbar"
                    aria-expanded="false" aria-label="Toggle navigation">
                    <i className="far fa-bars fa-fw"></i>
                </button>
            <div className="collapse mx-auto w-100 navbar-collapse" id="navbar">
                <ul className="navbar-nav mx-auto mb-4 my-lg-0">
                   <NavItems current={current.page}/>
                </ul>
                <div className='mb-xs-4 d-none d-lg-block'>
                    <button className="btn brand-btn btn-outline mr-2" type="submit">Book a demo</button>
                    <button className="btn brand-btn brand-btn-primary" type="submit">Try now</button>
                </div>
            </div>
        </nav>
        </>
    );
};