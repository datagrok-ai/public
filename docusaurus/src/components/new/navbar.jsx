import Markdown from 'markdown-to-jsx';
import React, { useEffect, useState, Component } from 'react';
import Link from '@docusaurus/Link';

const navItems = [
    {
        name: "Solutions",
        href: "#",
        items: [
            {
                name: "Cheminformatics",
                href: "/cheminformatics",
                iconFA: 'hexagon',
            },
            {
                name: "Bioinformatics",
                href: "/macromolecules",
                iconFA: 'dna',
            },
            {
                name: "Data Science",
                href: "/data-science",
                iconFA: 'chart-bar'
            },
            {
                name: "Enterprise IT",
                href: "/enterprise",
                iconFA: "building"
            },
            {
                name: "Academia",
                href: "/solutions/academia",
                iconFA: "university"
            },

        ]
    },
    {
        name: "Company",
        href: "#",
        items: [
            {
                name: "Careers",
                href: "/company/careers",
                iconFA: "briefcase"
            }, 
            {
                name: "Team",
                href: "/company/team",
                iconFA: "users"
            },
            {
                name: "Customer Stories",
                href: "/solutions",
                iconFA: "comment-alt-smile"
            }]
    },
    {
        name: "Resources",
        href: "#",
        items: [
            {
                name: "Wiki",
                href: "/help",
                iconFA: "book"
            }, 
            {
                name: "JS-API",
                href: "/js-api",
                iconFA: "cog"
            }, 
            {
                name: "Community",
                href: "https://community.datagrok.ai",
                iconFA: "comments-alt"
            }, 
            {
                name: "Youtube",
                href: "https://www.youtube.com/channel/UCXPHEjOd4gyZ6m6Ji-iOBYg",
                iconFA: "play"
            }, 
            {
                name: "User Meetings",
                href: "https://us02web.zoom.us/meeting/register/up0vfuCrpjgpHNzi371YEJIQ4GkMpTm4disW#/registration",
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
                    <Link to="mailto:sales@datagrok.ai" target="_blank" rel="noopener noreferrer">
                      <button className="btn brand-btn btn-outline mr-2" type="submit">Book a demo</button>
                    </Link>
                    <Link to="https://public.datagrok.ai" target="_blank" rel="noopener noreferrer">
                      <button className="btn brand-btn brand-btn-primary" type="submit">Try now</button>
                    </Link>
                </div>
            </div>
        </nav>
        </>
    );
};