#name: Python Anchors Count
#description: anchors count in html
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.7, pip=22.1.2, {pip: [beautifulsoup4,]}]
#input: string html
#output: int count

from bs4 import BeautifulSoup

soup = BeautifulSoup(html, 'html.parser')
count = len(soup.findAll('a'))
