# Simple Python Web App Demo Project #

This project consists of a web-form for entering a DNA sequence and selecting frames which the backend will translate and present in a HTML table.

The current version is running on Python Flask as the backend which uses Jinja2 as templating language. A Python-CGI version will follow.

## Installation ##

The project depends on the following Python packages, installable using Python *pip*:

* Flask
* Biopython
* werkzeug
* json

what is this? magic?
nah
just a readme

Install with `pip install flask biopython werkzeug json` or to force Python3 `pip3 install flask biopython` (preferably in a *virtualenv*)

## Running ##

Execute the **webapp.py** file: `python webapp.py` and visit http://127.0.0.1:5000/