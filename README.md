# Molecular simulation concepts in Python

Introduction to coding up some components of molecular simulations.

## Getting started

First, you need to get a copy of this repository. It has been configured as a
**template** so that you can make a copy for yourself. (In a collaborative
project, you might instead fork it so that you can contribute back using a pull
request.) To make a copy of the template, click the **Use this template** button
and add it to your GitHub account. Then, clone your repository to your computer:

```
git clone https://github.com/<myaccount>/warmup-py.git
cd warmup-py
```

All code should be run with Python 3. You can get a good working installation
for scientific computing using a distribution like Anaconda, or you can use the
one supplied by your system. There are different ways to manage packages, but
for generality, I will use Python's virtual environments with pip.

From the project directory, run:

```
python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
```

Don't forget to run `activate` the next time you want to use this environment.
You can configure a comparable environment with conda if you prefer!

## First step

With your environment set up, you should build the documentation that comes with
this project. From the project directory:

```
cd doc
make html
```

You can then open `build/html/index.html` with your browser.

## Moving forward

Read the introduction in the generated documentation, run the included unit tests
from the project directory, and then start to code up the missing components!
