#!/usr/bin/env python
# -*- coding:utf-8 -*-
###
# Author: xizhihui <zhihui_xi@qq.com>
# Date: 2019-06-16 12:56:01
# LastEditTime: 2019-06-16 16:10:20
# Description: entrance for the flask
###


from app import create_app
from flask import render_template

app = create_app()

@app.route("/")
def hello():
    return render_template("index.html")
