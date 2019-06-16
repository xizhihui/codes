#!/usr/bin/env python
# -*- coding:utf-8 -*-
###
# Author: xizhihui <zhihui_xi@qq.com>
# Date: 2019-06-16 13:16:14
# LastEditTime: 2019-06-16 14:07:18
# Description: environmental config
###


import os


class Config(object):
    DATABASE = os.path.realpath("app/db.json")
    SECRET_KEY = os.getenv("SECRET_KEY") or "You-can-guess-what"