#!/usr/bin/env python
# -*- coding:utf-8 -*-
###
# Author: xizhihui <zhihui_xi@qq.com>
# Date: 2019-06-16 12:57:11
# LastEditTime: 2019-06-16 14:07:54
# Description: 
###


from flask import Flask
from flask_jwt_extended import JWTManager
from app.api import api_bp
from app.extension import cors
from app.config import Config


def configure_extension(app):
    cors.init_app(app)
    JWTManager(app)
    app.register_blueprint(api_bp, url_prefix="/api")


def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)
    configure_extension(app)
    return app
