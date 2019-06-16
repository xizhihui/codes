#!/usr/bin/env python
# -*- coding:utf-8 -*-
###
# Author: xizhihui <zhihui_xi@qq.com>
# Date: 2019-06-16 12:57:28
# LastEditTime: 2019-06-16 15:50:31
# Description: Description here
###


"""
response code:
0: ok
1: some of email, password, username are empty
2: email is not registered
3: email has beed used when registering
4: password is wrong
"""

from flask import Blueprint, request, jsonify
from flask_jwt_extended import create_access_token, jwt_required, get_jwt_identity
from app.model import User


api_bp = Blueprint("api", __name__)  # url_prefix="/api"


@api_bp.route("/login", methods=["POST"])
def login():
    req_data = request.get_json()
    password = req_data.get("password", None)
    email = req_data.get("email", None)
    if not password or not email:
        return jsonify({
            "code": 1
        })
    user = User.get(email)
    if user:
        if User.validate_password(user["password"], password):
            return jsonify(dict(
                code=0,
                username=user["username"],
                email=user["email"],
                avatar=user["avatar"],
                token=create_access_token(identity=email)
            ))
        else:
            code = 4
    else:
        code = 2
    return jsonify({"code": code})


@api_bp.route("/register", methods=["POST"])
def register():
    req_data = request.get_json()
    username = req_data.get("username", None)
    email = req_data.get("email", None)
    password = req_data.get("password", None)
    if not username or not email or not password:
        code = 1
    else:
        user = User(username, email, password)
        if user.insert():
            code = 0
        else:
            code = 2
    return jsonify({"code": code})


@api_bp.route("/index", methods=["GET"])
@jwt_required
def index():
    email = get_jwt_identity()
    user = User.get(email)
    return jsonify(dict(
        code=0,
        email=email,
        username=user["username"],
        avatar=user["avatar"]
    ))
