#!/usr/bin/env python
# -*- coding:utf-8 -*-
###
# Author: xizhihui <zhihui_xi@qq.com>
# Date: 2019-06-16 13:12:07
# LastEditTime: 2019-06-16 15:41:20
# Description: model
###


from werkzeug.security import generate_password_hash, check_password_hash
from flask import jsonify
from app.util import avatar
from app.config import Config
import json


class User(object):
    """新建一个用户
    Args:
        username
        email
        password
    Returns:
        an User instance
    """
    def __init__(self, username, email, password):
        self.username = username
        self.email = email
        self.password = generate_password_hash(password)
        self.avatar = avatar(email, 128)
    
    def insert(self):
        with open(Config.DATABASE, "r") as f:
            db = json.load(f)
        if "users" not in db.keys():
            db["users"] = []
        if any(map(lambda user: self.email == user["email"], db["users"])):
            return False
        else:
            db["users"].append(self.to_dict())
            with open(Config.DATABASE, "w") as f:
                json.dump(db, f, indent=4)
            return True

    def to_dict(self):
        user = dict(
            email=self.email,
            password=self.password,
            username=self.username,
            avatar=self.avatar
        )
        return user

    @staticmethod
    def get(email):
        with open(Config.DATABASE, "r") as f:
            db = json.load(f)
            if not db or "users" not in db.keys():
                return None
            users = db["users"]
            user = list(filter(lambda u: u["email"] == email, users))
            if user:
                user = user[0]
                return dict(
                    email=user["email"],
                    password=user["password"],
                    username=user["username"],
                    avatar=user["avatar"]
                )
            else:
                return None

    @staticmethod
    def validate_password(pwdhash, pwd):
        return check_password_hash(pwdhash, pwd)
