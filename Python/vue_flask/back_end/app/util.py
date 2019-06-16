#!/usr/bin/env python
# -*- coding:utf-8 -*-
###
# Author: xizhihui <zhihui_xi@qq.com>
# Date: 2019-06-16 13:13:52
# LastEditTime: 2019-06-16 13:14:25
# Description: Description here
###


from hashlib import md5


def avatar(email, size):
    """ 根据邮箱生成全球统一的头像 """
    digest = md5(email.lower().encode('utf-8')).hexdigest()
    return 'https://www.gravatar.com/avatar/{}?d=identicon&s={}'.format(
        digest,
        size
    )