<!--
 * @Author: xizhihui <zhihui_xi@qq.com>
 * @Date: 2019-06-16 10:43:45
 * @LastEditTime: 2019-06-16 16:13:23
 * @Description: Login/Register page
 -->

<template>
<div class="wrapper">
<mu-container>
        <mu-tabs :value.sync="active" center color="teal">
            <mu-tab><mu-icon value="phone"></mu-icon>登陆</mu-tab>
            <mu-tab><mu-icon value="phone"></mu-icon>注册</mu-tab>
        </mu-tabs>
        <div class="login" v-if="active == 0">
            <mu-form ref="loginForm" :model="loginForm">
                <mu-form-item label="邮箱" prop="email" :rules="rules.email">
                    <mu-text-field v-model="loginForm.email" prop="email"></mu-text-field>
                </mu-form-item>
                <mu-form-item label="密码" prop="password" :rules="rules.password">
                    <mu-text-field type="password" v-model="loginForm.password" prop="password"></mu-text-field>
                </mu-form-item>
                <mu-button class="submitBtn" color="primary" @click="loginSubmit">登陆</mu-button>
            </mu-form>
        </div>
        <div class="register" v-if="active == 1">
            <mu-form ref="registerForm" :model="registerForm">
                <mu-form-item label="用户名" prop="username" :rules="rules.username">
                    <mu-text-field v-model="registerForm.username" prop="username"></mu-text-field>
                </mu-form-item>
                <mu-form-item label="邮箱" prop="email" :rules="rules.email">
                    <mu-text-field v-model="registerForm.email" prop="email"></mu-text-field>
                </mu-form-item>
                <mu-form-item label="密码" prop="password" :rules="rules.password">
                    <mu-text-field type="password" v-model="registerForm.password" prop="password"></mu-text-field>
                </mu-form-item>
                <mu-form-item label="重复密码" prop="confirm" :rules="rules.confirm">
                    <mu-text-field type="password" v-model="registerForm.confirm" prop="confirm"></mu-text-field>
                </mu-form-item>
                <mu-button class="submitBtn" color="primary" @click="registerSubmit">注册</mu-button>
            </mu-form>
        </div>
    </mu-container>
</div>
</template>


<script>
import config from "../config"
export default {
    name: "Login",
    data () {
        return {
            active: 0,
            domains: ["qq.com", "163.com", "vip.163.com", "263.net", "yeah.net", "sohu.com",
                      "sina.cn", "sina.com", "eyou.com", "gmail.com", "hotmail.com"],
            loginForm: {
                email: "",
                password: "",
            },
            registerForm: {
                username: "",
                email: "",
                password: "",
                confirm: ""
            },
            rules: {
                username: [
                    { validate: (val) => !!val, message: "请填写用户名" },
                    { validate: (val) => val.length >= 3, message: "用户名长度必需大于或等于3" }
                ],
                password: [
                    { validate: (val) => !!val, message: "请填写密码" },
                    { validate: (val) => val.length >= 3 && val.length <= 12, message: "密码长度在3-12之间" }
                ],
                email: [
                    { validate: (val) => !!val, message: "请填写邮箱" },
                    { validate: this.isEmail, message: "请填写合法的邮箱" }
                ],
                confirm: [
                    { validate: (val) => !!val, message: "请重复密码" },
                    { validate: (val) => this.registerForm.password == val, message: "填写的两次密码不一致" }
                ]
            }
        }
    },
    methods: {
        isEmail (val) {
            var pattern = /^([A-Za-z0-9_\-\.])+\@([A-Za-z0-9_\-\.])+\.([A-Za-z]{2,4})$/;
            if (pattern.test(val)) {
                var domain = val.substring(val.indexOf("@")+1);
                for(var i = 0; i< this.domains.length; i++) {
                    if(domain == this.domains[i]) {
                        return true;
                    }
                }
            }
            return false;
        },
        loginSubmit () {
            let that = this
            that.$refs.loginForm.validate().then(pass => {
                if (!pass) {
                    that.$toasted.show("请保证各项填写正确")
                    return null
                }
                that.axios
                    .post("/api/login", that.loginForm)
                    .then(response => {
                        if (response.data.code == 0) {
                            sessionStorage.setItem("token", response.data.token)
                            that.$toasted.show("登陆成功", {
                                theme: "outline",
                                position: "top-center",
                                duration: 2000
                            })
                            that.$router.push("/index")
                        } else {
                            that.$toasted.show(config.status[response.data.code])
                        }
                    }).catch(error => {
                        that.$toasted.show(error)
                    })
            })
        },
        registerSubmit () {
            let that = this
            // 提交注册成功后，就切换到登陆按钮， active = 0
            this.$refs.registerForm.validate().then(pass => {
                if (!pass) {
                    that.$toasted.show("请保证各项填写正确")
                    return null
                }
                that.axios
                    .post("/api/register", that.registerForm)
                    .then(response => {
                        if (response.data.code == 0) {
                            that.$toasted.show("注册成功")
                            // 转到登陆，重置密码为空
                            this.active = 0
                            this.$refs.registerForm.clear()
                            this.registerForm.password = ""
                            this.registerForm.confirm = ""
                            this.loginForm.email = this.registerForm.email
                        } else {
                            that.$toasted.show(config.status[response.data.code])
                        }
                    }).catch(error => {
                        that.$toasted.show(error)
                    })

            })
        }
    }
}
</script>

<style scoped>
.wrapper {
    position: absolute;
    left: 0;
    top: 0;
    right: 0;
    bottom: 0;
    background-image: url("../assets/biotech.jpg");
    background-size: cover;
}
.container {
    width: 30%;
    background: snow;
    margin-top: 60px;
    padding: 0 0 30px;
}
.login, .register {
    width: 80%;
    margin: 15px auto;
}
.submitBtn {
    display: block;
    max-width: 120px;
    margin: 15px auto;
}
</style>
