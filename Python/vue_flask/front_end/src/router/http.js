/*
 * @Author: xizhihui <zhihui_xi@qq.com>
 * @Date: 2019-06-16 11:05:29
 * @LastEditTime: 2019-06-16 14:09:05
 * @Description: axios, 请求拦截
**/

import Vue from 'vue'
import axios from "axios"
import router from "./index"
import config from "../config"

// axios
axios.defaults.timeout = config.axios_timeout
axios.defaults.baseURL = config.axios_baseurl
// 请求拦截器, 设置头部
axios.interceptors.request
    .use( config => {
        const token = sessionStorage.getItem("token")
        config.data = JSON.stringify(config.data)
        config.headers = {
            "Content-Type": "application/json"
            // 'Content-Type': 'application/x-www-form-urlencoded'
        }
        if (token) {
            config.headers["Authorization"] = "Bearer " + token
        }
        return config
    }, error => {
        return Promise.reject(error)
    })
// 响应拦截器
axios.interceptors.response
    .use( response => {
        // additional codes that do with reponse
        return response
    }, error => {
        let res = error
        console.log(error)
        if (error.response) {
            res = error.response
            switch (res.status) {
                case 401:
                    router.replace("/login")
                    break
                case 403:
                    // do something
                    break
                case 405:
                    // do something
                    break
                case 500:
                    // router.push("/login")
                    break
            }
        }
        return Promise.reject(res)
})

export default axios
