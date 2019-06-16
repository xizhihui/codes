/*
 * @Author: xizhihui <zhihui_xi@qq.com>
 * @Date: 2019-06-16 10:18:09
 * @LastEditTime: 2019-06-16 14:50:42
 * @Description:
 */
// The Vue build version to load with the `import` command
// (runtime-only or standalone) has been set in webpack.base.conf with an alias.
import Vue from 'vue'
import App from './App'
import Router from "vue-router"
import router from './router'
import axios from "./router/http"
import VueAxios from "vue-axios"
import MuseUI from 'muse-ui'
import 'muse-ui/dist/muse-ui.css'
import Toasted from "vue-toasted"


Vue.config.productionTip = false


Vue.use(MuseUI)
Vue.use(Toasted, {
    theme: "outline",
    position: "top-center",
    duration: 2000
})


Vue.use(VueAxios, axios)
Vue.use(Router)


/* eslint-disable no-new */
new Vue({
    el: '#app',
    router,
    components: { App },
    template: '<App/>'
})
