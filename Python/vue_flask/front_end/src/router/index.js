/*
 * @Author: xizhihui <zhihui_xi@qq.com>
 * @Date: 2019-06-16 10:18:09
 * @LastEditTime: 2019-06-16 14:09:37
 * @Description: Description here
 */
import Vue from 'vue'
import Router from 'vue-router'
import HelloWorld from "@/components/HelloWorld"
import Login from "@/components/Login"

// router
const router = new Router({
    routes: [
        {
            path: '/index',
            name: 'HelloWorld',
            component: HelloWorld
        },
        {
            path: "/login",
            name: "Login",
            component: Login
        }
    ]
})

// 如果检测不到 token, 势必要进入 login 页面
// 如果能检测到 token, 那么就不应该到 login 页面, 而是到主界面，这里是 index
router.beforeEach((to, from, next) => {
    const token = sessionStorage.getItem("token")
    if (token) {
        if (to.name == "Login") next("/index")
        else next()
    } else {
        if (to.name == "Login") next()
        else {
            next("/login")
        }
    }
})

export default router
