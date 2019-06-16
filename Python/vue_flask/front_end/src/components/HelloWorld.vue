<!--
 * @Author: xizhihui <zhihui_xi@qq.com>
 * @Date: 2019-06-16 10:18:09
 * @LastEditTime: 2019-06-16 16:02:09
 * @Description: index
 -->

<template>
<div>
    <Header :username="userinfos.username" :avatar="userinfos.avatar"></Header>
    <div class="main">
        <ul>
            <li v-for="(value,key) in userinfos" :key="key">{{key}} : {{value}}</li>
        </ul>
    </div>
</div>
</template>

<script>
import Header from "./Header"
import config from "../config"

export default {
    name: "HelloWorld",
    components: {
        Header
    },
    data () {
        return {
            userinfos: {
                username: "",
                email: "",
                avatar: ""
            }
        }
    },
    created () {
        let that = this
        that.axios
            .get("/api/index")
            .then(response => {
                if (response.data.code == 0) {
                    for (let key of ["username", "email", "avatar"]) {
                        that.userinfos[key] = response.data[key]
                    }
                } else {
                    that.$toasted.show(config.status[response.data.code])
                }
            }).catch( error => {
                that.$toasted.show(error)
            })
    }
}
</script>

<style scoped>
.main {
    max-width: 1104px;
    width: 100%;
    margin: 30px auto;
}
</style>

