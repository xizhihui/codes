<!--
 * @Author: xizhihui <zhihui_xi@qq.com>
 * @Date: 2019-06-10 14:46:50
 * @LastEditTime: 2019-06-16 16:13:41
 * @Description: Description here
 -->
<template>
    <div class="header">
        <header>
            <img src="../assets/logo.png" alt="vue-flask"><h2>Vue-Flask</h2>
            <div class="user-infos">
                <img :src="avatar" :alt="username" :title="username">
                <mu-button @click="logout" size="mini">退出</mu-button>
            </div>
        </header>
    </div>
</template>

<script>
export default {
    name: "Header",
    props: ["username", "avatar"],
    data () {
        return {}
    },
    methods: {
        logout () {
            let that = this
            //去掉 sessionStorage
            sessionStorage.removeItem("token")
            sessionStorage.removeItem("redirect")
            this.$router.push("/login")
        }
    }
}
</script>

<style scoped>
.header {
    width: 100%;
    height: 80px;
    padding-top: 10px;
    box-sizing: border-box;
    background-color: #3f51b5;
}
header {
    max-width: 1104px;
    margin: auto;
    line-height: 60px;
    height: 60px;
    text-align: left;
    overflow: hidden;
    color: white;
}
header > img {
    width: 60px;
    height: 60px;
    float: left;
}
header h2 {
    height: 60px;
    margin-left: 30px;
    float: left;
    color: white;
}
.user-infos {
    float: right;
}
.user-infos img {
    width: 30px;
    height: 30px;
    vertical-align: middle;
}
.user-infos button {
    height: 30px;
    width: 60px;
    margin-left: 15px;
}
</style>
