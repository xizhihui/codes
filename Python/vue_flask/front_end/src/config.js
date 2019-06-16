/*
 * @Author: xizhihui <zhihui_xi@qq.com>
 * @Date: 2019-06-16 10:37:39
 * @LastEditTime: 2019-06-16 14:21:13
 * @Description: Description here
**/


const config = {
    env: "development",
    axios_timeout: 5000,
    axios_baseurl: "http://localhost:5000",
    status: {
        "1": "some of email, password, username are empty",
        "2": "email is not registered",
        "3": "email has beed used when registering",
        "4": "password is wrong"
    }
}


export default config
