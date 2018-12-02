import curses
import random
import time
import copy

class Snake:
	def __init__(self, board=None):
		self.board = board or [40, 40]	# [width, height]
		self.snake = [
			[int(self.board[0] / 2), int(self.board[1] / 2)],
			[int(self.board[0] / 2)-1, int(self.board[1] / 2)]
		]
		self.directions = {
			"up": curses.KEY_UP,
			"down": curses.KEY_DOWN,
			"left": curses.KEY_LEFT,
			"right": curses.KEY_RIGHT
		}
		
	def begin(self, interval=100, useAI=False):
		self.useAI = useAI
		self.window = self.create(interval=interval)
		self.render()

	def create(self, interval, arrays=None):
		stdscr = curses.initscr()
		# 去掉光标
		curses.curs_set(0)
		""" create the space for snake moving """
		curses.start_color()
		#文字和背景色设置，设置了两个color pair，分别为1和2
		curses.init_pair(1, curses.COLOR_GREEN, curses.COLOR_BLACK)
		curses.init_pair(2, curses.COLOR_WHITE, curses.COLOR_BLACK)
		curses.init_pair(3, curses.COLOR_RED, curses.COLOR_BLACK)
		# 划定窗口区域
		height, width = stdscr.getmaxyx() # 获取最大窗口
		if height > self.board[1]: # 包含外围的栅栏
			height = self.board[1]
		if width > self.board[0]:
			width = self.board[0]
		window = curses.newwin(height, width, 0, 0)
		window.border(0)
		# 允许响应键盘输入
		window.keypad(True)
		# 如果 interval ms 没有用户输入，就直接进行下一步
		window.timeout(interval)

		# 用于 debug ，根据 arrays 画出当时的状态，arrays 最后一个元素是食物
		if arrays:	
			for arr in arrays:
				x, y = arr
				window.addch(y, x, "o", curses.color_pair(1))	# 蛇身，绿色
			window.addch(y, x, "$", curses.color_pair(2)) # 食物，白色
			window.addch(arrays[0][1], arrays[0][0], "o", curses.color_pair(3)) # 蛇头，红色
			window.addch(arrays[-2][1], arrays[-2][0], "#", curses.color_pair(3)) # 蛇尾，红色
			while True:
				next_key = window.getch()
				if next_key == ord("q"):
					curses.endwin()
					quit()
		return window

	def render(self):
		""" refresh the snake status """
		food = self.food()
		self.window.addch(food[1], food[0], "$", curses.color_pair(2))
		tail = None
		direct = curses.KEY_RIGHT
		while True:
			# 添加得分
			self.window.addstr(0, 2, str(len(self.snake)), curses.color_pair(2))
			new_direct = self.window.getch()

			if new_direct == ord("q"):
				self.clearall()

			# 启用 AI
			if self.useAI:
				direct = self.navigate(food)
			else:
				if not self.invalid_direct(direct, new_direct):
					direct = new_direct
			# 设定历史 direct，在 wander 函数走 S 形用到
			self.prev_direct = direct
			tail = self.move(direct)

			# 如果吃到食物
			if self.snake[0] == food:
				# 延长蛇
				self.snake.append(tail)
				self.window.addch(self.snake[-1][1], self.snake[-1][0], "#", curses.color_pair(3))
				# 添加新的食物
				food = self.food()
				self.window.addch(food[1], food[0], "$", curses.color_pair(2))
			
			# 检测停止边界：例如达到边界，吃到自己，那就会导致退出
			if self.terminate():
				self.clearall()

	def move(self, direct):
		""" 移动, 这里以添加头部去掉尾部实现 """
		# "up down left right"
		if direct in self.directions:
			direct = self.directions[direct]

		if direct == curses.KEY_LEFT:
			head = [self.snake[0][0]-1, self.snake[0][1]]
		elif direct == curses.KEY_RIGHT:
			head = [self.snake[0][0]+1, self.snake[0][1]]
		elif direct == curses.KEY_UP:
			head = [self.snake[0][0], self.snake[0][1]-1]
		else:
			head = [self.snake[0][0], self.snake[0][1]+1]
		self.window.addch(self.snake[0][1], self.snake[0][0], "o", curses.color_pair(2))	# 旧的头改成白色
		# 走到边界了，打印出 head 坐标
		print("line 111: ", head)
		self.window.addch(head[1], head[0], "o", curses.color_pair(3))			# 新的头是红色
		# 添加头部，去掉尾部, 在添加头部之后进行
		self.snake.insert(0, head)
		tail = self.snake.pop()
		self.window.addch(tail[1], tail[0], ' ')
		# 新的尾巴改成 “#”
		self.window.addch(self.snake[-1][1], self.snake[-1][0], "#", curses.color_pair(1))
		return tail

	def navigate(self, food):
		"""
		蛇的状态
		1、能吃食物，吃完食物后还能找到尾巴，这时候就直接去吃食物。
		   能吃食物，并且吃完之后能找到尾巴。这种情况只在蛇身较短的时候出现，而且这里吃完之后能找到尾巴，指的是沿最短路径去吃。
		2. 蛇不能吃食物，或者吃食物后找不到尾巴（找不到尾巴还去吃是不允许的，这回导致吃完后无路可走），这时候就跟着蛇尾巴走，走到能去吃食物为止
		   追着尾巴走，也是又讲究的，首先，要保证走完这一步后还能继续找到尾巴，这是基本前提，
		   其次，如果两个方向都能找到尾巴，应该选择离食物较远的方向（注意不是离食物较近或者离尾巴较近），只有这样，才能保证留出了足够空间给蛇去吃食物
		3. wander 一下
			在BFS无解后， 告诉蛇一个步数step(随机产生step)，让它在空白区域以S形运动step步。

		BFS -> food 如果有路径，虚拟蛇跑过去吃食物后，检查是否安全。(安全是指头和尾巴有路径)
			如果安全。走一步
			如果不安全，检测当前头部与尾巴是否有路径，有的话，朝尾巴走一步
		BFS -> food 如果没有路径，检测头到尾巴是否有路径(是否安全)
			如果有，朝尾巴走一步
		如果没有路径
			挑一步可以走的走
		核心原则：
		1. 目标是食物时，走最短路径 
		2. 目标是蛇尾时，走最长路径(即从蛇头邻接的三个点最长路径到蛇尾的那个点)
		3. 与食物和蛇尾都没路径存在的情况下，挑一步可行的步子来走，最短最长关系都不大了
		下一步如此循环

		策略1.如果吃完苹果还可以找到到自己尾巴路线的话，才去吃苹果。
		策略2.如果找不到吃苹果路线或者吃完苹果后，会发生找不到自己尾巴路线的话，
			那么就在头部的周围找一个格子，这个格子要满足两个条件，
			条件1，走完这个格子要能找到到尾巴的路线，条件2，这个格子到苹果的距离是最远的。
		策略3.如果前两个都失败了，在附近以 S 形走
		"""
		# 蛇到食物有路径，吃了食物后还能找到尾巴
		direct = None
		if self.can_eat_food(food):
			path = BFS(self.board, self.snake, food)
			direct = path[0]
		# 蛇能找到尾巴，走一步也能找到尾巴
		else:
			direct = self.follow_tail()

		if direct is None:
			direct = self.wander()

		return direct

	def can_eat_food(self, food):
		""" 蛇到食物有路径，吃了食物后还能找到尾巴 """
		path = BFS(self.board, self.snake, food)
		# 蛇到食物的路都没有，肯定不是吃东西了
		if len(path) == 0:
			return False

		dumpy_snake = copy.deepcopy(self.snake)
		length = len(path)
		# 生成吃了食物的虚拟蛇
		for i in range(length):
			direct = path[i]
			x, y = dumpy_snake[0]
			if direct == "left":
				dumpy_snake.insert(0, [x-1, y])
			elif direct == "right":
				dumpy_snake.insert(0, [x+1, y])
			elif direct == "up":
				dumpy_snake.insert(0, [x, y-1])
			else:
				dumpy_snake.insert(0, [x, y+1])
			# 如果还在途中，则是添加了新头要把原来的尾巴去掉；吃了食物就不用去掉尾巴（延长了）
			if i != length - 1:
				dumpy_snake.pop()
		# 此时的虚拟蛇找尾巴
		head_to_tail = BFS(self.board, dumpy_snake, dumpy_snake[-1], tail=True)
		return len(head_to_tail) > 0

	def follow_tail(self):
		# 蛇找尾巴
		follow_path = BFS(self.board, self.snake, self.snake[-1], tail=True)
		# 找不到尾巴
		if len(follow_path) == 0:
			return None

		# 找得到尾巴
		# 然后在四周的点中，虚拟蛇走到这几个点，看能否找到尾巴
		# 多个点找的到，就取最长的
		points = self.neighbors()
		directs = self.neighbors(direct=True)
		print("snake: ", self.snake)
		print(points)
		if len(points) == 0:
			return None

		# 虚拟蛇走到这几个点的路径
		paths = []
		for i in range(len(points)):
			point = points[i]
			snake = [point]
			snake.extend(self.snake)
			snake.pop()
			path = BFS(self.board, snake, snake[-1], tail=True)
			# 添加真实蛇到虚拟蛇要走的路径
			path.insert(0, directs[i])
			paths.append(path)
		longest_path = paths[0]
		for path in paths:
			if len(path) > len(longest_path):
				longest_path = path

		return longest_path[0]

	def wander(self):
		points = self.neighbors()
		directs = self.neighbors(direct=True)
		# 不是"只能走上一个方向"
		if len(directs) > 1 and (self.prev_direct in directs):
			idx = directs.index(self.prev_direct)
			points.pop(idx)
			directs.pop(idx)
		if len(directs) > 0:
			return directs[random.randrange(len(directs))]
		# 没有其他方向，只能朝同一个方向走
		else:
			return self.prev_direct

	def neighbors(self, direct=False):
		""" 获得下一步可以走的点 """
		x, y = self.snake[0]
		points = []
		directions = []
		direct_condidate = ["right", "left", "down", "up"]
		condidate = [
			[x+1, y], [x-1, y], [x, y+1], [x, y-1]
		]
		for i in range(len(condidate)):
			con = condidate[i]
			if ((0 < con[0] < self.board[0]-1) and
			   (0 < con[1] < self.board[1]-1) and
			   (con not in self.snake)):
			   points.append(con)
			   directions.append(direct_condidate[i])
		if direct:
			return directions
		return points

	def food(self, count=None):
		""" create food randomly """
		food = [0, 0]
		while True:
			food[0] = random.randrange(1, self.board[0]-1)
			food[1] = random.randrange(1, self.board[1]-1)
			# 不能出现在 snake 自己身上
			if food not in self.snake:
				break
		return food

	def invalid_direct(self, prev_key, next_key):
		# 不是方向键
		if next_key not in [curses.KEY_RIGHT, curses.KEY_LEFT, curses.KEY_UP, curses.KEY_DOWN]:
			return True
		# 虽然是，但是是反向的
		elif prev_key == curses.KEY_RIGHT:
			return next_key == curses.KEY_LEFT
		elif prev_key == curses.KEY_LEFT:
			return next_key == curses.KEY_RIGHT
		elif prev_key == curses.KEY_UP:
			return next_key == curses.KEY_DOWN
		else:
			return next_key == curses.KEY_UP

	def terminate(self):
		# 吃到自己了
		if self.snake in self.snake[1:]:
			return True
		# 吃到墙了
		head = self.snake[0]
		if head[0] == 0 or head[0] == self.board[0] or head[1] == 0 or head[1] == self.board[1]:
			return True

		return False
	
	def clearall(self):
		print("clearall")
		curses.endwin()
		quit()


def BFS(board, snake, food, tail=False):
	""" 寻找 snake 到 food 的路，tail 为真，表明是寻找到尾巴的路 """
	visited = [[0 for i in range(board[0])] for j in range(board[1])]
	directions = [[0 for i in range(board[0])] for j in range(board[1])]
	# 设置已经访问的
	for s in snake:
		try:
			visited[s[1]][s[0]] = 1
			directions[s[1]][s[0]] = 'right'
		except Exception as e:
			print(s)
			raise(e)
	if tail:
		s = snake[-1]
		visited[s[1]][s[0]] = 0
		directions[s[1]][s[0]] = 0

	start = copy.deepcopy(snake[0])
	stack = [start]
	target = None

	# 寻路
	while stack:
		current = stack.pop(0)
		# 标记已经访问
		visited[current[1]][current[0]] = 1

		# 是否符合条件
		if current[0] == food[0] and current[1] == food[1]:
			target = current
			break

		# 寻找相邻的, 开区间是因为在 board 四周有栅栏, 也是不能访问的
		x = current[0]
		y = current[1]
		if x+1 < board[0]-1 and not visited[y][x+1]:
			visited[y][x+1] = 1
			directions[y][x+1] = "right"
			stack.append([x+1, y])
		if 0 < x-1 and not visited[y][x-1]:
			visited[y][x-1] = 1
			directions[y][x-1] = "left"
			stack.append([x-1, y])
		if y+1 < board[1]-1 and not visited[y+1][x]:
			visited[y+1][x] = 1
			directions[y+1][x] = "down"
			stack.append([x, y+1])
		if 0 < y-1 and not visited[y-1][x]:
			visited[y-1][x] = 1
			directions[y-1][x] = "up"
			stack.append([x, y-1])
	

	# 回溯
	path = []
	while target:
		x, y = target
		direct = directions[y][x]
		path.append(direct)

		if target == start:
			path.pop()
			break

		if direct == "left":
			target = [x+1, y]
		elif direct == "right":
			target = [x-1, y]
		elif direct == "up":
			target = [x, y+1]
		else:
			target = [x, y-1]
	path.reverse()

	return path

def main():
	mysnake = Snake([30, 15])
	mysnake.begin(useAI=True, interval=10)
	# board = [20, 10]
	# snake = [[19, 8], [18, 8], [17, 8], [16, 8], 
	# 		[15, 8], [14, 8], [13, 8], [12, 8],
	# 		[11, 8], [11, 7], [10, 7], [10, 6],
	# 		[10, 5], [10, 4], [9, 4], [8, 4], 
	# 		[7, 4], [6, 4], [5, 4], [4, 4], 
	# 		[3, 4], [2, 4], [2, 3], [3, 3], 
	# 		[4, 3], [5, 3], [5, 2], [4, 2], 
	# 		[3, 2], [2, 2], [2, 1], [3, 1], 
	# 		[4, 1], [5, 1], [6, 1], [7, 1], 
	# 		[8, 1], [9, 1], [10, 1], [11, 1], 
	# 		[12, 1], [13, 1], [14, 1], [15, 1], 
	# 		[16, 1], [16, 2], [15, 2], [14, 2], 
	# 		[13, 2], [12, 2], [11, 2], [10, 2], 
	# 		[9, 2], [8, 2], [7, 2], [6, 2], 
	# 		[6, 3], [7, 3], [8, 3], [9, 3], 
	# 		[10, 3], [11, 3], [12, 3], [13, 3],
	# 		[14, 3], [15, 3], [16, 3], [17, 3], 
	# 		[17, 2], [17, 1], [18, 1], [18, 2], 
	# 		[18, 3], [18, 4], [17, 4], [16, 4], 
	# 		[15, 4], [14, 4], [13, 4], [13, 5], 
	# 		[14, 5], [15, 5], [16, 5], [17, 5], 
	# 		[18, 5], [18, 6], [17, 6], [16, 6], 
	# 		[15, 6], [15, 7], [16, 7], [17, 7], [18, 7], [18, 8]]
	# food = [19,8]
	# path = BFS(board, snake, food)
	# snake.append(food)
	# mysnake.create(interval=100, arrays=snake)
	

if __name__ == '__main__':
	main()