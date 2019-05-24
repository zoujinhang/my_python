import queue
import threading

class My_thread(threading.Thread):
	def __init__(self,function,qu,num):
		super().__init__()
		self.qu = qu
		self.num = num
		self.function = function
	def run(self):
		with self.num:
			targ = self.qu.get()
			self.function(targ)
			self.qu.task_done()


def my_thread(function,pa,num):
	num = threading.Semaphore(num)
	qu = queue.Queue()
	for i in pa:
		qu.put(i)
	n = len(pa)
	thread = []
	for i in range(n):
		t = My_thread(function,qu,num)
		t.start()
		thread.append(t)
	for i in thread:
		i.join()
	qu.join()



















