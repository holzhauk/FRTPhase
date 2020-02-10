from multiprocessing import Pool
from time import sleep
import random

def sum(task, a, b):
    sleepTime = random.randint(1, 4)
    print(task, " requires ", sleepTime, " seconds to finish")
    sleep(sleepTime)
    return a+b

def printResult(result):
    print(result)

if __name__=="__main__":
    myPool = Pool(5)

    for i in range(10):

        result =  myPool.apply_async(sum, args=("task"+str(i), 10, 20,), callback = printResult)

    #result1 = myPool.apply_async(sum, args=("task1", 10, 20,), callback = printResult)
    #result2 = myPool.apply_async(sum, args=("task2", 20, 30,), callback = printResult)
    #result3 = myPool.apply_async(sum, args=("task3", 30, 40,), callback = printResult)
    #result4 = myPool.apply_async(sum, args=("task4", 40, 50,), callback = printResult)
    #result5 = myPool.apply_async(sum, args=("task5", 50, 60,), callback = printResult)

    print("Submitted tasks to pool")

    myPool.close()
    myPool.join()
