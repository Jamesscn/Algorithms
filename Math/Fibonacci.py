fibonacciList = {
	0: 0,
	1: 1,
	2: 1
}

def fibonacci(n):
	n = round(n)
	if n not in fibonacciList:
		if n % 2 == 0:
			k = n // 2
			fibonacciList[n] = fibonacci(k) * (2 * fibonacci(k + 1) - fibonacci(k))
		else:
			k = (n - 1) // 2
			fibonacciList[n] = fibonacci(k + 1) * fibonacci(k + 1) + fibonacci(k) * fibonacci(k)
	return fibonacciList[n]