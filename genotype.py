def cal_CIPOS(std, num):
	pos = int(1.96 * std / num ** 0.5)
	return "-%d,%d"%(pos,pos)
