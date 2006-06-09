subsets <- function(n, k, v = 1:n, allow.repeat = F)
{
	#Reference:  Statistics and Computing, Venebles and Ripley, p49)
	#select all subsets of size k from the integers 1:n
	#either allow repeats or don't
	v.all <- if(allow.repeat) v else v[-1]
	n.all <- if(allow.repeat) n else n - 1
	if(k <= 0 | n <= 0)
		NULL
	else if(k == 1)
		matrix(v[1:n], ncol = 1)
	else if((k >= n) & !allow.repeat)
		v[1:n]
	else if((k >= n) & allow.repeat)
		rbind(cbind(v[1], Recall(n, k - 1, v, allow.repeat)), Recall(
			n - 1, k, v[-1], allow.repeat))
	else rbind(cbind(v[1], Recall(n.all, k - 1, v.all, allow.repeat)),
			Recall(n - 1, k, v[-1], allow.repeat))
}
