import os
import math
import random

init = 0
stop = 10
for x in range (init,stop):
	var='gnome-terminal -x sh -c "'+' ./bayesOpt '+str(x)+'"'
	print(var)
	os.system(var)
