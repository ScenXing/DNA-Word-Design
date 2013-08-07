# Author: Nguyen Truong Duy
# Contact: truongduy134@gmai.com

import algo_routine

mapConstraintSetToFunc = {'1-4': algo_routine.genDnaWord14, '1-2-3-4-5-6': algo_routine.genDnaWord1To6, '1-2-3-4-5-6-7': algo_routine.genDnaWord1To7, '1-2-3-7-8': algo_routine.genDnaWord12378, '1-2-3-4-5-6-7-8': algo_routine.genDnaWord1To8, '1-2-3-4-5-6-9': algo_routine.genDnaWord1To6And9}

def genDnaWords(n, mapTypeToParam):
	key = getCodedConstraintSet(mapTypeToParam)

	global mapConstraintSetToFunc
	if not (key in mapConstraintSetToFunc):
		raise RuntimeError("The program does not support this constraint set")
	
	func = mapConstraintSetToFunc[key]
	return func(n, mapTypeToParam)

def getCodedConstraintSet(mapTypeToParam):
	keyArr = [str(key) for key in mapTypeToParam]
	keyArr.sort()
	return '-'.join(keyArr)
