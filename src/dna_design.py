# Author: Nguyen Truong Duy
# Contact: truongduy134@gmail.com

import algo_routine

MAP_CONSTRAINT_SET_TO_FUNC = {'1-4': algo_routine.gendnaword14, '1-2-3-4-5-6': algo_routine.gendnaword1to6, '1-2-3-4-5-6-7': algo_routine.gendnaword1to7, '1-2-3-7-8': algo_routine.gendnaword12378, '1-2-3-4-5-6-7-8': algo_routine.gendnaword1to8, '1-2-3-4-5-6-9': algo_routine.gendnaword1to6and9}

def gen_dna_words(n, maptypetoparam):
    key = get_coded_constraint_set(maptypetoparam)

    if not (key in MAP_CONSTRAINT_SET_TO_FUNC):
        raise RuntimeError("The program does not support this constraint set")
    
    func = MAP_CONSTRAINT_SET_TO_FUNC[key]
    return func(n, maptypetoparam)
    
def get_coded_constraint_set(maptypetoparam):
    keyarr = [str(key) for key in maptypetoparam]
    keyarr.sort()
    return '-'.join(keyarr)

#def preprocess_input_map(maptypetoparam):
    # Remove invalid keys
    # Keys are valid if they are integers from 1 to 9
#    for key in maptypetoparam:
#        if (not isinstance(key, int)) or key < 0 or key > 9:
#            del maptypetoparam[key]

#    maxkey = max([key for key in maptypetoparam])
