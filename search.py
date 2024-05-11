from re import *

def get_sim_num(fname:str):
    r1 = compile(r'sim\d+\.gro') # check the the file is a .gro-file of a simulation part
    if (r1.match(fname)):        # otherwise return -1
        r2 = compile('[^0-9]')    
        num = r2.sub('', fname)  # remove all letters from .gro-file and return a num of a simulation part
        return int(num)
    else: return -1
