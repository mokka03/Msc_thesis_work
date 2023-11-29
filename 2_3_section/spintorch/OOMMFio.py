import numpy as np
import torch
from torch import nn
from matplotlib import pyplot as plt 


def OOMMF2torch(filename):

    f = open(filename, "r")

    line = ''
    while ('# Begin: Data Text' in line) == False:
        line = f.readline()
        if ('# xnodes:' in line) == True:
            words = line.split()
            nx = int(words[-1])
            
        if ('# ynodes:' in line) == True:
            words = line.split()
            ny = int(words[-1])

    datax = np.zeros((ny,nx))
    datay = np.zeros((ny,nx))
    dataz = np.zeros((ny,nx))

    data_np = np.loadtxt(filename)

    data_tensor = torch.zeros((3,ny,nx))
    data_tensor[0,] = torch.tensor(data_np[:,0].reshape((ny,nx)), dtype=torch.float32)
    data_tensor[1,] = torch.tensor(data_np[:,1].reshape((ny,nx)), dtype=torch.float32)
    data_tensor[2,] = torch.tensor(data_np[:,2].reshape((ny,nx)), dtype=torch.float32)
    
    return data_tensor