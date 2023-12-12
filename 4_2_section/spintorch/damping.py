"""Damping with absorbing boundaries"""
from torch import nn, ones
from skimage.draw import rectangle_perimeter

class Damping(nn.Module):
    # alpha = 0.0016    ##
    alpha = 13e-4    # damping coefficient ()
    alpha_max = 0.5 # maximum damping used on boundaries and for relax ()
    # region_width = 10   # width of absorbing region (cells)
    
    def __init__(self, dim: tuple, region_width: int):
        super().__init__()
        self.dim = dim
        self.region_width = region_width
        
        A = self.alpha*ones((1, 1,) + self.dim) # damping coefficient pointwise ()
        
        for i in range(self.region_width):
            x, y = rectangle_perimeter((i+1, i+1), (self.dim[0]-i-2, self.dim[1]-i-2))
            A[:, :, x, y] = (1-i/self.region_width)**2*(self.alpha_max-self.alpha) + self.alpha

        A[:,:,:,0:100] = A[:,:,:,int(dim[1]/2)-50:int(dim[1]/2)+50]
        A[:,:,:,dim[1]-100:dim[1]] = A[:,:,:,int(dim[1]/2)-50:int(dim[1]/2)+50]
            
        self.register_buffer("Alpha", A)

    def forward(self, relax=False):
        if relax:
            return self.alpha_max
        else:
            return self.Alpha



