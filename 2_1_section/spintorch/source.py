import skimage
import torch


class WaveSource(torch.nn.Module):
    def __init__(self, excitation_mask):
        super().__init__()

        self.register_buffer('excitation_mask', torch.tensor(excitation_mask))

    def forward(self, B, B_ext_0, Bt): ## B_ext_0
        B = B.clone()
        B = B_ext_0 + Bt*self.excitation_mask     ##
        return B

    def coordinates(self):
        return self.x.cpu().numpy(), self.y.cpu().numpy()


class WaveLineSource(WaveSource):
    def __init__(self, r0, c0, r1, c1, dim=0):
        x, y = skimage.draw.line(r0, c0, r1, c1)

        self.r0 = r0
        self.c0 = c0
        self.r1 = r1
        self.c1 = c1
        super().__init__(x, y, dim)

class CPWSource(WaveSource):
    def __init__(self,excitation_mask):
        super().__init__(excitation_mask)


