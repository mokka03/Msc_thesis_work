import torch
import time

usecolortext = True


def tic():
    global tic_t0, tic_t1
    tic_t0 = tic_t1 = time.perf_counter()


def toc():
    if usecolortext:
        green = '\x1b[1;32m'
        blue =  '\x1b[1;34m'
        end =   '\x1b[0m'
    else:
        green=blue=end=''
    t = time.perf_counter()
    global tic_t1
    if (t - tic_t1) > 0.01:
        print('Elapsed time %s%.2f s%s,' % (green, (t - tic_t0), end),
              '(dt = %s%.2f s%s)' % (blue,(t - tic_t1), end))
    else:
        print('Elapsed time %s%.2f s%s,' % (green, (t - tic_t0), end),
              '(dt = %s%.3f ms%s)' % (blue,((t - tic_t1)*1000), end))
    tic_t1 = t


def stat_cuda(msg=''):
    print('GPU memory usage ' + msg + ':')
    print('allocated: %dM (max %dM), cached: %dM (max %dM)'
          % (torch.cuda.memory_allocated() / 1024 / 1024,
             torch.cuda.max_memory_allocated() / 1024 / 1024,
             torch.cuda.memory_reserved() / 1024 / 1024,
             torch.cuda.max_memory_reserved() / 1024 / 1024))
    
def curl(M_datax,M_datay,M_dataz,dx,dy,dz):

    dFx_dy,dFx_dx,dFx_dz = torch.gradient(M_datax, spacing = dx)
    dFy_dy,dFy_dx,dFy_dz = torch.gradient(M_datay, spacing = dy)
    dFz_dy,dFz_dx,dFz_dz = torch.gradient(M_dataz, spacing = dz)

    curly = -(dFz_dy - dFy_dz)/2
    curlx = -(dFx_dz - dFz_dx)/2
    curlz = dFy_dx - dFx_dy

    return curlx, curly, curlz
