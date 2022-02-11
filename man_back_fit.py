import numpy as np
import scipy.signal
import scipy
import os
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
import matplotlib.backend_bases
plt.rc('text', usetex=True)
plt.rcParams['text.usetex'] = True
import sys
lambda_max = 1500
width_g_max = 3
width_l_max = 3
obs = np.loadtxt("MICAlibs.txt")
x = obs[:,0]
sample_index = input("Input sample index (number of COLUMN):   ")
obs[:,1]= obs[:,int(sample_index)]
data1 = np.loadtxt("MgI.txt")
mask1 = (data1[:,0]>200) & (data1[:,0]<1000)
data = data1[mask1]
T = 1
a0 = 0.56534
a1 = 0.0245
a2 = 0.0157
Stark_w = .0410e-02
l = data[:,0]
aki = data[:,1]
ei = data[:,2]
ek = data[:,3]
gi = data[:,4]
gk = data[:,5]
eq1 = r"\begin{eqnarray*}" + \
      r"p\cross\frac{a_g}{\sqrt{2\pi\sigma^2}\exp{-\frac{(\lambda-\lambda_0)^2}{2\sigma^2}}}\\" + \
      r"(1-p)\cross\frac{a_l}{\pi}\frac{(\flatfrac{\Gamma}{2})^2}{(\lambda - \lambda_0)^2-(\flatfrac{|gamma}{2})^2} " + \
      r"\end{eqnarray*}"
I =((aki*gk)/(l))*np.exp(-(ek)/(T*100))
peak,_ = scipy.signal.find_peaks(I, prominence=10000)
larray = [390.5523]
from matplotlib.pyplot import figure, show
import numpy

class ZoomPan:
    def __init__(self):
        self.press = None
        self.cur_xlim = None
        self.cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None
    def zoom_factory(self, ax1, base_scale = 2.):
        def zoom(event):
            cur_xlim = ax1.get_xlim()
            cur_ylim = ax1.get_ylim()
            xdata = event.xdata 
            ydata = event.ydata 
            if event.button == 'down':
                scale_factor = 1 / base_scale
            elif event.button == 'up':
                scale_factor = base_scale
            else:
                scale_factor = 1
                print(event.button)
            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor
            relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata)/(cur_ylim[1] - cur_ylim[0])
            ax1.set_xlim([xdata - new_width * (1-relx), xdata + new_width * (relx)])
            ax1.set_ylim([ydata - new_height * (1-rely), ydata + new_height * (rely)])
            ax1.figure.canvas.draw()
        fig = ax1.get_figure() # get the figure of interest
        fig.canvas.mpl_connect('scroll_event', zoom)
        return zoom
    def pan_factory(self, ax1):
        def onPress(event):
            if event.inaxes != ax1: return
            self.cur_xlim = ax1.get_xlim()
            self.cur_ylim = ax1.get_ylim()
            self.press = self.x0, self.y0, event.xdata, event.ydata
            self.x0, self.y0, self.xpress, self.ypress = self.press
        def onRelease(event):
            self.press = None
            ax1.figure.canvas.draw()
        def onMotion(event):
            if self.press is None: return
            if event.inaxes != ax1: return
            dx = event.xdata - self.xpress
            dy = event.ydata - self.ypress
            self.cur_xlim -= dx
            self.cur_ylim -= dy
            ax1.set_xlim(self.cur_xlim)
            ax1.set_ylim(self.cur_ylim)
            ax1.figure.canvas.draw()
        fig1 = ax1.get_figure()
        fig1.canvas.mpl_connect('button_press_event',onPress)
        fig1.canvas.mpl_connect('button_release_event',onRelease)
        fig1.canvas.mpl_connect('motion_notify_event',onMotion)
        return onMotion
fig1 = figure()
#for lam in larray:
for lam in l[peak]:
    mask = (x>lam-1) & (x<lam+1)
    c_o = []
    b_l = []
    A_ki = aki[l == lam]
    area = 0
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(obs[mask][:,0],obs[mask][:,1], label = 'obs')
    ax1.axvline(lam, label = 'Standerd Peak Position')
    ax1.set_ylabel(r'$I \rightarrow$')
    ax1.set_xlabel(r'$\lambda \longrightarrow$')
    ax1.set_title("Peak at $\lambda = $"+str(lam)+"  $A_{ki} = $"+str(A_ki))
    ax1.grid()
    ax1.legend(loc = 'upper left')
    scale = 1.1
    zp = ZoomPan()
    figZoom = zp.zoom_factory(ax1, base_scale = scale)
    figPan = zp.pan_factory(ax1)
    plt.title("Peak at $\lambda = $"+str(lam)+"  $A_{ki} = $"+str(A_ki))
    fig1.show()
    coords_wl = []
    coords_int = []
    co = []
    def onclick(event):
        global ix, iy
        ix, iy = event.xdata, event.ydata
        global coords_wl
        coords_wl.append(ix)
        global coords_int
        coords_int.append(iy)
        return coords_wl, coords_int
    cid = fig1.canvas.mpl_connect('button_press_event', onclick)
    while True:
        sele = input("Remove background" )
        if sele == 'n':
            break
        else:
            x0 = coords_wl[0]
            x1 = coords_wl[1]
            y0 = coords_int[0]
            y1 = coords_int[1]
            x_values = obs[mask][:,0]
            def backgroung(x_values, x0,x1,y0,y1):
                BG = (((y1-y0)/(x1-x0))*(x_values - x0)) + y0
                return BG
            correct_out = obs[mask][:,0],obs[mask][:,1] - background(x_values, x0,x1,y0,y1)
            plt.plot(obs[mask][:,0],obs[mask][:,1], label = 'obs')
            plt.plot(correct_out[0],correct_out[1], label = 'corrected')
            #plt.plot(obs[mask][:,0],b_l, '--',  label = 'background')
            plt.axvline(lam, label = 'Standerd Peak Position')
            maskl = correct_out[0] <= lam
            maskr = correct_out[0] >= lam
            dipl,_ = scipy.signal.find_peaks(-correct_out[1][maskl])
            dipr,_ = scipy.signal.find_peaks(-correct_out[1][maskr])
            plt.legend(loc = 'upper left')
            lim_left = dipl[-1]
            lim_right = dipr[0]
            mask_int = (correct_out[0][maskl][lim_left] <= correct_out[0])&(correct_out[0] <= correct_out[0][maskr][lim_right])
            x_int = correct_out[0][mask_int]
            y_int = correct_out[1][mask_int]
            for x1 in range(np.size(correct_out[0][mask_int]) - 1):
                area1 = ((x_int[x1+1]-x_int[x1])*y_int[x1+1]) - ((1/2)*(x_int[x1+1]-x_int[x1])*(y_int[x1+1]-y_int[x1]))
                area += area1
            plt.fill_between(correct_out[0][mask_int],correct_out[1][mask_int], color = 'gray')
            plt.axvline(correct_out[0][maskl][lim_left], color = 'r')
            plt.axvline(correct_out[0][maskr][lim_right], color = 'r')
            plt.title("Peak at $\lambda = $"+str(lam)+"  $A_{ki} = $"+str(A_ki))
            plt.ylabel(r'$I \rightarrow$')
            plt.xlabel(r'$\lambda \longrightarrow$')
            plt.grid()
            plt.show()
            while True:
                sel = input("Remove peak at "+str(lam)+ ":           [Input y]           " )
                if sel == 'y':
                    print("This Line is removed")
                    break
                else:
                    with open('Out/area'+str(lam)+'_'+str(sample_index)+'.txt', 'w') as f1:
                        f1.write(str(area))
                    with open('Out/correct_out'+str(lam+1)+'_'+str(sample_index)+'.txt', 'w') as f2:
                        f2.write(str(correct_out))
                    param_bounds=([0,0,0,0,0,0.0],[np.inf,np.inf,lambda_max,3,3,1])
                    p = 0.5
                    popt_pv = []
                    import numpy as np
                    [amp_g, amp_l, cen, width_g, width_l, p] = [np.max(y_int),np.max(y_int), lam, x_int[-1]-x_int[0] , x_int[-1]-x_int[0], .5]
                    def pv(x_int, amp_g, amp_l, cen, width_g, width_l, p ):
                        return p*((amp_g/np.sqrt(2*np.pi*width_g**2))*np.exp(-((x_int-cen)**2)/(2*width_g**2)))+(1-p)*((amp_l/np.pi)*((width_l/2)**2/((x_int-cen)**2+(width_l/2)**2)))
                    popt_pv,_ = scipy.optimize.curve_fit(pv, x_int, y_int, p0 = [np.max(y_int),np.max(y_int), lam, x_int[-1]-x_int[0] , x_int[-1]-x_int[0], p] ,bounds = param_bounds, maxfev=5000000000)
                    amp_g = popt_pv[0]
                    amp_l = popt_pv[1]
                    cen = popt_pv[2]
                    width_g = popt_pv[3]
                    width_l = popt_pv[4]
                    p = popt_pv[5]
                    def Gauss(x_int, amp_g, cen, width_g, p ):
                        return p*((amp_g/np.sqrt(2*np.pi*width_g**2))*np.exp(-((x_int-cen)**2)/(2*width_g**2)))
                    def Lorentz(x_int, amp_l, cen, width_l, p ):
                        return (1-p)*((amp_l/np.pi)*((width_l/2)**2/((x_int-cen)**2+(width_l/2)**2)))
                    x_int_smooth = np.arange(x_int[0],x_int[-1],0.0001)
                    PV = pv(x_int_smooth, *popt_pv)
                    C = 4.700e+21
                    FWHM_g=popt_pv[3]*np.sqrt(4*np.log(2))
                    #N_e = (C)/((popt_pv[2]/10)**2*(10**(-8))*(Ek-Ei))
                    gaussian = Gauss(x_int_smooth, amp_g, cen, width_g, p )
                    lorentzian = Lorentz(x_int_smooth, amp_l, cen, width_l, p )
                    def T_e(FWHM_l):
                        return np.exp((-a1+np.sqrt(np.abs((a1**2)-(4*a2*(a0-np.log(FWHM_l))))))/(2*a2))
                    plt.plot(x_int, y_int, label = 'Corrected observation')
                    plt.plot(x_int_smooth, gaussian, label = 'Gaussian component')
                    plt.plot(x_int_smooth, lorentzian, label = 'lorentzian component')
                    plt.plot(x_int_smooth, PV, label = 'PV Fit')
                    #plt.text(eq1, {'color': 'C2', 'fontsize': 18}, va="top", ha="right")
                    plt.title("Pesudo Voigt Fitting at $\lambda = $"+str(lam))
                    plt.xlabel(r'$\lambda \longrightarrow$')
                    plt.ylabel(r'$I \rightarrow$')
                    plt.legend(loc = 'upper left')
                    plt.show()
                    while True:
                        sel1 = input("do you want to change following parameters    :    [Input y] \n amp_g = "+str(popt_pv[0])+" \n amp_l = "+str(popt_pv[1])+"\n cen = "+str(popt_pv[2])+" \n width_g = "+str(popt_pv[3])+" \n width_l = "+str(popt_pv[4])+" \n p = "+str(popt_pv[5])+"\n \n 'y' for yes, 's' for save             Input response: ")
                        amp_g = popt_pv[0]
                        amp_l = popt_pv[1]
                        cen = popt_pv[2]
                        width_g = popt_pv[3]
                        width_l = popt_pv[4]
                        p = popt_pv[5]
                        FWHM_l=popt_pv[4]
                        FWHM_g=popt_pv[3]*np.sqrt(4*np.log(2))
                        IV,_ = quad(pv, x_int[0],x_int[-1], args = ( amp_g, amp_l, cen, width_g, width_l, p ))
                        IG,_ = quad(Gauss, x_int[0],x_int[-1], args = (amp_g, cen, width_g, p))
                        IL,_ = quad(Lorentz, x_int[0],x_int[-1], args = (amp_l, cen, width_l, p ))
                        N_e = (FWHM_l*(10**16))/(2*Stark_w)
                        output_data = lam, cen, p, amp_g, amp_l, width_g, width_l, N_e, IG, IL, IV
                        if sel1 == 'y':
                            k = input("press the index of parameter to be changed(1 for amp_g, 2 for amp_l and so on... ...):")
                            if k == '1':
                                param = input("input the new value of parameter:  ")
                                opt_param = popt_pv
                                opt_param[0] = param
                                param_bounds1=([float(param) - 1.00000000,0,0,0,0,0],[float(param) + 1.00000000,np.inf,np.inf,np.inf,np.inf,1])
                                popt_pv1,_ = scipy.optimize.curve_fit(pv, x_int, y_int, p0 = [float(param),amp_l, cen, width_g , width_l, p] ,bounds =param_bounds1, maxfev=5000000000)
                                PV = pv(x_int_smooth, *popt_pv1)
                                amp_g = popt_pv1[0]
                                amp_l = popt_pv1[1]
                                cen = popt_pv1[2]
                                width_g = popt_pv1[3]
                                width_l = popt_pv1[4]
                                p = popt_pv1[5]
                                gaussian = Gauss(x_int_smooth, amp_g, cen, width_g, p)
                                lorentzian = Lorentz(x_int_smooth, amp_l, cen, width_l, p)
                                plt.plot(x_int, y_int, label = 'Corrected observation')
                                plt.plot(x_int_smooth, gaussian, label = 'Gaussian component')
                                plt.plot(x_int_smooth, lorentzian, label = 'lorentzian component')
                                plt.plot(x_int_smooth, PV, label = 'PV Fit')
                                plt.title("Pesudo Voigt Fitting Compoments at $\lambda = $"+str(lam))
                                plt.xlabel(r'$\lambda \longrightarrow$')
                                plt.ylabel(r'$I \rightarrow$')
                                plt.legend(loc = 'upper left')
                                plt.show()
                                popt_pv = popt_pv1
                                continue
                            elif k == '2':
                                param = input("input the new value of parameter:  ")
                                opt_param = popt_pv
                                opt_param[1] = param
                                param_bounds1=([0, float(param) - 1.00000000,0,0,0,0],[np.inf, float(param) + 1.00000000,np.inf,np.inf,np.inf,1])
                                popt_pv1,_ = scipy.optimize.curve_fit(pv, x_int, y_int, p0 = [amp_g,float(param), cen, width_g , width_l, p] ,bounds =param_bounds1, maxfev=5000000000)
                                PV = pv(x_int_smooth, *popt_pv1)
                                amp_g = popt_pv1[0]
                                amp_l = popt_pv1[1]
                                cen = popt_pv1[2]
                                width_g = popt_pv1[3]
                                width_l = popt_pv1[4]
                                p = popt_pv1[5]
                                gaussian = Gauss(x_int_smooth, amp_g, cen, width_g, p)
                                lorentzian = Lorentz(x_int_smooth, amp_l, cen, width_l, p)
                                plt.plot(x_int, y_int, label = 'Corrected observation')
                                plt.plot(x_int_smooth, gaussian, label = 'Gaussian component')
                                plt.plot(x_int_smooth, lorentzian, label = 'lorentzian component')
                                plt.plot(x_int_smooth, PV, label = 'PV Fit')
                                plt.title("Pesudo Voigt Fitting Compoments at $\lambda = $"+str(lam))
                                plt.xlabel(r'$\lambda \longrightarrow$')
                                plt.ylabel(r'$I \rightarrow$')
                                plt.legend(loc = 'upper left')
                                plt.show()
                                popt_pv = popt_pv1
                                continue
                            elif k == '3':
                                param = input("input the new value of parameter:  ")
                                opt_param = popt_pv
                                opt_param[2] = param
                                param_bounds1=([0,0,float(param) - .300000000,0,0,0],[np.inf,np.inf , float(param) + .300000000, np.inf,np.inf,1])
                                popt_pv1,_ = scipy.optimize.curve_fit(pv, x_int, y_int, p0 = [amp_g,amp_l, float(param), width_g , width_l, p] ,bounds =param_bounds1, maxfev=5000000000)
                                PV = pv(x_int_smooth, *popt_pv1)
                                amp_g = popt_pv1[0]
                                amp_l = popt_pv1[1]
                                cen = popt_pv1[2]
                                width_g = popt_pv1[3]
                                width_l = popt_pv1[4]
                                p = popt_pv1[5]
                                gaussian = Gauss(x_int_smooth, amp_g, cen, width_g, p)
                                lorentzian = Lorentz(x_int_smooth, amp_l, cen, width_l, p)
                                plt.plot(x_int, y_int, label = 'Corrected observation')
                                plt.plot(x_int_smooth, gaussian, label = 'Gaussian component')
                                plt.plot(x_int_smooth, lorentzian, label = 'lorentzian component')
                                plt.plot(x_int_smooth, PV, label = 'PV Fit')
                                plt.title("Pesudo Voigt Fitting Compoments at $\lambda = $"+str(lam))
                                plt.xlabel(r'$\lambda \longrightarrow$')
                                plt.ylabel(r'$I \rightarrow$')
                                plt.legend(loc = 'upper left')
                                plt.show()
                                popt_pv = popt_pv1
                                continue
                            elif k == '4':
                                param = input("input the new value of parameter:  ")
                                opt_param = popt_pv
                                opt_param[0] = param
                                param_bounds1=([0,0,0,float(param) - .100000000,0,0],[np.inf,np.inf,np.inf,float(param) + .100000000,np.inf,1])
                                popt_pv1,_ = scipy.optimize.curve_fit(pv, x_int, y_int, p0 = [amp_g,amp_l, cen, float(param), width_l, p] ,bounds =param_bounds1, maxfev=5000000000)
                                PV = pv(x_int_smooth, *popt_pv1)
                                amp_g = popt_pv1[0]
                                amp_l = popt_pv1[1]
                                cen = popt_pv1[2]
                                width_g = popt_pv1[3]
                                width_l = popt_pv1[4]
                                p = popt_pv1[5]
                                gaussian = Gauss(x_int_smooth, amp_g, cen, width_g, p)
                                lorentzian = Lorentz(x_int_smooth, amp_l, cen, width_l, p)
                                plt.plot(x_int, y_int, label = 'Corrected observation')
                                plt.plot(x_int_smooth, gaussian, label = 'Gaussian component')
                                plt.plot(x_int_smooth, lorentzian, label = 'lorentzian component')
                                plt.plot(x_int_smooth, PV, label = 'PV Fit')
                                plt.title("Pesudo Voigt Fitting Compoments at $\lambda = $"+str(lam))
                                plt.xlabel(r'$\lambda \longrightarrow$')
                                plt.ylabel(r'$I \rightarrow$')
                                plt.legend(loc = 'upper left')
                                plt.show()
                                popt_pv = popt_pv1
                                continue
                            elif k == '5':
                                param = input("input the new value of parameter:  ")
                                opt_param[0] = param
                                param_bounds1=([0,0,0,0,float(param) - .100000000,0],[np.inf,np.inf,np.inf,np.inf,float(param) + .100000000,1])
                                popt_pv1,_ = scipy.optimize.curve_fit(pv, x_int, y_int, p0 = [amp_g,amp_l, cen, width_g , float(param), p] ,bounds =param_bounds1, maxfev=5000000000)
                                PV = pv(x_int_smooth, *popt_pv1)
                                amp_g = popt_pv1[0]
                                amp_l = popt_pv1[1]
                                cen = popt_pv1[2]
                                width_g = popt_pv1[3]
                                width_l = popt_pv1[4]
                                p = popt_pv1[5]
                                gaussian = Gauss(x_int_smooth, amp_g, cen, width_g, p)
                                lorentzian = Lorentz(x_int_smooth, amp_l, cen, width_l, p)
                                plt.plot(x_int, y_int, label = 'Corrected observation')
                                plt.plot(x_int_smooth, gaussian, label = 'Gaussian component')
                                plt.plot(x_int_smooth, lorentzian, label = 'lorentzian component')
                                plt.plot(x_int_smooth, PV, label = 'PV Fit')
                                plt.title("Pesudo Voigt Fitting Compoments at $\lambda = $"+str(lam))
                                plt.xlabel(r'$\lambda \longrightarrow$')
                                plt.ylabel(r'$I \rightarrow$')
                                plt.legend(loc = 'upper left')
                                plt.show()
                                popt_pv = popt_pv1
                                continue
                            elif k == '6':
                                param = input("input the new value of parameter:  ")
                                opt_param = popt_pv
                                opt_param[0] = param
                                param_bounds1=([0,0,0,0,0,float(param) - 0.100000000],[np.inf,np.inf,np.inf,np.inf, np.inf,float(param) + .100000000])
                                popt_pv1,_ = scipy.optimize.curve_fit(pv, x_int, y_int, p0 = [amp_g,amp_l, cen, width_g , width_l, float(param)] ,bounds =param_bounds1, maxfev=5000000000)
                                PV = pv(x_int_smooth, *popt_pv1)
                                amp_g = popt_pv1[0]
                                amp_l = popt_pv1[1]
                                cen = popt_pv1[2]
                                width_g = popt_pv1[3]
                                width_l = popt_pv1[4]
                                p = popt_pv1[5]
                                gaussian = Gauss(x_int_smooth, amp_g, cen, width_g, p)
                                lorentzian = Lorentz(x_int_smooth, amp_l, cen, width_l, p)
                                plt.plot(x_int, y_int, label = 'Corrected observation')
                                plt.plot(x_int_smooth, gaussian, label = 'Gaussian component')
                                plt.plot(x_int_smooth, lorentzian, label = 'lorentzian component')
                                plt.plot(x_int_smooth, PV, label = 'PV Fit')
                                plt.title("Pesudo Voigt Fitting Compoments at $\lambda = $"+str(lam))
                                plt.xlabel(r'$\lambda \longrightarrow$')
                                plt.ylabel(r'$I \rightarrow$')
                                plt.legend(loc = 'upper left')
                                plt.show()
                                popt_pv = popt_pv1
                                continue
                            else:
                                break
                        elif sel1 == 's':
                            with open('Out/output_data'+str(lam)+'_'+str(sample_index)+'.txt', 'w') as f3:
                                f3.write(str(output_data)) 
                            break          
        #data = lam, cen, p, amp_g, amp_l, width_g, width_l, N_e
        #p0 = [np.max(y_int),np.max(y_int), lam, x_int[-1]-x_int[0]] , x_int[-1]-x_int[0], p] ,
        
    
