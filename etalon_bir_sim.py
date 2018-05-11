#conf python
import numpy as np
import matplotlib.pyplot as plt
import dos_lib as dl #MY LIBRARY
from matplotlib.widgets import Slider, Button, RadioButtons

print('######################### Etalon widget tool ##########################')
print('## David Orozco (oroco@iaa.es) and Francisco Bailén (fbailen@iaa.es) ##')
print('## Instituto de Astrofísica de Granada (IAA-CSIC)                    ##')
print('## Current: Version 1.0 released on 6th April, 2018                  ##')
print('## Based on procedures in dos_lib.py Version 1.0                     ##')
print('#######################################################################')

#History
#added rotation slider April 2018.

#TO BE DONE
# Put the posibility of two etalons
# Introduce posibility of changing no.
# Introduce posibility of changing the frequency nad frequency range.

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 8}
plt.rc('font', **font)
plt.rc('xtick',labelsize=6)
plt.rc('ytick',labelsize=6)

lcentral = 525.6 #in nm
rangel = 0.2
lan = np.arange(lcentral-rangel, lcentral+rangel, 0.001)*1e-9

dl.Reflectance = 0.92 #Reflectance
dl.no=2.293 ##Ordinary refractive index
dl.n3=2.2 #n3
dl.h=251.63e-6 #Etalon width
dl.Absorptance=0 #Absorptance

s = dl.FPM(np.radians(0.50), lan , np.radians(0), np.radians(0.),doprint='d')

fig = plt.figure(figsize=(8, 7))
grid = plt.GridSpec(2, 4, hspace=0.2, wspace=0.3)
plt.subplots_adjust(left=0.25, bottom=0.45, top=0.98)

main_abx = fig.add_subplot(grid[0,:2])
main_cdx = fig.add_subplot(grid[0,2:])
a_prof = fig.add_subplot(grid[1, 0])
b_prof = fig.add_subplot(grid[1, 1])
c_prof = fig.add_subplot(grid[1, 2])
d_prof = fig.add_subplot(grid[1, 3])

#axis limits
ssmax = np.maximum(s[0,0,:]-s[1,0,:],s[0,0,:]+s[1,0,:])
ssmin = np.minimum(s[0,0,:]-s[1,0,:],s[0,0,:]+s[1,0,:])
main_abx_max = np.amax(ssmax)
main_abx_min = np.amin(ssmin)
main_abx.axis([-rangel, rangel, main_abx_min, main_abx_max])
ssmax = np.maximum(s[2,2,:]-s[2,3,:],s[2,2,:]+s[2,3,:])
ssmin = np.minimum(s[2,2,:]-s[2,3,:],s[2,2,:]+s[2,3,:])
main_cdx_max = np.amax(ssmax)
main_cdx_min = np.amin(ssmin)
main_cdx.axis([-rangel, rangel, main_cdx_min, main_cdx_max])
a_max = np.amax(s[0,0,:])
a_min = np.amin(s[0,0,:])
a_prof.axis([-rangel, rangel, a_min, a_max])
b_max = np.amax(s[1,0,:])
b_min = np.amin(s[1,0,:])
b_prof.axis([-rangel, rangel, b_min, b_max])
c_max = np.amax(s[2,2,:])
c_min = np.amin(s[2,2,:])
c_prof.axis([-rangel, rangel, c_min, c_max])
d_max = np.amax(s[2,3,:])
d_min = np.amin(s[2,3,:])
d_prof.axis([-rangel, rangel, d_min, d_max])

l11, = main_abx.semilogy(lan*1e9-lcentral, s[0,0,:]+s[1,0,:], lw=1, color='red',label='a+b')
l21, = main_abx.semilogy(lan*1e9-lcentral, s[0,0,:]-s[1,0,:], lw=1, color='blue',label='a-b')
main_abx.legend(bbox_to_anchor=(-0.3, 1), loc=1, borderaxespad=0.)
l12, = main_cdx.plot(lan*1e9-lcentral, s[2,2,:]+s[2,3,:], lw=1, color='red',label='c+d')
l22, = main_cdx.plot(lan*1e9-lcentral, s[2,2,:]-s[2,3,:], lw=1, color='blue',label='c-d')
main_cdx.legend(bbox_to_anchor=(1.3, 1), loc=1, borderaxespad=0.)
l3, = a_prof.plot(lan*1e9-lcentral, s[0,0,:], lw=1, color='green',label='a')
l4, = b_prof.plot(lan*1e9-lcentral, s[1,0,:], lw=1, color='black',label='b')
l5, = c_prof.plot(lan*1e9-lcentral, s[2,2,:], lw=1, color='green',label='c')
l6, = d_prof.plot(lan*1e9-lcentral, s[2,3,:], lw=1, color='black',label='d')

#fig, ax = plt.subplots()
#plt.subplots_adjust(left=0.25, bottom=0.35, top=0.98)
#l1, = plt.semilogy(lan-617.234, s[0][0][:]+s[1][0][:], lw=1, color='red',label='a+b')
#plt.axis([-0.5, 0.5, 1e-3, 1.+1.e-5])
#l2, = plt.semilogy(lan-617.234, s[0][0][:]-s[1][0][:], lw=1, color='blue',label='a-b')
#l3, = plt.semilogy(lan-617.234, s[0][0][:], lw=1, color='green',linestyle=':',label='a')
#l4, = plt.semilogy(lan-617.234, s[1][0][:], lw=1, color='black',linestyle=':',label='b')

axcolor = 'lightgoldenrodyellow'
ax_volt =    plt.axes([0.25, 0.28, 0.65, 0.03], facecolor=axcolor)
ax_rot =     plt.axes([0.25, 0.24, 0.65, 0.03], facecolor=axcolor)
ax_theta_i = plt.axes([0.25, 0.20, 0.65, 0.03], facecolor=axcolor)
ax_h =       plt.axes([0.25, 0.16, 0.65, 0.03], facecolor=axcolor)
ax_theta_3 = plt.axes([0.25, 0.12, 0.65, 0.03], facecolor=axcolor)
ax_r =       plt.axes([0.25, 0.08, 0.65, 0.03], facecolor=axcolor)
ax_ne =      plt.axes([0.25, 0.04, 0.65, 0.03], facecolor=axcolor)

s_volt =    Slider(ax_volt   , 'Voltage (volts).', -5000.0, 5000.0, valinit=0.,valfmt='%0.4f')
s_rot =     Slider(ax_rot    , 'Rotation angle (deg).', 0, 90.0, valinit=0.,valfmt='%0.4f')
s_theta_i = Slider(ax_theta_i, 'Theta_i (degree)', 0., 1.0, valinit=0.55,valfmt='%0.4f')
s_h =       Slider(ax_h      , 'H (mm)', 0.1, 1.0, valinit=0.251,valfmt='%0.4f')
s_theta_3 = Slider(ax_theta_3, 'Theta_3 (degree)', 0., 90.0, valinit=0.0,valfmt='%0.4f')
s_r =       Slider(ax_r      , 'Reflt.', 0.05, 1, valinit=0.92,valfmt='%0.3f')
s_ne =      Slider(ax_ne     , 'Refraction index (e).', 2.0, 3.0, valinit=2.2,valfmt='%0.4f')


def update(val):
    dl.Volt = s_volt.val #Etalon voltage
    dl.h = s_h.val*1e-3 #Etalon width
    dl.Reflectance = s_r.val #Reflectance
    dl.n3 = s_ne.val #Reflectance
    s = dl.FPM(np.radians(s_theta_i.val), lan , np.radians(s_theta_3.val), np.radians(s_rot.val),doprint='false')
    ssmax = np.maximum(s[0,0,:]-s[1,0,:],s[0,0,:]+s[1,0,:])
    ssmin = np.minimum(s[0,0,:]-s[1,0,:],s[0,0,:]+s[1,0,:])
    main_abx_max = np.amax(ssmax)
    main_abx_min = np.amin(ssmin)
    main_abx.axis([-rangel, rangel, main_abx_min, main_abx_max])
    ssmax = np.maximum(s[2,2,:]-s[2,3,:],s[2,2,:]+s[2,3,:])
    ssmin = np.minimum(s[2,2,:]-s[2,3,:],s[2,2,:]+s[2,3,:])
    main_cdx_max = np.amax(ssmax)
    main_cdx_min = np.amin(ssmin)
    main_cdx.set_ylim([main_cdx_min, main_cdx_max])
    a_max = np.amax(s[0,0,:])
    a_min = np.amin(s[0,0,:])
    a_prof.set_ylim([a_min, a_max])
    b_max = np.amax(s[1,0,:])
    b_min = np.amin(s[1,0,:])
    b_prof.set_ylim([b_min, b_max])
    c_max = np.amax(s[2,2,:])
    c_min = np.amin(s[2,2,:])
    c_prof.set_ylim([c_min, c_max])
    d_max = np.amax(s[2,3,:])
    d_min = np.amin(s[2,3,:])
    d_prof.set_ylim([d_min, d_max])
    l11.set_ydata(s[0,0,:]+s[1,0,:])
    l21.set_ydata(s[0,0,:]-s[1,0,:])
    l12.set_ydata(s[2,2,:]+s[2,3,:])
    l22.set_ydata(s[2,2,:]-s[2,3,:])
    l3.set_ydata(s[0,0,:])
    l4.set_ydata(s[1,0,:])
    l5.set_ydata(s[2,2,:])
    l6.set_ydata(s[2,3,:])
    fig.canvas.draw_idle()
s_volt.on_changed(update)
s_rot.on_changed(update)
s_theta_i.on_changed(update)
s_h.on_changed(update)
s_theta_3.on_changed(update)
s_r.on_changed(update)
s_ne.on_changed(update)

resetax = plt.axes([0.03, 0.03, 0.05, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975',)

def reset(event):
    s_volt.reset()
    s_rot.reset()
    s_theta_i.reset()
    s_h.reset()
    s_theta_3.reset()
    s_r.reset()
    s_ne.reset()
button.on_clicked(reset)

# rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
# radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)
# def colorfunc(label):
#     l1.set_color(label)
#     fig.canvas.draw_idle()
# radio.on_clicked(colorfunc)

plt.show()

raise SystemExit
