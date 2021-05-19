'''
Constant Parameters
'''

'''
1) Parameters of short time fourier analysis
'''
Fs_ref = 16e3 #1.1
M_ref = 512 #1.2
Mo_ref = 0.75 * M_ref #1.3


'''
2) Parameters of noise spectrum estimate
'''
w = 1
alpha_s_ref = 0.9
Nwin = 8
Vwin = 15
delta_s = 1.67
Bmin = 1.66
delta_y = 4.6
delta_yt = 3
alpha_d_ref = 0.85


'''
3) Parameters of a Priori Probability for Signal Absense Estimate
'''
alpha_xi_ref = 0.7
w_xi_local = 1
w_xi_global=15
f_u=10e3
f_l=50
P_min=0.005
xi_lu_dB=-5
xi_ll_dB=-10
xi_gu_dB=-5
xi_gl_dB=-10
xi_fu_dB=-5
xi_mu_dB=10
xi_ml_dB=0
q_max=0.998

'''
4) Parameters of 'Decision-Directed' a Priori SNR Estimate
'''
alpha_eta_ref=0.95
eta_min_dB=-18

'''
5) Flags
'''
broad_flag=1
tone_flag=1
nonstat='medium'


'''
'''
alpha_d_long = 0.99
eta_min = 10 ** (eta_min_dB / 10)
G_f=eta_min**0.5