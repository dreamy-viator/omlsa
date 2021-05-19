from scipy.io import wavfile
from params import *
from utils import lnshift, mat_hanning

import math
import numpy as np

import os

def omlsa(fin, fout):
    Fs, y_origin = wavfile.read(fin)
    N = len(y_origin)
    if Fs is not Fs_ref:
        M = 2 ** round(math.log2(Fs / Fs_ref * M_ref))
        Mo = Mo_ref / M_ref * M
        alpha_s = alpha_s_ref ** (M_ref / M * Fs / Fs_ref)
        alpha_d = alpha_d_ref ** (M_ref / M * Fs / Fs_ref)
        alpha_eta = alpha_eta_ref ** (M_ref / M * Fs / Fs_ref)
        alpha_xi = alpha_xi_ref ** (M_ref / M * Fs / Fs_ref)
    else:
        M = M_ref
        Mo = Mo_ref
        alpha_s = alpha_s_ref
        alpha_d = alpha_d_ref
        alpha_eta = alpha_eta_ref
        alpha_xi = alpha_xi_ref

    win = np.hamming(M)
    win2 = np.power(win, 2)
    Mno = int(M - Mo)
    W0 = win2[:Mno]
    for k in range(Mno, M - 1, Mno):
        swin2 = lnshift(win2, k)
        W0 = W0 + swin2[:Mno]

    W0 = np.mean(W0) ** 0.5
    win = win / W0
    Cwin = np.sum(np.power(win, 2)) ** 0.5
    win = win / Cwin

    Nframes = int((N - Mo) / Mno)
    out = np.zeros(M)
    b = mat_hanning(2 * w + 1)
    b = b / np.sum(b)
    b_xi_local = mat_hanning(2 * w_xi_local + 1)
    b_xi_local = b_xi_local / np.sum(b_xi_local)
    b_xi_global = mat_hanning(2 * w_xi_global + 1)
    b_xi_global = b_xi_global / np.sum(b_xi_global)

    l_mod_lswitch = 0
    M21 = int(M / 2 + 1)
    k_u = round(f_u / Fs * M + 1)
    k_l = round(f_l / Fs * M + 1)
    k_u = min(k_u, M21)
    k2_local = round(500 / Fs * M + 1)
    k3_local = round(3500 / Fs * M + 1)
    eta_2term = 1
    xi = 0
    xi_frame = 0
    Ncount = round(Nframes / 10)

    l_fnz = 0 # python index starts with 0
    fnz_flag = 0
    zero_thres = 1e-10

    for l in range(Nframes):
        with open(fin, 'rb') as finID:
            if l == 0:
                y = y_origin[:M] / (2 ** 15)

            else:
                y0 = y_origin[M + (Mno * (l - 1)): M + (Mno * l)] / (2 ** 15)
                y = np.hstack((y[Mno:M], y0))


            if (not fnz_flag and (abs(y[0]) > zero_thres)) or (fnz_flag and np.any(abs(y) > zero_thres)):
                fnz_flag = 1
                '''1. Short Time Fourier Analysis'''
                Y = np.fft.fft(win * y)
                Ya2 = abs(Y[:M21]) ** 2
                if l is l_fnz:
                    lambda_d = Ya2

                temp = np.full(lambda_d.shape, 1e-10)
                gamma = Ya2 / np.maximum(lambda_d, temp)
                eta = alpha_eta * eta_2term + (1 - alpha_eta) * np.maximum(gamma - 1, 0)
                eta = np.maximum(eta, eta_min)
                v = gamma * eta / (1 + eta)

                '''2.1 smooth over freqeuncy'''
                Sf = np.convolve(b, Ya2)
                Sf = Sf[w: M21 + w]

                if l is l_fnz:
                    Sy = Ya2
                    S = Sf
                    St = Sf
                    lambda_dav = Ya2
                else:
                    S = alpha_s * S + (1 - alpha_s) * Sf

                if l < 14 + l_fnz:
                    Smin = S
                    SMact = S
                else:
                    Smin = np.minimum(Smin, S)
                    SMact = np.minimum(SMact, S)

                '''Local Minima Search'''

                I_f = np.less(Ya2, delta_y * Bmin * Smin).astype(int) & np.less(S, delta_s * Bmin * Smin).astype(int)
                conv_I = np.convolve(b, I_f)
                conv_I = conv_I[w: M21 + w]
                Sft = St.copy()
                idx = conv_I.nonzero()
                if not np.all(idx==0):
                    if w:
                        conv_Y = np.convolve(b, I_f * Ya2)
                        conv_Y = conv_Y[w:M21+w]
                        Sft[idx] = conv_Y[idx] / conv_I[idx]
                    else:
                        Sft[idx] = Ya2[idx]

                if l < 14 + l_fnz:
                    St = S
                    Smint = St
                    SMactt = St
                else:
                    St = alpha_s * St + (1 - alpha_s) * Sft
                    Smint = np.minimum(Smint, St)
                    SMactt = np.minimum(SMactt, St)

                qhat = np.ones(M21)
                phat = np.zeros(M21)

                if nonstat == "low":
                    gamma_mint = Ya2 / Bmin / np.maximum(Smin, 1e-10)
                    zetat = S / Bmin / np.maximum(Smin, 1e-10)
                else:
                    gamma_mint = Ya2 / Bmin / np.maximum(Smint, 1e-10)
                    zetat = S / Bmin / np.maximum(Smint, 1e-10)

                cond = np.greater(gamma_mint, 1) & np.less(gamma_mint, delta_yt) & np.less(zetat, delta_s)
                idx = cond.nonzero()
                qhat[idx] = (delta_yt - gamma_mint[idx]) / (delta_yt - 1)
                phat[idx] = 1 / (1 + qhat[idx]) / (1 - qhat[idx]) * (1 + eta[idx]) * np.exp(-v[idx])
                phat[np.greater_equal(gamma_mint, delta_yt) | np.greater_equal(zetat, delta_s)] = 1
                alpha_dt = alpha_d + (1 - alpha_d) * phat
                lambda_dav = alpha_dt * lambda_dav + (1 - alpha_dt) * Ya2

                if l < 14 + l_fnz:
                    lambda_dav_long = lambda_dav
                else:
                    alpha_dt_long = alpha_d_long + (1 - alpha_d_long) * phat
                    lambda_dav_long = alpha_dt_long * lambda_dav_long + (1 - alpha_dt_long) * Ya2

                l_mod_lswitch = l_mod_lswitch + 1
                if l_mod_lswitch is Vwin:
                    l_mod_lswitch = 0
                    if l == (Vwin - 1 + l_fnz):
                        SW = np.transpose(np.tile(S, [Nwin, 1]))
                        SWt = np.transpose(np.tile(St, [Nwin, 1]))
                    else:
                        SW = np.column_stack((SW[:, 1:Nwin], SMact))
                        Smin = np.amin(SW, axis=1)
                        SMact = S.copy()
                        SWt = np.column_stack((SWt[:, 1:Nwin], SMactt))
                        Smint = np.amin(SWt, axis=1)
                        SMactt = St.copy()




                '''2.4 Noise Spectrum Estimate'''
                #TODO:

if __name__ == "__main__":

    omlsa("./clean_fileid_0.wav", "./clean_fileid_0(out).wav")