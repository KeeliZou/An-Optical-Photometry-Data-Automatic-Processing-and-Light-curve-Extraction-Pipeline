from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import os
import warnings
import glob

#super_bias:超级flat路径
#super_flat:超级bias路径
#band_info:波段信息
#raw_data_file_path:原始图像路径
#reduced_data_file_path:
#此函数用于对图像进行本底平场改正工作
def bfc_raw(super_bias,super_flat,band_info,raw_data_file_path,reduced_data_file_path):
    #处理图像输出路径
    outdir = reduced_data_file_path
    print("\033[1;31m Correctting Images!(data-bias)/(flat-bias)\033[0m")
    for j in range(len(band_info)):
        sfiles = glob.glob(raw_data_file_path[j])
        sfiles.sort()
        print("compiling "+band_info[j]+" band images!")
        for i,ifile in enumerate(sfiles):
            #print("reducing "+band_info[j]+" band (debias,flat-fielding, and flipping):",i + 1, len(sfiles), ifile)
            indir, infile = os.path.split(ifile)
            rootname,_ = os.path.splitext(infile)
            progress = (i+1) / len(sfiles) * 100
            print(f"working process: [{int(progress):3d}%],",end='\r')

            outfile = os.path.join(outdir,"Cor_" + rootname + '.fits')
            head = fits.getheader(ifile,output_verifystr="silentfix")
            data = fits.getdata(ifile)
            if i ==0:
                flat = fits.getdata(super_flat[j])
                bias = fits.getdata(super_bias)
            data = (data - bias) / flat
            fits.writeto(outfile,data,header=head,overwrite = True,output_verify = "silentfix")
        print(band_info[j]+" band completed!")




#此函数用于修剪图像边缘，去掉噪声较大的区域
def cor_cut(cor_data_file_path,cut_file_path,band_info,edge_cutx,edge_cuty):
    outdir = cut_file_path
    print("\033[1;31m Cutting Images!\033[0m")
    for j in range(len(band_info)):
        sfiles = glob.glob(cor_data_file_path[j])
        sfiles.sort()
        print("cutting "+band_info[j]+" band images!")
        for i,ifile in enumerate(sfiles):
            data = fits.getdata(ifile)
            head = fits.getheader(ifile)
            naxis1 = head['NAXIS1']
            naxis2 = head['NAXIS2']
            data = data[(edge_cuty-1):(naxis2-edge_cuty-1),(edge_cutx-1):(naxis1-edge_cutx-1)]
            head['NAXIS1'] = head['NAXIS1'] - 2*edge_cutx
            head['NAXIS2'] = head['NAXIS2'] - 2*edge_cuty

            progress = (i+1) / len(sfiles) * 100
            print(f"working process: [{int(progress):3d}%],",end='\r')

            indir, infile = os.path.split(ifile)
            rootname,_ = os.path.splitext(infile)
            outfile = os.path.join(outdir,"Cut_" + rootname + '.fits')

            fits.writeto(outfile,data,header=head,overwrite=True,output_verify = "silentfix")
        print(band_info[j]+" band completed!")
